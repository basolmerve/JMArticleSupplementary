library(parallel)
library(JM)
library(ggplot2)
library(dplyr)
library(magrittr)
library(readr)
source("helper_functions.R")

# Verilerin R'a okutulması
#
# longitudinal verileri.
pdData <- read_delim(file = "longitudinal60.txt", delim = "\t", col_names = TRUE, 
                     locale = locale(decimal_mark = "."))

# Ek bilgilerin yer aldığı veri seti
tmp <- read_delim("ekveriler.csv", delim = ",", locale = locale(decimal_mark = "."))
pdData <- pdData %>%
  full_join(tmp)

# String olan değişkenler "character" sınıfında okutuluyor. Bu değişkenler "factor" sınıfına dönüştürülür.
# JMbayes kütüphanesi "data.frame" sınıfına göre yazılmış. "tibble" sınıfını bu nedenle düzenledikten sonra 
# "data.frame" olarak kaydediyoruz.
dataLong <- pdData %>%
  mutate_if(is.character, as.factor) %>%
  mutate(preHD = factor(preHD, levels = c(0, 1), labels = c("No", "Yes")),
         highlow = factor(highlow, levels = c(1, 2), labels = c("Low", "High"))) %>%
  as.data.frame

dataLong <- dataLong %>%
  dplyr::select(id, surv, event, time, ALB, age, highlow, peritonitrate, preHD, ill_count, gender) %>%
  filter(complete.cases(.)) %>%
  filter(!(id %in% c(42, 160, 245, 355)))  # 4 gözlem veri setinden çıkartıldı. Test amaçlı kullanılacak.

# Fit Joint Models
# Seçilen değişkenlerde eksik olan gözlemler silindikten sonra "lme" modeli kuruluyor. "longitudinal" parça için.
lmeFit <- dataLong %>%
  lme(ALB ~ age + highlow + peritonitrate + time,
      data = .,
      random = ~ time | id)

# Seçilen değişkenlerde eksik olan gözlemler silindikten sonra "coxph" modeli kuruluyor. "survival" parça için.
# survFit <- train_data %>%
#   filter(!duplicated(id)) %>%
#   coxph(Surv(surv, event) ~ age + gender + ill_count + BKI + highlow + peritonitrate, data = ., x = TRUE)

survFit <- dataLong %>%
  filter(!duplicated(id)) %>%
  coxph(Surv(surv, event) ~ age + preHD + ill_count + peritonitrate, 
        data = ., x = TRUE)

# jmControl <- list(tol1 = 1e-7, tol2 = 1e-7)
# MODEL
jointFit <- jointModel(lmeObject = lmeFit, survObject = survFit, timeVar = "time", 
                       method = "weibull-PH-aGH")

n <- c(100, 400, 800)  # number of samples
fupTimeMax <- c(60, 120)  # Maximum followup time. Longitudinal measurements are taken at 6-months intervals within 0 and maximum of followup time.
B = 200  # Number of simulated data sets in each scenario
scenarios <- expand.grid(nSamples = n, maxFUpTime = fupTimeMax, Repeat = 1:B)  

set.seed(23)
seeds <- sample(1:10000, nrow(scenarios), replace = FALSE)
scenarios[["seed"]] <- seeds
scenarios[["SimID"]] <- 1:nrow(scenarios)

runSim <- function(simID = NULL, .model, .data = NULL, plist = NULL, save = TRUE, 
                   path = file.path(getwd(), "results/JM"), landmarkTimes = NULL, landMarkStep = 12,
                   deltaT = NULL, horizonTime = NULL, nIterMax = 5, ...){
  
  library(JM)
  #library("splines")
  library(dplyr)
  library(magrittr)
  
  # dt ve horiz ikisi birden NULL olamaz. Bu durumda deltaT değeri 1 olarak alınıyor.
  if (all(is.null(horizonTime), is.null(deltaT))){
    deltaT <- 1 + 1e-7
  } else if (is.null(horizonTime) & !is.null(deltaT)){
    deltaT <- deltaT + 1e-7
  }
  
  # dt ve horiz ikisi birden tanımlanmış ise dt dikkate alınır.
  if (all(!is.null(horizonTime), !is.null(deltaT))){
    horizonTime <- NULL
    deltaT <- deltaT + 1e-7
  }
  
  plist <- plist[plist[["SimID"]] == simID, ]
  if (is.null(landmarkTimes)){
    landmarkTimes = seq(landMarkStep, plist[["maxFUpTime"]], by = landMarkStep)
  }
  
  scen <- suppressWarnings(data.frame(plist, Tstart = landmarkTimes))
  
  measures <- data.frame(scen, Iteration = 1, CPUtimeData = NA, 
                         AUC = NA, BS = NA)
  
  iter <- 1
  while (iter <= nIterMax){
    cputime_simData <- system.time({
      # simulatedData <- try({
      #   simulateJM(nsim = 1, nsub = plist[["nSamples"]] * 2, thetas = .model$coefficients,
      #              times = landmarkTimes, Data = .data,
      #              formulas = list(Yfixed = .model$termsYx, Yrandom = .model$termsYz,
      #                              Tfixed = .model$formT, timeVar = "time"), max.FUtime = 144)[[1]]
      # })
      
      simulatedData <- try({
        simulateJM(nsim = 1, nsub = plist[["nSamples"]] * 2, thetas = .model$coefficients,
                   times = landmarkTimes, Data = .data, seed = plist[["seed"]],
                   formulas = list(Yfixed = .model$termsYx, Yrandom = .model$termsYz,
                                   Tfixed = .model$formT, timeVar = "time"), ...)[[1]]
      })
    })[3]
    
    # Veri başarılı şekilde üretildi ise
    if (!inherits(simulatedData, "try-error")){
      # Train/Test
      IDs <- sort(unique(simulatedData[["id"]]))
      idx.train <- sort(sample(IDs, size = floor(length(IDs) / 2), replace = FALSE))
      
      train_data <- simulatedData %>%
        filter(id %in% idx.train)
      
      test_data <- simulatedData %>%
        filter(!(id %in% idx.train))
      
      # Model fit
      lmeFitTRUE <- train_data %>%
        lme(y ~ age + highlow + peritonitrate + time,
            data = .,
            random = ~ time | id)
      
      lmeFitFALSE <- train_data %>%
        lme(y ~ 1 + time,
            data = .,
            random = ~ 1 | id)
      
      survFitTRUE <- train_data %>%
        filter(!duplicated(id)) %>%
        coxph(Surv(Time, Event) ~ age + preHD + ill_count + peritonitrate,
              data = ., x = TRUE)
      
      survFitFALSE <- train_data %>%
        filter(!duplicated(id)) %>%
        coxph(Surv(Time, Event) ~ gender, data = ., x = TRUE)
      
      # Bütün modeller hatasız bir şekilde kurulabildi ise jointModel aşamasına geçiliyor.
      if (!all(inherits(lmeFitTRUE, "try-error"), inherits(lmeFitFALSE, "try-error"),
               inherits(survFitTRUE, "try-error"), inherits(survFitFALSE, "try-error"))){
        jointFitTRUE <- try(jointModel(lmeObject = lmeFitTRUE, survObject = survFitTRUE, timeVar = "time"))
        jointFitFALSE <- try(jointModel(lmeObject = lmeFitFALSE, survObject = survFitFALSE, timeVar = "time"))
      }
      
      # Joint modeller hatasız şekilde kuruldu ise performans hesaplamalarına geçiliyor.
      if (!all(inherits(jointFitTRUE, "try-error"), inherits(jointFitFALSE, "try-error"))){
        # Calculate Performance Measures
        # AUC
        auc_fun <- function (.trueModel, .falseModel, time, ...) {
          aucTRUE <- try({
            aucJM(.trueModel, Tstart = time, ...)[["auc"]]
          })
          
          if (inherits(aucTRUE, "try-error")){
            return(NA)
          } else {
            aucFALSE <- try({
              aucJM(.falseModel, Tstart = time, ...)[["auc"]]
            })
          }
          
          if (inherits(aucFALSE, "try-error")){
            return(NA)
          }
          
          if (any(is.na(aucTRUE), is.nan(aucTRUE), is.na(aucFALSE), is.nan(aucFALSE))){
            return(NA)
          }
          return(aucTRUE - aucFALSE)
        }
        JM_AUCs <- sapply(scen[["Tstart"]], auc_fun, .trueModel = jointFitTRUE, .falseModel = jointFitFALSE, 
                          Dt = deltaT, Thoriz = horizonTime, newdata = test_data)
        
        # PE
        pe_fun <- function(.trueModel, .falseModel, time, horiz = NULL, dt = NULL, ...) {
          if (is.null(horiz)){
            horiz <- time + dt
          }
          
          bsTRUE <- try({
            prederrJM(.trueModel, Tstart = time, Thoriz = horiz, ...)[["prederr"]]
          })
          
          if (inherits(bsTRUE, "try-error")){
            return(NA)
          } else {
            bsFALSE <- try({
              prederrJM(.falseModel, Tstart = time, Thoriz = horiz, ...)[["prederr"]]
            })
          }
          
          if (inherits(bsFALSE, "try-error")){
            return(NA)
          }
          if (any(is.na(bsTRUE), is.nan(bsTRUE), is.na(bsFALSE), is.nan(bsFALSE))){
            return(NA)
          }
          
          return(bsTRUE - bsFALSE)
        }
        
        JM_PEs <- sapply(scen[["Tstart"]], pe_fun, .trueModel = jointFitTRUE, .falseModel = jointFitFALSE, 
                         dt = deltaT, horiz = horizonTime, newdata = test_data, lossFun = "square")
        
        
        measures <- data.frame(scen, Iteration = iter, CPUtimeData = round(as.numeric(cputime_simData) / 60, 1), 
                               AUC = JM_AUCs, BS = JM_PEs)
        break
      }
    }
    
    iter <- iter + 1
  }
  
  if (save){
    if (!dir.exists(path)){
      dir.create(path, recursive = TRUE)
    }
    
    fname <- paste0("SimID_", unique(scen[["SimID"]]), ".Rda")
    save(measures, file = file.path(path, fname))
  } else {
    return(measures)
  }
}

# Parallel hesaplamalar.
cl <- makeCluster(30)
clusterExport(cl, c("dataLong"))

system.time({
  results <- parSapply(cl, X = scenarios[["SimID"]], runSim, .model = jointFit, .data = dataLong, 
                       plist = scenarios, path = "/home/opencpu/externalvolume/Containers/merve/Tez/Analizler/results/JM/Simulation", 
                       deltaT = 12, max.FUtime = 144, simplify = FALSE) %>%
    bind_rows
})

# Stop running clusters.
stopCluster(cl)

results <- data.frame()
for (file in dir("results/JM/Simulation/")){
  load(file.path("results/JM/Simulation/", file))
  results <- rbind(results, measures)
}

load("results/JM/Simulation/results.Rda")
results <- results %>%
  mutate(nSamples = factor(nSamples, levels = c(100, 400, 800), labels = c("Small (n = 100)", "Medium (n = 400)", "Large (n = 800)")),
         maxFUpTime = factor(maxFUpTime, levels = c(60, 120), labels = c("Max(t) = 60", "Max(t) = 120")),
         Tstart = as.factor(Tstart)) %>%
  as_tibble

summary(results)

ids <- unique(results$SimID)
missing <- c(1:1200)[!(1:1200 %in% ids)]
filter(scenarios, SimID %in% missing)

pp <- ggplot(results, aes(x = Tstart, y = BS)) + 
  theme_bw(base_size = 14) + 
  geom_boxplot(fill = "gray90") +
  geom_hline(yintercept = 0, lty = 2) + 
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(margin = margin(t = 5, b = 5)),
        axis.text.y = element_text(margin = margin(r = 5, l = 5))) + 
  facet_wrap(maxFUpTime ~ nSamples, nrow = 2, ncol = 3, scales = "free") + 
  labs(x = "Time (months)", y = expression(paste(Delta, "BS")))

ggsave(filename = "simulationBS-EN.png", plot = pp, device = "png", path = "figures/", width = 10, height = 6, units = "in", dpi = 450)


