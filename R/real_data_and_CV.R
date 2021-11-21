# Required R libraries
library(readr)
library(dplyr)
library(magrittr)
library(JM)
library(ggplot2)


#----- Internal functions ------#
map <- function(x, min.y = 0, max.y = 1, xmin = NULL, xmax = NULL, ..){
  if (is.null(xmin)){
    xmin <- min(x)
  }
  
  if (is.null(xmax)){
    xmax <- max(x)
  }
  
  ratio <- (x - xmin) / (xmax - xmin)
  min.y + ratio * (max.y - min.y)
}

# Read data into R
# Dataset is available in two separate files. Complete dataset is the combinatin of two different data files.
part1 <- read_delim(file = "data/longitudinal60.txt", delim = "\t", col_names = TRUE, 
                    locale = locale(decimal_mark = "."))
part2 <- read_delim("data/ekveriler.csv", delim = ",", locale = locale(decimal_mark = "."))

pdData <- part1 %>%
  full_join(part2) %>%
  mutate_if(is.character, as.factor) %>%
  mutate(preHD = factor(preHD, levels = c(0, 1), labels = c("No", "Yes")),
         highlow = factor(highlow, levels = c(1, 2), labels = c("Low", "High")))

# Remove rows with missing data in the variables of interest.
trainData <- pdData %>%
  filter(complete.cases(age, time, ill_count, event, surv, peritonitrate, highlow, preHD, ALB)) 

#----- Fit Joint Models --------#
# Part 1: Linear mixed effects model (LME)
lmeFit <- trainData %>%
  lme(ALB ~ age + highlow + peritonitrate + time, 
      data = ., random = ~ time | id)

# Part 2: Cox proportional hazard model (Cox PH)
survFit <- trainData %>%
  filter(!duplicated(id)) %>%
  coxph(Surv(surv, event) ~ age + preHD + ill_count + peritonitrate, data = ., x = TRUE)


set.seed(2128)
#  MODEL1 - current value  
jointFit1 <- jointModel(lmeFit, survFit, timeVar = "time", verbose=T)
summary(jointFit1)

# MODEL 2
# dForm <- list(fixed = ~  -1 + I(rep(1/0.0105, length(time))), indFixed = 5,  random = ~ -1 + I(rep(1/0.0105, length(time))) , indRandom = 2)
dForm <- list(fixed = ~  1, indFixed = 5,  random = ~ 1, indRandom = 2)
jointFit2_tmp <- update(jointFit1, interFact = NULL)
jointFit2 <- update(jointFit2_tmp, parameterization = "both", derivForm = dForm)
summary(jointFit2)

# MODEL3 
g <- function(u, pow = 0){
  f <- function(t){
    integrate(function(s) s ^ pow * dnorm(t - s), 0, t)$value
  }
  sapply(u, f)
}

iformW <- list(fixed = ~ -1 + g(time) + I(g(time)* age) + I(g(time) * as.numeric(highlow)) + I(g(time) * as.numeric(peritonitrate)) + I(g(time, 1)), 
               indFixed = 1:5, random = ~ -1 + g(time) + I(g(time, 1)), indRandom = 1:2)

jointFit3 <- update(jointFit2, parameterization = "slope", derivForm = iformW)

summary(jointFit3)

# MODEL4  
#dform2 <- list(fixed = ~  -1 + I(rep(1/0.0105, length(time))), indFixed = 5,  random = ~ -1 + I(rep(1/0.0105, length(time))) , indRandom = 2)
dform2 <- list(fixed = ~ 1, indFixed = 5, random = ~ 1, indRandom = 2)
jointFit4 <- update(jointFit2, parameterization = "slope", derivForm = dform2)
summary(jointFit4)



####### CROSS-VALIDATION ##########
n <- length(unique(dataLong$id)) # number of subjects
V <- 5 # number of folds in the cross-validation
M <- 10 # number of times to replicate the cross-validation procedure
set.seed(1)
seedList <- sample(1:10000, V * M, FALSE)

folds <- lapply(seq_len(M), function(x, nFolds = V, nSamples = n, ...){
  set.seed(100 + x)
  splits <- split(seq_len(nSamples), sample(rep(seq_len(nFolds), length.out = nSamples)))
  return(splits)
})

allFolds <- list()
for (i in 1:length(folds)){
  for (j in 1:length(folds[[i]])){
    tmp <- list(list(seed = seedList[i * j], rep = i, fold = j, id = folds[[i]][[j]]))
    allFolds <- c(allFolds, tmp)
  }
}

runFold <- function(data, foldList = NULL, nboot = 10000, save = TRUE, 
                    path = file.path(getwd(), "results/JM"), landmarkTimes = seq(6, 60, by = 6),
                    deltaT = 6.1, horizonTime = NULL, ...){
  
  
  set.seed(foldList[["seed"]])
  
  # load data
  foldIDs <- foldList[["id"]]
  train <- data[!data$id %in% foldIDs, ]
  test <- data[data$id %in% foldIDs, ]
  train$id <- match(train$id, unique(train$id))
  test$id <- match(test$id, unique(test$id))
  
  calculate_REs <- function (lmeObject, newdata) {
    data <- lmeObject$data
    formYx <- formula(lmeObject)
    mfX <- model.frame(terms(formYx), data = data)
    TermsX <- attr(mfX, "terms")
    mfX_new <- model.frame(TermsX, data = newdata)
    X_new <- model.matrix(formYx, mfX_new)
    formYz <- formula(lmeObject$modelStruct$reStruct[[1]])
    mfZ <- model.frame(terms(formYz), data = data)
    TermsZ <- attr(mfZ, "terms")
    mfZ_new <- model.frame(TermsZ, data = newdata)
    Z_new <- model.matrix(formYz, mfZ_new)
    y_new <- model.response(mfX_new, "numeric")
    idVar <- names(lmeObject$modelStruct$reStruct)
    if (length(idVar) > 1)
      stop("the current version of the function only works with a single grouping variable.\n")
    if (is.null(newdata[[idVar]]))
      stop("subject id variable not in newdata.")
    id <- match(newdata[[idVar]], unique(newdata[[idVar]]))
    n <- length(unique(id))
    betas <- fixef(lmeObject)
    D <- lapply(pdMatrix(lmeObject$modelStruct$reStruct), "*",
                lmeObject$sigma^2)[[1]]
    modes <- matrix(0.0, n, ncol(Z_new))
    for (i in seq_len(n)) {
      id_i <- id == i
      X_new_id <- X_new[id_i, , drop = FALSE]
      Z_new_id <- Z_new[id_i, , drop = FALSE]
      Vi_inv <- solve(Z_new_id %*% tcrossprod(D, Z_new_id) + 
                        lmeObject$sigma^2 * diag(sum(id_i)))
      DZtVinv <- tcrossprod(D, Z_new_id) %*% Vi_inv
      modes[i, ] <- c(DZtVinv %*% (y_new[id_i] - X_new_id %*% betas))
    }
    modes
  }
  
  fitJointModels <- function(train_data, test_data, landmark = NULL, dt = NULL,
                             horiz = NULL, B = nboot) {
    
    # dt ve horiz ikisi birden NULL olamaz. Bu durumda deltaT değeri 1 olarak alınıyor.
    if (all(is.null(horiz), is.null(dt))){
      dt <- 1
    }
    
    # dt ve horiz ikisi birden tanımlanmış ise dt dikkate alınır.
    if (all(!is.null(horiz), !is.null(dt))){
      horiz <- NULL
    }
    
    
    modelFitSuccess <- try({
      ##################
      # Joint Modeling #
      ##################
      
      # train_data <- train_data %>%
      #   filter(complete.cases(ALB, age, time, ill_count, surv, event))
      
      # Fit Joint Models
      # Seçilen değişkenlerde eksik olan gözlemler silindikten sonra "lme" modeli kuruluyor. "longitudinal" parça için.
      lmeFit <- train_data %>%
        lme(ALB ~ age + highlow + peritonitrate + time,
            data = .,
            random = ~ time | id)
      
      # Seçilen değişkenlerde eksik olan gözlemler silindikten sonra "coxph" modeli kuruluyor. "survival" parça için.
      # survFit <- train_data %>%
      #   filter(!duplicated(id)) %>%
      #   coxph(Surv(surv, event) ~ age + gender + ill_count + BKI + highlow + peritonitrate, data = ., x = TRUE)
      
      survFit <- train_data %>%
        filter(!duplicated(id)) %>%
        coxph(Surv(surv, event) ~ age + preHD + ill_count + peritonitrate, 
              data = ., x = TRUE)
      
      # jmControl <- list(tol1 = 1e-7, tol2 = 1e-7)
      # MODEL1
      jointFit1 <- jointModel(lmeObject = lmeFit, survObject = survFit, timeVar = "time", 
                              method = "weibull-PH-aGH")
      
      # MODEL 2
      dForm <- list(fixed = ~ 1, indFixed = 5,  random = ~ 1 , indRandom = 2)
      jointFit2_tmp <- update(jointFit1, interFact = NULL)
      jointFit2 <- update(jointFit2_tmp, parameterization = "both", derivForm = dForm)
      
      # MODEL3 
      g <- function(u, pow = 0){
        f <- function(t){
          integrate(function(s) s ^ pow * dnorm(t - s), 0, t)$value
        }
        sapply(u, f)
      }
      
      iformW <- list(fixed = ~ -1 + g(time) + I(g(time)* age) + I(g(time) * as.numeric(highlow)) + I(g(time) * as.numeric(peritonitrate)) + I(g(time, 1)), 
                     indFixed = 1:5, random = ~ -1 + g(time) + I(g(time, 1)), indRandom = 1:2)
      jointFit3 <- update(jointFit2, parameterization = "slope", derivForm = iformW)
      
      # MODEL4  
      dform2 <- list(fixed = ~ 1, indFixed = 5, random = ~ 1, indRandom = 2)
      jointFit4 <- update(jointFit2, parameterization = "slope", derivForm = dform2)
      
      TRUE   # Model kurma sürecinde hata ile karşılaşılmaz ise "modelFitSuccess" içerisine TRUE yazdırmak için bu satır eklendi.
    }) 
    
    modelNames <- c("Model1", "Model2", "Model3", "Model4")
    combos <- expand.grid("model_name" = modelNames, 
                          "time" = landmark, stringsAsFactors = FALSE)
    
    if (!inherits(modelFitSuccess, "try-error")){
      # Calculate Performance Measures
      JM_Models <- list(jointFit1, jointFit2, jointFit3, jointFit4)
      names(JM_Models) <- modelNames
      
      # AUC
      auc_fun <- function (model_name, time, ...) {
        auc <- try({
          aucJM(JM_Models[[model_name]], Tstart = time, ...)[["auc"]]
        })
        if (inherits(auc, "try-error")){
          return(NA)
        }
        
        if (any(is.na(auc), is.nan(auc))){
          return(NA)
        }
        return(auc)
      }
      JM_AUCs <- mapply(auc_fun, combos[["model_name"]], combos[["time"]], 
                        MoreArgs = list(Dt = dt, newdata = test_data, Thoriz = horiz))
      
      # PE
      pe_fun <- function(model_name, time, horiz = NULL, dt = NULL, ...) {
        bs <- try({
          if (is.null(horiz)){
            horiz <- time + dt
          }
          prederrJM(JM_Models[[model_name]], Tstart = time,  
                    Thoriz = horiz, ...)[["prederr"]]
        })
        if (inherits(bs, "try-error")){
          return(NA)
        }
        if (any(is.na(bs), is.nan(bs))){
          return(NA)
        }
        return(bs)
      }
      JM_PEs <- mapply(pe_fun, combos[["model_name"]], combos[["time"]], 
                       MoreArgs = list(lossFun = "square", newdata = test_data, 
                                       horiz = horiz, dt = dt))
      
      
      measures <- data.frame(Fold = foldList[["fold"]], Repeat = foldList[["rep"]], 
                             Model = combos[["model_name"]], Tstart = combos[["time"]], 
                             AUC = JM_AUCs, BS = JM_PEs)
      
    } else {
      measures <- data.frame(Fold = foldList[["fold"]], Repeat = foldList[["rep"]], 
                             Model = combos[["model_name"]], Tstart = combos[["time"]], 
                             AUC = NA, BS = NA)
    }
    
    return(measures)
  }
  
  res <- fitJointModels(train, test, landmark = landmarkTimes, 
                        dt = deltaT, horiz = horizonTime)
  
  if (save){
    if (!dir.exists(path)){
      dir.create(path, recursive = TRUE)
    }
    
    fname <- paste0("Rep", foldList[["rep"]], "_", "Fold", foldList[["fold"]], ".Rda")
    save(res, file = file.path(path, fname))
  }
  return(res)
}

# Parallel hesaplamalar.
cl <- makeCluster(30)
clusterExport(cl, c("dataLong", "allFolds"))

system.time({
  results <- parLapply(cl, X = allFolds, runFold, data = dataLong, nboot = 10000,
                       path = "/home/opencpu/externalvolume/Containers/merve/Tez/Analizler/results/JM/RealData_CV/DeltaT") %>%
    bind_rows
})

# Stop running clusters.
stopCluster(cl)

# "Results" klasöründe saklanan tüm sonuçların tek bir data.frame içerisinde birleştirilmesi.
results <- data.frame()
for (file in dir("results/JM/RealData_CV/DeltaT/")){
  load(file.path("results/JM/RealData_CV/DeltaT/", file))
  results <- rbind(results, res)
}

###### GRAFİKLER ######
# Line plot
thm <- theme(
  axis.text.x = element_text(margin = margin(t = 5, b = 5)),
  axis.text.y = element_text(margin = margin(r = 5, l = 5))
)
bs <- 18

lineplotAUC <- results %>%
  group_by(Model, Tstart) %>%
  summarise(AUC = mean(AUC, na.rm = TRUE)) %>%
  ungroup(Model) %>%
  ggplot(aes(x = Tstart, y = AUC, shape = Model, linetype = Model)) +
  geom_point(size = 3) +
  geom_line() + 
  theme_bw(base_size = bs) + 
  thm +
  ylim(c(0.50, 1.00)) + 
  scale_x_continuous(breaks = seq(6, 60, by = 6)) + 
  labs(x = "Time (months)", y = "AUC") + 
  guides(linetype = guide_legend(title = ""), shape = guide_legend(title = ""))

lineplotBS <- results %>%
  filter(!is.na(BS), !is.infinite(BS)) %>%
  group_by(Model, Tstart) %>%
  summarise(BS = mean(BS, na.rm = TRUE)) %>%
  ungroup(Model) %>%
  ggplot(aes(x = Tstart, y = BS, shape = Model, linetype = Model)) +
  geom_point(size = 3) +
  geom_line() + 
  theme_bw(base_size = bs) + 
  thm + 
  scale_x_continuous(breaks = seq(6, 60, by = 6)) + 
  labs(x = "Time (months)", y="BS") + 
  guides(linetype = guide_legend(title = ""), shape = guide_legend(title = ""))


ggsave(filename = "lineplotAUC.png", plot = lineplotAUC, device = "png", 
       width = 7, height = 5, units = "in", dpi = 600, path = "figures/")

ggsave(filename = "lineplotBS.png", plot = lineplotBS, device = "png", 
       width = 7, height = 5, units = "in", dpi = 600, path = "figures/")
