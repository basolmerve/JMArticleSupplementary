library(ggplot2)
library(magrittr)
library(readr)
library(ggpubr)

plotData <- read_delim(file = "data/iteration_results.txt", delim = "\t", col_names = TRUE)

# Plots for model (and model parameters) convergence
# Log-likelihood (Model)

axisSpacing <- theme(
  axis.text.x = element_text(margin = margin(5, 0, 5, 0)),
  axis.text.y = element_text(margin = margin(0, 5, 0, 5))
)

addExtraSpaceY <- theme(axis.title.y = element_text(margin = margin(0, 5, 0, 20)))
  
fontSize <- 13

logL <- ggplot(plotData, aes(x = No, y = logL / 1000)) +
  geom_line() + 
  theme_classic(base_size = fontSize) +
  axisSpacing + 
  addExtraSpaceY + 
  labs(x = "Iteration", y = expression('log-likelihood of the model (x1000)')) + 
  ylim(c(-2, -1.8))


# Model parameters
p.age <- ggplot(plotData, aes(x = No, y = beta2)) +
  geom_line() + 
  theme_classic(base_size = fontSize) +
  axisSpacing + 
  labs(x = "Iteration", y = expression('Model parameter'~~('age, '~beta[1])))

p.tppm <- ggplot(plotData, aes(x = No, y = beta3)) +
  geom_line() + 
  theme_classic(base_size = fontSize) +
  axisSpacing + 
  addExtraSpaceY +
  labs(x = "Iteration", y = expression('Model parameter'~~('TPPM, '~beta[2])))

p.pr <- ggplot(plotData, aes(x = No, y = beta4)) +
  geom_line() + 
  theme_classic(base_size = fontSize) +
  axisSpacing + 
  labs(x = "Iteration", y = expression('Model parameter'~~('PR, '~beta[3])))

p.time <- ggplot(plotData, aes(x = No, y = beta5)) +
  geom_line() + 
  theme_classic(base_size = fontSize) +
  axisSpacing + 
  addExtraSpaceY +
  labs(x = "Iteration", y = expression('Model parameter'~~('time, '~beta[4])))

p.assoc <- ggplot(plotData, aes(x = No, y = alpha)) +
  geom_line() + 
  theme_classic(base_size = fontSize) +
  axisSpacing + 
  labs(x = "Iteration", y = expression('Model parameter'~~('Assoc., '~alpha[1])))

p.age.surv <- ggplot(plotData, aes(x = No, y = gamma2)) +
  geom_line() + 
  theme_classic(base_size = fontSize) +
  axisSpacing + 
  labs(x = "Iteration", y = expression('Model parameter'~~('age, '~gamma[1])))

p.pd.surv <- ggplot(plotData, aes(x = No, y = gamma3)) +
  geom_line() + 
  theme_classic(base_size = fontSize) +
  axisSpacing + 
  addExtraSpaceY +
  labs(x = "Iteration", y = expression('Model parameter'~~('PD history, '~gamma[2])))

p.ni.surv <- ggplot(plotData, aes(x = No, y = gamma4)) +
  geom_line() + 
  theme_classic(base_size = fontSize) +
  axisSpacing + 
  labs(x = "Iteration", y = expression('Model parameter'~~('NI, '~gamma[3])))

p.pr.surv <- ggplot(plotData, aes(x = No, y = gamma5)) +
  geom_line() + 
  theme_classic(base_size = fontSize) +
  axisSpacing + 
  addExtraSpaceY +
  labs(x = "Iteration", y = expression('Model parameter'~~('PR, '~gamma[4])))


# Combine convergerce plots for longitidunal part of Joint Model.
converge_lme <- ggarrange(p.age, p.tppm, p.pr, p.time, ncol = 2, nrow = 2)
converge_surv <- ggarrange(p.age.surv, p.pd.surv, p.ni.surv, p.pr.surv, p.assoc, logL, 
                           ncol = 2, nrow = 3)

ggsave(filename = "figures/converge_lme.png", plot = converge_lme, 
       width = 24, height = 20, units = "cm", dpi = 450, device = "png")

ggsave(filename = "figures/converge_surv.png", plot = converge_surv, 
       width = 24, height = 30, units = "cm", dpi = 450, device = "png")

