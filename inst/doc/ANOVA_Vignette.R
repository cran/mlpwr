## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----eval=FALSE---------------------------------------------------------------
#  # Example code for performing ANOVA in R
#  result <- aov(outcome ~ group, data = dataset)

## ----echo=FALSE---------------------------------------------------------------
groupA <- rnorm(2, mean = 5, sd = 1)
groupB <- rnorm(2, mean = 7, sd = 1)
groupC <- rnorm(2, mean = 6, sd = 1)
dat <- data.frame("score" = c(groupA, groupB, groupC), "group" = rep(c("A","B", "C"), each=2))
dat

## -----------------------------------------------------------------------------
res <- aov(score ~ group, data = dat)
summary(res)

## -----------------------------------------------------------------------------
TukeyHSD(res)

## ----eval=FALSE---------------------------------------------------------------
#  install.packages("mlpwr")
#  install.packages("tidyr")

## ----message=FALSE, warning=FALSE, results='hide'-----------------------------
library(mlpwr)
library(tidyr)

## -----------------------------------------------------------------------------
simfun_anova <- function(n, n.groups) {
  groupmeans <- rnorm(n.groups, sd = 0.2)
  dat <- sapply(groupmeans, function(x) rnorm(n, mean = x, sd = 1))
  dat <- gather(as.data.frame(dat))
  res <- aov(value ~ key, data = dat)
  summary(res)[[1]][1, 5] < 0.01
}

## ----eval=FALSE---------------------------------------------------------------
#  simfun_ANOVA_example <- example.simfun("anova")

## -----------------------------------------------------------------------------
costfun_anova <- function(n, n.groups) {n + n.groups*10}

## ----echo=FALSE, results='hide'-----------------------------------------------
# The following loads the precomputed results of the next chunk to reduce the vignette creation time
ver <- as.character(packageVersion("mlpwr"))
file = paste0("/extdata/ANOVA_Vignette_results_", ver, ".RData")
file_path <- paste0(system.file(package="mlpwr"),file)
if (!file.exists(file_path)) {
set.seed(112)
res <- find.design(simfun = simfun_anova, boundaries = list(n = c(5,200), n.groups = c(5,30)),
                   power = .8, costfun = costfun_anova, evaluations = 2000)
save(res, file = paste0("../inst",file))
} else {
  load(file_path) 
}

## ----results='hide', eval=FALSE-----------------------------------------------
#  set.seed(112)
#  res <- find.design(simfun = simfun_anova, boundaries = list(n = c(5,200), n.groups = c(5,30)),
#                     power = .8, costfun = costfun_anova, evaluations = 2000)

## ----echo=TRUE----------------------------------------------------------------
summary(res)

## -----------------------------------------------------------------------------
plot(res)

