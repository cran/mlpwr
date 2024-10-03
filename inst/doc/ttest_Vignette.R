## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----eval=FALSE---------------------------------------------------------------
#  # Example code for performing Welch's t-test in R
#  t.test(group_a, group_b, var.equal = FALSE)
#  

## ----echo=FALSE---------------------------------------------------------------
groupA <- rnorm(5, mean = 5, sd = 2)
groupB <- rnorm(5, mean = 7, sd = 1)
dat <- data.frame("blood_pressure" = c(groupA, groupB), "group" = rep(c("A","B"), each=5))
dat

## -----------------------------------------------------------------------------
t.test(dat[dat$group=="A",]$blood_pressure, dat[dat$group=="B",]$blood_pressure, 
       alternative = "less", var.equal = FALSE)


## ----eval=FALSE---------------------------------------------------------------
#  install.packages("mlpwr")

## ----message=FALSE, warning=FALSE, results='hide'-----------------------------
library(mlpwr)

## -----------------------------------------------------------------------------
simfun_ttest <- function(nA, nB) {
  groupA <- rnorm(nA, mean = 5, sd = 2)
  groupB <- rnorm(nB, mean = 6, sd = 1)
  res <- t.test(groupA, groupB, alternative = "less", var.equal = FALSE)
  res$p.value < 0.01
}

## ----eval=FALSE---------------------------------------------------------------
#  simfun_ttest_example <- example.simfun("ttest")

## -----------------------------------------------------------------------------
costfun_ttest <- function(nA, nB) {1.5*nA + 1*nB}

## ----echo=FALSE, results='hide'-----------------------------------------------
# The following loads the precomputed results of the next chunk to reduce the vignette creation time
ver <- as.character(packageVersion("mlpwr"))
file = paste0("/extdata/ttest_Vignette_results_", ver, ".RData")
file_path <- paste0(system.file(package="mlpwr"),file)
if (!file.exists(file_path)) {
set.seed(111)
res <- find.design(simfun = simfun_ttest, boundaries = list(nA = c(5,200), nB = c(5,200)),
                   power = .8, costfun = costfun_ttest, evaluations = 1000)
save(res, file = paste0("../inst",file))
} else {
  load(file_path) 
}

## ----warning=FALSE, eval = FALSE----------------------------------------------
#  set.seed(111)
#  res <- find.design(simfun = simfun_ttest, boundaries = list(nA = c(5,200), nB = c(5,200)),
#                     power = .8, costfun = costfun_ttest, evaluations = 1000)

## ----echo=TRUE----------------------------------------------------------------
summary(res)

## -----------------------------------------------------------------------------
plot(res)

