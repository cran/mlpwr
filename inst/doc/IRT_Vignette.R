## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----eval=FALSE---------------------------------------------------------------
#  install.packages("mlpwr")
#  install.packages("mirt")

## ----message=FALSE, warning=FALSE, results='hide'-----------------------------
library(mlpwr)
library(mirt)

## ----message=FALSE, warning=FALSE, results='hide'-----------------------------
# Defining intercepts and slopes
a <- c(1.04, 1.2, 1.19, 0.61, 1.31, 0.83, 1.46, 1.27, 0.51, 0.81)
d <- c(0.06, -1.79, -1.15, 0.88, -0.2, -1.87, 1.23, -0.08, -0.71, 0.6)

# Setting number of observations 
N <- 100

# Itemtype 
itemtype <- "2PL"

## ----echo=TRUE----------------------------------------------------------------
# Simulate Data
sim_data <- simdata(a = a, d = d, N = N, itemtype = itemtype)

# First 5 rows if simulated data
sim_data[1:5,]

## ----message=FALSE, warning=FALSE, results='hide'-----------------------------
# Fit 2PL model
mod <- mirt(sim_data)

# Rasch contsraint for items 1-4
constrained <- "F = 1-4
          CONSTRAIN = (1-4, a1)"
# Fit constrained model
mod_constrained <- mirt(sim_data, constrained)  # Fit 2PL with equal item discrimination

# Compare model fit
res <- anova(mod_constrained, mod)  

## ----echo=TRUE----------------------------------------------------------------
res$p[2] < 0.01  # extract significance

## ----message=FALSE, warning=FALSE, results='hide'-----------------------------
simfun_irt1 <- function(N) {

    # generate data
    dat <- simdata(a = c(1.04, 1.2, 1.19, 0.61, 1.31,
        0.83, 1.46, 1.27, 0.51, 0.81), d = c(0.06,
        -1.79, -1.15, 0.88, -0.2, -1.87, 1.23, -0.08,
        -0.71, 0.6), N = N, itemtype = "2PL")

    # test hypothesis
    mod <- mirt(dat)  # Fit 2PL Model
    constrained <- "F = 1-4
          CONSTRAIN = (1-4, a1)"
    mod_constrained <- mirt(dat, constrained)  # Fit 2PL with equal slopes

    res <- anova(mod_constrained, mod)  # perform model comparison
    res$p[2] < 0.01  # extract significance
}

## ----eval = FALSE-------------------------------------------------------------
#  example.simfun("irt1")

## ----echo=FALSE, results='hide'-----------------------------------------------
# The following loads the precomputed results of the next chunk to reduce the vignette creation time
ver <- as.character(packageVersion("mlpwr"))
file = paste0("/extdata/IRT_Vignette_results1_", ver, ".RData")
file_path <- paste0(system.file(package="mlpwr"),file)
if (!file.exists(file_path)) {
set.seed(123)
res <- find.design(simfun = simfun_irt1, boundaries = c(40,
    100), power = .95, evaluations = 2000)
save(res, file = paste0("../inst",file))
} else {
  load(file_path) 
}

## ----results='hide', eval=FALSE-----------------------------------------------
#  set.seed(123)
#  res <- find.design(simfun = simfun_irt1, boundaries = c(40,
#      100), power = .95, evaluations = 2000)

## ----echo=TRUE----------------------------------------------------------------
summary(res)

## -----------------------------------------------------------------------------
plot(res)

## ----message=FALSE, warning=FALSE, results='hide'-----------------------------
library(mlpwr)
library(mirt)

## ----message=FALSE, warning=FALSE, results='hide'-----------------------------
simfun_irt2 <- function(N1, N2) {
  # generate data
  a1 <- a2 <- c(1.04, 1.2, 1.19, 0.61, 1.31, 0.83,
                1.46, 1.27, 0.51, 0.81)
  d1 <- d2 <- c(0.06, -1.79, -1.15, 0.88, -0.2, -1.87,
                1.23, -0.08, -0.71, 0.6)
  a2[1] <- a2[1] + 0.3
  d2[1] <- d2[1] + 0.5
  dat1 <- simdata(a = a1, d = d1, N = N1, itemtype = "2PL")
  dat2 <- simdata(a = a2, d = d2, N = N2, itemtype = "2PL")
  dat <- as.data.frame(rbind(dat1, dat2))
  group <- c(rep("1", N1), rep("2", N2))
  
  # Fit model
  model <- multipleGroup(dat, 1, group = group)
  # Perform DIF for item 1
  dif <- DIF(model, which.par = c("a1", "d"), items2test = c(1))
  #Extract significance
  dif$p[1] < 0.05
}

## ----echo=FALSE---------------------------------------------------------------
a1 <- a2 <- c(1.04, 1.2, 1.19, 0.61, 1.31, 0.83,
    1.46, 1.27, 0.51, 0.81)
d1 <- d2 <- c(0.06, -1.79, -1.15, 0.88, -0.2, -1.87,
    1.23, -0.08, -0.71, 0.6)
a2[1] <- a2[1] + 0.3
d2[1] <- d2[1] + 0.5
dat1 <- simdata(a = a1, d = d1, N = 5, itemtype = "2PL")
dat2 <- simdata(a = a2, d = d2, N = 5, itemtype = "2PL")
group <- c(rep("1", 5), rep("2", 5))
dat <- as.data.frame(rbind(dat1, dat2))
dat$group <- group
dat[, c("group", names(dat)[names(dat) != "group"])]

## -----------------------------------------------------------------------------
costfun_irt2 <- function(N1, N2) 5 * N1 + 7 * N2

## ----echo=FALSE, results='hide'-----------------------------------------------
# The following loads the precomputed results of the next chunk to reduce the vignette creation time
ver <- as.character(packageVersion("mlpwr"))
file = paste0("/extdata/IRT_Vignette_results2_", ver, ".RData")
file_path <- paste0(system.file(package="mlpwr"),file)
if (!file.exists(file_path)) {
set.seed(111)
res <- find.design(simfun = simfun_irt2, boundaries = list(N1 = c(100,
    700), N2 = c(100, 700)), costfun = costfun_irt2,
    cost = 4000, evaluations = 1000)
save(res, file = paste0("../inst",file))
} else {
  load(file_path) 
}

## ----eval=FALSE---------------------------------------------------------------
#  set.seed(111)
#  res <- find.design(simfun = simfun_irt2, boundaries = list(N1 = c(100,
#      700), N2 = c(100, 700)), costfun = costfun_irt2,
#      cost = 4000, evaluations = 1000)

## ----echo=TRUE----------------------------------------------------------------
summary(res)

## -----------------------------------------------------------------------------
plot(res)

