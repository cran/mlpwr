## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----eval = FALSE, message=FALSE, warning=FALSE, results='hide'---------------
#  install.packages("mlpwr")
#  install.packages("lme4")
#  install.packages("lmerTest")

## ----message=FALSE, warning=FALSE, results='hide'-----------------------------
library(mlpwr)
library(lme4)
library(lmerTest)

## ----echo=FALSE, message=FALSE, warning=FALSE---------------------------------
# generate data
N = 10
# generate data
  params <- list(theta = 0.5, beta = c(2, 0.2))
  num_classes <- ceiling(N/4)
  class_id <- rep(1:num_classes, length.out = N)
  group <- rep(1:2, times=c(floor(N/2), ceiling(N/2)))
  
  dat <- data.frame(class_id = class_id, group = group)
  dat$x <- simulate(~group + (1 | class_id), newdata = dat,
      family = poisson, newparams = params)[[1]]
dat[1:5,]

## ----warning=FALSE------------------------------------------------------------
simfun_multi1 <- function(N) {
  
  # generate data
  params <- list(theta = 0.5, beta = c(2, 0.2))
  num_classes <- ceiling(N/4) # We can recruit a max of 4 people per class
  class_id <- rep(1:num_classes, length.out = N)
  group <- rep(1:2, times=c(floor(N/2), ceiling(N/2)))
  
  dat <- data.frame(class_id = class_id, group = group)
  dat$x <- simulate(~group + (1 | class_id), newdata = dat,
      family = poisson, newparams = params)[[1]]

  # model
  mod <- glmer(x ~ group + (1 | class_id), data = dat,
      family = poisson)  # fit model
  
  # Extract P-Value and coefficient
  p_value <- summary(mod)$coefficient["group", "Pr(>|z|)"]
  group_coef <- summary(mod)$coefficient["group", "Estimate"]

  # Check if coefficient is significantly positive
  p_value < 0.01 & group_coef > 0
}


## ----echo=FALSE, results='hide'-----------------------------------------------
# The following loads the precomputed results of the next chunk to reduce the vignette creation time
ver <- as.character(packageVersion("mlpwr"))
file = paste0("/extdata/MLM_Vignette_results1_", ver, ".RData")
file_path <- paste0(system.file(package="mlpwr"),file)
if (!file.exists(file_path)) {
set.seed(111)
res <- find.design(simfun = simfun_multi1, boundaries = c(20,
    200), power = .8, evaluations = 2000)
save(res, file = paste0("../inst",file))
} else {
  load(file_path) 
}

## ----warning=FALSE, eval=FALSE------------------------------------------------
#  set.seed(111)
#  res <- find.design(simfun = simfun_multi1, boundaries = c(20,
#      200), power = .8, evaluations = 2000)

## ----echo=TRUE----------------------------------------------------------------
summary(res)

## -----------------------------------------------------------------------------
plot(res)

## ----message=FALSE, warning=FALSE, results='hide'-----------------------------
library(mlpwr)
library(lme4)
library(lmerTest)

## -----------------------------------------------------------------------------
logistic <- function(x) 1/(1 + exp(-x))
set.seed(109)

# 300 participants from 20 countries
N.original <- 300
n.countries.original <- 20

# generate original data
dat.original <- data.frame(country = rep(1:n.countries.original,
    length.out = N.original), iq = rnorm(N.original),
    cortisol = rnorm(N.original))
country.intercepts <- rnorm(n.countries.original, sd = 0.5)
dat.original$intercepts <- country.intercepts[dat.original$country]
beta <- c(1, 0.4, -0.3)  # parameter weights
prob <- logistic(as.matrix(dat.original[c("intercepts",
    "iq", "cortisol")]) %*% as.matrix(beta))  # get probability
dat.original$criterion <- rbinom(N.original, 1, prob)  # draw according to probability
dat.original <- dat.original[,names(dat.original)!="intercepts"]

# fit original model to obtain parameters
mod.original <- glmer(criterion ~ iq + cortisol + 0 +
    (1 | country), data = dat.original, family = binomial)
dat.original[1:5,]

## ----warning=FALSE------------------------------------------------------------
simfun_multi2 <- function(n, n.countries) {

    # generate data
    dat <- data.frame(country = rep(1:n.countries,
        length.out = n * n.countries), iq = rnorm(n *
        n.countries), cortisol = rnorm(n * n.countries))
    dat$criterion <- simulate(mod.original, nsim = 1,
        newdata = dat, allow.new.levels = TRUE, use.u = FALSE) |>
        unlist()  # criterion data from the fitted model

    # test hypothesis
    mod <- glmer(criterion ~ iq + cortisol + 0 + (1 |
        country), data = dat, family = binomial)
    summary(mod)[["coefficients"]]["cortisol", "Pr(>|z|)"] <
        0.01
}

## ----eval = FALSE-------------------------------------------------------------
#  example.simfun("multi2")

## -----------------------------------------------------------------------------
costfun_multi2 <- function(n, n.countries) 5 * n +
    100 * n.countries

## ----echo=FALSE, results='hide'-----------------------------------------------
# The following loads the precomputed results of the next chunk to reduce the vignette creation time
ver <- as.character(packageVersion("mlpwr"))
file = paste0("/extdata/MLM_Vignette_results2_", ver, ".RData")
file_path <- paste0(system.file(package="mlpwr"),file)
if (!file.exists(file_path)) {
set.seed(112)
res <- find.design(simfun = simfun_multi2, boundaries = list(n = c(10,
    300), n.countries = c(5, 20)), costfun = costfun_multi2,
    cost = 2000, evaluations = 2000)
save(res, file = paste0("../inst",file))
} else {
  load(file_path) 
}

## ----warning=FALSE,eval=FALSE-------------------------------------------------
#  set.seed(112)
#  res <- find.design(simfun = simfun_multi2, boundaries = list(n = c(10,
#      300), n.countries = c(5, 20)), costfun = costfun_multi2,
#      cost = 2000, evaluations = 2000)

## ----echo=TRUE----------------------------------------------------------------
summary(res)

## -----------------------------------------------------------------------------
plot(res)

