## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  tidy = TRUE,
  tidy.opts=list(arrow=TRUE,width.cutoff = 50),
  eval=F
)

## -----------------------------------------------------------------------------
#  simfun <- function(N) {
#      # Generate a data set
#      # Test the hypothesis
#  }

## -----------------------------------------------------------------------------
#  library(mlpwr)

## -----------------------------------------------------------------------------
#  simfun_ttest <- function(N) {
#      # Generate a data set
#      dat <- rnorm(n = N, mean = 0.3)
#      # Test the hypothesis
#      res <- t.test(dat)
#      res$p.value < 0.01
#  }

## ----eval = F-----------------------------------------------------------------
#   res <- find.design(simfun = simfun_ttest,
#       boundaries = c(100,300), power = .95)

## -----------------------------------------------------------------------------
#  library(mlpwr)

## -----------------------------------------------------------------------------
#  simfun_anova <- function(n, n.groups) {
#  
#      # Generate a data set
#      groupmeans <- rnorm(n.groups, sd = 0.2)  # generate groupmeans using cohen's f=.2
#      dat <- sapply(groupmeans, function(x) rnorm(n,
#          mean = x, sd = 1))  # generate data
#      dat <- dat |>
#          as.data.frame() |>
#          gather()  # format
#  
#      # Test the hypothesis
#      res <- aov(value ~ key, data = dat)  # perform ANOVA
#      summary(res)[[1]][1, 5] < 0.01  # extract significance
#  }
#  

## ----eval = F-----------------------------------------------------------------
#  res <- find.design(simfun = simfun_anova,
#     costfun = function(n,n.groups) 5*n+20*n.groups,
#     boundaries = list(n = c(10, 150), n.groups = c(5, 30)),
#     power = .95)

## -----------------------------------------------------------------------------
#  library(mlpwr)

## -----------------------------------------------------------------------------
#  dat.original <- data.frame(
#    counts = c(18, 17, 15, 20,
#      10, 20, 25, 13, 12),
#    treatment = gl(3, 1, 9),
#    outcome = gl(3, 3))
#  mod.original <- glm(counts ~ outcome + treatment, data = dat.original,
#      family = poisson) # setting up the generalized linear model
#  summary(mod.original)

## -----------------------------------------------------------------------------
#  simfun_glm1 <- function(N) {
#  
#      # generate data
#      dat <- data.frame(outcome = gl(3, 1, ceiling(N/3)),
#          treatment = gl(3, ceiling(N/3)))[1:N, ] # predictors
#      a <- predict(mod.original, newdata = dat, type = "response") # criterion 'raw'
#      dat$counts <- rpois(N, a) # criterion applying poisson distribution
#  
#      # test hypothesis
#      mod <- glm(counts ~ outcome + treatment, data = dat,
#          family = poisson) # fit a glm
#      summary(mod)$coefficients["treatment2", "Pr(>|z|)"] <
#          0.01 # test the coefficient of the treatment
#  }

## ----eval = F-----------------------------------------------------------------
#  res <- find.design(simfun = simfun_glm1,
#       boundaries = c(20,100), power = .95)

## -----------------------------------------------------------------------------
#  logistic <- function(x) 1/(1 + exp(-x)) # logistic function

## -----------------------------------------------------------------------------
#  simfun_glm2 <- function(N) {
#  
#      # generate data
#      dat <- data.frame(pred1 = rnorm(N), pred2 = rnorm(N))
#      beta <- c(1.2, 0.8)  # parameter weights
#      prob <- logistic(as.matrix(dat) %*% beta)  # get probability
#      dat$criterion <- runif(N) < prob  # draw according to probability
#  
#      # test hypothesis
#      mod <- glm(criterion ~ pred1 + pred2, data = dat,
#          family = binomial) # fit a glm
#      summary(mod)$coefficients["pred2", "Pr(>|z|)"] <
#          0.01 # test the coefficient of the predictor
#  }

## ----eval = F-----------------------------------------------------------------
#  res <- find.design(simfun = simfun_glm2,
#       boundaries = c(90,200), power = .95)

## -----------------------------------------------------------------------------
#  library(mlpwr)
#  library(mirt)

## -----------------------------------------------------------------------------
#  simfun_irt1 <- function(N) {
#  
#      # generate data
#      dat <- simdata(a = c(1.04, 1.2, 1.19, 0.61, 1.31,
#          0.83, 1.46, 1.27, 0.51, 0.81), d = c(0.06,
#          -1.79, -1.15, 0.88, -0.2, -1.87, 1.23, -0.08,
#          -0.71, 0.6), N = N, itemtype = "2PL") # uses a 2PL model with a and d parameters
#  
#      # test hypothesis
#      mod <- mirt(dat)  # Fit 2PL Model
#      constrained <- "F = 1-4
#            CONSTRAIN = (1-4, a1)" # specifying that the slopes should be kept equal for items 1 to 4
#      mod_constrained <- mirt(dat, constrained)  # Fit 2PL with equal slopes
#  
#      res <- anova(mod_constrained, mod)  # perform model comparison
#      res$p[2] < 0.01  # extract significance
#  }

## ----eval = F-----------------------------------------------------------------
#  res <- find.design(simfun = simfun_irt1,
#       boundaries = c(100,500), power = .95,evaluations =500)

## -----------------------------------------------------------------------------
#  costfun_irt2 <- function(N1, N2) 5 * N1 + 7 * N2 # specifying a cost function

## -----------------------------------------------------------------------------
#  simfun_irt2 <- function(N1, N2) {
#  
#      # generate data
#      a1 <- a2 <- c(1.04, 1.2, 1.19, 0.61, 1.31, 0.83,
#          1.46, 1.27, 0.51, 0.81) # specifying the slope
#      d1 <- d2 <- c(0.06, -1.79, -1.15, 0.88, -0.2, -1.87,
#          1.23, -0.08, -0.71, 0.6) # specifying the difficulty
#      a2[1] <- a2[1] + 0.3 # the slope is different for the first item in group2
#      d2[1] <- d2[1] + 0.5 # the difficulty is different for the first item in group2
#      dat1 <- simdata(a = a1, d = d1, N = N1, itemtype = "2PL") # creating artificial data for both groups
#      dat2 <- simdata(a = a2, d = d2, N = N2, itemtype = "2PL")
#      dat <- as.data.frame(rbind(dat1, dat2)) # combining the data sets into one object
#      group <- c(rep("1", N1), rep("2", N2)) # create a variable that indicates the group membership
#  
#      # fit models
#      mod1 <- multipleGroup(dat, 1, group = group) # fit model with different parameters for each group
#      mod2 <- multipleGroup(dat, 1, group = group, invariance = c("slopes",
#          "intercepts")) # fit model with the same parameters for each group
#  
#      # test hypothesis
#      res <- anova(mod2, mod1)
#  
#      # extract significance
#      res$p[2] < 0.01
#  }

## ----eval = F-----------------------------------------------------------------
#  res <- find.design(simfun = simfun_irt2,
#       boundaries = list(N1 = c(100,700), N2 = c(100,700)),
#       costfun = costfun_irt2,
#       power = .95)

## -----------------------------------------------------------------------------
#  library(mlpwr)
#  library(lme4)
#  library(lmerTest)

## -----------------------------------------------------------------------------
#  simfun_multi1 <- function(N) {
#  
#      # generate data
#      params <- list(theta = 0.5, beta = c(2, -0.2, -0.4,
#          -0.6)) # specifying the standard deviation of the random effects and parameter weights
#      dat <- expand.grid(herd = 1:ceiling(N/4), period = factor(1:4))[1:N,
#          ] # creating predictors
#      dat$x <- simulate(~period + (1 | herd), newdata = dat,
#          family = poisson, newparams = params)[[1]] # creating criterion
#  
#      # test hypothesis
#      mod <- glmer(x ~ period + (1 | herd), data = dat,
#          family = poisson)  # fit model
#      pvalues <- summary(mod)[["coefficients"]][2:4,
#          "Pr(>|z|)"] # extract p-values
#      any(pvalues < 0.01) # test hypothes that any is significant
#  }

## ----eval = F-----------------------------------------------------------------
#  res <- find.design(simfun = simfun_multi1,
#       boundaries = c(100, 500),
#       power = .95)

## -----------------------------------------------------------------------------
#  logistic <- function(x) 1/(1 + exp(-x))
#  
#  N.original <- 300
#  n.countries.original <- 20
#  
#  # generate original data
#  dat.original <- data.frame(country = rep(1:n.countries.original,
#      length.out = N.original), pred1 = rnorm(N.original),
#      pred2 = rnorm(N.original)) # creating predictors
#  country.intercepts <- rnorm(n.countries.original, sd = 0.5) # creating random intercepts
#  dat.original$intercepts <- country.intercepts[dat.original$country] # add interecepts to data
#  beta <- c(1, 0.4, -0.3)  # parameter weights
#  prob <- logistic(as.matrix(dat.original[c("intercepts",
#      "pred1", "pred2")]) %*% as.matrix(beta))  # get probability
#  dat.original$criterion <- runif(N.original) < prob  # draw according to probability
#  
#  # fit original model to obtain parameters
#  mod.original <- glmer(criterion ~ pred1 + pred2 + 0 +
#      (1 | country), data = dat.original, family = binomial)

## -----------------------------------------------------------------------------
#  simfun_multi2 <- function(n, n.countries) {
#  
#      # generate data
#      dat <- data.frame(country = rep(1:n.countries,
#          length.out = n * n.countries), pred1 = rnorm(n *
#          n.countries), pred2 = rnorm(n * n.countries))
#      dat$criterion <- simulate(mod.original, nsim = 1,
#          newdata = dat, allow.new.levels = TRUE, use.u = FALSE) |>
#          unlist()  # criterion data from the fitted model
#  
#      # test hypothesis
#      mod <- glmer(criterion ~ pred1 + pred2 + 0 + (1 |
#          country), data = dat, family = binomial)
#      summary(mod)[["coefficients"]]["pred2", "Pr(>|z|)"] <
#          0.01 # check if significant
#  }

## -----------------------------------------------------------------------------
#  costfun_multi2 <- function(n, n.countries) 5 * m +
#      100 * n.countries
#  

## ----eval = F-----------------------------------------------------------------
#  res <- find.design(simfun = simfun_multi2,
#       boundaries = list(n=c(10,40),n.countries=c(5,20)),
#       costfun = costfun_multi2,
#       power = .95)

