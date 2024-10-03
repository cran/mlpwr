## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  tidy = TRUE,
  tidy.opts=list(arrow=TRUE,width.cutoff = 50),
  eval=T
)

## ----setup, echo = F----------------------------------------------------------
library(mlpwr)
set.seed(1)

## -----------------------------------------------------------------------------
data(extensions_results)

## -----------------------------------------------------------------------------
N = 100
alpha = .01
goal_power = .95

## -----------------------------------------------------------------------------
simfun_sensitivity = function(esize) {

  # Generate a data set
  dat <- rnorm(n = N, mean = esize)
  # Test the hypothesis
  res <- t.test(dat)
  res$p.value < alpha
}

## -----------------------------------------------------------------------------
# Effect Size Finding
# res1 <- find.design(simfun = simfun_sensitivity,
#                    boundaries = c(0,1),integer=FALSE, power = goal_power,surrogate = "gpr",evaluations=8000)
res1 = extensions_results[[1]]
summary(res1)

## -----------------------------------------------------------------------------
library(pwr)
esize_correct = pwr.t.test(n=N,sig.level=alpha,power=goal_power,type="one.sample")$d
esize_correct

## -----------------------------------------------------------------------------
a = replicate(10000,simfun_sensitivity(esize_correct))
mean(a)

## -----------------------------------------------------------------------------
N = 100
esize = .3
desired_ratio = 1

## -----------------------------------------------------------------------------
simfun_compromise = function(crit) {

  # Generate a data set
  dat <- rnorm(n = N, mean = 0)
 
  # Test the hypothesis
  res <- t.test(dat)
  a = res$statistic>crit
  
  # Generate a data set
  dat <- rnorm(n = N, mean = esize)
  
  # Test the hypothesis
  res <- t.test(dat)
  b = res$statistic<crit

  return(c(a,b))
}

## -----------------------------------------------------------------------------
simfun_compromise(.1)

## -----------------------------------------------------------------------------
aggregate_fun = function(x) {y=rowMeans(matrix(x,nrow=2));y[2]/y[1]}

## -----------------------------------------------------------------------------
# res2 <- find.design(simfun = simfun_compromise, boundaries = c(1,2), power = desired_ratio,integer = FALSE, aggregate_fun = aggregate_fun,surrogate="svr",use_noise=FALSE,evaluations=8000)
res2 = extensions_results[[2]]
summary(res2)

## -----------------------------------------------------------------------------
fun = \(alpha) {
  beta = alpha * desired_ratio
  abs(
    qt(1-alpha,ncp=0, df=N-1) # crit for alpha
    -qt(beta,ncp=esize*sqrt(N), df=N-1) # crit for beta
  )
}
alpha = optim(0.1,fun,lower=10e-10,upper=1-10e-10,method = "L-BFGS-B")$par
crit = qt(1-alpha,ncp=0, df=N-1)
crit

## -----------------------------------------------------------------------------
a = replicate(10000,simfun_compromise(crit))
aggregate_fun(a)

## -----------------------------------------------------------------------------
library(tidyr)
simfun_anova <- function(n, n.groups) {

  # Generate a data set
  groupmeans <- rnorm(n.groups, sd = 0.2)  # generate groupmeans using cohen's f=.2
  dat <- sapply(groupmeans, function(x) rnorm(n,
                                              mean = x, sd = 1))  # generate data
  dat <- dat |>
    as.data.frame() |>
    gather()  # format

  # Test the hypothesis
  res <- aov(value ~ key, data = dat)  # perform ANOVA
  summary(res)[[1]][1, 5] < 0.01  # extract significance
}


## -----------------------------------------------------------------------------
prices = c(rep(10,5),rep(15,5),rep(20,5),rep(25,5),rep(30,5),rep(35,5),rep(40,5))

costfun = function(n, n.groups) {
  5 * n + n.groups * sum(prices[1:n.groups])
}

## -----------------------------------------------------------------------------
# res3 <- find.design(
#   simfun = simfun_anova,
#   costfun = costfun,
#   boundaries = list(n = c(10, 150), n.groups = c(5, 30)),
#   power = .95
# )
res3 = extensions_results[[3]]
summary(res3)

## -----------------------------------------------------------------------------
library(lme4)
library(lmerTest)

## -----------------------------------------------------------------------------
simfun_3d <- function(n.per.school,n.schools, n.obs) {

  # generate data
  school = rep(1:n.schools,each=n.per.school*n.obs)
  student = rep(1:(n.schools*n.per.school),each=n.obs)
  pred = factor(rep(c("old","new"),n.per.school*n.schools*n.obs,each=n.obs),levels=c("old","new"))
  dat = data.frame(school = school, student = student, pred = pred)
  

  params <- list(theta = c(.5,0,.5,.5), beta = c(0,1),sigma = 1.5)
  names(params$theta) = c("school.(Intercept)","school.prednew.(Intercept)","school.prednew","student.(Intercept)")
  names(params$beta) = c("(Intercept)","prednew")
  dat$y <- simulate.formula(~pred + (1 + pred | school) + (1 | student), newdata = dat, newparams = params)[[1]]
  
  # test hypothesis
  mod <- lmer(y ~ pred + (1 + pred | school) + (1 | student), data = dat)
  pvalue <- summary(mod)[["coefficients"]][2,"Pr(>|t|)"]
  pvalue < .01
}

costfun_3d <- function(n.per.school, n.schools,n.obs) {
  100 * n.per.school + 200 * n.schools + .1 * n.obs * n.per.school * n.schools
  }

## -----------------------------------------------------------------------------
# res4 = find.design(simfun = simfun_3d, costfun = costfun_3d, boundaries = list(n.per.school = c(5, 25), n.schools = c(10, 30), n.obs = c(3,10)), power = .95)
res4 = extensions_results[[4]]
summary(res4)

