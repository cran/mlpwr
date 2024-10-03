## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----echo=FALSE---------------------------------------------------------------
set.seed(111)
dat.original <- data.frame(accidents = rep(1:3, 20) + sample(1:30, 20, replace = TRUE), road = gl(3, 1, 20), weather = gl(3,20))
dat.original[1:5,]

## -----------------------------------------------------------------------------
mod.original <- glm(accidents ~  road + weather, data = dat.original,
    family = poisson)
summary(mod.original)

## ----eval=FALSE---------------------------------------------------------------
#  install.packages("mlpwr")

## ----message=FALSE, warning=FALSE, results='hide'-----------------------------
library(mlpwr)

## -----------------------------------------------------------------------------
simfun_glm1 <- function(N) {

    # generate data
    dat <- data.frame(road = gl(3, 1, ceiling(N/3)),
        weather = gl(3, ceiling(N/3)))[1:N, ]
    a <- predict(mod.original, newdata = dat, type = "response")
    dat$accidents <- rpois(N, a)

    # test hypothesis
    mod <- glm(accidents ~ road + weather, data = dat,
        family = poisson)
    summary(mod)$coefficients["road3", "Pr(>|z|)"] <
        0.05
}

## ----echo=FALSE, results='hide'-----------------------------------------------
# The following loads the precomputed results of the next chunk to reduce the vignette creation time
ver <- as.character(packageVersion("mlpwr"))
file = paste0("/extdata/GLM_Vignette_results_", ver, ".RData")
file_path <- paste0(system.file(package="mlpwr"),file)
if (!file.exists(file_path)) {
set.seed(111)
res <- find.design(simfun = simfun_glm1, boundaries = c(10,1000),
                   power = .8, evaluations = 2000)
save(res, file = paste0("../inst",file))
} else {
  load(file_path) 
}

## ----warning=FALSE, eval=FALSE------------------------------------------------
#  set.seed(111)
#  res <- find.design(simfun = simfun_glm1, boundaries = c(10,1000),
#                     power = .8, evaluations = 2000)

## ----echo=TRUE----------------------------------------------------------------
summary(res)

## -----------------------------------------------------------------------------
plot(res)

