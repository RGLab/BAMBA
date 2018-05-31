context("BAMBA Model Data Building")

source("simulate_data.r")

## Make data up
nTp <- 2
nGrp <- 2
nSubj <- 2
nAg <- 2
nRe <- 2

agAlpha <- 1
agBeta <- 1
reAlpha <- 5
reBeta <- 2
sigmaAlpha <- 20
sigmaBeta <- 20

mu0 <- 0.5
groupProbs <- c(0.93, 0.97)
agProbs <- rbeta(nAg, 0.34, 0.16)
reProbs <- rbeta(nRe, 1.04, 0.18)

testdata <- simulate_fc_data(nTp, nGrp, nSubj, nAg, nRe,
                             agAlpha, agBeta, reAlpha, reBeta,
                             sigmaAlpha, sigmaBeta, mu0,
                             groupProbs, agProbs, reProbs)

simulated_data <- testdata$simData %>%
    select(subjectId, group, re, ag, val, tp)

prepared_data <- prepare_data(simulated_data, "fc")

## ## Test aspects of the list returned by build_model_data

model_data <- build_model_data(prepared_data, "fc")

test_that("The build_model_data function counts things correctly", {
    expect_equal(nTp * nGrp * nSubj * nAg * nRe, model_data$N)
    expect_equal(nTp, model_data$T)
    expect_equal(nGrp, model_data$N_grp)
    expect_equal(nAg, model_data$N_ag)
    expect_equal(nRe, model_data$N_re)
})
