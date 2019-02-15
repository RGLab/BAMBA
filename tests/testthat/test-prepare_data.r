context("BAMBA Data Preparation")

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
    select(subjectId, group, re, ag, val, tp) %>%
    tbl_df()


## test column name requirement
bad_name_data <- simulated_data %>% rename(subject = subjectId)
missing_column_data <- simulated_data %>% select(-ag)
extra_column_data <- simulated_data %>% mutate(extraCol = 1:n())

test_that("The columns are checked and selected correctly", {
    expect_error( prepare_data(bad_name_data) )

    expect_error( prepare_data(missing_column_data) )
})

## test non-numeric/NA values
single_na_data <- simulated_data
single_na_data$val[10] <- NA

multiple_na_data <- simulated_data
multiple_na_data$val[c(10, 14, 23:28)] <- NA

non_numeric_data <- simulated_data
non_numeric_data$val[22] <- "bad"
    
na_non_numeric_data <- non_numeric_data
na_non_numeric_data$val[11] <- NA

test_that("The val column is properly checked and filtered", {
    nrow_sim <- nrow(simulated_data)

    prepData <- NULL

    expect_warning({prepData <- prepare_data(single_na_data)})
    expect_equal(nrow_sim - 1, nrow(prepData))

    expect_warning({prepData <- prepare_data(multiple_na_data)})
    expect_equal(nrow_sim - 8, nrow(prepData))

    expect_warning({prepData <- prepare_data(non_numeric_data)})
    expect_equal(nrow_sim - 1, nrow(prepData))

    expect_warning({prepData <- prepare_data(na_non_numeric_data)})
    expect_equal(nrow_sim - 2, nrow(prepData))
})

## test column addition
nogroup_data <- simulated_data %>% select(-group)
fixed_data <- NULL
test_that("The group column is added if not present", {
    expect_message({fixed_data = prepare_data(nogroup_data)})
    expect("group" %in% colnames(fixed_data))
})
    
nore_data <- simulated_data %>% select(-re)
test_that("The re column is added if not present", {
  expect_message({fixed_data = prepare_data(nore_data)})
    expect("re" %in% colnames(fixed_data))
})

notp_data <- simulated_data %>% select(-tp)
test_that("The tp column is added if not present", {
  expect_message({fixed_data = prepare_data(notp_data)})
  expect("tp" %in% colnames(fixed_data))
})

## test non-standard tp
non_numeric_tp_data <- simulated_data
non_numeric_tp_data$tp[3] <- "tp"

non_dense_data <- simulated_data
non_dense_data$tp[c(1,4,7)] <- "1.5"

test_that("The tp column is properly checked and adjusted", {
    expect_error( prepare_data(non_numeric_tp_data) )
    dense_data <- NULL
    expect_warning({dense_data <- prepare_data(non_dense_data)})
    dense_tps <- non_dense_data$tp
    dense_tps[dense_tps == 2] <- 3
    dense_tps[c(1,4,7)] <- 2
    non_dense_data$tp <- as.numeric(dense_tps)
    comp_data <- dense_data %>%
        left_join(non_dense_data, by = c("subjectId", "group", "re", "ag", "val"))
    expect_equal(comp_data$tp.x, comp_data$tp.y)
})

