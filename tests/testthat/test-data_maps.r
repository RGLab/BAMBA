context("BAMBA Data Maps")

## Test antigen map
ags <- paste0("ag", 1:5)
agData <- data_frame(ag = sample(rep(ags, 5), 25))
agMap <- data_frame(ag = ags, agId = 1:5)
test_that("The antigen_map function works appropriately", {
    expect_equal(antigen_map(agData), agMap)
})

## Test reagent map
res <- paste0("re", 1:5)
reData <- data_frame(re = sample(rep(res, 5), 25))
reMap <- data_frame(re = res, reId = 1:5)
test_that("The reagent_map function works appropriately", {
    expect_equal(reagent_map(reData), reMap)
})

## Test group map
groups <- 1:5
groupData <- data_frame(group = sample(rep(groups, 5), 25))
groupMap <- data_frame(group = as.character(groups), groupId = 1:5)
test_that("The group_map function works appropriately", {
    expect_equal(group_map(groupData), groupMap)
})

