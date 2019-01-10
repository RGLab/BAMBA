bamaData <- prepare_bama_data()

fit <- BAMBA(bamaData,
              nChains = 1,
              nIter = 200, ## set higher for a real analysis
              outFolder = NULL)
