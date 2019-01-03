#' Adds parameter summary slots to a BAMBAResult object
#'
#' @param result The BAMBAResult object to fill with parameter summaries.
#'
#' @return The filled BAMBAResult object.
#'
#' @noRd
add_parameter_summaries <- function(result) {

    agtbl <- result$data %>% dplyr::select(ag, agId) %>% distinct()
    
    result$mu0 <- extract_params(result, "mu0") %>%
        mutate(tp = as.numeric(str_sub(var, 5, -2))) %>%
        rename(mu0 = mean) %>%
        dplyr::select(tp, mu0, sd, q025, q975, n_eff, Rhat)

    result$mu_ag <- extract_params(result, "mu_ag") %>%
        mutate(indices = str_sub(var, 7, -2)) %>%
        separate(indices, c("agId", "tp"), ",") %>%
        mutate(agId = as.numeric(agId), tp = as.numeric(tp)) %>%
        left_join(agtbl, by = "agId") %>%
        rename(mu_ag = mean) %>%
        dplyr::select(ag, agId, tp, mu_ag, sd, q025, q975, n_eff, Rhat)

    result$omega_t <- extract_params(result, "omega_t") %>%
        mutate(tp = as.numeric(str_sub(var, 9, -2))) %>%
        rename(omega_t = mean) %>%
        dplyr::select(tp, omega_t, sd, q025, q975, n_eff, Rhat)

    result$omega_ag <- extract_params(result, "omega_ag") %>%
        mutate(agId = as.numeric(str_sub(var, 10, -2))) %>%
        rename(omega_ag = mean) %>%
        left_join(agtbl, by = "agId") %>%
        dplyr::select(ag, agId, omega_ag, sd, q025, q975, n_eff, Rhat)

    result$hyperparameters <- result$parameters %>%
        dplyr::filter(var %in% c("alpha_ag", "beta_ag", "alpha_re", "beta_re",
                                 "phi_grp", "lambda_grp", "phi_ag", "lambda_ag",
                                 "phi_re", "lambda_re", "phi_t", "lambda_t",
                                 "a_grp", "b_grp", "a_ag", "b_ag",
                                 "a_re", "b_re", "a_t", "b_t")) %>%
        rename(hyperparameter = var) %>%
        dplyr::select(hyperparameter, mean, sd, q025, q975, n_eff, Rhat)

    if ( result$dataType == "fc" ) {
        retbl <- result$data %>% dplyr::select(re, reId) %>% distinct()
        grptbl <- result$data %>% dplyr::select(group, groupId) %>% distinct()
        
        result$mu_re <- extract_params(result, "mu_re") %>%
            mutate(indices = str_sub(var, 7, -2)) %>%
            separate(indices, c("reId", "tp"), ",") %>%
            mutate(reId = as.numeric(reId), tp = as.numeric(tp)) %>%
            left_join(retbl, by = "reId") %>%
            rename(mu_re = mean) %>%
            dplyr::select(re, reId, tp, mu_re, sd, q025, q975, n_eff, Rhat)

        result$mu_ar <- extract_params(result, "mu_ar") %>%
            mutate(indices = str_sub(var, 7, -2)) %>%
            separate(indices, c("tp", "agId", "reId"), ",") %>%
            mutate(agId = as.numeric(agId),
                   reId = as.numeric(reId),
                   tp = as.numeric(tp)) %>%
            left_join(agtbl, by = "agId") %>%
            left_join(retbl, by = "reId") %>%
            rename(mu_ar = mean) %>%
            dplyr::select(ag, agId, re, reId, tp,
                   mu_ar, sd, q025, q975, n_eff, Rhat)

        result$omega_re <- extract_params(result, "omega_re") %>%
            mutate(reId = as.numeric(str_sub(var, 10, -2))) %>%
            rename(omega_re = mean) %>%
            left_join(retbl, by = "reId") %>%
            dplyr::select(re, reId, omega_re, sd, q025, q975, n_eff, Rhat)

        result$omega_grp <- extract_params(result, "omega_grp") %>%
            mutate(groupId = as.numeric(str_sub(var, 11, -2))) %>%
            rename(omega_grp = mean) %>%
            left_join(grptbl, by = "groupId") %>%
            dplyr::select(group, groupId, omega_grp, sd, q025, q975, n_eff, Rhat)


    }

    result
}
