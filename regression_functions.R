
library(Hmisc)
library(jtools)
library(car)
library(psych)
library(labelled)
library(knitr)
library(kableExtra)
library(tidyverse)
library(qqplotr)
source("../statcan-data-input/data_from_txt_and_sas.R")


# 
# # raw CTNS PUMF
# ctns <- load_statcan_data(path_txt_data,
#                           path_i,
#                           path_lb,
#                           path_pf,
#                           path_fmt)
# 
# 
# 
# dat <- ctns %>%
#   mutate(across(TBC_30AR:TBC_30GR,
#                 ~ case_when(.x == 0 ~ 0,
#                             .x > 95 ~ NA,
#                             TRUE ~ .x*5-2)), # tbc from ordinal to continuous (mean)
#          agerange = case_when(AGEGROUP < 3 ~ 1,
#                               TRUE ~ AGEGROUP - 1),
#          otp_any = case_when(if_any(OTP_05A:OTP_05ER, ~ .x < 4) ~ 1,
#                              if_any(OTP_05A:OTP_05ER, ~ .x == 4) ~ 0)) %>%
#   rowwise() %>%
#   mutate(cigavg = mean(c_across(TBC_30AR:TBC_30GR),
#                             na.rm = TRUE)) %>%
#   ungroup() %>%
#   mutate(across(all_of(c("DV_AL30R",
#                          "DV_CN30R",
#                          "DV_VC30R",
#                          "DV_VP30R",
#                          "FIRSTTRR",
#                          "GENDER",
#                          "PROV_C")),
#                 as_factor),
#          cigavg = ifelse(is.nan(cigavg), NA, cigavg), # NAs
#          across(all_of(c("DV_AL30R",
#                          "DV_CN30R",
#                          "DV_VC30R",
#                          "DV_VP30R",# exclude missings
#                          "FIRSTTRR",
#                          "HHLDSIZE",
#                          "GENDER")),
#                 ~ case_when(.x %in% c("Valid skip",
#                                      "Don't know",
#                                      "Refusal",
#                                      "Not stated") ~ NA,
#                             TRUE ~ .x)),
#          otp_any = fct_recode(factor(otp_any),
#                               Yes = "1",
#                               No = "0"),
#          across(where(is.factor), fct_drop)) %>%
#   select(cigavg,
#          FIRSTTRR,
#          DV_AL30R,
#          DV_CN30R,
#          DV_VC30R,
#          DV_VP30R,
#          otp_any,
#          HHLDSIZE,
#          GENDER,
#          PROV_C,
#          agerange)
# 
# #across(where(is.factor), fct_infreq),

freq_table <- function(dat, v) {
  freq_tab <- table(dat[[v]], useNA = "always")
  freqs <- as.data.frame(freq_tab)$Freq
  valid_freq_tab <- table(dat[[v]])
  valid_freqs <- as.data.frame(valid_freq_tab)$Freq
  listwise_freq_tab <- table(drop_na(dat)[[v]])
  listwise_freqs <- as.data.frame(listwise_freq_tab)$Freq
  out <- tibble(
    value = names(freq_tab),
    frequency = freqs,
    proportion = prop.table(freqs) * 100,
    cumulative = cumsum(freqs),
    cumulative_proportion = (cumsum(freqs) / sum(freqs)) * 100,
    valid_frequency = c(valid_freqs, NA),
    valid_proportion = c(prop.table(valid_freqs) * 100, NA),
    valid_cumulative = c(cumsum(valid_freqs), NA),
    valid_cumulative_proportion = 
      c((cumsum(valid_freqs) / sum(valid_freqs)) * 100, NA),
    listwise_frequency = c(listwise_freqs, NA),
    listwise_loss = c(listwise_freqs - valid_freqs, NA),
    listwise_carryover = c((listwise_freqs / valid_freqs) * 100, NA),
    listwise_proportion = c(prop.table(listwise_freqs) * 100, NA),
    listwise_cumulative = c(cumsum(listwise_freqs), NA),
    listwise_cumulative_proportion = 
      c((cumsum(listwise_freqs) / sum(listwise_freqs)) * 100, NA)
  ) %>%
    bind_rows(summarise(
      .,
      across(c(-value, -contains("cumulative")), \(x) round(sum(x, na.rm = TRUE))),
      value = "Total",
      listwise_loss = sum(listwise_freqs) - sum(valid_freqs),
      listwise_carryover = (sum(listwise_freqs) / sum(valid_freqs)) * 100,
      across(contains("cumulative"), \(x) max(x, na.rm = TRUE))
    )) 
  return(out)
}

freq_tables <- function(dat) {
  cat <-
    bind_rows(sapply(names(select(dat, where(is.factor))), 
                     \(x) freq_table(dat, x), 
                     simplify = FALSE)) %>% 
    mutate(variable = unlist(
      sapply(
        names(select(dat, where(is.factor))), 
        \(x) rep(x, length(levels(dat[[x]])) + 2))
    )) %>% 
    relocate(variable)
  
  return(list(cat = cat, dat = dat))
}

print_freq_table <- function(freq_tables) {
  cat <-  freq_tables$cat
  dat <-  freq_tables$dat
  # create group index of variable names: # of lines
  index <- setNames(unlist(sapply(unique(cat$variable),
                                  \(x) nrow(filter(cat, variable == x)))), 
                    
                    unlist(sapply(unique(cat$variable),
                                  \(x) if (!is.null(var_label(dat[[x]])))
                                    paste(x, var_label(dat[[x]]), sep = " - ")
                                  else
                                    x)))
  
  cat_tbl <- cat[-1] %>% # drop variable from output
    kbl(col.names = str_to_title(gsub("[_]",
                                      " ",
                                      names(.))),
        booktabs = TRUE,
        digits = 2, 
        format.args = list(big.mark = ",",
                           scientific = FALSE),
        caption = "Descriptive Statistics for Categorical Variables") %>% 
    pack_rows(index = index) %>%  # group lines
    add_header_above(c(" ", 
                       "All Cases" = 4, 
                       "Valid Cases" = 4, 
                       "Listwise Valid Cases" = 6)) %>%
    row_spec(cumsum((index)), bold = T) %>% 
    row_spec(cumsum((index)) - 1, italic = T, bold = F) %>%
    row_spec(0, align = "c", bold = T) %>% 
    kable_classic(full_width = F, html_font = "Cambria", fixed_thead = T) 
  
  return(cat_tbl)
}

quant_describe <- function(dat) {
  out <-
    describe(
      dat,
      omit = TRUE,
      quant = c(.25, .5, .75),
      IQR = TRUE,
      skew = FALSE
    ) %>% as_tibble() %>% 
    rename(
      valid_frequency = n,
      variable = vars,
      Q1 = Q0.25,
      median = Q0.5,
      Q3 = Q0.75
    ) %>% mutate(
      frequency = nrow(dat),
      pairwise_valid_frequency = as_tibble(describe(dat, na.rm = FALSE, omit = TRUE))$n,
      variable = names(dat)[variable]
    ) %>% select(
      variable,
      frequency,
      valid_frequency,
      pairwise_valid_frequency,
      mean,
      se,
      sd,
      min,
      Q1,
      median,
      Q3,
      max,
      IQR
    )
  return(out)
}

print_quant <- function(quant) {
  quant_tbl <- quant %>% 
    kbl(
      col.names = gsub("^([SI][edq]r?)$", 
                       "\\U\\1", 
                       str_to_title(gsub("[_]",
                                         " ",
                                         names(quant))), 
                       perl = TRUE),
      booktabs = TRUE,
      digits = 2, 
      format.args = list(big.mark = ",",
                         scientific = FALSE),
      caption = "Descriptive Statistics for Numerical Variables") %>% 
    row_spec(0, align = "c", bold = T) %>% 
    kable_classic(full_width = F, html_font = "Cambria", fixed_thead = T)
  
  return(quant_tbl)
}

descriptives <- function(dat) {
  quant <- quant_describe(dat)
  cat <- freq_tables(dat)$cat
  
  return(list(cat = cat, quant = quant, dat = dat))
  
}

print_descriptives <- function(desc_in) {
  cat_tbl <- print_freq_table(list(cat = desc_in$cat, dat = desc_in$dat))
  quant_tbl <- print_quant(desc_in$quant)
  
  print(cat_tbl)  
  cat('\r\n\r\n\r\n\r\n')
  print(quant_tbl)
}

dummies <- function(dat){
  out <- model.matrix(
    ~ .,
    data = dat,
    contrasts.arg = lapply(dat[, 
                               sapply(dat, is.factor), 
                               drop = FALSE],
                           contrasts, contrasts = FALSE))
  
  return(out)
}

dummy_names <- function(dat) {
  detail = sapply(names(dat), 
                  \(x) sapply(levels(dat[[x]]), 
                              \(y) paste(x, y, sep = " - ")))
  detail[sapply(names(dat), 
                \(x) length(levels(dat[[x]]))) == 0] <-
    names(dat)[sapply(names(dat), 
                      \(x) length(levels(dat[[x]]))) == 0]
  detail <- unlist(detail)
  names(detail) <- NULL
  
  return(detail)
}

correlation_table <- function(dat) {
  model_all_factors <- dummies(dat)
  correlations <- rcorr(model_all_factors)
  cor_coef <- correlations$r[-1,-1]
  cor_p <- correlations$P[-1,-1]
  mycor <- data.frame(matrix(nrow = nrow(cor_coef), ncol = ncol(cor_coef)))
 
   for (i in 1:nrow(cor_coef)) {
    for (j in 1:ncol(cor_coef)) {
      if (!is.na(cor_p[i,j]) & cor_p[i,j] < 0.001) p <- "***" else
        if (!is.na(cor_p[i,j]) & cor_p[i,j] < 0.01) p <- "**" else
          if (!is.na(cor_p[i,j]) & cor_p[i,j] < 0.05) p <- "*" else 
            if (!is.na(cor_p[i,j]) & cor_p[i,j] < 0.1) p <- "'" else
              p <- ""
            mycor[i,j] <- paste0(round(cor_coef[i,j], 3), p)
    }
  }
  
  names(mycor) <- names(cor_coef[,1])
  mycor$var <- dummy_names(dat)
  mycor <- relocate(mycor, var)
  
  return(mycor)
}

print_correlations <- function(cor) {
  cor_tbl <- cor %>% 
    kbl(
      col.names = c("Variable", cor$var),
      booktabs = TRUE,
      caption = "Pearson Correlation Coefficients") %>% 
    row_spec(0, align = "l", bold = F, font_size = 10) %>% 
    kable_classic(full_width = F, html_font = "Cambria", fixed_thead = T, font_size = 10) %>%
    footnote(general = "\\*\\*\\*p < 0.001; \\*\\*p<0.01; \\*p<0.05; 'p<0.1") %>% 
    column_spec(1, border_right = T)
  
  print(cor_tbl)
}

regress_and_check <- function(formula, dat) {

  dat_w_dummies <- as.data.frame(dummies(dat)[-1,-1]) 
  
  descr_stats <- quant_describe(dat_w_dummies)
  
  corr_tbl <- correlation_table(dat_w_dummies)
  
  model <- lm(formula, dat)
  
  model_anova <- anova(model)
  
  model_summary <- summary(model, correlation = TRUE)
  
  model_summary2 <- summ(model, 
                         vifs = TRUE, 
                         part.corr = TRUE, 
                         confint = TRUE)
  
  model_summary3 <- summ(model, robust = TRUE)
  
  model_zsummary <- summ(model, 
                         scale = TRUE,
                         confint = TRUE)
  model_zsummary2 <- summ(model,
                          scale = TRUE)
  model_zsummary3 <- summ(model,
                          scale = TRUE,
                          robust = TRUE)
  
  fstat_p <- pf(model_summary$fstatistic[1], 
                model_summary$fstatistic[2], 
                model_summary$fstatistic[3], 
                lower.tail = FALSE)
  
  if (fstat_p < 0.001) fstat_sig <- "***" else
    if (fstat_p < 0.01) fstat_sig <- "**" else
      if (fstat_p < 0.05) fstat_sig <- "*" else
        if (fstat_p < 0.1) fstat_sig <- "'" else
          fstat_sig <- ""
  
  RSE <- sqrt(deviance(model)/df.residual(model))
  
  coef <- data.frame(
    coefficient = model$coefficients,
    sig = case_when(model_summary$coefficients[,4] < 0.001 ~ "***",
                    model_summary$coefficients[,4] < 0.01 ~ "**",
                    model_summary$coefficients[,4] < 0.05 ~ "*",
                    model_summary$coefficients[,4] < 0.1 ~ "'",
                    TRUE ~ ""),
    SE = model_summary$coefficients[,2],
    robust_SE = model_summary3$coeftable[,2],
    CILB = model_summary2$coeftable[,2],
    CIUB = model_summary2$coeftable[,3],
    t = model_summary2$coeftable[,4],
    scaled_coef = model_zsummary$coeftable[,1],
    scaled_SE = model_zsummary2$coeftable[,2],
    scaled_robust_SE = model_zsummary3$coeftable[,2],
    scaled_CILB = model_zsummary2$coeftable[,2],
    scaled_CIUB = model_zsummary2$coeftable[,3],
    p = model_summary$coefficients[,4],
    correlation = model_summary$correlation[,1],
    partial_correlation = model_summary2$coeftable[,7],
    semipartial_correlation = model_summary2$coeftable[,8],
    vif = model_summary2$coeftable[,6],
    tolerance = 1/model_summary2$coeftable[,6],
    covariance_intercept = model_summary$cov.unscaled[,1]
  )
  
  fit_stats <- c(r2 = model_summary$r.squared,
                 adj.r2 = model_summary$adj.r.squared,
                 f = model_summary$fstatistic[1],
                 sig = fstat_sig,
                 p = fstat_p,
                 residual_standard_error = RSE)
  
  
  dw = durbinWatsonTest(model)
  
  influence = influence.measures(model)
  
  zvals = zvals(model)
  
  return(list(
    model = model,
    descriptives = descr_stats,
    anova = model_anova,
    coefficients = coef,
    fit = fit_stats))
}

zvals <- function(model){
  zvals = data.frame(zresid = rstandard(model),
                     zpred = scale(model$fitted.values)) %>% 
    mutate(pdist = pnorm(zresid))
  return(zvals)
}

residual_plots <- function(model){
  
  zpred_zresid = zvals(model)
  
  zpred_plot <- ggplot(data = zpred_zresid, aes(zpred, zresid)) +
    geom_point() +
    geom_abline(aes(intercept = 0, slope = 0)) + 
    labs(x = "Standardized Predicted Values", 
         y = "Standardized Residuals",
         title = "Residuals against Predicted Values")
  
  zresid_plot <- ggplot(data = zpred_zresid, aes(zresid)) +
    geom_histogram(binwidth = 0.3) +
    stat_function(
      fun = function(x)
        dnorm(x, mean = 0, sd = 1) * length(zpred_zresid$zresid) * 0.3
    ) +
    labs(x = "Standardized Residuals",
         y = NULL,
         title = "Distribution of Residuals") +
    scale_y_discrete(labels = NULL, breaks = NULL)
  
  pp_plot <- ggplot(data = zpred_zresid, mapping = aes(sample = zresid)) +
    stat_pp_band() +
    stat_pp_line() +
    stat_pp_point() +
    labs(x = "Observed Probability", 
         y = "Expected Probability",
         title = "P-P Plot of Standardized Residuals")
  
  print(zpred_plot)
  print(zresid_plot)
  print(pp_plot)
  

}

partial_plots <- function(model, layout = NA){
  crPlots(model, 
          layout = layout, 
          main = "Partial Plots", 
          ask = FALSE)
}

allplots <- function(model, partial_layout = NA) {
  residual_plots(model)
  partial_plots(model, partial_layout)
}
