
library(Hmisc)
library(jtools)
library(car)
library(psych)
library(labelled)
library(knitr)
library(kableExtra)
library(tidyverse)
source("../statcan-data-input/data_from_txt_and_sas.R")



# raw CTNS PUMF
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

#across(where(is.factor), fct_infreq),

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



descriptives <- function(dat) {
  quant <- quant_describe(dat)
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
  
  return(list(cat = cat, quant = quant))
  
}

print_descriptives <- function(desc_in, dat) {

  # create group index of variable names: # of lines
  index <- setNames(unlist(sapply(unique(desc_in$cat$variable),
                                  \(x) nrow(filter(desc_in$cat, variable == x)))), 
                    
                    unlist(sapply(unique(desc_in$cat$variable),
                      \(x) if (!is.null(var_label(dat[[x]])))
                        paste(x, var_label(dat[[x]]), sep = " - ")
                      else
                        x)))

  cat_tbl <- desc_in$cat[-1] %>% # drop variable from output
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
  
  quant_tbl <- desc_in$quant %>% 
    kbl(
      col.names = gsub("^([SI][edq]r?)$", 
                       "\\U\\1", 
                       str_to_title(gsub("[_]",
                                         " ",
                                         names(desc_in$quant))), 
                       perl = TRUE),
      booktabs = TRUE,
      digits = 2, 
      format.args = list(big.mark = ",",
                         scientific = FALSE),
      caption = "Descriptive Statistics for Numerical Variables") %>% 
    row_spec(0, align = "c", bold = T) %>% 
    kable_classic(full_width = F, html_font = "Cambria", fixed_thead = T)
  
  print(cat_tbl)  
  cat('\r\n\r\n\r\n\r\n')
  print(quant_tbl)
}

regress_and_check <- function(data, formula) {
  
  model_all_factors <- 
    model.matrix(
      ~ .,
      data = dat,
      contrasts.arg = lapply(dat[, 
                                 sapply(dat, is.factor), 
                                 drop = FALSE],
                             contrasts, contrasts = FALSE))
  
  descriptives <- psych::describe(as.data.frame(model_all_factors))
  
  correlations <- rcorr(model_all_factors)
  
  fit <- lm(formula, data)
  
  fit_anova <- anova(fit)
  
  fit_summary <- summary(fit, correlation = TRUE)
  
  fit_summary2 <- summ(fit, 
                       vifs = TRUE, 
                       part.corr = TRUE, 
                       confint = TRUE)
  
  fit_zsummary <- summ(fit, 
                       scale = TRUE,
                       vifs = TRUE, 
                       part.corr = TRUE, 
                       confint = TRUE)
  
  
  
}
