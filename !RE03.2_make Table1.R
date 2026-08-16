source('functions.R')
library(table1)
library(dplyr)

setwd('../')

mg = read.csv('data/Tables/input_Table1.csv', check.names=F, row.names=1)

# clean mg
mg <- mg %>%
  mutate(
    `Race` = factor(`Race/Ethnicity`, levels = c("NHW", "NHB")),
    `Disease State` = factor(`Disease State`, levels = c("HRL", "DCIS", "BC")),
    Diabetes = factor(Diabetes, levels = c("N", "Y")),
    Obesity = factor(Obesity, levels = c("N", "Y")),
    Subtype = factor(Subtype),
    ER = factor(ER),
    PR = factor(PR),
    HER2 = factor(HER2),
    T = factor(T),
    N = factor(N),
    Grade = factor(Grade)
  )

pvalue <- function(x, ...) {
    # Construct vectors of data y, and groups (strata) g
    y <- unlist(x)
    g <- factor(rep(1:length(x), times=sapply(x, length)))
    if (is.numeric(y)) {
        # For numeric variables, perform a standard 2-sample t-test
        p <- t.test(y ~ g)$p.value
    } else {
        # For categorical variables, perform a chi-squared test of independence
        p <- chisq.test(table(y, g))$p.value
    }
    # Format the p-value, using an HTML entity for the less-than sign.
    # The initial empty string places the output on the line below the variable label.
    c("", sub("<", "&lt;", format.pval(p, digits=3, eps=0.001)))
}

tbl_patient <- table1(
  ~ `Age at Diagnosis` +
    `Disease State` +
    Diabetes +
    Obesity 
  | `Race`,
  data = mg,
  overall = F,
  extra.col = list(`P value` = pvalue)
)

tbl_patient