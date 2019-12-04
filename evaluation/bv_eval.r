#!/usr/bin/env Rscript

library(Hmisc)
library(readr)

script.dir <- dirname( sys.frame(1)$ofile )
base.dir <- file.path( getwd(), ".." )
# base.dir <- file.path( script.dir, ".." )
scrapp_measures <- read_csv( file.path( base.dir, "workdir", "scrapp", "result.csv" ) )
place_measures <- read_csv( file.path( base.dir, "workdir", "bwpd", "result.csv" ) )
guppy_measures <- read_csv( file.path( base.dir, "workdir", "guppy-fpd", "result.csv" ) )
colnames( guppy_measures ) <- paste( colnames( guppy_measures ), "guppy", sep=".")

meta <- read_delim( file.path( base.dir, "datasets", "bv", "00_data", "meta_simple.csv" ), "\t" )

measures <- merge( scrapp_measures, place_measures, by="sample", suffixes=c(".scrapp",".place") )
measures <- merge( measures, guppy_measures, by.x="sample", by.y="placerun.guppy"  )

# get some order in the df
library(dplyr)
measures <- measures[,-1] %>% select(sort(names(.))) %>% cbind(sample=measures[,1], .)

# correlation dendrogram
library(ggplot2)
library(ggdendro)
correl <- rcorr( as.matrix(measures[c(-1)]), type=c("pearson"))
correl$r <- 1 - correl$r
clust <- hclust( as.dist(correl$r) )
ggdendrogram( clust, rotate = TRUE )

# leave-one-out cross validation of linear regression of per-measure prediction of amsel criterion
measures_and_meta <- merge( measures, meta, by.x="sample", by.y="specimen" )

require(ISLR)
library(boot)

log_reg_results <- data.frame(measure=character(0), accuracy.amsel=numeric(0), r_squared.nugent=numeric(0), stringsAsFactors = FALSE)
for (m in colnames(measures[c(-1)])) {
  # log amsel model
  formula.log <- as.formula(sprintf("amsel ~ %s", m))
  print(formula.log)
  model.log <- glm( formula.log, data = measures_and_meta, family = binomial)
  validation <- cv.glm( measures_and_meta, model.log )
  accuracy.amsel <- 1 - validation$delta[1]
  # lin nugent model
  formula.lin <- as.formula(sprintf("nugent ~ %s", m))
  print(formula.lin)
  model.lin <- lm( formula.lin, data = measures_and_meta )
  r_squared.nugent <- summary(model.lin)$r.squared
  log_reg_results[nrow(log_reg_results) + 1,] = list(m, accuracy.amsel, r_squared.nugent)
}

glm.probs <- predict(model.log, measures_and_meta, type="response")
glm.pred <- ifelse(glm.probs > 0.5, 1, 0)
table(glm.pred,measures_and_meta$amsel)
mean(glm.pred == measures_and_meta$amsel)
