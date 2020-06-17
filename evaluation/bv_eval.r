#!/usr/bin/env Rscript

library(Hmisc)
library(readr)

# script.dir <- dirname( sys.frame(1)$ofile )
base.dir <- file.path( getwd(), ".." )
# base.dir <- file.path( script.dir, ".." )
scrapp_measures <- read_csv( file.path( base.dir, "workdir", "bv", "scrapp", "result.csv" ) )
names(scrapp_measures) <- c( names(scrapp_measures[,1]), paste(names(scrapp_measures[,-1]), "scrapp", sep = ".") )

place_measures <- read_csv( file.path( base.dir, "workdir", "bv", "bwpd", "result.csv" ) )

guppy_measures <- read_csv( file.path( base.dir, "workdir", "bv", "guppy-fpd_pendant", "result.csv" ) )
colnames( guppy_measures ) <- paste( colnames( guppy_measures ), "guppy", sep=".")

otu_measures <- read_csv( file.path( base.dir, "workdir", "bv", "otu", "result.csv" ) )
otu_measures <- otu_measures[, c("Shannon", "Simpson", "ACE", "Chao1", "sample")]

meta <- read_delim( file.path( base.dir, "datasets", "bv", "data", "meta_simple.csv" ), "\t" )

measures <- scrapp_measures
# measures <- merge( measures, place_measures, by="sample", suffixes=c(".scrapp",".place") )
measures <- merge( measures, guppy_measures, by.x="sample", by.y="placerun.guppy")
measures <- merge( measures, otu_measures, by="sample")

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
library(data.table)

log_reg_results <- data.frame(measure=character(0), accuracy.amsel=numeric(0), r_squared.nugent=numeric(0), amsel.p.value=numeric(0), stringsAsFactors = FALSE)
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

  # two sample t-test
  formula.t.test <- as.formula(sprintf("%s ~ amsel", m))
  tt <- t.test(formula.t.test, data = measures_and_meta, var.equal = TRUE)

  # log results
  log_reg_results[nrow(log_reg_results) + 1,] = list(m, accuracy.amsel, r_squared.nugent, tt$p.value)
}
res_for_print <- log_reg_results[!(log_reg_results$measure %in% c("bwpd.guppy", "bwpd.scrapp")),]
res_for_print$rank <- rowMeans( cbind(frank(-res_for_print$accuracy.amsel), frank(-res_for_print$r_squared.nugent), frank(res_for_print$amsel.p.value)) )
res_for_print <- res_for_print[order(res_for_print$rank),]
res_for_print <- format(res_for_print, digits=3)
names(res_for_print) <- c("Measure", "Amsel accuracy", "Nugent $R^{2}$", "Amsel $p$-value", "Mean rank")


out.path <- file.path( base.dir, "evaluation", "bv_pend" )
if ( !dir.exists(out.path) ) {dir.create(out.path, recursive=TRUE)}
write.csv(res_for_print, file.path( out.path, "result.csv" ), row.names = FALSE)


# glm.probs <- predict(model.log, measures_and_meta, type="response")
# glm.pred <- ifelse(glm.probs > 0.5, 1, 0)
# table(glm.pred,measures_and_meta$amsel)
# mean(glm.pred == measures_and_meta$amsel, na.rm = TRUE)
