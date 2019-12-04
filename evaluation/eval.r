#!/usr/bin/env Rscript

library(Hmisc)
library(readr)

# ds <- "bv"
ds <- "hmp2"

script.dir <- dirname( sys.frame(1)$ofile )
base.dir <- file.path( getwd(), ".." )
# base.dir <- file.path( script.dir, ".." )
scrapp_measures <- read_csv( file.path( base.dir, "workdir", ds, "scrapp", "result.csv" ) )
place_measures <- read_csv( file.path( base.dir, "workdir", ds, "bwpd", "result.csv" ) )
guppy_measures <- read_csv( file.path( base.dir, "workdir", ds, "guppy-fpd", "result.csv" ) )
colnames( guppy_measures ) <- paste( colnames( guppy_measures ), "guppy", sep=".")

# read, clean up the metadata table a bit to make it less unwieldy
keep_names <- c("External ID", "diagnosis", "Subject has an acute gastrointestinal illness:", "Subject has an acute gastrointestinal infection or IBD:"
                , "is_inflamed", "CRP (mg/L)", "VSL #3", "Blood in the stool")
meta <- read_delim( file.path( base.dir, "datasets", ds, "data", "hmp2_metadata.csv" ), "," )[ ,keep_names]
names(meta) <- make.names( names(meta) )

measures <- merge( scrapp_measures, place_measures, by="sample", suffixes=c(".scrapp",".place") )
measures <- merge( measures, guppy_measures, by.x="sample", by.y="placerun.guppy"  )

measures$sample <- sub( ".fasta", "", measures$sample, fixed=TRUE )

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
measures_and_meta <- merge( measures, meta, by.x="sample", by.y="External.ID" )
# add a binary IBD column
measures_and_meta <- transform(measures_and_meta, ibd= ifelse(diagnosis=="nonIBD", as.integer(0), as.integer(1)))
measures_and_meta <- transform(measures_and_meta, is_inflamed_bin= ifelse(is_inflamed=="No", as.integer(0), as.integer(1)))

require(ISLR)
library(boot)

results <- data.frame(measure=character(0),  accuracy.log=numeric(0), stringsAsFactors = FALSE)
response <- "ibd"
for (m in colnames(measures[c(-1)])) {
  m<-"bwpd_0.place"
  # log model
  formula.log <- as.formula(sprintf("%s ~ %s", response, m))
  print(formula.log)
  # rarefy
  # measures_rarefied <- measures_and_meta %>% group_by( ibd ) %>% sample_n( min( table(unlist(measures_and_meta$ibd) ) ) )
  
  measures_rarefied <- measures_and_meta %>% 
    group_by( .dots=response ) %>% 
    sample_n( min( table(unlist(measures_and_meta[[response]]) ) ) )
    # glm( formula.log,  family=binomial, . )
  model.log <- glm( formula.log, data = measures_rarefied, family = binomial)
  validation <- cv.glm( measures_rarefied, model.log, )
  accuracy.log <- 1 - validation$delta[1]
  results[nrow(results) + 1,] = list(m, accuracy.log)
  break
  # lin nugent model
  # formula.lin <- as.formula(sprintf("CRP..mg.L. ~ %s", m))
  # print(formula.lin)
  # model.lin <- lm( formula.lin, data = measures_and_meta )
  # r_squared <- summary(model.lin)$r.squared
  # log_reg_results[nrow(log_reg_results) + 1,] = list(m, r_squared)
}

glm.probs <- predict(model.log, measures_and_meta, type="response")
glm.pred <- ifelse(glm.probs > 0.5, 1, 0)
table(glm.pred,measures_and_meta$ibd)
mean(glm.pred == measures_and_meta$ibd)


measures_diag <- measures_and_meta[, c("bwpd_0.scrapp", "diagnosis", "ibd")]
ggplot(aes(ibd, bwpd_0.scrapp, color=diagnosis), data = measures_diag)+
  # geom_point()+
  geom_jitter(position = position_jitter(width = 0.05))+
  scale_x_continuous(breaks=c(0,1),labels=c("No", "Yes"))


ggplot(aes(is_inflamed_bin, bwpd_0.scrapp, color=is_inflamed), data = measures_and_meta[, c("bwpd_0.scrapp", "is_inflamed", "is_inflamed_bin")])+
  geom_jitter(position = position_jitter(width = 0.05))+
  scale_x_continuous(breaks=c(0,1),labels=c("No", "Yes"))

