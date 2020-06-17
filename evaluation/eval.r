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
# otu_measures <- read_csv( file.path( base.dir, "workdir", ds, "otu", "result.csv" ) )
colnames( guppy_measures ) <- paste( colnames( guppy_measures ), "guppy", sep=".")

# read, clean up the metadata table a bit to make it less unwieldy
keep_names <- c("External ID", "diagnosis", "Subject has an acute gastrointestinal illness:", "Subject has an acute gastrointestinal infection or IBD:"
                , "is_inflamed", "CRP (mg/L)", "VSL #3", "Blood in the stool")
meta <- read_delim( file.path( base.dir, "datasets", ds, "data", "hmp2_metadata.csv" ), "," )
# [ ,keep_names]
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
measures_and_meta <- Filter(function(x)!all(is.na(x)), measures_and_meta)
# add a binary IBD column
measures_and_meta <- transform(measures_and_meta, ibd= ifelse(diagnosis=="nonIBD", as.integer(0), as.integer(1)))
measures_and_meta <- transform(measures_and_meta, is_inflamed_bin= ifelse(is_inflamed=="No", as.integer(0), as.integer(1)))
# measures_and_meta <- transform(measures_and_meta, uvsc= ifelse(diagnosis=="nonIBD", as.integer(0), as.integer(1)))

require(ISLR)
library(boot)
library(stats)
library(nnet)
library(data.table)
library(caret)


results <- data.frame(measure=character(0),  accuracy.log=numeric(0), ttest.pvalue=numeric(0), mn.acc=numeric(0), mn.pvalue=numeric(0), stringsAsFactors = FALSE)
for (m in colnames(measures[c(-1)])) {
  # m<-"bwpd_0.place"
  # log model
  response <- "diagnosis"
  formula.log <- as.formula(sprintf("%s ~ %s", response, m))
  print(formula.log)
  measures_rarefied <- measures_and_meta %>% 
    filter( diagnosis %in% c("UC", "CD", "nonIBD") ) %>%
    mutate( diagnosis = ifelse( diagnosis=="nonIBD",as.integer(0), as.integer(1)) ) %>%
    group_by( .dots=response ) %>% 
    sample_n( min( table(unlist(.[[response]]) ) ) )
    # glm( formula.log,  family=binomial, . )
  model.log <- glm( formula.log, data = measures_rarefied, family = binomial)
  validation <- cv.glm( measures_rarefied, model.log )
  accuracy.log <- 1 - validation$delta[1]
  
  # ttest
  formula.t.test <- as.formula(sprintf("%s ~ %s", m, response))
  tt <-  measures_and_meta %>% 
    mutate( diagnosis = ifelse( diagnosis=="nonIBD",as.integer(0), as.integer(1)) ) %>%
    t.test(formula.t.test, data=., var.equal = TRUE)

  # linear regression
  # response <- "biopsy_location"
  # formula.lin <- as.formula(sprintf("%s ~ %s", response, m))
  # print(formula.lin)
  # model.lin <- measures_and_meta %>%
  #   mutate( biopsy_location := case_when(
  #     biopsy_location=="Rectum" ~ 0,
  #     biopsy_location=="Illeum" ~ 4,
  #     biopsy_location=="Sigmoid Colon" ~ 1,
  #     biopsy_location=="Descending (left-sided) colon" ~ 2,
  #     biopsy_location=="Transverse colon" ~ 3,
  #     TRUE ~ NA_real_
  #   ) )%>%
  #   lm( formula.lin, data=., na.action=na.exclude )
  # r_squared <- summary(model.lin)$r.squared
  
  
  # kmeans clustering
  # cl <- kmeans(measures_and_meta[m], 5)
  # r_squared <- cl$betweenss / cl$totss
  #
  
  # multinomial logistic regression
  response <- "diagnosis"
  mn.formula <- as.formula(sprintf("%s ~ %s", response, m))
  print(mn.formula)
  
  set.seed(123)
  train.smpls <- measures_and_meta[[response]] %>% 
    createDataPartition(p = 0.8, list = FALSE)
  train.data  <- measures_and_meta[train.smpls, ]
  test.data <- measures_and_meta[-train.smpls, ]
  
  mn <- multinom(mn.formula, data = train.data)
  mn.z <- summary(mn)$coefficients/summary(mn)$standard.errors
  mn.pvalue <- (1 - pnorm(abs(mn.z), 0, 1)) * 2
  pred <- mn %>% predict(test.data)
  mn.acc <- mean(pred == test.data[[response]], na.rm = TRUE )
  
  results[nrow(results) + 1,] = list(m, accuracy.log, tt$p.value, mn.acc, mean(mn.pvalue[,m]) )
  # log_reg_results[nrow(log_reg_results) + 1,] = list(m, r_squared)
}
results$rank <- rowMeans( frank(-results$accuracy.log), frank(results$ttest.pvalue), frank(-mn.acc) )



