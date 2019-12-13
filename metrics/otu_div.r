#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
if ( length(args) != 1 ) {stop("usage: otu_div.r <dataset>")}

library(readr)
#source('http://bioconductor.org/biocLite.R')
#biocLite('phyloseq')
library(phyloseq)

initial.options <- commandArgs(trailingOnly = FALSE)
file.arg.name <- "--file="
script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
script.dir <- dirname(script.name)
base.dir <- normalizePath( file.path( script.dir, ".." ) )

ds <- args[1]

data.dir <- normalizePath( file.path( base.dir, "datasets", ds ) )
if ( !dir.exists(data.dir) ) {stop("No such dataset")}

# read and prepare the otu table
p <-read_delim(file.path( base.dir, "datasets", ds, "otu", "otu.table"),"\t")
p <- subset(p, select=-c(OTU, total, cloud, abundance, spread))
pm <- as.matrix( p[, -1] )
rownames(pm) <- p$amplicon

otu.table = otu_table(pm, taxa_are_rows = TRUE)
res <- estimate_richness(otu.table)
library(data.table)
res <- setDT(res, keep.rownames = "sample")[]

out.path <- file.path( base.dir, "workdir", ds, "otu")
if ( !dir.exists(out.path) ) {dir.create(out.path, recursive=TRUE)}
write.csv(res, file.path( out.path, "result.csv" ), quote=FALSE, row.names = FALSE)
