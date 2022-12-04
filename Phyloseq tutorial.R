install.packages("phyloseq")
source('http://bioconductor.org/biocLite.R')
biocLite('phyloseq')
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()
#double colon as phyloseq not in cram but in biomanager

yes

library(phyloseq)

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("phyloseq")

yes

library("phyloseq")

install.packages("vctrs")
yes

packageVersion("phyloseq")

library("ggplot2")
packageVersion("ggplot2")

#below you are adjusting the default theme for ggplot
theme_set(theme_bw())

options(url.method = "libcurl")

zipftp = "ftp://ftp.microbio.me/pub"
zipfile = tempfile("RestroomBiogeography")
download.file(zipftp,zipfile)




import_dir <- tempdir()
unzip(zipfile, exdir = import_dir)
