

setwd(this.path::this.dir())
#Generate help files
roxygen2::roxygenize(package.dir = ".", clean = T)
print(Rcpp::compileAttributes(pkgdir = ".", verbose = T))

#Check package
unlink("./inst/doc",recursive=T)
devtools::check(args="--as-cran")

#Build manual
devtools::build_manual(path=".")
#devtools::check_win_devel()
#build package
unlink("./inst/doc",recursive=T)
devtools::build()

# build vignettes
tools::buildVignettes(dir = ".", tangle=TRUE)
dir.create("./inst/doc")
file.copy(dir("vignettes", full.names=TRUE), "inst/doc", overwrite=TRUE)

remotes::install_local("./../rapidsplithalf_0.3.tar.gz")

# browseVignettes("rapidsplithalf")
# vignette("rapidsplithalf")
