


#Generate help files
roxygen2::roxygenize(package.dir = ".", clean = T)
print(Rcpp::compileAttributes(pkgdir = ".", verbose = T))
#Check package
devtools::check(args="--as-cran")
#Build manual
devtools::build_manual(path=".")
#devtools::check_win_devel()
#build package
devtools::build()


