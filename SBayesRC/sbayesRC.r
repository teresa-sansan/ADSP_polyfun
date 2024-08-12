
# Run in R to install dependencies
##install.packages(c("Rcpp", "data.table", "stringi", "BH",  "RcppEigen"))
install.packages(c("Rcpp", "data.table", "stringi", "BH",  "RcppEigen"))
# Install SBayesRC package
install.packages("https://github.com/zhilizheng/SBayesRC/releases/download/v0.2.6/SBayesRC_0.2.6.tar.gz",
                 repos=NULL, type="source")
# If R report problem when installing, try alternative version (worse performance and an old version)
install.packages("https://github.com/zhilizheng/SBayesRC/releases/download/v0.2.0/SBayesRC_0.2.0_comp.tar.gz", repos=NULL, type="source")

library(SBayesRC)
packageVersion(SBayesRC)


for (package_name in sort(loadedNamespaces())) {
    print(paste(package_name, packageVersion(package_name)))
}

