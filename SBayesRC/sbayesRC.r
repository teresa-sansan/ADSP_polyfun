
# Run in R to install dependencies
##install.packages(c("Rcpp", "data.table", "stringi", "BH",  "RcppEigen"))

# Install SBayesRC package
install.packages("https://github.com/zhilizheng/SBayesRC/releases/download/v0.2.5/SBayesRC_0.2.5.tar.gz",
                 repos=NULL, type="source")

# If R report problem when installing, try alternative version (worse performance and an old version)
install.packages("https://github.com/zhilizheng/SBayesRC/releases/download/v0.2.0/SBayesRC_0.2.0_comp.tar.gz", repos=NULL, type="source")
