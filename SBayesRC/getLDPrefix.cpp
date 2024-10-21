#include <Rcpp.h>
using namespace Rcpp;

// Define the C++ function that you want to call from R
// [[Rcpp::export]]
List getLDPrefix(std::string mldm) {
    // Function implementation goes here
    List result;
    result["prefix"] = mldm + "_processed"; // Example processing
    return result;
}
