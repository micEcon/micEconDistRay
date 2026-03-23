.onAttach <- function( libname, pkgname ) {
   packageStartupMessage(
      paste0( "\nPlease cite the 'micEconDistRay' package in publications as:\n",
         "Price, J.J. and Henningsen, A. (2023): ",
         "A Ray-Based Input Distance Function to Model Zero-Valued Output Quantities: ", 
         "Derivation and an Empirical Application. ",
         "Journal of Productivity Analysis 60, p. 179-188. ",
         "DOI: 10.1007/s11123-023-00684-1. ",
         "URL: https://doi.org/10.1007/s11123-023-00684-1.\n\n",
         "If you have questions, suggestions, or comments ",
         "regarding the 'micEconDistRay' package, ",
         "please use 'Issues' at the package's GitHub site:\n",
         "https://github.com/micEcon/micEconDistRay/issues" ),
      domain = NULL,  appendLF = TRUE )
}
