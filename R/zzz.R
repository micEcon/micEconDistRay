.onAttach <- function( libname, pkgname ) {
   packageStartupMessage(
      paste0( "\nPlease cite the 'micEconDistRay' package in publications as:\n",
         "Price, J.J. and Henningsen, A. (2022): ",
         "A Ray-Based Input Distance Function to Model Zero-Valued Output Quantities: ", 
         "Derivation and an Empirical Application. ",
         "IFRO Working Paper 2022/03. ",
         "Department of Food and Resource Economics, University of Copenhagen. ", 
         "url: https://ideas.repec.org/p/foi/wpaper/2022_03.html.\n\n",
         "If you have questions, suggestions, or comments ",
         "regarding the 'micEconDistRay' package, ",
         "please use 'Issues' at the package's GitHub site:\n",
         "https://github.com/micEcon/micEconDistRay/issues" ),
      domain = NULL,  appendLF = TRUE )
}
