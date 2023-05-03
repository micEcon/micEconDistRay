.onAttach <- function( libname, pkgname ) {
   packageStartupMessage(
      paste0( "\nPlease cite the 'micEconDistRay' package as:\n",
         "Price and Henningsen (forthcoming): ", 
         "A Ray-Based Input Distance Function to Model Zero-Valued Output Quantities: ", 
         "Derivation and an Empirical Application. ",
         "Journal of Productivity Analysis.\n\n",
         "If you have questions, suggestions, or comments ",
         "regarding the 'micEconDistRay' package, ",
         "please use 'Issues' at the package's GitHub site:\n",
         "https://github.com/micEcon/micEconDistRay/issues" ),
      domain = NULL,  appendLF = TRUE )
}
