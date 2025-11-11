#' Setting up the package dartR.popgenomics
#'
#' Setting up dartR.captive
#' @importFrom utils packageVersion read.csv read.delim read.table write.csv write.table
#' @importFrom methods getPackageName is new
#' @importFrom grDevices rainbow hcl
#' @importFrom graphics lines par
#' @importFrom stats ave pchisq var variable.names complete.cases pnorm quantile sd
#' @import adegenet
#' @import dartR.base
#' @import ggplot2
#' @importFrom crayon red yellow green blue cyan
#' @keywords internal


# needed to avoid error
zzz <- NULL

build <- "Jody"
error <- crayon::red
warn <- crayon::yellow
report <- crayon::green
important <- crayon::blue
code <- crayon::cyan

# WELCOME MESSAGE
.onAttach <- function(...) {
  pn <- getPackageName()
  packageStartupMessage(important(
    paste(
      "**** Welcome to", pn, "[Version",
      packageVersion(pn),
      "] ****\n"
    )
  ))
}

.onLoad <- function(libname, pkgname) {
  # Only set a default if user hasn’t set it already
  if (is.null(getOption("dartR_fbm"))) {
    val <- Sys.getenv("dartR_fbm", "")
    # Accept a few truthy values: 1, true, yes, on (case-insensitive)
    if (val=="TRUE") options(dartR_fbm = TRUE) else options(dartR_fbm=FALSE)
  }
}

# Fix error using . as placeholder 
if(getRversion() >= "2.15.1") utils::globalVariables(".")
