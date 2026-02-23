#' Setting up dartR.captive
#' @importFrom utils packageVersion read.csv read.delim read.table write.csv write.table
#' @importFrom methods getPackageName is new
#' @importFrom grDevices rainbow hcl
#' @importFrom graphics lines par
#' @importFrom stats ave pchisq var variable.names complete.cases pnorm quantile sd
#' @rawNamespace import(adegenet, except = c(glMean, glSum))
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

# Fix error using . as placeholder 
if(getRversion() >= "2.15.1") utils::globalVariables(".")

## returns NULL if the 'fbm' slot is missing OR is NULL
.fbm_or_null <- function(x) {
  if (methods::.hasSlot(x, "fbm")) {
    val <- methods::slot(x, "fbm")
    return(if (is.null(val)) NULL else val)
  }
  NULL
}

.onLoad <- function(libname, pkgname) {
  # Only set a default if user hasn’t set it already
  if (is.null(getOption("dartR_fbm"))) {
    val <- Sys.getenv("dartR_fbm", "")
    # Accept a few truthy values: 1, true, yes, on (case-insensitive)
    if (val=="TRUE") options(dartR_fbm = TRUE) else options(dartR_fbm=FALSE)
  }
}

