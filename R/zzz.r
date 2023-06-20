#' Setting up the package dartR.popgenomics
#'
#' Setting up dartR.captive
#' @importFrom utils packageVersion read.csv read.delim read.table write.csv write.table
#' @importFrom methods getPackageName is new
#' @importFrom grDevices rainbow hcl
#' @importFrom graphics lines par
#' @importFrom stats ave pchisq var variable.names complete.cases
#' @import adegenet
#' @import dartR.base
#' @import ggplot2
#' @import stringr
#' @import data.table
#' @importFrom dplyr "%>%" arrange bind_rows everything group_by mutate mutate_all mutate_at n n_distinct rename select starts_with ungroup vars

#' @keywords internal


#needed to avoid error
zzz<-NULL


build = "Jody"
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
      "**** Welcome to",pn,"[Version",
      packageVersion(pn),
      "] ****\n"
    )
  ))
}

