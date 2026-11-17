# gl.relatedness runs the closed-source program 'dartR_coancestry', which is
# downloaded separately. Tests that need it run only where the folder holding it
# is given in the environment variable DARTR_COANCESTRY_PATH (so never on CRAN).
coancestry_path <- function() Sys.getenv("DARTR_COANCESTRY_PATH", "")

skip_if_no_coancestry <- function() {
  p <- coancestry_path()
  exe <- if (.Platform$OS.type == "windows") "dartR_coancestry.exe" else "dartR_coancestry"
  testthat::skip_if(!nzchar(p) || !file.exists(file.path(p, exe)),
                    "the dartR_coancestry program is not available")
}
