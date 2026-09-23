#' @name gl.diagnostics.relatedness
#'
#' @title Run simulations and relatedness analyses on genlight objects
#'
#' @description
#' This function wraps a variety of methods for estimating relatedness, such
#' that they can be directly compared for accuracy and precision. It also
#' provides the ability to run the gl.sim function for a minimum of 3 generations,
#' providing further functionality with regards to estimating gene flow and population
#' dynamics. It supports multiple simulation back ends, correlation
#' output, error checking, RMSE/variance summaries, and optional plotting.
#'
#' @param x A genlight object containing SNP data [required].
#' @param cleanup Logical. Apply callrate, heterozygosity and all-NA filters
#'   before simulation [default = FALSE].
#' @param ref_variables Path to the reference variable file for
#'   \code{dartR.sim::gl.sim.WF.table} [default NULL, which uses
#'   ref_variables.csv shipped with dartR.sim].
#' @param sim_variables Path to the simulation variable file for
#'   \code{dartR.sim::gl.sim.WF.run} [default NULL, which uses
#'   sim_variables.csv shipped with dartR.sim].
#' @param which_tests Character vector of relatedness tests to apply
#'   [default = "wang"].
#' @param run_sim Logical. If TRUE, run simulations [default = FALSE].
#' @param IncludePlots Logical. If TRUE, generate and return plots
#'   [default = FALSE].
#' @param plotOut Logical. If TRUE, prints the generated plots to the graphics
#'   device (requires \code{IncludePlots = TRUE}) [default = FALSE].
#' @param varOut Logical. If TRUE, return variance results [default = FALSE].
#' @param rmseOut Logical. If TRUE, return the root mean square error of each
#'   estimator against the pedigree kinship, by relationship class
#'   [default = FALSE].
#' @param numberIterations Integer. Number of simulation iterations
#'   [default = 1].
#' @param numberGenerations Integer. Number of generations to simulate;
#'   minimum 3 [default = 3].
#' @param genToSave Either "all" or a numeric vector of generations to save
#'   [default = "all"].
#' @param run.e9 Logical. If TRUE, include EMIBD9 analysis [default = FALSE].
#' @param E9Inbreed Logical. If TRUE, then runs EMIBD9 twice - once with inbreeding once w/out
#'   [default = FALSE].
#' @param e9Path Path to external EMIBD9 binary [optional].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#'   progress log; 3, progress and results summary; 5, full report
#'   [default 2, unless specified using gl.set.verbosity].
#' @param e9parallel Logical. Run EMIBD9 in parallel [default = FALSE].
#' @param nCores Integer. Number of cores if running EMIBD9 in parallel
#'   [default = 1].
#' @param includedPed Logical. If TRUE, the input has a pedigree attached in
#'   \code{x@other$ind.metrics} (columns id, dad, mom; missing parents 0 or
#'   NA) [default = FALSE]
#'
#' @details
#' The function manages filtering, simulation setup, correlation
#' and relatedness outputs, and optional plotting. It handles quality
#' control checks on input objects and file paths before analysis.
#'
#' The relatedness estimators in \code{which_tests} are computed by
#' \code{coancestry()} from the package \code{related}, which is not on CRAN;
#' install it with
#' \code{devtools::install_github("timothyfrasier/related")}. Its Fortran code
#' prints progress to the console at every verbosity.
#'
#' All estimates are on the kinship scale. When a pedigree is available (from
#' the simulation or attached with includedPed), every pair of individuals
#' gets its exact pedigree kinship (column rel), computed recursively from the
#' pedigree with founders unrelated and not inbred, and its closest
#' relationship class (column RelDegree: parent_offspring, full_sibs,
#' half_sibs, full_first_cousins, half_first_cousins, second_cousins or
#' unrelated). All pairs are kept, whatever their estimated kinship. RMSE is
#' the root mean square difference between each estimator and rel, by
#' relationship class.
#'
#' @return Returns an S4 object containing simulation and/or relatedness
#'   outputs. The slots for the output class are as follows:
#'   \itemize{
#'     \item @InputDf: The genlight input (after filtering when
#'       \code{cleanup = TRUE})
#'     \item @SimOutput: Genlight object of simulation outputs
#'     \item @MergedDf: Kinship estimates per iteration, one row per pair;
#'       with a pedigree, also the columns RelDegree and rel
#'     \item @corOutList: Results of correlation analysis
#'     \item @corVals: Output of correlation results between methods
#'     \item @plotList: List of plots
#'   }
#'
#' @author Author(s): Ethan, Luis Mijangos. Custodian: Luis Mijangos -- Post
#'   to \url{https://groups.google.com/d/forum/dartr}
#'
#' @examples
#' \dontrun{
#' # requires the package 'related'
#' res <- gl.diagnostics.relatedness(possums.gl, run_sim = TRUE,
#'                                   rmseOut = TRUE, IncludePlots = TRUE)
#' }
#'
#' @seealso \code{\link[dartR.base]{gl.filter.callrate}},
#'   \code{\link[dartR.base]{gl.filter.heterozygosity}}
#'
#' @family captive management
#' @export
#' @import methods
#' @import stats
#' @import dartR.sim
#' @importFrom gridExtra tableGrob
#' @importFrom gridExtra ttheme_default
#' @importFrom magrittr %>%
#' @importFrom tidyr pivot_wider
#' @importFrom reshape2 acast
gl.diagnostics.relatedness <- function(
    x,
    cleanup = FALSE,
    ref_variables = NULL,
    sim_variables = NULL,
    which_tests = "wang",
    run_sim = FALSE,
    IncludePlots = FALSE,
    plotOut = FALSE,
    varOut = FALSE,
    rmseOut = FALSE,
    numberIterations = 1,
    numberGenerations = 3,
    genToSave = "all",
    run.e9 = FALSE,
    E9Inbreed = FALSE,
    e9Path = NULL,
    verbose = NULL,
    e9parallel = FALSE,
    nCores = 1,
    includedPed = FALSE
) {

  # SET VERBOSITY ----
  verbose <- gl.check.verbosity(verbose)

  # FLAG SCRIPT START ----
  funname <- match.call()[[1]]
  utils.flag.start(func = funname, verbose = verbose)

  # CHECK DATATYPE ----
  datatype <- utils.check.datatype(x, verbose = verbose)

  # FUNCTION SPECIFIC ERROR CHECKING ----
  # coancestry() from 'related' computes the which_tests estimators; everything
  # else this function needs is a declared dependency. 'related' is not on CRAN
  # so it cannot be declared, hence the variable indirection to keep R CMD
  # check's dependency scan quiet.
  pkg <- "related"
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(error(
      "  gl.diagnostics.relatedness needs the package 'related' (not on CRAN).\n",
      "  Install it with:\n",
      "    devtools::install_github('timothyfrasier/related')\n"))
  }

  # related and gl.grm need SNP genotypes
  if (datatype == "SilicoDArT") {
    stop(error(
      "  Only SNP data are supported; x contains SilicoDArT data\n"
    ))
  }

  # the simulation needs both variable files; default to those shipped with
  # dartR.sim
  if (is.null(ref_variables)) {
    ref_variables <- system.file("extdata", "ref_variables.csv",
                                 package = "dartR.sim")
  }
  if (is.null(sim_variables)) {
    sim_variables <- system.file("extdata", "sim_variables.csv",
                                 package = "dartR.sim")
  }

  # Validate input parameters
  if (!is.logical(cleanup) || length(cleanup) != 1) {
    stop(error("cleanup must be TRUE or FALSE\n"))
  }

  # Check file paths
  if (!is.null(ref_variables) && !is.null(sim_variables)) {
    if (!file.exists(ref_variables)) {
      stop(error("ref_variables file not found: ", ref_variables, "\n"))
    }
    if (!file.exists(sim_variables)) {
      stop(error("sim_variables file not found: ", sim_variables, "\n"))
    }
  }

  if (!is.character(which_tests)) {
    stop(error("which_tests must be a character vector\n"))
  }
  if (!is.logical(run_sim) || length(run_sim) != 1) {
    stop(error("run_sim must be TRUE or FALSE\n"))
  }
  if (!is.logical(IncludePlots) || length(IncludePlots) != 1) {
    stop(error("IncludePlots must be TRUE or FALSE\n"))
  }
  if (!is.numeric(numberIterations) || numberIterations <= 0) {
    stop(error("numberIterations must be a positive number\n"))
  }
  if (!is.numeric(numberGenerations) || numberGenerations < 3) {
    stop(error("numberGenerations must be at least 3\n"))
  }
  if (!(identical(genToSave, "all") || is.numeric(genToSave))) {
    stop(error("genToSave must be 'all' or a numeric vector\n"))
  }
  if (includedPed && is.null(x@other$ind.metrics)) {
    stop(error(
      "includedPed set to TRUE but input does not contain a pedigree; set",
      " includedPed to FALSE, or attach a pedigree\n"))
  }
  if (run.e9 && is.null(e9Path)) {
    stop(error(
      "Cannot run EMIBD9 without the binary; please provide e9Path, the path",
      " to the folder containing the EMIBD9 binaries\n"))
  }
  if (includedPed && run_sim && verbose >= 1) {
    cat(warn(
      "  run_sim and includedPed both TRUE - the attached pedigree will be",
      " ignored and the pedigree from the simulation used instead. To use the",
      " attached pedigree, run separately with run_sim = FALSE\n"))
  }
  if (!(run_sim || includedPed) && (varOut || rmseOut)) {
    stop(error(
      "Cannot calculate variance or RMSE without a pedigree, either from",
      " simulation or attached to the input. Set varOut and rmseOut to FALSE",
      " or run_sim to TRUE\n"))
  }
  if (nLoc(x) == 0 && verbose >= 1) {
    cat(warn(
      "  Input genlight object has no loci - results may be meaningless\n"))
  }

  # DO THE JOB ----

  if (cleanup) {
    x <- gl.filter.callrate(x, threshold = 1, verbose = 0, mono.rm = FALSE)
    x <- gl.filter.heterozygosity(x, verbose = 0)
    x <- gl.filter.allna(x, verbose = 0)
  }

  pedOrSim <- run_sim || includedPed

  finalClassValues <- list(InputDf = x)

  # GUARD AGAINST related's FORTRAN ABORTS ----
  # coancestry()'s Fortran reads a space-delimited genotype file and calls STOP
  # on malformed input; on Windows that kills the whole R session rather than
  # raising an R error. Two known triggers are guarded here: whitespace in
  # individual names (breaks the file parsing) and loci scored NA across all
  # retained individuals (a common state after subsetting a filtered object by
  # population).
  orig.names <- indNames(x)
  safe.names <- gsub("[[:space:]]+", "_", orig.names)
  names.changed <- !identical(safe.names, orig.names)
  if (names.changed) {
    if (anyDuplicated(safe.names)) {
      stop(error(
        "Individual names differ only by whitespace and collide once spaces",
        "are replaced; please make individual names unique\n"))
    }
    indNames(x) <- safe.names
    if (includedPed) {
      for (cn in intersect(c("id", "dad", "mom"),
                           colnames(x@other$ind.metrics))) {
        x@other$ind.metrics[[cn]] <-
          gsub("[[:space:]]+", "_", as.character(x@other$ind.metrics[[cn]]))
      }
    }
    if (verbose >= 1) {
      cat(warn(
        "  Individual names contain spaces; sanitised names are used",
        "internally and the original names restored in the output\n"))
    }
  }
  name.map <- stats::setNames(orig.names, safe.names)

  am <- as.matrix(x)
  keep.loc <- colSums(!is.na(am)) > 0
  keep.ind <- rowSums(!is.na(am)) > 0
  if (!all(keep.loc)) x <- x[, keep.loc]
  if (!all(keep.ind)) x <- x[keep.ind, ]
  if ((!all(keep.loc) || !all(keep.ind)) && verbose >= 1) {
    cat(warn(sprintf(
      "  Dropped %d all-NA loci and %d all-NA individuals before the relatedness analysis\n",
      sum(!keep.loc), sum(!keep.ind))))
  }

  defaultAnalysisDf <- list(x)

  pedigreeDfFinal <- NULL
  corOutList <- NULL
  finalOutputPlots <- NULL

  # 1. Run sim and store output
  if (run_sim) {
    sim_new <- new("DartSim",
                   input_data = x,
                   table_input = ref_variables,
                   sim_input = sim_variables,
                   gen_number = numberGenerations,
                   number_iterations = numberIterations)

    dartSim <- do_sim(sim_new)

    if (is.numeric(genToSave)) {
      for (i in seq_along(dartSim)) {
        dartSim[[i]] <- dartSim[[i]][genToSave]
      }
    }

    # Combine simulation outputs
    finalSimOutput <- lapply(dartSim, function(sim) do.call(rbind, sim))
    finalClassValues[["SimOutput"]] <- finalSimOutput
    defaultAnalysisDf <- finalSimOutput

    # Pedigree (id, dad, mom) of each iteration
    simPedigrees <- lapply(seq_along(dartSim), function(i) {
      ExtractParents(dartSim, iteration = i)
    })

  }

  # 2. Run analysis
  analysisOutputDf <- lapply(defaultAnalysisDf, cleanup_rel,
                             testSelect = which_tests)
  analysisOutputDf <- lapply(analysisOutputDf, na.omit)
  which_tests <- c(which_tests, "rrBLUP")

  if (isTRUE(run.e9)) {
    which_tests <- c(which_tests, "E9")
    for (i in seq_along(defaultAnalysisDf)) {
      analysisOutputDf[[i]] <- runE9(defaultAnalysisDf[[i]],
                                     e9Path,
                                     e9parallel = e9parallel,
                                     numCores = nCores) %>%
        mergeE9Related(analysisOutputDf[[i]], test_select = which_tests)
    }

    if (isTRUE(E9Inbreed)) {
      which_tests <- c(which_tests, "E9_Inbred")
      for (i in seq_along(defaultAnalysisDf)) {
        analysisOutputDf[[i]] <- runE9(defaultAnalysisDf[[i]],
                                       e9Path,
                                       e9parallel = e9parallel,
                                       numCores = nCores,
                                       E9Inbreed = TRUE) %>%
          mergeE9Related(analysisOutputDf[[i]], test_select = which_tests)
      }
    }

  }

  finalClassValues[["MergedDf"]] <- analysisOutputDf

  # 3. Pedigree calculation - either after sim/added pedigree
  # When both run_sim and includedPed are TRUE, the simulated pedigree wins
  # (warned above).
  if (run_sim) {
    pedigreeDfFinal <- mapply(mergePedigreeTruth,
                              relatedDf = analysisOutputDf,
                              ped = simPedigrees,
                              SIMPLIFY = FALSE)
    finalClassValues[["MergedDf"]] <- pedigreeDfFinal
  } else if (includedPed) {
    pedigreeDfFinal <- list(mergePedigreeTruth(analysisOutputDf[[1]],
                                               attachedPedigree(x)))
    finalClassValues[["MergedDf"]] <- pedigreeDfFinal
  }

  # Restore original individual names in the output tables
  if (names.changed) {
    finalClassValues[["MergedDf"]] <- lapply(
      finalClassValues[["MergedDf"]], function(df) {
        for (cn in intersect(c("ind1", "ind2", "id1", "id2"), colnames(df))) {
          v <- as.character(df[[cn]])
          hit <- v %in% names(name.map)
          v[hit] <- name.map[v[hit]]
          df[[cn]] <- v
        }
        df
      })
    if (!is.null(pedigreeDfFinal)) {
      pedigreeDfFinal <- finalClassValues[["MergedDf"]]
    }
  }

  # 4. Construct correlation values
  if (rmseOut || varOut) {
    if (rmseOut) {
      corOutList[["rmseDf"]] <- calcRMSE(pedigreeDfFinal, which_tests) %>%
        tableOut()
    }

    if (varOut) {
      corOutList[["varDf"]] <- calcVar(pedigreeDfFinal, which_tests) %>%
        tableOut()
    }

    # Select which columns are doubles to then calculate r^2
    numTrue <- vapply(pedigreeDfFinal[[1]], function(col) {
      typeof(col[[1]]) == "double"
    }, logical(1))

    corVals <- lapply(pedigreeDfFinal, function(df) cor(df[, numTrue]))

    finalClassValues[["corOutList"]] <- new("corOutList",
                                            rmsePlot = corOutList[["rmseDf"]],
                                            varPlot = corOutList[["varDf"]])
    finalClassValues[["corVals"]] <- new("corVals",
                                         corVals = corVals)

  }

  # 5. Construct plots
  if (IncludePlots) {
    for (i in seq_along(finalClassValues[["MergedDf"]])) {
      finalOutputPlots[[paste0("Iteration", i)]] <-
        relatedLevelPlots(finalClassValues[["MergedDf"]][[i]],
                          which_tests = which_tests,
                          pedSim = pedOrSim)
    }

    finalClassValues[["plotList"]] <- finalOutputPlots

    if (plotOut) {
      # relatedLevelPlots returns a single ggplot, or a list of them when a
      # pedigree is available
      for (p in finalOutputPlots) {
        if (inherits(p, "ggplot")) print(p) else lapply(p, print)
      }
    }
  }

  # 6. Construct final class to store everything
  templateClass <- new("finalOutputClass")
  templateClass@InputDf <- finalClassValues[["InputDf"]]
  templateClass@SimOutput <- finalClassValues[["SimOutput"]]
  templateClass@MergedDf <- finalClassValues[["MergedDf"]]
  templateClass@corVals <- finalClassValues[["corVals"]]
  templateClass@corOutList <- finalClassValues[["corOutList"]]
  templateClass@plotList <- finalClassValues[["plotList"]]

  return(templateClass)

}
