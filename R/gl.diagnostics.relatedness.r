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
#' @param x A genlight object containing SNP or SilicoDArT data [required].
#' @param cleanup Logical. Apply callrate, heterozygosity and all-NA filters
#'   before simulation [default = FALSE].
#' @param ref_variables Path to reference variable file [optional].
#' @param sim_variables Path to simulation variable file [optional].
#' @param which_tests Character vector of relatedness tests to apply
#'   [default = "wang"].
#' @param run_sim Logical. If TRUE, run simulations [default = FALSE].
#' @param IncludePlots Logical. If TRUE, generate and return plots
#'   [default = FALSE].
#' @param plotOut Logical. If TRUE, prints plots [default = FALSE]. 
#' @param varOut Logical. If TRUE, return variance results [default = FALSE].
#' @param rmseOut Logical. If TRUE, return RMSE results [default = FALSE].
#' @param numberIterations Integer. Number of simulation iterations
#'   [default = 1].
#' @param numberGenerations Integer. Number of generations to simulate
#'   [default = 3].
#' @param genToSave Either "all" or a numeric vector of generations to save
#'   [default = "all"].
#' @param runE9 Logical. If TRUE, include E9 analysis [default = FALSE].
#' @param E9Inbreed Logical. If TRUE, then runs EMIBD9 twice - once with inbreeding once w/out 
#'   [default = FALSE].
#' @param e9Path Path to external E9 binary [optional].
#' @param verbose Verbosity level: 0–5. If NULL, set by
#'   \code{gl.set.verbosity()} [default = NULL].
#' @param e9parallel Logical. Run E9 in parallel [default = FALSE].
#' @param nCores Integer. Number of cores if running E9 in parallel
#'   [default = 1].
#' @param includedPed Logical. If TRUE then input file has attache pedigree 
#'   [default = FALSE]
#'
#' @details
#' The function manages filtering, simulation setup, correlation
#' and relatedness outputs, and optional plotting. It handles quality
#' control checks on input objects and file paths before analysis.
#'
#' @return Returns an S4 object containing simulation and/or relatedness
#'   outputs. The slots for the output class are as follows: 
#'   \itemize{
#'     \item @InputDf: Original genlight input
#'     \item @SimOutput: Genlight object of simulation outputs
#'     \item @corOutList: Results of correlation analysis
#'     \item @corVals: Output of correlation results between methods
#'     \item @plotList: List of plots
#'   }
#'
#' @author Ethan, Luis (Post to
#'   \url{https://groups.google.com/d/forum/dartr})
#'
#' @examples
#' \dontrun{
#' gl.diagnostics.relatedness(possums.gl, run_sim = TRUE, IncludePlots = TRUE)
#' }
#'
#' @seealso \code{\link{gl.filter.callrate}},
#'   \code{\link{gl.filter.heterozygosity}}
#'
#' @export
#' @import methods
#' @import stats
#' @import dartR.sim
#' @import digest
#' @importFrom gridExtra tableGrob
#' @importFrom gridExtra ttheme_default
#' @importFrom magrittr %>%
#' @importFrom tidyr pivot_wider
#' @importFrom reshape2 acast
#' @importFrom related coancestry
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
    runE9 = FALSE,
    E9Inbreed = FALSE, 
    e9Path = NULL,
    verbose = NULL,
    e9parallel = FALSE,
    nCores = 1, 
    includedPed = FALSE
) {
  
  ID1 <- NA
  ID2 <- NA
  RelDegree <- NA
  child1 <- NA
  child2 <- NA
  dad <- NA
  id <- NA
  id1 <- NA
  id2 <- NA
  ind1 <- NA
  ind2 <- NA
  mom <- NA
  relationship <- NA
  value <- NA
  variable <- NA
  yintercept<- NA
  
  # SET VERBOSITY ----
  verbose <- gl.check.verbosity(verbose)
  
  # FLAG SCRIPT START ----
  funname <- match.call()[[1]]
  utils.flag.start(func = funname, build = "Jody", verbose = verbose)
  
  # CHECK DATATYPE ----
  datatype <- utils.check.datatype(x, verbose = verbose)
  
  # FUNCTION SPECIFIC ERROR CHECKING ----
  # Check required packages
  needed_pkgs <- c(
    "dartRverse", "related", "Rcpp",
    "ggplot2", "tidyverse", "reshape2", "data.table"
  )
  for (pkg in needed_pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("Package ", pkg,
           " needed for this function to work. Please install it."
      )
    }
  }
  
  # Validate input parameters
  if (!inherits(x, "genlight")) {
    stop("Input x must be a genlight object\n")
  }
  if (!is.logical(cleanup) || length(cleanup) != 1) {
    stop("cleanup must be TRUE or FALSE\n")
  }
  
  # Check file paths
  if (!is.null(ref_variables) && !is.null(sim_variables)) {
    if (!file.exists(ref_variables)) {
      stop("ref_variables file not found: ", ref_variables, "\n")
    }
    if (!file.exists(sim_variables)) {
      stop("sim_variables file not found: ", sim_variables, "\n")
    }
  }
  
  if (!is.character(which_tests) && !is.vector(which_tests)) {
    stop("which_tests must be a character vector\n")
  }
  if (!is.logical(run_sim) || length(run_sim) != 1) {
    stop("run_sim must be TRUE or FALSE\n")
  }
  if (!is.logical(IncludePlots) || length(IncludePlots) != 1) {
    stop("IncludePlots must be TRUE or FALSE\n")
  }
  if (!is.numeric(numberIterations) || numberIterations <= 0) {
    stop("numberIterations must be a positive number\n")
  }
  if (!is.numeric(numberGenerations) || numberGenerations <= 0) {
    stop("numberGenerations must be a positive number\n")
  }
  if (!(identical(genToSave, "all") || is.numeric(genToSave))) {
    stop("genToSave must be 'all' or a numeric vector\n")
  }
  if(includedPed && is.null(x$other$ind.metrics)){
    stop("includedPed set to True but input does not contain pedigree, set includedPed
         to False, or attach pedigree")
  }
  if(runE9 && is.null(e9Path)){
    stop("Cannot run EMIBD9 without necessary files, please provide path to file containing EMIBD9
         binaries.")
  }
  if(includedPed && run_sim){
    warning("run_sim and includedPed both True - attached pedigree will be overwritten, 
            pedigree from simulation to be used instead. To use attached pedigree,
            run separately without simulation (ie run_sim=F)")
  }
  if (!(run_sim | includedPed) & (varOut | rmseOut)) {
    stop("Cannot calculate variance or RMSE without pedigree, either from simulation or 
         appended to original dataset. Please set varOut and rmseOut to F or run_sim=T")
  }
  
  # Warn if no loci present
  if (nLoc(x) == 0) {
    warning(
      "Input genlight object has no loci - results may be meaningless\n"
    )
  }
  
  if (numberGenerations < 3) {
    stop("numberGenerations must be at least 3\n")
  }
  
  
  # do the job 
  
  if (!is.null(ref_variables) && !is.null(sim_variables)) {
    if (!file.exists(ref_variables) || !file.exists(sim_variables)) {
      stop("Input must have valid file path")
    }
  }
  
  if (!is.logical(cleanup)) stop("Cleanup value must be TRUE or FALSE")
  if (cleanup) {
    x <- x %>% 
      {. <- gl.filter.callrate(., threshold = 1, verbose = 0, mono.rm = F);.} %>%
      {. <- gl.filter.heterozygosity(.);.} %>%
      {. <- gl.filter.allna(.);.}
  }
  
  datatype <- utils.check.datatype(x)
  corOutList <- NULL
  pedOrSim <-FALSE
  if(run_sim ==T | includedPed==T){
    pedOrSim <- TRUE
  }
  
  finalClassValues <- NULL
  finalClassValues[["InputDf"]] <- x
  defaultAnalysisDf <- NULL
  defaultAnalysisDf[[1]] <- x
  
  pedigreeDfFinal <- NULL
  corOutList <- NULL
  corValStore <- NULL
  finalOutputPlots <- NULL
  finalSimOutput <- NULL
  slotNeedsClass <- NULL
  
  
  # 1. Run sim and store output 
  if(run_sim){
    defaultAnalysisDf <- NULL
    # --- Define class for storing Simulation output --- 
    sim_new <- new("DartSim",
                   input_data = x,
                   table_input = ref_variables,
                   sim_input = sim_variables, 
                   gen_number = numberGenerations, 
                   number_iterations = numberIterations)
    
    dartSim <- do_sim(sim_new)
    
    
    if(genToSave == "all" | genToSave == "All"){
      dartSim <- dartSim
    }else if(typeof(genToSave) == "double"){
      for(i in 1:length(dartSim)){
        dartSim[[i]] <- dartSim[[i]][genToSave]
      }
    }else{
      stop("Must be double or 'all'; Come on BRO!!")
    }
    
    # Combine simulation outputs
    finalSimOutput <- lapply(dartSim, function(sim) do.call(rbind, sim))
    finalClassValues[["SimOutput"]] <- finalSimOutput
    defaultAnalysisDf <- finalSimOutput
    
    # Extract pedigree from simulation 
    RelatedDataTable <- lapply(seq_along(dartSim), function(i) {
      ExtractParents(dartSim, iteration = numberIterations) %>%
        CleanupExtractParents() %>%
        as.matrix()
    })
    
    # Apply manual recoding, add relationship level, and fix column names
    RelatedManualRecode <- NULL
    for(i in 1:length(RelatedDataTable)){
      RelatedManualRecode[[i]] <- RelatedDataTable[[i]] %>%
        {colnames(.) <- c("id1", "id2", "RelDegree", "relationship");.} 
    }
    
  }
  
  # 2. Run analysis 
  analysisOutputDf <- lapply(defaultAnalysisDf, cleanup_rel, testSelect = which_tests)
  for(i in 1:length(analysisOutputDf)){
    analysisOutputDf[[i]] <- na.omit(analysisOutputDf[[i]])
  }
  which_tests <- c(which_tests, "rrBLUP")
  
  if (isTRUE(runE9)) {
    which_tests <- c(which_tests, "E9")
    for(i in 1:length(defaultAnalysisDf)){
      e9Run <- runE9(defaultAnalysisDf[[i]], 
                     e9Path, 
                     e9parallel=e9parallel, 
                     numCores=nCores) %>%
        {. <- mergeE9Related(., 
                             analysisOutputDf[[i]], 
                             test_select = which_tests);.}
      analysisOutputDf[[i]] <- e9Run
    }
    
    if(isTRUE(E9Inbreed)){
      which_tests <- c(which_tests, "E9_Inbred")
      for(i in 1:length(defaultAnalysisDf)){
        e9Run <- runE9(defaultAnalysisDf[[i]], 
                       e9Path, 
                       e9parallel=e9parallel, 
                       numCores=nCores, 
                       E9Inbreed = T) %>%
          {. <- mergeE9Related(., 
                               analysisOutputDf[[i]], 
                               test_select = which_tests);.}
        analysisOutputDf[[i]] <- e9Run
      }
    }
    
  }
  
  finalClassValues[["MergedDf"]] <- analysisOutputDf
  
  
  # 3. Pedigree calculation - either after sim/added pedigree
  if((!(includedPed) && run_sim) || (includedPed && run_sim)){
    pedigreeDfFinal <- mapply(mergeRelatedManual, 
                              relatedDf = analysisOutputDf, 
                              RecodeDf = RelatedManualRecode, 
                              SIMPLIFY = FALSE)
    finalClassValues[["MergedDf"]] <- pedigreeDfFinal
  }else if (includedPed && !(run_sim)){
    pedigreeDfFinal <- generateRelatedTableBaseInput(x, analysisOutputDf[[1]])
    finalClassValues[["MergedDf"]] <- pedigreeDfFinal
  }
  
  
  # 4. Construct correlation values 
  if(rmseOut || varOut){
    if(rmseOut == T){
      rmseDf <- calcRMSE(pedigreeDfFinal, which_tests) %>%
        tableOut()
      corOutList[["rmseDf"]] <- rmseDf
    }
    
    if(varOut==T){
      varDf <- calcVar(pedigreeDfFinal, which_tests) %>%
        tableOut()
      corOutList[["varDf"]] <- varDf
    }
    
    # Select which columns are doubles to then calculate r^2
    numTrue <- vapply(pedigreeDfFinal[[1]], function(col) {
      typeof(col[[1]]) == "double"
    }, logical(1))
    
    corVals <- lapply(pedigreeDfFinal, function(df) cor(df[, numTrue]))
    
    #slotNeedsClass[["corOutList"]] <- createTemplateClass
    #slotNeedsClass[["corVals"]] <- createCorOutput 
    finalClassValues[["corOutList"]] <- new("corOutList", 
                                            rmsePlot = corOutList[["rmseDf"]], 
                                            varPlot = corOutList[["varDf"]])
    finalClassValues[["corVals"]] <- new("corVals", 
                                         corVals = corVals)
    
  }
  
  
  # 5. Construct plots 
  if(IncludePlots){
    for(i in 1:length(finalClassValues[["MergedDf"]])){
      finalOutputPlots[[paste("Iteration", i,sep="")]] <- 
        relatedLevelPlots(finalClassValues[["MergedDf"]][[i]], 
                          which_tests=which_tests,
                          pedSim=pedOrSim)
    }
    
    # If run_sim True - includes function to create layered list of sim plots 
    # for different iterations
    if(run_sim){
      #slotNeedsClass[["plotList"]] <- createTemplateClass
      finalClassValues[["plotList"]] <- finalOutputPlots
    }
    
    finalClassValues[["plotList"]] <- finalOutputPlots
    
    templateClass <- finalOutputPlots
    
  }
  
  finalClassValues <- finalClassValues
  
  # 6. Construct final class to store everything 
  #finalOutput <- createOutputClass(finalClassValues, 
  #                                 slotNeedsClass, 
  #                                 slotDescriptions)
  
  templateClass <- new("finalOutputClass")
  templateClass@InputDf <- x
  templateClass@SimOutput <- finalClassValues[["SimOutput"]]
  templateClass@corVals <- finalClassValues[["corVals"]]
  templateClass@corOutList <- finalClassValues[["corOutList"]]
  templateClass@plotList <- finalClassValues[["plotList"]]
  
  if(plotOut){
    finalClassValues[["plotList"]]
  }
  
  return(templateClass)
  
  
}