#' @name gl.run.EMIBD9
#' @title Run program EMIBD9
#' @description
#' Runs the program EMIBD9 (Wang 2022) on a genlight object with SNP data and
#' returns pairwise kinship, the IBD mode probabilities (delta1 to delta9) and
#' individual inbreeding.
#' @param x Name of the genlight object containing the SNP data [required].
#' @param outfile File name (without a path) of the EMIBD9 output file; it is
#' copied to outpath [default "EMIBD9_Res.ibd9"].
#' @param outpath Path where to save the output file. Use outpath=getwd() or
#' outpath='.' when calling this function to direct output files to your working
#' or current directory [default tempdir(), mandated by CRAN].
#' @param emibd9.path Path to the folder with the EMIBD9 files.
#'  Please note there are 2 different executables depending on your OS:
#'  EM_IBD_P.exe (=Windows) EM_IBD_P (=Mac, Linux).
#'  You only need to point to the folder (the function will recognise which OS 
#'  you are running) [default getwd()].
#' @param OutAlleleFre Whether to output allele frequencies (TRUE/FALSE or
#'  1/0) [default FALSE].
#' @param EM_Method An integer that indicates the method to use for the expectation
#'  maximization (EM) algorithm. 1, the standard EM method;
#'  2, the EM method with a quasi-Newton acceleration; 3, the EM method with a
#'  SQUAREM acceleration [default 1].
#' @param Inbreed A boolean that indicates whether to compute inbreeding (i.e. delta1 to delta6) [default FALSE].
#' @param palette_convergent A character vector of colours to use for the heatmap plot.
#'  If NULL, the default palette from gl.colors("div") will be used [default NULL].
#' @param parallel A boolean that indicates whether to run the parallel version of EM
#' IBD9 (EM_IBD_P_mpi) [default FALSE].
#' @param ncores An integer specifying the number of cores to use when parallel is TRUE
#' [default 1].
#' @param ISeed An integer specifying the random seed to use for the EM algorithm
#' [default 42].
#' @param plot.out A boolean that indicates whether to plot the results
#'  [default TRUE].
#' @param plot.dir Directory to save the plot RDS files [default as specified
#' by the global working directory or tempdir()]
#' @param plot.file Name for the RDS binary file to save (base name only, 
#' exclude extension) [default NULL]
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log, including the EMIBD9 console output; 3, progress and results
#' summary; 5, full report [default 2, unless specified using
#' gl.set.verbosity].
#' @details
#' The results of EMIBD9 include the identical in state (IIS) values for each mode 
#'(S1 - 9) and nine condensed identical by descent (IBD) modes (delta1 - delta9) as well as 
#' the relatedness coefficient (r). Alleles are IIS if they are the same. Similarly,
#' IBD describes a matching allele between two individuals that has been inherited 
#' from a common ancestor or common gene. In a pairwise comparison, delta1 to delta9 are the
#'  probabilities associated with each IBD mode. delta1 to delta6 take value > 0 in presence
#'  of inbreeding and hence are only computed when this option is selected. 
#'  
#'EMIBD9 uses an expectation maximization (EM) algorithm based on the maximum
#' likelihood expectations (MLE) of \eqn{\delta} to estimate both allele frequencies (p) 
#' and \eqn{\delta} jointly from genotype data. By iteratively calculating p and \eqn{\delta}, 
#' relatedness can be modified to reduce biases due to small sample sizes. 
#' Wang J. (2022) suggest the resulting r coefficient is therefore more robust 
#' compared to previous methods.
#'
#'The kinship coefficient is the probability that two alleles at a random locus
#'  drawn from two individuals are IBD.
#'
#'Below is a table modified from Speed & Balding (2015) showing kinship values,
#'and their confidence intervals (CI), for different relationships.
#'
#' \tabular{lll}{
#'   \strong{Relationship} \tab \strong{Kinship} \tab \strong{95\% CI} \cr
#'   Identical twins / clones / same individual \tab 0.5   \tab -              \cr
#'   Sibling / Parent-Offspring                \tab 0.25  \tab (0.204, 0.296)\cr
#'   Half-sibling                              \tab 0.125 \tab (0.092, 0.158)\cr
#'   First cousin                              \tab 0.062 \tab (0.038, 0.089)\cr
#'   Half-cousin                               \tab 0.031 \tab (0.012, 0.055)\cr
#'   Second cousin                             \tab 0.016 \tab (0.004, 0.031)\cr
#'   Half-second cousin                        \tab 0.008 \tab (0.001, 0.020)\cr
#'   Third cousin                              \tab 0.004 \tab (0.000, 0.012)\cr
#'   Unrelated                                 \tab 0     \tab -              \cr
#' }
#'
#'For greater detail on the methods employed by EMIBD9, we encourage you to 
#'read Wang, J. (2022).
#'
#' Download the program from here:
#'
#' https://www.zsl.org/about-zsl/resources/software/emibd9
#'
#' For Windows, Mac and Linux install the program then point to the folder
#'  where you find: EM_IBD_P.exe (=Windows) and EM_IBD_P (=Mac, Linux). If 
#'  running really slow you may want to create the files using the function 
#'  and then run in parallel using the documentation provided by the authors
#'   [you need to have mpiexec installed].
#'   
#'  Please note individual names must have a maximal length of 20 characters. 
#'  The IDs must NOT contain blank space and other illegal characters 
#'  (such as /), and must be unique among all sampled individuals (i.e. NO 
#'  duplications). Any string longer than 20 characters for individual ID will 
#'  be truncated to have 20 characters.
#'
#' Each call runs EMIBD9 in its own temporary folder. EMIBD9 does not return
#' an error status when it fails, so the function stops with the last lines
#' of the EMIBD9 console output when no output file is written.
#'
#' @return A list with four elements:
#' \itemize{
#' \item rel -- a square matrix of pairwise kinship coefficients (theta; the
#' EMIBD9 column r(1,2)), with self-comparisons on the diagonal
#' (0.5 x (1 + F)), tagged attr(rel, "scale") = "kinship".
#' \item raw -- the raw EMIBD9 table, all pairs including self-comparisons,
#' with numeric columns.
#' \item processed -- the table without self-comparisons and redundant pairs
#' (e.g. the second individual with the first, when the first with the
#' second is already present).
#' \item inbreeding -- a table of individual inbreeding values; EMIBD9 writes
#' it whether or not Inbreed is TRUE.
#' }
#'
#' @author Author(s): Luis Mijangos. Custodian: Luis Mijangos -- Post to
#' \url{https://groups.google.com/d/forum/dartr}
#' @examples
#' \dontrun{
#' #To run this function needs EMIBD9 installed in your computer
#' if (isTRUE(getOption("dartR_fbm"))) platypus.gl <- gl.gen2fbm(platypus.gl)
#' t1 <- gl.filter.allna(platypus.gl)
#' res_rel <- gl.run.EMIBD9(t1)
#' }
#'
#' @references
#' \itemize{
#' \item Wang, J. (2022). A joint likelihood estimator of relatedness and allele
#'  frequencies from a small sample of individuals. Methods in Ecology and
#'  Evolution, 13(11), 2443-2462.
#' }
#' 
#' @family captive management
#' @importFrom utils combn
#' @importFrom stringr str_split
#' @rawNamespace import(data.table)
#' @export

gl.run.EMIBD9 <- function(x,
                          outfile = "EMIBD9_Res.ibd9",
                          outpath = tempdir(),
                          emibd9.path = getwd(),
                          OutAlleleFre = 0,
                          EM_Method = 1,
                          Inbreed = FALSE,
                          palette_convergent = NULL,
                          parallel = FALSE,
                          ncores = 1,
                          ISeed = 42,
                          plot.out = TRUE,
                          plot.dir = NULL,
                          plot.file = NULL,
                          verbose = NULL) {
  # SET VERBOSITy
  verbose <- gl.check.verbosity(verbose)
  
  # SET WORKING DIRECTORY
  plot.dir <- gl.check.wd(plot.dir, verbose = 0)
  
  # FLAG SCRIPT START
  funname <- match.call()[[1]]
  utils.flag.start(func = funname, verbose = verbose)
  
  # CHECK DATATYPE
  datatype <- utils.check.datatype(x, verbose = verbose)
  
  # FUNCTION SPECIFIC ERROR CHECKING
  # EMIBD9 reads 0/1/2 genotypes; presence/absence would be read as SNPs
  if (datatype == "SilicoDArT") {
    stop(error(
      "  Only SNP data are supported; x contains SilicoDArT data\n"
    ))
  }
  
  # each call runs in its own folder, so output left by an earlier call can
  # never be read as the result of this one
  rundir <- tempfile("EMIBD9_")
  dir.create(rundir)
  
  if (Sys.info()["sysname"] == "Windows") {
    prog <- c("EM_IBD_P.exe", "impi.dll", "libiomp5md.dll")
    cmd <- "EM_IBD_P.exe"
    cmd.args <- "INP:MyData.par"
  } 
  
  if (Sys.info()["sysname"] %in% c("Linux", "Darwin")) {
    if(parallel){
      prog <- "EM_IBD_P_mpi"
      cmd <- "mpirun"
      cmd.args <- c("-np", ncores, "--use-hwthread-cpus", "./EM_IBD_P_mpi",
                    "INP:MyData.par")
    }else{
      prog <- "EM_IBD_P"
      cmd <- "./EM_IBD_P"
      cmd.args <- "INP:MyData.par"
    }
  }
  
  # check if file program can be found
  if (all(file.exists(file.path(emibd9.path, prog)))) {
    file.copy(file.path(emibd9.path, prog),
              to = rundir,
              overwrite = TRUE)
    if (verbose >= 2) {
      cat(report("  Found necessary files to run EMIBD9.\n"))
    }
    
  } else {
    stop(error(
      "  Cannot find", paste(prog, collapse = ", "),
      "in the folder given by emibd9.path:", emibd9.path, "\n"
    ))
  }
  
  # Resolve no visible global function definition 
  J <- NULL
  
  # individual IDs must have a maximal length of 20 characters. The IDs must NOT
  # contain blank space and other illegal characters (such as /), and must be
  # unique among all sampled individuals (i.e. NO duplications). Any string longer
  # than 20 characters for individual ID will be truncated to have 20 characters.

  x2 <- x  #copy to work only on the copied data set

  # EMIBD9 requires individual IDs to be unique, to contain no spaces or
  # other illegal characters, and truncates IDs longer than 20 characters
  # (which can create duplicates). Sanitise a copy of the names for the
  # EMIBD9 input files and map results back to the original names afterwards.
  hold_names <- indNames(x)
  safe_names <- gsub("[^A-Za-z0-9_.-]", "_", hold_names)
  safe_names[is.na(safe_names) | safe_names == ""] <- "ind"
  safe_names <- make.unique(strtrim(safe_names, 15), sep = "_")
  indNames(x2) <- safe_names
  if (!identical(safe_names, hold_names) && verbose > 0) {
    cat(warn(
      "  Individual names were adjusted to meet EMIBD9 requirements (unique,",
      "no spaces or special characters, maximum 20 characters). Original",
      "names are restored in the results.\n"
    ))
  }

  NumIndiv <- nInd(x2)
  NumLoci <- nLoc(x2)
  DataForm <- 2
  if (Inbreed) {
    Inbreed <- 1
  } else{
    Inbreed <- 0
  }
  
  GtypeFile <- "EMIBD9_Gen.dat"
  OutFileName <- outfile
  RndDelta0 <- 1
  EM_Method <- EM_Method
  # EMIBD9 reads 0/1; a logical would be written as TRUE/FALSE and crash it
  OutAlleleFre <- as.integer(isTRUE(as.logical(OutAlleleFre)))

  param <- paste(NumIndiv,
    NumLoci,
    DataForm,
    Inbreed,
    GtypeFile,
    OutFileName,
    ISeed,
    RndDelta0,
    EM_Method,
    OutAlleleFre,
    sep = "\n"
  )
  
  IndivID <- paste(indNames(x2))
  
  gl_mat <- as.matrix(x2)
  gl_mat[is.na(gl_mat)] <- 3
  
  tmp <- cbind(apply(gl_mat, 1, paste, collapse = ""))
  
  tmp <- rbind(paste(indNames(x2), collapse = " "), tmp)

  # run EMIBD9
  # change into tempdir (run it there)
  old.path <- getwd()
  on.exit(setwd(old.path))
  setwd(rundir)
  write.table(tmp,
              file = GtypeFile,
              quote = FALSE,
              row.names = FALSE,
              col.names = FALSE
  )
  write.table(
    param,
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE,
    file = "MyData.par")
  # EMIBD9 exits with status 0 when it fails, so the output file is the test
  console.file <- file.path(rundir, "EMIBD9_console.txt")
  console.to <- if (verbose >= 2) "" else console.file
  status <- system2(cmd, args = cmd.args, stdout = console.to,
                    stderr = console.to)
  if (!file.exists(outfile) ||
      !any(grepl("^IBD", readLines(outfile, warn = FALSE)))) {
    console.tail <- if (file.exists(console.file)) {
      utils::tail(readLines(console.file, warn = FALSE), 10)
    } else {
      "(EMIBD9 console output shown above)"
    }
    stop(error(paste0(
      "  EMIBD9 did not write its results (exit status ", status, "):\n",
      paste(console.tail, collapse = "\n"), "\n"
    )))
  }
  
  ### get output  

  x_lines <- readLines(outfile)
  strt <- which(grepl("^IBD", x_lines)) + 2
  stp <- which(grepl("Indiv genotypes", x_lines)) - 4
  linez_headings <- x_lines[strt]
  linez_data <- x_lines[(strt + 1):stp]
  tmp_headings <- unlist(stringr::str_split(linez_headings, " "))
  tmp_data <- stringr::str_split(linez_data, " ")
  
  # Raw data
  tmp_data_raw_1 <- lapply(tmp_data, "[", c(2:22))
  tmp_data_raw_2 <- do.call("rbind", tmp_data_raw_1)
  tmp_data_raw_3 <- as.data.frame(tmp_data_raw_2)
  colnames(tmp_data_raw_3) <- tmp_headings[2:22]
  # every column except the two individual IDs is numeric
  num.cols <- setdiff(colnames(tmp_data_raw_3), c("Indiv1", "Indiv2"))
  tmp_data_raw_3[num.cols] <- lapply(tmp_data_raw_3[num.cols], as.numeric)
  
  # Kick out self & redundant comparisons
  # the parsed Indiv1/Indiv2 columns hold the sanitised names, so the pairs
  # to keep are built from safe_names (unique by construction)
  unq_pairs <- data.table(t(combn(safe_names, 2)))
  setnames(unq_pairs, new = c("Indiv1", "Indiv2"))
  
  # table_output <- data.table(apply(tmp_data_raw_3, 2, as.numeric))
  table_output <- cbind(Ind1=rep(indNames(x), each=nInd(x)), 
                        Ind2=rep(indNames(x), nInd(x)),
                        tmp_data_raw_3)
  table_output <- as.data.table(table_output)
  setkeyv(table_output, c("Indiv1", "Indiv2"))
  
  table_output <- table_output[J(unq_pairs), c(1, 2, 14:23), with=FALSE]
  
  #Relatedness
  # work on individual indices (unique by construction) and restore the
  # original names on the finished matrix
  df <- data.frame(ind1 = match(tmp_data_raw_3$Indiv1, safe_names),
                   ind2 = match(tmp_data_raw_3$Indiv2, safe_names),
                   rel = as.numeric(unlist(tmp_data_raw_3$`r(1,2)`)))

  res <- reshape2::acast(df, ind1 ~ ind2, value.var = "rel")
  res <- res[order(as.integer(rownames(res))),
             order(as.integer(colnames(res))), drop = FALSE]
  dimnames(res) <- list(hold_names, hold_names)
  # EMIBD9's r(1,2) is the kinship coefficient; gl.kin, gl.grm.network and
  # utils.kin.as.kinship read this tag
  attr(res, "scale") <- "kinship"

  # restore original individual names in the raw table
  tmp_data_raw_3$Indiv1 <- hold_names[match(tmp_data_raw_3$Indiv1, safe_names)]
  tmp_data_raw_3$Indiv2 <- hold_names[match(tmp_data_raw_3$Indiv2, safe_names)]

# Inbreeding 
 inbreedStart <- which(grepl("^Indiv genotypes at polymorphic loci", x_lines)) + 1
 if(length(inbreedStart)>0) {
   if (verbose >= 2){
     cat(
       report("  Exporting individual diversity and inbreeding values\n"))
   }
   
   inbTable <- fread(file = OutFileName, nrows = nInd(x), skip = inbreedStart)
   # restore original individual names
   if ("Indiv" %in% colnames(inbTable)) {
     inbTable$Indiv <- hold_names[match(inbTable$Indiv, safe_names)]
   }
 }
 
  #return to old path
  setwd(old.path)
  
  #compile the two dataframes into on list for output

  if (verbose >= 3){
    cat(
      report(
        "  Returning a list with the EMIBD9 results:\n",
        "          $rel -- a square matrix of pairwise kinship\n",
        "          $raw -- raw EMIBD9 results table\n",
        "          $processed -- EMIBD9 results without self and redundant comparisons\n",
        "          $inbreeding -- individual diversity and inbreeding\n"
      )
    )
  }
  
  # PRINTING OUTPUTS
  
  if (plot.out || !is.null(plot.file)) {
    
    if (is.null(palette_convergent)) {
      palette_convergent <- gl.colors("div")
    } 
    
    # gl.plot.heatmap returns the plot only while drawing it; to save a plot
    # that is not displayed, draw it on a null device
    if (!plot.out) {
      grDevices::pdf(NULL)
    }
    p1 <- gl.plot.heatmap(res,
                          palette.divergent = palette_convergent,
                          plot.out = TRUE,
                          verbose = 0)
    if (!plot.out) {
      grDevices::dev.off()
    }
  }
  
  # copy the EMIBD9 output file to outpath
  file.copy(from = file.path(rundir, outfile),
            to = file.path(outpath, outfile),
            overwrite = TRUE)
  
  # Optionally save the plot ---------------------
  if (!is.null(plot.file)) {
    tmp <- utils.plot.save(p1,
                           dir = plot.dir,
                           file = plot.file,
                           verbose = 0)
  }
  
  # FLAG SCRIPT END
  
  if (verbose >= 1) {
    cat(report("Completed:", funname, "\n"))
  }
  
  # Make a list
  results <-
    list(
      rel = res,
      raw = tmp_data_raw_3,
      processed = table_output)
  
      if(length(inbreedStart) > 0) {
        results[["inbreeding"]] <- inbTable
      }
  
  return(results)
}
