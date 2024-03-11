#' @name gl.run.EMIBD9
#' @title Run program EMIBD9
#' @description
#' Run program EMIBD9
#' @param x Name of the genlight object containing the SNP data [required].
#' @param outfile A string, giving the path and name of the output file
#' [default "EMIBD9_Res.ibd9"].
#' @param outpath Path where to save the output file. Use outpath=getwd() or
#' outpath='.' when calling this function to direct output files to your working 
#' or current directory [default tempdir(), mandated by CRAN].
#' @param emibd9.path Path to the folder emidb files.
#'  Please note there are 2 different executables depending on your OS:
#'  EM_IBD_P.exe (=Windows) EM_IBD_P (=Mac, Linux). 
#'  You only need to point to the folder (the function will recognise which OS you
#'  are running) [default getwd()].
#' @param Inbreed A Boolean, taking values 0 or 1 to indicate inbreeding is
#' not and is allowed in estimating IBD coefficients [default 1].
#' @param nclusters The number of cores or processor threads used during 
#' calculations [default 1].
#' @param ISeed An integer used to seed the random number generator [default 42].
#' @param plot.out A boolean that indicates whether to plot the results [default TRUE].
#' @param plot.dir Directory to save the plot RDS files [default as specified 
#' by the global working directory or tempdir()]
#' @param plot.file Name for the RDS binary file to save (base name only, exclude extension) [default NULL]
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#'  [default NULL, unless specified using gl.set.verbosity]
#' @details
#' 'The results of EMIBD9 include the identical in state (IIS) values for each mode 
#'(S1 - 9) and nine condensed identical by descent (IBD) modes (∆1 - ∆9) as well as #'the relatedness coefficient (r). Alleles are IIS if they are the same. Similarly,
#' IBD describes a matching allele between two individuals that has been inherited from a common ancestor or common gene. In a pairwise comparison, ∆1 to ∆9 are the
#'  probabilities associated with each IBD mode. In inbreeding populations, only  ∆1
#'   to  ∆6 can can occur. In contrast, ∆7 to ∆9 can only occur in large, panmictic
#'    outbred populations. 
#' 
#'EMIBD9 uses an expectation maximization (EM) algorithm based on the maximum
#' likelihood expectations (MLE) of ∆ to estimate both allele frequencies (p) and ∆
#'  jointly from genotype data. By iteratively calculating p and ∆, relatedness 
#'  can be modified to reduce biases due to small sample sizes. Wang J. (2022) 
#'  suggest the resulting r coefficient is therefore more robust compared to 
#'  previous methods.
#'
#'The kinship coefficient is the probability that two alleles at a random locus
#'  drawn from two individuals are IBD.
#'
#'Below is a table modified from Speed & Balding (2015) showing kinship values,
#'and their confidence intervals (CI), for different relationships that could 
#'be used to guide the choosing of the relatedness threshold in the function.
#'
#'|Relationship                               |Kinship  |     95% CI       |
#'
#'|Identical twins/clones/same individual     | 0.5     |        -         |
#'
#'|Sibling/Parent-Offspring                   | 0.25    |    (0.204, 0.296)|
#'
#'|Half-sibling                               | 0.125   |    (0.092, 0.158)|
#'
#'|First cousin                               | 0.062   |    (0.038, 0.089)|
#'
#'|Half-cousin                                | 0.031   |    (0.012, 0.055)|
#'
#'|Second cousin                              | 0.016   |    (0.004, 0.031)|
#'
#'|Half-second cousin                         | 0.008   |    (0.001, 0.020)|
#' 
#'|Third cousin                               | 0.004   |    (0.000, 0.012)|
#'
#'|Unrelated                                  | 0       |        -         | 
#'
#'For greater detail on the methods employed by EMIBD9, we encourage you to read Wang, J. (2022). 
#'
#' Download the program from here:
#'
#' https://www.zsl.org/about-zsl/resources/software/emibd9
#'
#' For Windows, Mac and Linux install the program then point to the folder where you find:
#' EM_IBD_P.exe (=Windows) and EM_IBD_P (=Mac, Linux). If running really slow you may 
#' want to create the files using the function and then run in parallel using the
#' documentation provided by the authors [you need to have mpiexec installed].
#' 
#'
#' @return A matrix with pairwise relatedness
#' @author Custodian: Luis Mijangos -- Post to
#' \url{https://groups.google.com/d/forum/dartr}
#' @examples
#' \dontrun{
#' #To run this function needs EMIBD9 installed in your computer
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
#' @importFrom stringr str_split
#' @export

gl.run.EMIBD9 <- function(x,
                          outfile = "EMIBD9_Res.ibd9",
                          outpath = tempdir(),
                          emibd9.path = getwd(),
                          Inbreed = TRUE,
                          nclusters = 1,
                          ISeed = 42,
                          plot.out = TRUE,
                          plot.dir=NULL,
                          plot.file = NULL,
                          verbose = NULL) {
  
  
  # SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)
  
  # SET WORKING DIRECTORY
  plot.dir <- gl.check.wd(plot.dir, verbose = 0)
  
  # FLAG SCRIPT START
  funname <- match.call()[[1]]
  utils.flag.start(
    func = funname,
    build = "Jody",
    verbose = verbose
  )
  
  # CHECK DATATYPE
  datatype <- utils.check.datatype(x, verbose = verbose)
 
  #check if embid9 is available
  
  os <- Sys.info()["sysname"]
  
  if (Sys.info()["sysname"] == "Windows") {
    prog <- c("EM_IBD_P.exe", "impi.dll", "libiomp5md.dll", "mpiexec.exe")
    nclusters = nclusters
    cmd <- "EM_IBD_P.exe INP:MyData.par"
    #cmd <- paste("mpiexec -n ", nclusters , " EM_IBD_P INP:MyData.par", sep = "")
  }
  
  #note: I am working on the windows function first, then we can deal with Mac and Linux...
  if (Sys.info()["sysname"] == "Linux") {
    prog <- "EM_IBD_P"
    cmd <- "mpirun -register M EM_IBD_P_mpi INP:MyData.par"
  }
  
  if (Sys.info()["sysname"] == "Darwin") {
    prog <- "EM_IBD_P"
    cmd <- "mpirun -n M EM_IBD_P_impi INP:MyData.par"
  }
  
  # check if file program can be found
  if (all(file.exists(file.path(emibd9.path, prog)))) {
    file.copy(file.path(emibd9.path, prog),
              to = tempdir(),
              overwrite = TRUE
    )
    if (verbose > 0) cat(report("Found necessary files to run EMIBD9."))
  } else {
    message(
      error(
        "  Cannot find",
        prog,
        "in the specified folder given by emibd9.path:",
        emibd9.path,
        "\n"
      )
    )
    stop()
  }
  
  
  rundir <- tempdir()
  
   
  
  # individual IDs must have a maximal length of 20 characters. The IDs must NOT
  # contain blank space and other illegal characters (such as /), and must be
  # unique among all sampled individuals (i.e. NO duplications). Any string longer
  # than 20 characters for individual ID will be truncated to have 20 characters.

  
  
  x2 <- x  #copy to work only on the copied data set
  hold_names <- indNames(x)
  indNames(x2) <- 1:nInd(x2)
  
  NumIndiv <- nInd(x2)
  NumLoci <- nLoc(x2)
  DataForm <- 2
  if (Inbreed) Inbreed <- 1 else Inbreed <- 0
  # Inbreed <- Inbreed
  GtypeFile <- "EMIBD9_Gen.dat"
  OutFileName <-  outfile
  # ISeed <- ISeed
  RndDelta0 <- 1
  EM_Method <- 1
  OutAlleleFre <- 0

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

  write.table(param,
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE,
    file = file.path(rundir, "MyData.par")
  )

  IndivID <- paste(indNames(x2))

  gl_mat <- as.matrix(x2)
  gl_mat[is.na(gl_mat)] <- 3

  tmp <- cbind(apply(gl_mat, 1, function(y) {
    Reduce(paste0, y)
  }))

  tmp <- rbind(paste(indNames(x2), collapse = " "), tmp)

  write.table(tmp,
    file = file.path(rundir, "EMIBD9_Gen.dat"),
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )
  

  # run EMIBD9
  # change into tempdir (run it there)
  old.path <- getwd()
  on.exit(setwd(old.path))
  setwd(rundir)
  system(cmd)
  
  ### get output  
  
  x_lines <- readLines("EMIBD9_Res.ibd9")
  strt <- which(grepl("^IBD", x_lines)) + 2
  stp <- which(grepl("Indiv genotypes", x_lines)) - 4
  linez_headings <- x_lines[strt]
  linez_data <- x_lines[(strt + 1):stp]
  tmp_headings <- unlist(stringr::str_split(linez_headings, " "))
  tmp_data <- stringr::str_split(linez_data, " ")
  #Raw data 
  tmp_data_raw_1 <- lapply(tmp_data, "[", c(2:22))
  tmp_data_raw_2 <- do.call("rbind", tmp_data_raw_1)
  tmp_data_raw_3 <- as.data.frame(tmp_data_raw_2)
  tmp_data_raw_3$V3 <- lapply(tmp_data_raw_3$V3, as.numeric)
  colnames(tmp_data_raw_3) <- tmp_headings[2:22]
  
  df <- data.frame(ind1=tmp_data_raw_3$Indiv1, ind2=tmp_data_raw_3$Indiv2,rel= tmp_data_raw_3$`r(1,2)`)
  df<- apply(df, 2, as.numeric)
  #Relatedness
  res <- matrix(NA, nrow = nInd(x), ncol = nInd(x))
  
  for (i in 1:nrow(df)) {
    res[df[i, 1], df[i, 2]] <- df[i, 3]

rel.list <- data.frame(as.numeric(tmp_data_raw_3[, 21]))

inbreeding <- data.frame(as.numeric(tmp_data_raw_3[, 21]))
  }

  colnames(res) <- indNames(x)
  rownames(res) <- indNames(x)
  
  colnames(rel.list) <- c("Relatedness")
  
  #return to old path
  setwd(old.path)
  
  #compile the two dataframes into on list for output
  if (verbose>0)
  {
  cat(
    report(
      "Returning a list containing the input gl object, a square matrix  of pairwise kinship, and the raw EMIBD9 results table as follows:\n",
      "          $rel -- a square matrix of relatedness \n",
      "          $raw -- raw EMIBD9 results table \n")
  )
  }

  # PRINTING OUTPUT
  
 p1 <- heatmap.2(res, 
                 dendrogram = c("none"),
                 key = T, 
                 main = "Heatmap of Relatedness", 
                 col = viridis(60), 
                 trace = c("none"), 
                 scale = "column")
 
                 
 p2 <- ggplot(rel.list, aes(x = Relatedness)) +
    geom_histogram(binwidth = 0.01) +
    ggtitle("Histogram of relatedness") +
    theme_classic()
 
 p3 <- list(p1, p2)
 
 if (plot.out)  print(p3)

  # Optionally save the plot ---------------------
  if(!is.null(plot.file)){
    tmp <- utils.plot.save(p2,
                           dir=plot.dir,
                           file=plot.file,
                           verbose=verbose)
  }
  
  #Make a list
  results <-
    list(
      rel = res,
      raw = tmp_data_raw_3
    )
  
  return(results)
}
