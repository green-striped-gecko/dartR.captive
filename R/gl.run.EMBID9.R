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
#'  Please note there are 3 different executables depending on your OS:
#'  EM_IBD_P.exe (=Windows) [and the two dlls impi.dll, libiomp5md.dll], 
#'  EM_IBD_P (=Mac, Linux). You only need to point
#'  to the folder (the function will recognise which OS you are running)
#'  [default getwd()].
#' @param Inbreed A Boolean, taking values 0 or 1 to indicate inbreeding is not
#'  and is allowed in estimating IBD coefficients [default 1].
#' @param GtypeFile A string, giving the path and name of the genotype file
#' [default "EMIBD9_Gen.dat"].
#' @param OutFileName_par A string, giving the path and name of the parameter
#'  file [default "MyData.par"].
#' @param OutFileName A string, giving the path and name of the output file
#' [default "EMIBD9_Res.ibd9"].
#' @param ISeed An integer used to seed the random number generator [default 52].
#' @param plot.dir Directory to save the plot RDS files [default as specified 
#' by the global working directory or tempdir()]
#' @param plot.file Name for the RDS binary file to save (base name only, exclude extension) [default NULL]
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#'  [default NULL, unless specified using gl.set.verbosity]
#' @details
#' Download the program from here:
#'
#' https://www.zsl.org/about-zsl/resources/software/emibd9
#'
#' For Windows, install the program and then move the following files to your
#'  working directory: "EM_IBD_P.exe", "impi.dll" and "libiomp5md.dll".
#'
#' For Mac move the file "EM_IBD_P" to the working directory.
#'
#' @return A matrix with pairwise relatedness
#' @author Custodian: Luis Mijangos -- Post to
#' \url{https://groups.google.com/d/forum/dartr}
#' @examples
#' \dontrun{
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
                          Inbreed = TRUE,
                          outfile=
                          GtypeFile = "EMIBD9_Gen.dat",
                          OutFileName_par = "MyData.par",
                          ISeed = 52,
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
  # GtypeFile <- GtypeFile
  # OutFileName <- OutFileName
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
    file = file.path(outpath, OutFileName_par)
  )

  IndivID <- paste(indNames(x2))

  gl_mat <- as.matrix(x2)
  gl_mat[is.na(gl_mat)] <- 3

  tmp <- cbind(apply(gl_mat, 1, function(y) {
    Reduce(paste0, y)
  }))

  tmp <- rbind(paste(indNames(x), collapse = " "), tmp)

  write.table(tmp,
    file = file.path(outpath,GtypeFile),
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )
  
  # Find executable makeblastdb if unix
  if (grepl("unix", .Platform$OS.type, ignore.case = TRUE)) {
    system("./EM_IBD_P INP:MyData.par")
  }
  ## if windows
  if (!grepl("unix", .Platform$OS.type, ignore.case = TRUE)) {
    system("EM_IBD_P.exe INP:MyData.par")
  }

  x_lines <- readLines("EMIBD9_Res.ibd9")
  strt <- which(grepl("^IBD", x_lines)) + 2
  stp <- which(grepl("Indiv genotypes", x_lines)) - 4
  linez_headings <- x_lines[strt]
  linez_data <- x_lines[(strt + 1):stp]
  tmp_headings <- unlist(stringr::str_split(linez_headings, " "))
  tmp_data <- stringr::str_split(linez_data, " ")
  tmp_data_2 <- lapply(tmp_data, "[", c(2, 3, 22))
  tmp_data_3 <- do.call("rbind", tmp_data_2)
  tmp_data_4 <- as.data.frame(tmp_data_3)
  tmp_data_4$V3 <- lapply(tmp_data_4$V3, as.numeric)
  colnames(tmp_data_4) <- c("ind1", "ind2", "rel")

  res <- as.matrix(reshape2::acast(tmp_data_4, ind1 ~ ind2, value.var = "rel"))

  restore_names_2 <- merge(data.frame(id2 = rownames(res)), restore_names, by = "id2")

  res <- apply(res, 2, as.numeric)

  colnames(res) <- restore_names_2$id
  rownames(res) <- restore_names_2$id

  order_mat <- colnames(res)[order(colnames(res))]

  out <- list()

  out$res <- res[order_mat, order_mat]
  out$table <- 1:3
  
  if(!is.null(plot.file)){
    tmp <- utils.plot.save(p3,
                           dir=plot.dir,
                           file=plot.file,
                           verbose=verbose)
  }
  
  
  return(res)
}
