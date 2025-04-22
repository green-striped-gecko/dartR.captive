#' @name gl2colony
#' @title Export a COLONY2 input file from a genlight object
#' @description
#' Export a formatted text file compatible with the COLONY2 software from a
#'  \code{genlight} object containing parental and offspring information
#'  stored in the individual metadata.
#'
#' @param x A \code{genlight} object with individual metadata columns
#' 'offspring', 'mother', and 'father' indicating 'yes'/'no' for each sample
#' [required].
#' @param outfile File name of the output file (including extension)
#' [default "colony2.dat"].
#' @param outpath Path where to save the output file [default global working 
#' directory or if not specified, tempdir()].
#' @param project.name Project name to include in the file header
#'  [default 'my_project'].
#' @param output.name Output name to include in the file header
#' [default 'my_project'].
#' @param probability.father Probability that the father of an offspring is
#' included among candidates [default 0.5].
#' @param probability.mother Probability that the mother of an offspring is
#' included among candidates [default 0.5].
#' @param seed Seed for the random number generator [default NULL].
#' @param update.allele.freq 0 = do not update allele frequencies; 1 = update
#' [default 0].
#' @param di.mono.ecious 2 = dioecious species; 1 = monoecious species
#' [default 2].
#' @param inbreed 0 = no inbreeding; 1 = inbreeding allowed [default 0].
#' @param haplodiploid 0 = diploid species; 1 = haplodiploid species
#'  [default 0].
#' @param polygamy.male 0 = polygamy; 1 = monogamy for males [default 0].
#' @param polygamy.female 0 = polygamy; 1 = monogamy for females [default 0].
#' @param clone.inference 0 = no clone inference; 1 = infer clones [default 1].
#' @param scale.shibship 0 = do not scale full sibship; 1 = scale [default 1].
#' @param sibship.prior 0–4 specifying sibship prior strength (No, Weak,
#' Medium, Strong, Optimal) [default 0].
#' @param known.allele.freq 0 = unknown allele frequencies; 1 = known
#' [default 0].
#' @param num.runs Number of runs [default 1].
#' @param length.run 1–4 specifying run length (short, medium, long, very
#' long) [default 2].
#' @param monitor.method 0 = monitor by iteration number; 1 = monitor by time
#'  (seconds) [default 0].
#' @param monitor.interval Interval for monitoring (either iteration count or
#'  seconds) [default 10000].
#' @param windows.gui 0 = no Windows GUI; 1 = use Windows GUI [default 0].
#' @param likelihood 0–2 specifying likelihood scoring (PairLikelihood,
#' FullLikelihood, FPLS) [default 0].
#' @param precision.fl 0–3 specifying precision level for full-likelihood (Low,
#'  Medium, High, VeryHigh) [default 2].
#' @param marker.id Marker IDs string for all loci [default 'mk@'].
#' @param marker.type Marker types string for all loci (0@ for codominant, 1@
#' for dominant) [default '0@'].
#' @param allelic.dropout Allelic dropout rate string per locus
#' [default '0.000@'].
#' @param other.typ.err Other typing error rate string per locus
#' [default '0.05@'].
#' @param paternity.exclusion.threshold Threshold for paternity exclusion
#' ("0 0") [default '0 0'].
#' @param maternity.exclusion.threshold Threshold for maternity exclusion
#' ("0 0") [default '0 0'].
#' @param paternal.sibship Number of known paternal sibships [default 0].
#' @param maternal.sibship Number of known maternal sibships [default 0].
#' @param excluded.paternity Number of offspring with excluded paternity
#' [default 0].
#' @param excluded.maternity Number of offspring with excluded maternity
#'  [default 0].
#' @param excluded.paternal.sibships Number of excluded paternal sibships
#'  [default 0].
#' @param excluded.maternity.sibships Number of excluded maternal sibships
#' [default 0].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#'  progress log ; 3, progress and results summary; 5, full report
#'  [default 2 or as specified using gl.set.verbosity].
#'
#' @details
#' This function formats and writes a COLONY2-compatible text file, including
#' header, offspring genotypes, parental candidate probabilities, and
#' candidate genotypes, based on the \code{genlight} object's individual
#' metadata and genotype matrix.
#'
#' @return
#' Invisibly returns the output filename.
#'
#' @author
#' Jesús Castrejón-Figueroa, Diana A. Robledo-Ruiz -- Post to
#' \url{https://groups.google.com/d/forum/dartr}
#'
#' @examples
#' # gl2colony(x = platypus.gl,
#' #            project.name = "parentage_fish_2022",
#' #            output.name = "parentage_fish_jul_2022",
#' #            seed = 1234,
#' #            probability.father = 0.6,
#' #            probability.mother = 0.4,
#' #            update.allele.freq = 1,
#' #            allelic.dropout = '0.01',
#' #            other.typ.err = '0.001')
#'
#' @references
#' Wang, J. (2011). COLONY: a program for parentage and sibship inference
#' from multilocus genotype data. Molecular Ecology Resources 10: 551–555.
#'
#' @importFrom utils write.table
#' @export

gl2colony <- function(x,
                      outfile = "colony2.dat",
                      outpath = NULL,
                      project.name = 'my_project',
                      output.name = 'my_project',
                      probability.father = 0.5,
                      probability.mother = 0.5,
                      seed = NULL,
                      update.allele.freq = 0,
                      di.mono.ecious = 2,
                      inbreed = 0,
                      haplodiploid = 0,
                      polygamy.male = 0,
                      polygamy.female = 0,
                      clone.inference = 1,
                      scale.shibship = 1,
                      sibship.prior = 0,
                      known.allele.freq = 0,
                      num.runs = 1,
                      length.run = 2,
                      monitor.method = 0,
                      monitor.interval = 10000,
                      windows.gui = 0,
                      likelihood = 0,
                      precision.fl = 2,
                      marker.id = 'mk@',
                      marker.type = '0@',
                      allelic.dropout = '0.000@',
                      other.typ.err = '0.05@',
                      paternity.exclusion.threshold = '0 0',
                      maternity.exclusion.threshold = '0 0',
                      paternal.sibship = 0,
                      maternal.sibship = 0,
                      excluded.paternity = 0,
                      excluded.maternity = 0,
                      excluded.paternal.sibships = 0,
                      excluded.maternity.sibships = 0,
                      verbose = NULL) {
  # SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)
  
  # SET WORKING DIRECTORY
  outpath <- gl.check.wd(outpath,verbose=0)
  outfilespec <- file.path(outpath, outfile)
  
  # FLAG SCRIPT START
  funname <- match.call()[[1]]
  utils.flag.start(func = funname,
                   build = "Jody",
                   verbose = verbose)
  
  # CHECK DATATYPE
  datatype <- utils.check.datatype(x, verbose = verbose)
  
  # SET RANDOM SEED
  if (is.null(seed)) {
    seed <- sample.int(65535, 1)
  }
  cat(code(sprintf('Random seed set to %d', seed)), "\n")
  
  if(any(!c("offspring","mother","father") %in% 
         colnames(x$other$ind.metrics))){
  x$other$ind.metrics$offspring <- "yes"
  x$other$ind.metrics$mother <- "no"
  x$other$ind.metrics$father <- "no"
  
  cat(warn(
    "  The colums offspring, mother and father were not found in the genligth object. Setting all the individuals as offspring.\n"
    ))
  }
  
  x$other$ind.metrics$id <- indNames(x)
  
  # EXTRACT PARENTAL IDS
  ids <- parental.ids(x)
  offspring.ids <- ids$offs
  dad.ids       <- ids$dad
  mum.ids       <- ids$mum
  
  # COUNTS
  n.offspring <- length(offspring.ids)
  n.dads      <- length(dad.ids)
  n.mums      <- length(mum.ids)
  loci        <- nLoc(x)
  n.total     <- n.offspring + n.dads + n.mums
  
  cat(report(
    sprintf(
      '%d offspring, %d fathers, %d mothers detected.',
      n.offspring,
      n.dads,
      n.mums
    )
  ), "\n")
  
  # WARN IF OFFSPRING MISSING
  if (n.offspring == 0) {
    stop(error('No offspring IDs found in metadata.'))
  }
  
  # CONVERT TO STRUCTURE FORMAT
  cat(report('Exporting genlight object to COLONY2 format...'), "\n")
  struct.mat <- gl2structure(x)
  
  # SUBSET GENOTYPES
  offspring.gen <- struct.mat[offspring.ids, , drop = FALSE]
  
  mum.gen       <- if (n.mums > 0){
    struct.mat[mum.ids, , drop = FALSE]
  }else{
    NULL
  }
  
  dad.gen       <- if (n.dads > 0){
    struct.mat[dad.ids, , drop = FALSE]
  }else{
    NULL
  }
  
  if (n.mums == 0){
    probability.mother <- 0
  }
  
  if (n.dads == 0){
    probability.father <- 0
  }
  
  # PREPARE HEADER
  head.comments <- c(
    '! No. offspring',
    '! No. of loci',
    '! Seed for RNG',
    '! 0/1 = update allele freq',
    '! 2/1 = dioecious/monoecious',
    '! 0/1 = no inbreeding/inbreeding',
    '! 0/1 = diploid/haplodiploid',
    '! polygamy male female',
    '! clone inference',
    '! scale sibship',
    '! sibship prior',
    '! known allele freq',
    '! num runs',
    '! run length',
    '! monitor method',
    '! monitor interval',
    '! windows GUI',
    '! likelihood',
    '! precision FL',
    '',
    '! Marker Ids',
    '! Marker types',
    '! Allelic dropout rate',
    '! Other typing error rate'
  )
  head.values <- list(
    n.offspring,
    loci,
    seed,
    update.allele.freq,
    di.mono.ecious,
    inbreed,
    haplodiploid,
    paste(polygamy.male, polygamy.female),
    clone.inference,
    scale.shibship,
    sibship.prior,
    known.allele.freq,
    num.runs,
    length.run,
    monitor.method,
    monitor.interval,
    windows.gui,
    likelihood,
    precision.fl,
    '',
    marker.id,
    marker.type,
    allelic.dropout,
    other.typ.err
  )
  
  # WRITE HEADER
  sink(outfilespec)
  cat(project.name, '\n')
  cat(output.name, '\n')
  for (i in seq_along(head.values)) {
    cat(head.values[[i]], '\t', head.comments[i], '\n')
  }
  sink()
  
  # WRITE OFFSPRING
  write.table(
    offspring.gen,
    file = outfilespec,
    append = TRUE,
    quote = FALSE,
    col.names = FALSE
  )
  
  # WRITE CANDIDATE PROBABILITIES
  sink(outfilespec, append = TRUE)
  cat('\n')
  cat(
    paste(probability.father, probability.mother),
    '\t',
    '! Parental inclusion probabilities',
    '\n'
  )
  cat(paste(n.dads, n.mums), '\t', '! Number of candidates', '\n')
  cat('\n')
  sink()
  
  # WRITE DADS
  if (n.dads > 0) {
    cat(report('Writing paternal genotypes...'), "\n")
    write.table(
      dad.gen,
      file = outfilespec,
      append = TRUE,
      quote = FALSE,
      col.names = FALSE
    )
  }
  
  # WRITE MUMS
  if (n.mums > 0) {
    cat(report('Writing maternal genotypes...'), "\n")
    write.table(
      mum.gen,
      file = outfilespec,
      append = TRUE,
      quote = FALSE,
      col.names = FALSE
    )
  }
  
  # WRITE EXCLUSION & SIBSHIP PARAMETERS
  last.comments <- c(
    '! Offspring known paternity threshold',
    '! known maternity threshold',
    '',
    '! known paternal sibship',
    '! known maternal sibship',
    '',
    '! excluded paternity',
    '! excluded maternity',
    '',
    '! excluded paternal sibships',
    '! excluded maternal sibships'
  )
  last.values <- list(
    paternity.exclusion.threshold,
    maternity.exclusion.threshold,
    '',
    paternal.sibship,
    maternal.sibship,
    '',
    excluded.paternity,
    excluded.maternity,
    '',
    excluded.paternal.sibships,
    excluded.maternity.sibships
  )
  sink(outfilespec, append = TRUE)
  cat('\n')
  for (i in seq_along(last.values)) {
    cat(last.values[[i]], '\t', last.comments[i], '\n')
  }
  sink()
  
  # COMPLETION MESSAGE
  cat(report('(100%) COLONY2 file successfully exported!'), "\n")
  
  if (verbose >= 3) {
    cat(report(paste(
      "Records written to", outfilespec, "\n"
    )))
  }
  
  # FLAG SCRIPT END
  if (verbose >= 1) {
    cat(report('\nCompleted:, ', funname, '\n'))
  }
  
  return(invisible(outfilespec))
}


######################### Define function parental.ids #########################
## This function extracts parental information in a list of 3 elements (vectors
## with offspring, dads and mums IDs, respectively)
parental.ids <- function(gen.data) {
  # Read metadata and convert to lowercase
  indv.metadata <- gen.data@other$ind.metrics
  names(indv.metadata) <- tolower(names(indv.metadata))
  
  # Remove leading/trailing white spaces
  indv.metadata$mother    <- tolower(indv.metadata$mother)
  indv.metadata$father    <- tolower(indv.metadata$father)
  indv.metadata$offspring <- tolower(indv.metadata$offspring)
  
  # Subset metadata
  mum.ids  <- indv.metadata[indv.metadata$mother    %in% c("yes", " yes", "yes "), 'id']
  dad.ids  <- indv.metadata[indv.metadata$father    %in% c("yes", " yes", "yes "), 'id']
  offs.ids <- indv.metadata[indv.metadata$offspring %in% c("yes", " yes", "yes "), 'id']
  
  # Make them vectors
  mum.ids  <- as.vector(na.omit(mum.ids))
  dad.ids  <- as.vector(na.omit(dad.ids))
  offs.ids <- as.vector(na.omit(offs.ids))
  
  # Make a list with the 3 vectors
  x = list(offs = offs.ids, dad = dad.ids, mum = mum.ids)
  return(x)
}
################################################################################


######################### Define function gl2structure #########################
## This function converts gl matrix to Structure format and from 2-row-per-ind
## to 1-row-per-ind
gl2structure <- function(x,
                         addtlColumns = NULL,
                         ploidy = 2,
                         exportMarkerNames = FALSE) {
  genmat <- as.matrix(x)
  indNames <- dimnames(genmat)[[1]]
  nInd <- dim(genmat)[1] # number of individuals
  
  # Make sets of possible genotypes
  G <- list()
  for (i in 0:ploidy) {
    G[[i + 1]] <- c(rep(1, ploidy - i), rep(2, i))
  }
  #G[[ploidy + 2]] <- rep(-9, ploidy) # for missing data
  G[[ploidy + 2]] <- rep(0, ploidy) # for missing data
  
  # Set up data frame for Structure
  StructTab <- data.frame(ind = rep(indNames, each = ploidy))
  
  # Add any additional columns
  if (!is.null(addtlColumns)) {
    for (i in 1:dim(addtlColumns)[2]) {
      StructTab <- data.frame(StructTab, rep(addtlColumns[, i], each = ploidy))
      if (!is.null(dimnames(addtlColumns)[[2]])) {
        names(StructTab)[i + 1] <- dimnames(addtlColumns)[[2]][i]
      } else {
        names(StructTab)[i + 1] <- paste("X", i, sep = "")
      }
    }
  }
  
  # Add genetic data
  for (i in 1:dim(genmat)[2]) {
    thesegen <- genmat[, i] + 1
    thesegen[is.na(thesegen)] <- ploidy + 2
    StructTab[[dimnames(genmat)[[2]][i]]] <- unlist(G[thesegen])
  }
  
  # return(StructTab)  # Returning the value of gl2struct dartR function
  
  data <- StructTab
  # Define dimensions of the matrix (only genotypes, not Ids)
  out <- matrix(NA, nrow = (nrow(data) / 2), # no. of rows divided by 2
                ncol = (2 * (ncol(data) - 1)))  # no. of columns minus Ids column times 2
  
  # Select first row per ind, leaving behind first column (Ids), then assign as first column per ind
  out[, seq(1, ncol(out), by = 2)] <- as.matrix(data[seq(1, nrow(data), by = 2), -1])
  # Select second row per ind, leaving behind first column (Ids), then assign as second column per ind
  out[, seq(2, ncol(out), by = 2)] <- as.matrix(data[seq(2, nrow(data), by = 2), -1])
  
  # Select Id column (only first row per ind) and make it rownames for matrix
  rownames(out) <- data[seq(1, nrow(data), by = 2), 1]
  return(out)
}
