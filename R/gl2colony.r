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
#' @param filename_out Path to the output COLONY2 file [default "colony2.dat"].
#' @param project_name Project name to include in the file header
#'  [default 'my_project'].
#' @param output_name Output name to include in the file header
#' [default 'my_project'].
#' @param probability_father Probability that the father of an offspring is
#' included among candidates [default 0.5].
#' @param probability_mother Probability that the mother of an offspring is
#' included among candidates [default 0.5].
#' @param seed Seed for the random number generator [default NULL].
#' @param update_allele_freq 0 = do not update allele frequencies; 1 = update
#' [default 0].
#' @param di_mono_ecious 2 = dioecious species; 1 = monoecious species
#' [default 2].
#' @param inbreed 0 = no inbreeding; 1 = inbreeding allowed [default 0].
#' @param haplodiploid 0 = diploid species; 1 = haplodiploid species
#'  [default 0].
#' @param polygamy_male 0 = polygamy; 1 = monogamy for males [default 0].
#' @param polygamy_female 0 = polygamy; 1 = monogamy for females [default 0].
#' @param clone_inference 0 = no clone inference; 1 = infer clones [default 1].
#' @param scale_shibship 0 = do not scale full sibship; 1 = scale [default 1].
#' @param sibship_prior 0–4 specifying sibship prior strength (No, Weak,
#' Medium, Strong, Optimal) [default 0].
#' @param known_allele_freq 0 = unknown allele frequencies; 1 = known
#' [default 0].
#' @param num_runs Number of runs [default 1].
#' @param length_run 1–4 specifying run length (short, medium, long, very
#' long) [default 2].
#' @param monitor_method 0 = monitor by iteration number; 1 = monitor by time
#'  (seconds) [default 0].
#' @param monitor_interval Interval for monitoring (either iteration count or
#'  seconds) [default 10000].
#' @param windows_gui 0 = no Windows GUI; 1 = use Windows GUI [default 0].
#' @param likelihood 0–2 specifying likelihood scoring (PairLikelihood,
#' FullLikelihood, FPLS) [default 0].
#' @param precision_fl 0–3 specifying precision level for full-likelihood (Low,
#'  Medium, High, VeryHigh) [default 2].
#' @param marker_id Marker IDs string for all loci [default 'mk@'].
#' @param marker_type Marker types string for all loci (0@ for codominant, 1@
#' for dominant) [default '0@'].
#' @param allelic_dropout Allelic dropout rate string per locus
#' [default '0.000@'].
#' @param other_typ_err Other typing error rate string per locus
#' [default '0.05@'].
#' @param paternity_exclusion_threshold Threshold for paternity exclusion
#' ("0 0") [default '0 0'].
#' @param maternity_exclusion_threshold Threshold for maternity exclusion
#' ("0 0") [default '0 0'].
#' @param paternal_sibship Number of known paternal sibships [default 0].
#' @param maternal_sibship Number of known maternal sibships [default 0].
#' @param excluded_paternity Number of offspring with excluded paternity
#' [default 0].
#' @param excluded_maternity Number of offspring with excluded maternity
#'  [default 0].
#' @param excluded_paternal_sibships Number of excluded paternal sibships
#'  [default 0].
#' @param excluded_maternity_sibships Number of excluded maternal sibships
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
#' #            filename_out = "colony2.dat",
#' #            project_name = "parentage_fish_2022",
#' #            output_name = "parentage_fish_jul_2022",
#' #            seed = 1234,
#' #            probability_father = 0.6,
#' #            probability_mother = 0.4,
#' #            update_allele_freq = 1,
#' #            allelic_dropout = '0.01',
#' #            other_typ_err = '0.001')
#'
#' @references
#' Wang, J. (2011). COLONY: a program for parentage and sibship inference
#' from multilocus genotype data. Molecular Ecology Resources 10: 551–555.
#'
#' @importFrom utils write.table
#' @export

gl2colony <- function(x,
                      filename_out = "colony2.dat",
                      project_name = 'my_project',
                      output_name = 'my_project',
                      probability_father = 0.5,
                      probability_mother = 0.5,
                      seed = NULL,
                      update_allele_freq = 0,
                      di_mono_ecious = 2,
                      inbreed = 0,
                      haplodiploid = 0,
                      polygamy_male = 0,
                      polygamy_female = 0,
                      clone_inference = 1,
                      scale_shibship = 1,
                      sibship_prior = 0,
                      known_allele_freq = 0,
                      num_runs = 1,
                      length_run = 2,
                      monitor_method = 0,
                      monitor_interval = 10000,
                      windows_gui = 0,
                      likelihood = 0,
                      precision_fl = 2,
                      marker_id = 'mk@',
                      marker_type = '0@',
                      allelic_dropout = '0.000@',
                      other_typ_err = '0.05@',
                      paternity_exclusion_threshold = '0 0',
                      maternity_exclusion_threshold = '0 0',
                      paternal_sibship = 0,
                      maternal_sibship = 0,
                      excluded_paternity = 0,
                      excluded_maternity = 0,
                      excluded_paternal_sibships = 0,
                      excluded_maternity_sibships = 0,
                      verbose = NULL) {
  # SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)
  
  # FLAG SCRIPT START
  funname <- match.call()[[1]]
  utils.flag.start(func = funname,
                   build = "Jody",
                   verbose = verbose)
  
  # CHECK DATATYPE
  datatype <- utils.check.datatype(x, verbose = verbose)
  
  # PARAMETER VALIDATION
  if (missing(filename_out) || filename_out == "") {
    stop(error("'filename_out' is required."))
  }
  
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
  offspring_ids <- ids$offs
  dad_ids       <- ids$dad
  mum_ids       <- ids$mum
  
  # COUNTS
  n_offspring <- length(offspring_ids)
  n_dads      <- length(dad_ids)
  n_mums      <- length(mum_ids)
  loci        <- nLoc(x)
  n_total     <- n_offspring + n_dads + n_mums
  
  cat(report(
    sprintf(
      '%d offspring, %d fathers, %d mothers detected.',
      n_offspring,
      n_dads,
      n_mums
    )
  ), "\n")
  
  # WARN IF OFFSPRING MISSING
  if (n_offspring == 0) {
    stop(error('No offspring IDs found in metadata.'))
  }
  
  # CONVERT TO STRUCTURE FORMAT
  cat(report('Exporting genlight object to COLONY2 format...'), "\n")
  struct_mat <- gl2structure(x)
  
  # SUBSET GENOTYPES
  offspring_gen <- struct_mat[offspring_ids, , drop = FALSE]
  
  mum_gen       <- if (n_mums > 0){
    struct_mat[mum_ids, , drop = FALSE]
  }else{
    NULL
  }
  
  dad_gen       <- if (n_dads > 0){
    struct_mat[dad_ids, , drop = FALSE]
  }else{
    NULL
  }
  
  if (n_mums == 0){
    probability_mother <- 0
  }
  
  if (n_dads == 0){
    probability_father <- 0
  }
  
  # PREPARE HEADER
  head_comments <- c(
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
  head_values <- list(
    n_offspring,
    loci,
    seed,
    update_allele_freq,
    di_mono_ecious,
    inbreed,
    haplodiploid,
    paste(polygamy_male, polygamy_female),
    clone_inference,
    scale_shibship,
    sibship_prior,
    known_allele_freq,
    num_runs,
    length_run,
    monitor_method,
    monitor_interval,
    windows_gui,
    likelihood,
    precision_fl,
    '',
    marker_id,
    marker_type,
    allelic_dropout,
    other_typ_err
  )
  
  # WRITE HEADER
  sink(filename_out)
  cat(project_name, '\n')
  cat(output_name, '\n')
  for (i in seq_along(head_values)) {
    cat(head_values[[i]], '\t', head_comments[i], '\n')
  }
  sink()
  
  # WRITE OFFSPRING
  write.table(
    offspring_gen,
    file = filename_out,
    append = TRUE,
    quote = FALSE,
    col.names = FALSE
  )
  
  # WRITE CANDIDATE PROBABILITIES
  sink(filename_out, append = TRUE)
  cat('\n')
  cat(
    paste(probability_father, probability_mother),
    '\t',
    '! Parental inclusion probabilities',
    '\n'
  )
  cat(paste(n_dads, n_mums), '\t', '! Number of candidates', '\n')
  cat('\n')
  sink()
  
  # WRITE DADS
  if (n_dads > 0) {
    cat(report('Writing paternal genotypes...'), "\n")
    write.table(
      dad_gen,
      file = filename_out,
      append = TRUE,
      quote = FALSE,
      col.names = FALSE
    )
  }
  
  # WRITE MUMS
  if (n_mums > 0) {
    cat(report('Writing maternal genotypes...'), "\n")
    write.table(
      mum_gen,
      file = filename_out,
      append = TRUE,
      quote = FALSE,
      col.names = FALSE
    )
  }
  
  # WRITE EXCLUSION & SIBSHIP PARAMETERS
  last_comments <- c(
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
  last_values <- list(
    paternity_exclusion_threshold,
    maternity_exclusion_threshold,
    '',
    paternal_sibship,
    maternal_sibship,
    '',
    excluded_paternity,
    excluded_maternity,
    '',
    excluded_paternal_sibships,
    excluded_maternity_sibships
  )
  sink(filename_out, append = TRUE)
  cat('\n')
  for (i in seq_along(last_values)) {
    cat(last_values[[i]], '\t', last_comments[i], '\n')
  }
  sink()
  
  # COMPLETION MESSAGE
  cat(report('(100%) COLONY2 file successfully exported!'), "\n")
  
  # FLAG SCRIPT END
  if (verbose >= 1) {
    cat(report('\nCompleted:, ', funname, '\n'))
  }
  
  return(invisible(filename_out))
}


######################### Define function parental.ids #########################
## This function extracts parental information in a list of 3 elements (vectors
## with offspring, dads and mums IDs, respectively)
parental.ids <- function(gen_data) {
  # Read metadata and convert to lowercase
  indv.metadata <- gen_data@other$ind.metrics
  names(indv.metadata) <- tolower(names(indv.metadata))
  
  # Remove leading/trailing white spaces
  indv.metadata$mother    <- tolower(indv.metadata$mother)
  indv.metadata$father    <- tolower(indv.metadata$father)
  indv.metadata$offspring <- tolower(indv.metadata$offspring)
  
  # Subset metadata
  mum_ids  <- indv.metadata[indv.metadata$mother    %in% c("yes", " yes", "yes "), 'id']
  dad_ids  <- indv.metadata[indv.metadata$father    %in% c("yes", " yes", "yes "), 'id']
  offs_ids <- indv.metadata[indv.metadata$offspring %in% c("yes", " yes", "yes "), 'id']
  
  # Make them vectors
  mum_ids  <- as.vector(na.omit(mum_ids))
  dad_ids  <- as.vector(na.omit(dad_ids))
  offs_ids <- as.vector(na.omit(offs_ids))
  
  # Make a list with the 3 vectors
  x = list(offs = offs_ids, dad = dad_ids, mum = mum_ids)
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
