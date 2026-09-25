#' @name gl2colony
#' @title Export a COLONY2 input file from a genlight object
#' @description
#' Export a formatted text file compatible with the COLONY2 software from a
#'  \code{genlight} object containing parental and offspring information
#'  stored in the individual metadata.
#'
#' @param x A \code{genlight} object with SNP data and individual metadata
#' columns 'offspring', 'mother', and 'father' indicating 'yes'/'no' for each
#' sample. Column names and values are matched ignoring case. A missing column
#' is filled with 'yes' (offspring) or 'no' (mother, father). Individual names
#' must not contain whitespace [required].
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
#' @param sibship.prior Sibship prior; only 0 (no prior) is supported, because
#' other values require mean sibship sizes that this function does not write
#' [default 0].
#' @param known.allele.freq 0 = unknown allele frequencies. Known
#' frequencies (1) are not supported [default 0].
#' @param num.runs Number of runs [default 1].
#' @param length.run 1-4 specifying run length (short, medium, long, very
#' long) [default 2].
#' @param monitor.method 0 = monitor by iteration number; 1 = monitor by time
#'  (seconds) [default 0].
#' @param monitor.interval Interval for monitoring (either iteration count or
#'  seconds) [default 10000].
#' @param windows.gui 0 = no Windows GUI; 1 = use Windows GUI [default 0].
#' @param likelihood 0-2 specifying likelihood scoring (PairLikelihood,
#' FullLikelihood, FPLS) [default 0].
#' @param precision.fl 0-3 specifying precision level for full-likelihood (Low,
#'  Medium, High, VeryHigh) [default 2].
#' @param marker.id Marker IDs string; a trailing '@' applies one value to
#' all loci [default 'mk@'].
#' @param marker.type Marker types string (0@ for codominant) [default '0@'].
#' @param allelic.dropout Allelic dropout rate string; a single value without
#' '@' (e.g. '0.01') is applied to all loci [default '0.000@'].
#' @param other.typ.err Other typing error rate string; a single value without
#' '@' is applied to all loci [default '0.05@'].
#' @param paternity.exclusion.threshold Number of offspring with known father
#' and the exclusion threshold. Only a count of 0 is supported, because other
#' values require a list of offspring-father pairs [default '0 0'].
#' @param maternity.exclusion.threshold Number of offspring with known mother
#' and the exclusion threshold. Only a count of 0 is supported [default '0 0'].
#' @param paternal.sibship Number of known paternal sibships; only 0 is
#' supported [default 0].
#' @param maternal.sibship Number of known maternal sibships; only 0 is
#' supported [default 0].
#' @param excluded.paternity Number of offspring with excluded paternity; only
#' 0 is supported [default 0].
#' @param excluded.maternity Number of offspring with excluded maternity; only
#' 0 is supported [default 0].
#' @param excluded.paternal.sibships Number of excluded paternal sibships; only
#' 0 is supported [default 0].
#' @param excluded.maternity.sibships Number of excluded maternal sibships;
#' only 0 is supported [default 0].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#'  progress log; 3, progress and results summary; 5, full report
#'  [default 2, unless specified using gl.set.verbosity].
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
#' Author(s): Jesús Castrejón-Figueroa, Diana A. Robledo-Ruiz. Custodian: Luis
#' Mijangos -- Post to \url{https://groups.google.com/d/forum/dartr}
#'
#' @examples
#' t1 <- testset.gl[1:30, 1:50]
#' t1@other$ind.metrics$offspring <- rep(c("yes", "no", "no"), each = 10)
#' t1@other$ind.metrics$father <- rep(c("no", "yes", "no"), each = 10)
#' t1@other$ind.metrics$mother <- rep(c("no", "no", "yes"), each = 10)
#' gl2colony(x = t1,
#'           outpath = tempdir(),
#'           seed = 1234,
#'           probability.father = 0.6,
#'           probability.mother = 0.4,
#'           allelic.dropout = '0.01',
#'           other.typ.err = '0.001')
#'
#' @references
#' Jones, O. R., & Wang, J. (2010). COLONY: a program for parentage and
#' sibship inference from multilocus genotype data. Molecular Ecology
#' Resources, 10(3), 551-555.
#'
#' @family captive management
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
  utils.flag.start(func = funname, verbose = verbose)
  
  # CHECK DATATYPE
  datatype <- utils.check.datatype(x, verbose = verbose)
  
  # FUNCTION SPECIFIC ERROR CHECKING
  # COLONY codes codominant genotypes; presence/absence would be exported as
  # heterozygotes
  if (datatype == "SilicoDArT") {
    stop(error(
      "  Only SNP data are supported; x contains SilicoDArT data\n"
    ))
  }
  
  # COLONY reads whitespace-delimited records, so a space in a name shifts
  # every allele of that individual by one position
  bad.names <- indNames(x)[grepl("[[:space:]]", indNames(x))]
  if (length(bad.names) > 0) {
    stop(error(
      "  Individual names must not contain whitespace. Rename:",
      paste(shQuote(bad.names), collapse = ", "), "\n"
    ))
  }
  
  # COLONY reads the on/off switches as integers and stops on TRUE/FALSE, so
  # accept logicals and write them as 0/1
  flags <- list(update.allele.freq = update.allele.freq, inbreed = inbreed,
                haplodiploid = haplodiploid, polygamy.male = polygamy.male,
                polygamy.female = polygamy.female,
                clone.inference = clone.inference,
                scale.shibship = scale.shibship,
                known.allele.freq = known.allele.freq,
                monitor.method = monitor.method, windows.gui = windows.gui)
  bad.flags <- names(flags)[!vapply(flags, function(f) {
    length(f) == 1 && !is.na(f) && (is.logical(f) || f %in% c(0, 1))
  }, logical(1))]
  if (length(bad.flags) > 0) {
    stop(error(
      "  These settings must be 0 or 1 (or FALSE/TRUE):",
      paste(bad.flags, collapse = ", "), "\n"
    ))
  }
  flags <- lapply(flags, as.integer)
  update.allele.freq <- flags$update.allele.freq
  inbreed <- flags$inbreed
  haplodiploid <- flags$haplodiploid
  polygamy.male <- flags$polygamy.male
  polygamy.female <- flags$polygamy.female
  clone.inference <- flags$clone.inference
  scale.shibship <- flags$scale.shibship
  known.allele.freq <- flags$known.allele.freq
  monitor.method <- flags$monitor.method
  windows.gui <- flags$windows.gui

  # these settings need extra data blocks (sibship sizes, allele
  # frequencies, lists of known or excluded relatives) that this function
  # does not write; COLONY rejects the file without them
  known.counts <- c(
    paternity.exclusion.threshold = as.numeric(
      strsplit(trimws(paternity.exclusion.threshold), "[[:space:]]+")[[1]][1]),
    maternity.exclusion.threshold = as.numeric(
      strsplit(trimws(maternity.exclusion.threshold), "[[:space:]]+")[[1]][1]),
    paternal.sibship = paternal.sibship,
    maternal.sibship = maternal.sibship,
    excluded.paternity = excluded.paternity,
    excluded.maternity = excluded.maternity,
    excluded.paternal.sibships = excluded.paternal.sibships,
    excluded.maternity.sibships = excluded.maternity.sibships
  )
  unsupported <- c(
    if (sibship.prior != 0) "sibship.prior",
    if (known.allele.freq != 0) "known.allele.freq",
    names(known.counts)[is.na(known.counts) | known.counts != 0]
  )
  if (length(unsupported) > 0) {
    stop(error(
      "  Not supported by gl2colony (COLONY needs extra data for them); set",
      "to 0:", paste(unsupported, collapse = ", "), "\n"
    ))
  }
  
  # a single value without '@' would be read by COLONY as the value of the
  # first locus only
  at.all.loci <- function(s) {
    s <- trimws(s)
    if (!grepl("@", s) && !grepl("[[:space:]]", s)) paste0(s, "@") else s
  }
  marker.id <- at.all.loci(marker.id)
  marker.type <- at.all.loci(marker.type)
  allelic.dropout <- at.all.loci(allelic.dropout)
  other.typ.err <- at.all.loci(other.typ.err)
  
  # SET RANDOM SEED
  if (is.null(seed)) {
    seed <- sample.int(65535, 1)
  }
  if (verbose >= 2) {
    cat(report(sprintf("  Random seed set to %d\n", seed)))
  }
  
  # ROLE COLUMNS: add only the missing ones, matching names ignoring case
  if (is.null(x@other$ind.metrics)) {
    x@other$ind.metrics <- data.frame(id = indNames(x))
  }
  role.defaults <- c(offspring = "yes", mother = "no", father = "no")
  missing.roles <- setdiff(names(role.defaults),
                           tolower(colnames(x@other$ind.metrics)))
  for (role in missing.roles) {
    x@other$ind.metrics[[role]] <- role.defaults[[role]]
  }
  if (length(missing.roles) > 0 && verbose >= 1) {
    cat(warn(
      "  Warning: column(s)", paste(missing.roles, collapse = ", "),
      "not found in ind.metrics; set to",
      paste0(missing.roles, " = '", role.defaults[missing.roles], "'",
             collapse = ", "), "for all individuals\n"
    ))
  }
  
  # EXTRACT PARENTAL IDS
  ids <- utils.colony.parental.ids(x)
  offspring.ids <- ids$offs
  dad.ids       <- ids$dad
  mum.ids       <- ids$mum
  
  # COUNTS
  n.offspring <- length(offspring.ids)
  n.dads      <- length(dad.ids)
  n.mums      <- length(mum.ids)
  loci        <- nLoc(x)
  n.total     <- n.offspring + n.dads + n.mums
  
  if (verbose >= 2) {
    cat(report(sprintf("  %d offspring, %d fathers, %d mothers detected.\n",
                       n.offspring, n.dads, n.mums)))
  }
  
  # WARN IF OFFSPRING MISSING
  if (n.offspring == 0) {
    stop(error('No offspring IDs found in metadata.'))
  }
  
  # CONVERT TO STRUCTURE FORMAT
  if (verbose >= 2) {
    cat(report("  Exporting genlight object to COLONY2 format\n"))
  }
  struct.mat <- utils.colony.genotypes(x)
  
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
  cat(project.name, '\n', file = outfilespec)
  cat(output.name, '\n', file = outfilespec, append = TRUE)
  for (i in seq_along(head.values)) {
    cat(head.values[[i]], '\t', head.comments[i], '\n',
        file = outfilespec, append = TRUE)
  }
  
  # WRITE OFFSPRING
  write.table(
    offspring.gen,
    file = outfilespec,
    append = TRUE,
    quote = FALSE,
    col.names = FALSE
  )
  
  # WRITE CANDIDATE PROBABILITIES
  cat('\n', file = outfilespec, append = TRUE)
  cat(
    paste(probability.father, probability.mother),
    '\t',
    '! Parental inclusion probabilities',
    '\n',
    file = outfilespec, append = TRUE
  )
  cat(paste(n.dads, n.mums), '\t', '! Number of candidates', '\n',
      file = outfilespec, append = TRUE)
  cat('\n', file = outfilespec, append = TRUE)
  
  # WRITE DADS
  if (n.dads > 0) {
    if (verbose >= 2) cat(report("  Writing paternal genotypes\n"))
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
    if (verbose >= 2) cat(report("  Writing maternal genotypes\n"))
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
  cat('\n', file = outfilespec, append = TRUE)
  for (i in seq_along(last.values)) {
    cat(last.values[[i]], '\t', last.comments[i], '\n',
        file = outfilespec, append = TRUE)
  }
  
  if (verbose >= 3) {
    cat(report(paste(
      "Records written to", outfilespec, "\n"
    )))
  }
  
  # FLAG SCRIPT END
  if (verbose >= 1) {
    cat(report("Completed:", funname, "\n"))
  }
  
  return(invisible(outfilespec))
}


###################### Define function utils.colony.parental.ids ##################
## This function extracts parental information in a list of 3 elements (vectors
## with offspring, dads and mums IDs, respectively). IDs are taken from
## indNames(x), not from an ind.metrics id column, so they always match the
## genotype matrix row names.
utils.colony.parental.ids <- function(gen.data) {
  # Read metadata and convert column names to lowercase
  indv.metadata <- gen.data@other$ind.metrics
  names(indv.metadata) <- tolower(names(indv.metadata))
  
  # TRUE where a role column says "yes", ignoring case and surrounding spaces
  is.yes <- function(v) {
    v <- trimws(tolower(as.character(v)))
    !is.na(v) & v == "yes"
  }
  
  ind.ids <- indNames(gen.data)
  mum.ids  <- ind.ids[is.yes(indv.metadata$mother)]
  dad.ids  <- ind.ids[is.yes(indv.metadata$father)]
  offs.ids <- ind.ids[is.yes(indv.metadata$offspring)]
  
  # Make a list with the 3 vectors
  x = list(offs = offs.ids, dad = dad.ids, mum = mum.ids)
  return(x)
}
################################################################################


###################### Define function utils.colony.genotypes #####################
## This function converts gl matrix to Structure format and from 2-row-per-ind
## to 1-row-per-ind. Not dartR.base::gl2structure, which has a different
## signature and output.
utils.colony.genotypes <- function(x,
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
