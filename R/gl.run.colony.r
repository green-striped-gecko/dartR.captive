#' @name gl.run.colony
#' @title Run COLONY2
#' @description
#' An R wrapper for the COLONY pedigree-inference software (Jones & Wang
#' 2010), to run full-pedigree likelihood analyses of SNP genotypes from R.
#' The function writes the COLONY input file with \code{gl2colony}, runs the
#' COLONY executable and reads the best configuration back into R.
#'
#' @inheritParams gl2colony
#' @param colony.path Path to the folder that contains the COLONY executable:
#' Colony2p.exe (Windows), colony2s.ifort.out (Linux) or colony2s.out (macOS)
#' [default getwd()].
#' @param outpath Folder for the COLONY input file and for all the output
#' files COLONY writes [default global working directory or if not
#' specified, tempdir()].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#'  progress log, including the COLONY console output; 3, progress and
#'  results summary; 5, full report [default 2, unless specified using
#'  gl.set.verbosity].
#'
#' @details
#' COLONY implements a full-pedigree likelihood method that simultaneously
#' infers sibships and parentage by considering the likelihood of entire
#' pedigree configurations rather than pairwise comparisons.
#'
#' COLONY writes its output files (named after output.name) into outpath.
#'
#' COLONY truncates individual IDs to 20 characters and splits them at
#' spaces. If any name is longer than 20 characters or contains whitespace,
#' COLONY runs on the IDs ind1, ind2, ... instead; the original names are
#' restored in best.config, and the file output.name.IDmap.csv in outpath
#' maps the IDs used in COLONY's own output files to the original names.
#' COLONY does not return an error status when it rejects its input, so the
#' function checks for the COLONY error file and for a new BestConfig file
#' and stops with COLONY's message if the run failed.
#'
#' @return Invisibly, a list with two elements:
#' \itemize{
#' \item files -- full paths of the output files written by this run
#' (including output.name.IDmap.csv when IDs were replaced).
#' \item best.config -- the best configuration (output.name.BestConfig) as
#' a data frame with columns OffspringID, FatherID, MotherID, CloneIndex and
#' ClusterIndex. Parents not among the candidates are given COLONY's
#' inferred labels (e.g. *1 for fathers, #1 for mothers).
#' }
#'
#' @author
#' Author(s): Jesús Castrejón-Figueroa, Diana A. Robledo-Ruiz, Luis Mijangos.
#' Custodian: Luis Mijangos -- Post to
#' \url{https://groups.google.com/d/forum/dartr}
#'
#' @examples
#' \dontrun{
#' # requires the COLONY executable
#' t1 <- testset.gl[1:30, 1:50]
#' t1@other$ind.metrics$offspring <- rep(c("yes", "no", "no"), each = 10)
#' t1@other$ind.metrics$father <- rep(c("no", "yes", "no"), each = 10)
#' t1@other$ind.metrics$mother <- rep(c("no", "no", "yes"), each = 10)
#' res <- gl.run.colony(t1, colony.path = "path/to/colony", seed = 1234)
#' head(res$best.config)
#' }
#'
#' @references
#' Jones, O. R., & Wang, J. (2010). COLONY: a program for parentage and
#' sibship inference from multilocus genotype data. Molecular Ecology
#' Resources, 10(3), 551-555.
#'
#' @family captive management
#' @export

gl.run.colony <- function(x,
                          colony.path = getwd(),
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
  
  # FLAG SCRIPT START
  funname <- match.call()[[1]]
  utils.flag.start(func = funname, verbose = verbose)
  
  # CHECK DATATYPE
  datatype <- utils.check.datatype(x, verbose = verbose)
  
  # FUNCTION SPECIFIC ERROR CHECKING
  os <- Sys.info()[["sysname"]]
  exe.name <- switch(os,
                     Windows = "Colony2p.exe",
                     Linux = "colony2s.ifort.out",
                     Darwin = "colony2s.out",
                     NA_character_)
  if (is.na(exe.name)) {
    stop(error("  No COLONY executable is known for", os, "\n"))
  }
  exe <- file.path(path.expand(colony.path), exe.name)
  if (!file.exists(exe)) {
    stop(error(
      "  COLONY executable not found:", exe,
      "\n  Set colony.path to the folder that contains", exe.name, "\n"
    ))
  }
  
  # DO THE JOB
  # COLONY truncates IDs to 20 characters and splits them at whitespace, so
  # such names are replaced by short IDs for the run and restored afterwards
  original.names <- indNames(x)
  rename.ids <- any(nchar(original.names) > 20 |
                      grepl("[[:space:]]", original.names))
  if (rename.ids) {
    colony.ids <- paste0("ind", seq_along(original.names))
    indNames(x) <- colony.ids
    if (verbose >= 2) {
      cat(report(
        "  Some individual names are longer than 20 characters or contain",
        paste0("spaces; COLONY runs on IDs ind1 to ind", length(colony.ids)),
        "and the names are restored in best.config\n"
      ))
    }
  }
  
  # gl2colony resolves outpath (gl.check.wd) and returns the input file path
  outfilespec <- gl2colony(
    x =  x,
    outfile =  outfile,
    outpath =  outpath,
    project.name =  project.name,
    output.name =  output.name,
    probability.father =  probability.father,
    probability.mother =  probability.mother,
    seed =  seed,
    update.allele.freq =  update.allele.freq,
    di.mono.ecious =  di.mono.ecious,
    inbreed =  inbreed,
    haplodiploid =  haplodiploid,
    polygamy.male =  polygamy.male,
    polygamy.female =  polygamy.female,
    clone.inference =  clone.inference,
    scale.shibship =  scale.shibship,
    sibship.prior =  sibship.prior,
    known.allele.freq =  known.allele.freq,
    num.runs =  num.runs,
    length.run =  length.run,
    monitor.method =  monitor.method,
    monitor.interval =  monitor.interval,
    windows.gui =  windows.gui,
    likelihood =  likelihood,
    precision.fl =  precision.fl,
    marker.id =  marker.id,
    marker.type =  marker.type,
    allelic.dropout =  allelic.dropout,
    other.typ.err =  other.typ.err,
    paternity.exclusion.threshold =  paternity.exclusion.threshold,
    maternity.exclusion.threshold =  maternity.exclusion.threshold,
    paternal.sibship =  paternal.sibship,
    maternal.sibship =  maternal.sibship,
    excluded.paternity =  excluded.paternity,
    excluded.maternity =  excluded.maternity,
    excluded.paternal.sibships =  excluded.paternal.sibships,
    excluded.maternity.sibships =  excluded.maternity.sibships,
    verbose =  verbose
  )
  
  outpath <- dirname(outfilespec)
  
  # COLONY writes its output files into its working directory
  old.wd <- setwd(outpath)
  on.exit(setwd(old.wd), add = TRUE)
  
  # COLONY exits with status 0 when it rejects its input, so success is
  # judged by the files this run writes, not by the exit status alone
  # (OP 8821)
  start.time <- trunc(Sys.time(), "secs")
  colony.console <- if (verbose >= 2) "" else FALSE
  # system2() quotes the command itself; arguments must be quoted here
  status <- system2(exe,
                    args = shQuote(paste0("IFN:", outfilespec)),
                    stdout = colony.console,
                    stderr = colony.console)
  
  if (is.na(status) || status != 0) {
    stop(error(paste0(
      "  COLONY failed to run (exit status ", status, "): ", exe, "\n"
    )))
  }
  
  is.new <- function(f) file.exists(f) & file.mtime(f) >= start.time
  error.file <- file.path(outpath, "Colony2.ErrorMessage")
  if (is.new(error.file)) {
    stop(error(
      "  COLONY stopped with an error:\n",
      paste(readLines(error.file, warn = FALSE), collapse = "\n"), "\n"
    ))
  }
  best.config.file <- file.path(outpath, paste0(output.name, ".BestConfig"))
  if (!is.new(best.config.file)) {
    stop(error(
      "  COLONY finished without writing", best.config.file,
      "\n  Check the COLONY console output (verbose >= 2) for the cause.\n"
    ))
  }
  
  out.files <- list.files(outpath, full.names = TRUE)
  out.files <- out.files[startsWith(basename(out.files),
                                    paste0(output.name, ".")) &
                           is.new(out.files)]
  # '#' marks COLONY's inferred mothers, so comments must be off
  best.config <- utils::read.table(best.config.file, header = TRUE,
                                   comment.char = "",
                                   stringsAsFactors = FALSE)
  
  if (rename.ids) {
    # inferred parents (*1, #1) are not in the map and keep their labels
    restore <- function(v) {
      i <- match(v, colony.ids)
      ifelse(is.na(i), v, original.names[i])
    }
    for (col in c("OffspringID", "FatherID", "MotherID")) {
      best.config[[col]] <- restore(best.config[[col]])
    }
    # COLONY's own output files keep the short IDs; the map translates them
    map.file <- file.path(outpath, paste0(output.name, ".IDmap.csv"))
    utils::write.csv(data.frame(colony.id = colony.ids,
                                name = original.names),
                     map.file, row.names = FALSE)
    out.files <- c(out.files, map.file)
  }
  
  if (verbose >= 3) {
    cat(report(paste0(
      "  ", length(out.files), " COLONY output files written to ", outpath,
      "\n"
    )))
  }
  
  # FLAG SCRIPT END
  if (verbose >= 1) {
    cat(report("Completed:", funname, "\n"))
  }
  
    return(invisible(list(files = out.files, best.config = best.config)))
  
}