#' @name gl.run.colony
#' @title Run COLONY2
#' @description
#' A convenient R wrapper for the COLONY pedigree‐inference software 
#' (Jones & Wang 2010), allowing users to perform full‐pedigree likelihood 
#' analyses of multilocus genotype data directly from R. This function 
#' automates the creation of the required `Colony2.DAT` input file and runs 
#' the COLONY executable.
#'
#' @param x A \code{genlight} object with individual metadata columns
#' 'offspring', 'mother', and 'father' indicating 'yes'/'no' for each sample
#' [required].
#' @param colony.path Path to the colony executable [default getwd()].
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
#' COLONY implements a Bayesian full‐pedigree likelihood method that 
#' simultaneously infers sibships and parentage by considering the likelihood 
#' of entire pedigree configurations rather than pairwise comparisons. 
#' @return
#' Invisibly returns the output filename.
#'
#' @author
#' Jesús Castrejón-Figueroa, Diana A. Robledo-Ruiz & Luis Mijangos-- Post to
#' \url{https://groups.google.com/d/forum/dartr}
#'
#' @examples
#' # gl2colony(x = platypus.gl)
#'
#' @references
#' Wang, J. (2011). COLONY: a program for parentage and sibship inference
#' from multilocus genotype data. Molecular Ecology Resources 10: 551–555.
#'
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
  
 gl2colony(
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
  
  os <- Sys.info()["sysname"]
  
  if (Sys.info()["sysname"] == "Windows") {
    system(paste0("./Colony2p.exe IFN:",outfilespec))
  }
  
  if (Sys.info()["sysname"] == "Linux") {
    system(paste0("./colony2s.ifort.out IFN:",outfilespec))
  }
  
  if (Sys.info()["sysname"] == "Darwin") {
    system(paste0("./colony2s.out IFN:",outfilespec))
  }
  
  
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