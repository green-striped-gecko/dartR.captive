#' @name gl.assign.pa
#' @title Use private alleles to identify populations as possible source populations for an 
#' individual of unknown provenance.
#' @description
#' This script identifies as putative source populations,
#' those for which the individual has an expected number of private alleles. The
#' putative source populations are retained and returned in a genlight object.
#'
#' The algorithm calculates an expectation based on the number of private alleles each individual
#' in the putative source population has in comparison with the other members of that population. 
#' From the distribution of these values, an expectation
#' is established as a mean and standard deviation. The private alleles possessed by the unknown
#' individual in comparison with the putative source population is compared to this expectation.
#' Significant departures from expectation renders a population unlikely to be the source 
#' for the focal unknown individual.
#' 
#' An excessive count of private alleles is an indication that the unknown does
#' not belong to a target population (provided that the sample size is
#' adequate, say >=10).
#'  
#' WARNING: If a putative population is not in Hardy-Weinberg equilibrium, as might occur if it
#' includes F1 hybrids and backcrosses, then the standard deviation for the expectation will
#' be inflated. This inflation may result in false identification of the population
#' as a putative source for the focal unknown individual. For this reason, you may wish to 
#' remove populations that contain individuals likely
#' to be subject to contemporary hybridization or admixture.
#'
#' @param x Name of the input genlight object [required].
#' @param unknown SpecimenID label (indName) of the focal individual whose
#' provenance is unknown [required].
#' @param nmin Minimum sample size for a target population to be included in the
#' analysis [default 10].
#' @param alpha The critical value used to select populations for which the unknown individual
#' has a count of private alleles within expectation [default 0.001]
#' @param n.best If given a value, dictates the best n=n.best populations to
#' retain for consideration (or more if their are ties) based on private alleles. If not
#' specified, then the putative source populations identified as significant (p < alpha)
#' are retained. [default NULL].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#'  [default 2 or as specified using gl.set.verbosity].
#' @return A genlight object containing the focal individual (assigned to
#' population 'unknown') and populations for which the focal individual is not
#' distinctive. If no such populations, the genlight object contains only data
#' for the unknown individual with a warning.
#'
#' @importFrom stats dnorm qnorm
#' @export
#'
#' @author Script: Arthur Georges. Custodian: Arthur Georges -- Post to
#'   \url{https://groups.google.com/d/forum/dartr}
#'   
#' @examples
#' # Test run with a focal individual from the Macleay River (EmmacMaclGeor)
#' test <- gl.assign.pa(testset.gl,unknown='UC_00146',nmin=10,verbose=3)
#'
#' @seealso \code{\link{gl.assign.pca}}

gl.assign.pa <- function(x,
                         unknown,
                         nmin = 10,
                         #                         threshold = 0,
                         n.best = NULL,
                         alpha=0.01,
                         verbose = NULL) {
  
  # SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)
  
  # FLAG SCRIPT START
  funname <- match.call()[[1]]
  utils.flag.start(func = funname,
                   build = "v2025-05-29",
                   verbose = verbose)
  
  # CHECK DATATYPE
  datatype <- utils.check.datatype(x, verbose = verbose)
  
  if (!is(x, "dartR")) {
    class(x) <- "dartR"  
    if (verbose>2) {
      cat(warn("Warning: Standard adegenet genlight object encountered. Converted to compatible dartR genlight object\n"))
      cat(warn("                    Should you wish to convert it back to an adegenet genlight object for later use outside dartR, 
                 please use function dartR2gl\n"))
    }
  }
  
  # FUNCTION SPECIFIC ERROR CHECKING
  
  if (any(duplicated(indNames(x)))) {
    stop(error("Fatal Error: Duplicate individual names in genlight object.\n"))
  }
  
  if (any(is.na(indNames(x)))) {
    stop(error("Fatal Error: NA found in individual names.\n"))
  }
  
  test <- unknown %in% indNames(x)
  if (!all(test, na.rm = FALSE)) {
    stop(
      error(
        "Fatal Error: nominated focal individual (of unknown provenance) 
                is not present in the dataset!\n"
      )
    )
  }
  
  if (is.null(pop(x))) {
    stop(error("Fatal Error: Population assignments (pop(x)) are NULL.\n"))
  }
  
  if (any(is.na(pop(x)))) {
    stop(error("Fatal Error: NA values found in population assignments.\n"))
  }
  
  if (nmin <= 0) {
    if(verbose >= 1){cat(
      warn(
        "  Warning: the minimum size of the target population must be 
                greater than zero, set to 10\n"
      )
    )}
    nmin <- 10
  }
  
  if (!is.null(n.best)) {
    if (n.best < 1) {
      if(verbose >= 1){cat(
        warn(
          "  Warning: the n.best parameter for retention of best 
                    match populations must be a positive integer, set to NULL\n"
        )
      )}
      n.best <- NULL
    }
  }
  
  # DO THE JOB
  
  # Set a hard recommended minimum population size
  hard.min <- 10
  if (nmin < hard.min) {
    if(verbose >= 1){cat(warn(
      "  Warning: The specified minimum sample size is less than",hard.min, 
      "individuals\n"
    ))
      cat(warn(
        
        "    Risk of alleles present in the unknown being missed during
                sampling of populations with sample sizes less than 10\n"
        
      ))}
  }
  
  # Separate unknown individual from x
  unknown.ind <- gl.keep.ind(x,ind.list=unknown,verbose=0)
  pop(unknown.ind) <- "unknown"
  knowns <- gl.drop.ind(x,ind.list=unknown,verbose=0)
  if(any(popNames(knowns)=="unknowns")){
    knowns <- gl.drop.pop(x,pop.list="unknowns",verbose=0)
  }   
  
  # Remove all known populations with less than nmin individuals
  pop.keep <- levels(pop(knowns))[table(pop(knowns)) >= nmin]
  pop.toss <- levels(pop(knowns))[table(pop(knowns)) < nmin]
  if (verbose >= 2) {
    cat("  Discarding",length(pop.toss),"populations with sample size <",nmin,":\n")
    if (verbose >= 3) {
      cat(paste(pop.toss, collapse = ", "), "\n")
      knowns <- gl.keep.pop(knowns, pop.list = pop.keep, verbose = 0)      }
  }
  if (length(pop.keep) == 0) {
    stop(error("Fatal Error: All target populations excluded based on minimum sample size.\n"))
  }
  
  # Remove loci scored as NA for the unknown
  # Fuck, it changed the locus names, replaced hyphens with periods, use check.names=FALSE
  b <- data.frame(as.matrix(unknown.ind),check.names=FALSE)  
  c <- names(b)[is.na(b)]
  if (length(c) > 0) {
    knowns <- gl.drop.loc(knowns, loc.list = c, verbose = 0)
    unknown.ind <- gl.drop.loc(unknown.ind, loc.list=c,verbose=0)
  }
  
  # Split the genlight object into a list of matricies
  
  # Split into a list of genlight objects by population
  pop.list <- seppop(knowns)
  
  # Convert each to a matrix
  matrix.list <- lapply(pop.list, function(g) as.matrix(g))
  
  # Name the list elements by population
  names(matrix.list) <- names(pop.list)
  
  # Now we have each population represented by a matrix in the list
  # matrix.list plus an additional population containing the unknown.
  
  # Make a function that calculates the mean and standard deviation
  # of the number of private alleles for each individual in a population.
  # Note that the distribution of private alleles can be approximately
  # represented by a Poisson distribution (or a Negative Bionomial) whereby
  # a log10 transformation will render it approximately normal.
  
  # Function to count private alleles in a focal genotype vs a population matrix
  
  count.pa <- function(focal, genmat) {
    # Ensure focal is a vector
    focal <- as.vector(focal)
    # Identify loci where the focal individual has a non-zero genotype
    nonzero <- which(focal != 0 & !is.na(focal))
    # At each of those loci, check if all others are zero or NA
    private <- sapply(nonzero, function(i) {
      all(genmat[, i] == 0 | is.na(genmat[, i]))
    })
    return(sum(private))
  }
  
  # Prepare results container
  result <- data.frame(
    population = names(matrix.list),
    logmean = NA_real_,
    logsd = NA_real_,
    z = NA_real_,
    p = NA_real_,
    flag = ""
  )
  
  # Loop through each population 
  for (i in 1:length(matrix.list)) {
    
    pop.mat <- matrix.list[[i]]
    n <- nrow(pop.mat)
    
    # Calculate private allele counts for all individuals in the population i against the others
    n.pa <- sapply(1:n, function(i) count.pa(pop.mat[i, ], pop.mat[-i, , drop = FALSE]))
    log.n.pa <- log10(n.pa + 1)
    
    # Mean and SD
    mu <- mean(log.n.pa)
    sigma <- sd(log.n.pa)
    
    # Store results
    result$logmean[i] <- mu
    result$logsd[i] <- sigma
    
    # Calculate private allele counts for the focal individual against the others in each population i
    n.pa <- numeric(n)
    #names(n.pa) <- names(matrix.list[[i]])
    
    result$counts[i] <- count.pa(focal = as.matrix(unknown.ind), genmat = pop.mat)
    result$counts.log[i] <- log10(result$counts[i] +1)
    # Z-score and upper-tail probability of unknown count
    z.score <- (result$counts.log[i] - mu) / sigma
    p.value <- 1 - pnorm(z.score)
    result$z[i] <- z.score
    result$p[i] <- round(p.value,6)
  }
  
  # Apply significance
  result$flag[result$p < alpha] <- "no"
  result$flag[result$p >= alpha] <- "yes"
  
  result <- result[order(result$p,decreasing=TRUE),]
  
  # View result
  result <- result[, c("population", "counts", "z", "p", "flag")]
  names(result) <- c("pop","count","Z-score","p-value","assign")
  print(result)
  
  # Retain only those n.best populations, or if n.best not set, the non-significant
  # populations
  
  if(!is.null(n.best)){
    pop.keep <- result$pop[1:n.best]
  } else {
    pop.keep <- result$pop[result$assign == "yes"]
  }
  
  gl.out <- gl.keep.pop(knowns,pop.list=pop.keep,verbose=0)
  gl.out <- gl.join(gl.out,unknown.ind,method="join.by.loc",verbose=0)
  
  gl.out <- gl.filter.monomorphs(gl.out,verbose=0)
  
  if (nInd(gl.out) == 1) {
    cat(warn("Warning: Final genlight object contains only the unknown individual, no populations assigned.\n"))
  }
  
  # FLAG SCRIPT END
  
  if (verbose > 0) {
    cat(report("Completed:", funname, "\n"))
  }
  
  return(gl.out)
}
