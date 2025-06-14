#' @name gl.assign.on.genotype
#' @title Use genotype to identify populations as possible source populations for an 
#' individual of unknown provenance.
#' @description
#' This script identifies populations
#' for which the unknown individual has a reasonable expectation of having been drawn
#' from those populations given its genotype and the allele frequencies in the 
#' putative source populations. The putative source populations that survive
#' are retained and returned in a genlight object.
#'
#' The algorithm computes the log-likelihood of the focal genotype under Hardy-Weinberg (HWE), then computes
#' a Z-score and one-tailed p-value by comparing the unknown individual’s log-likelihood to those from 
#' individuals in each putative source population. Significant departures from expectation renders 
#' a population unlikely to be the source for the focal unknown individual.
#' 
#' A suitable estimate of the expectation for the log likelihoods requires that the sample size is
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
#' @param aic.threshold The critical value used to select populations for which their is considered
#' some support as a putative source based on AIC weights [default 0.05]
#' @param n.best If given a value, dictates the best n=n.best populations to
#' retain for consideration (or more if their are ties) based on AIC weight. If not
#' specified, then the putative source populations identified as possibilities (AIC.wt >= aic.threshold)
#' are retained. [default NULL].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#'  [default 2 or as specified using gl.set.verbosity].
#'  
#' @return A genlight object containing the focal individual (assigned to
#' population 'unknown') and putative source populations based on AIC weights
#' If no such populations, the genlight object contains only data
#' for the unknown individual with a warning.
#'
#' @export
#'
#' @author Script: Arthur Georges. Custodian: Arthur Georges -- Post to
#'   \url{https://groups.google.com/d/forum/dartr}
#'   
#' @examples
#' # Test run with a focal individual from the Macleay River (EmmacMaclGeor)
#' test <- gl.assign.on.genotype(testset.gl,unknown='UC_00146',nmin=10,verbose=3)
#'
#' @seealso \code{\link{gl.assign.pca}}, \code{\link{gl.assign.pa}}, \code{\link{gl.assign.mahalanobis}}

gl.assign.on.genotype <- function(x,
                         unknown,
                         nmin = 10,
                         n.best = NULL,
                         aic.threshold=0.05,
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
      if(verbose >=1){cat(
            warn(
                "  Warning: the minimum size of the target population must be 
                greater than zero, set to 10\n"
            )
        )}
        nmin <- 10
    }

    if (!is.null(n.best)) {
        if (n.best < 1) {
            if(verbose >=1){cat(
                warn(
                    "  Warning: the n.best parameter for retention of best 
                    match populations must be a positive integer, set to NULL\n"
                )
            )}
            n.best <- NULL
        }
    }
    
    if (aic.threshold < 0 || aic.threshold > 1) {
      if(verbose >= 1){cat(
        warn(
          "  Warning: the aic.threshold must be 
                greater than zero and less than 1, set to default 0.05\n"
        )
      )}
      aic.threshold <- 0.05
    }
    
# DO THE JOB
    
  # Set a hard recommended minimum population size
    hard.min <- 10
    if (nmin < hard.min) {
        if(verbose >= 1){cat(warn(
            "  Warning: The specified minimum sample size is less than",hard.min, 
            "individuals\n"
        ))
        cat(warn("    Risk of the allele profile of the putative source populations
                lacking sufficient representation of the populations from which they
                are drawn is quite high.\n"
        ))}
    }
    
    # Separate unknown individual from x
    unknown.ind <- gl.keep.ind(x,ind.list=unknown,verbose=0)
    pop(unknown.ind) <- "unknown"
    
    # Convert unknown to a vector
    unknown.vec <- as.vector(as.matrix(unknown.ind))
    
    # Knowns
    knowns <- gl.drop.ind(x,ind.list=unknown,verbose=0)
    if(any(popNames(knowns)=="unknowns")){
       knowns <- gl.drop.pop(x,pop.list="unknowns",verbose=0)
    }   
    
    # Remove all known populations with less than nmin individuals
    pop.keep <- levels(pop(knowns))[table(pop(knowns)) >= nmin]
    pop.toss <- levels(pop(knowns))[table(pop(knowns)) < nmin]
    if (verbose >= 3) {
      cat("  Discarding",length(pop.toss),"populations with sample size <",nmin,":\n")
      cat(paste(pop.toss, collapse = ", "), "\n")
     }
    if (length(pop.keep) == 0) {
      stop(error("Fatal Error: All target populations excluded based on minimum sample size.\n"))
    }
    knowns <- gl.keep.pop(knowns, pop.list = pop.keep, verbose = 0) 
 
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
    
    # Function to calculate log liklihood for a focal genotype vector vs a population matrix
    
    utils.gen.prob <- function(focal, popmat) {
      
      # Ensure dimensions are compatible
      stopifnot(length(focal) == ncol(popmat))
      
      # Calculate allele frequencies (0/1/2 encoding) and HWE genotype probabilities
      p <- colMeans(popmat, na.rm = TRUE) / 2  # frequency of 'B' allele
      pAA <- (1 - p)^2
      pAB <- 2 * p * (1 - p)
      pBB <- p^2
      
      # Likelihood of focal genotype under HWE
      logL_focal <- mapply(function(g, pa, pb, pc) {
        if (is.na(g)) return(NA)
        if (g == 0) return(log(pa + 1e-10))
        if (g == 1) return(log(pb + 1e-10))
        if (g == 2) return(log(pc + 1e-10))
        return(NA)
      }, focal, pAA, pAB, pBB)
      
      logL_focal_total <- sum(logL_focal, na.rm = TRUE)
      
      # Likelihoods for each individual in the population
      logLs_pop <- apply(popmat, 1, function(ind) {
        mapply(function(g, pa, pb, pc) {
          if (is.na(g)) return(NA)
          if (g == 0) return(log(pa + 1e-10))
          if (g == 1) return(log(pb + 1e-10))
          if (g == 2) return(log(pc + 1e-10))
          return(NA)
        }, ind, pAA, pAB, pBB)
      })
      
      logL_pop_totals <- colSums(logLs_pop, na.rm = TRUE)
      logL = logL_focal_total
      
      return(logL)
    }
    
    # Prepare results container
    result <- data.frame(
      population = names(matrix.list),
      logL = NA_real_,
      flag = ""
    )
    
    # Loop through each population 
    for (i in 1:length(matrix.list)) {
      
      pop.mat <- matrix.list[[i]]
      n <- nrow(pop.mat)
      
      # Calculate the stats
      
      logL <- utils.gen.prob(unknown.vec,pop.mat)
      # Store results
      result$logL[i] <- logL

    } # End loop
    
    #results$population <- names(matrix.list)
    # Calculate AIC values
    result$aic <- -2*result$logL
    result$delta.aic <- result$aic - min(result$aic)
    result$aic.wt <- exp(-0.5*result$delta.aic)/sum(exp(-0.5*result$delta.aic)) 
    
    # Apply significance
    result$flag[result$aic.wt < aic.threshold] <- "no"
    result$flag[result$aic.wt >= aic.threshold] <- "yes"

    result <- result[order(result$aic.wt,decreasing=TRUE),]
    
    # View result
    result <- result[, c("population", "logL", "aic", "delta.aic", "aic.wt", "flag")]
    names(result) <- c("population","Log Likelihood","AIC","dAIC","AIC.wt","assign")
    print(result)
    
    if(verbose >=3){
      cat("\n  Best prospect for source population is",result$population[1],"\n\n")
    }
    
    # Retain only those n.best populations, or if n.best not set, the non-significant
    # populations
    
    if(!is.null(n.best)){
      pop.keep <- result$pop[1:n.best]
    } else {
      pop.keep <- result$pop[result$assign == "yes"]
    }
    
    gl.out <- gl.keep.pop(knowns,pop.list=pop.keep,verbose=0)
    gl.out <- gl.join(gl.out,unknown.ind,method="end2end",verbose=0)
    
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
