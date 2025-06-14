#' @name gl.assign.pca
#' @title Eliminate from consideration putative source populations for a specified 
#' individual of unknown provenance using PCA
#' @description
#' This script eliminates from consideration putative source populations for a specified
#' individual of unknown provinence based on its proximity to each putative source
#' population defined by a confidence ellipse in ordinated space of two dimensions.
#'
#' The following process is followed:
#' \enumerate{
#' \item The space defined by the loci is ordinated to yield a series of
#' orthogonal axes (independent) and the top two dimensions are considered.
#' Populations for which the unknown individual lies outside the specified confidence
#' limits are set aside to allow further examination.
#' }
#' @details
#' There are three considerations to assignment. First, consider only those
#' populations for which the unknown has no private alleles. Substanial numbers of
#' private alleles are an indication that the unknown does not belong to a target population
#' (provided that the sample size is adequate, say >=10). This can be evaluated
#'  with gl.assign.pa().
#'
#' A next step is to consider the PCA plot for populations where no private
#' alleles have been detected and the position of the unknown in relation to 
#' confidence ellipses as produced by this script. Note, this plot is
#' considering only the top two dimensions of the ordination. This is justified because
#' an unknown lying outside the confidence ellipse in two dimensions cannot lie
#' within the confidence envelope incorporating deeper dimensions. It can be 
#' unambiguously interpreted as it lying outside the confidence envelope. 
#' However, if the unknown lies inside the confidence ellipse in two dimensions, 
#' then it may still lie outside the confidence envelope in deeper dimensions. 
#' 
#' As with the first step using gl.assign.pa(), this second step is good for 
#' eliminating populations from consideration, but does not provide confidence 
#' in assignment.
#'
#' The third step is to consider the assignment probabilities, using the script
#' gl.assign.mahalanobis(). This approach calculates the squared Generalised 
#' Linear Distance (Mahalanobis distance) of the unknown from the centroid 
#' for each remaining putative source population, and calculates the probability 
#' associated with its quantile under the zero truncated normal distribution. 
#' This index takes into account position of the unknown in relation to the 
#' confidence envelope in all selected dimensions of the ordination. 
#'
#' Each of these approaches provides evidence, none are 100% definitive. They
#'  need to be interpreted cautiously. They are best applied sequentially.
#'  
#' In deciding the assignment, the script considers an individual to be an
#' outlier with respect to a particular population at alpha = 0.001 as default.
#'
#' @param x Name of the input genlight object [required].
#' @param unknown Identity label of the focal individual whose provenance is
#' unknown [required].
#' @param nmin Minimum sample size for a target population to be included in the
#' analysis [default 10].
#' @param plevel Probability level for bounding ellipses in the PCoA plot
#' [default 0.999].
#' @param plot.out If TRUE, plot the 2D PCA showing the position 
#' of the unknown [default TRUE]
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default 2 or as specified using gl.set.verbosity].
#'
#' @return A genlight object containing only those populations that are
#' putative source populations for the unknown individual. 
#'
#' @importFrom stats dnorm qnorm
#' @export
#'
#' @author Script: Arthur Georges. Custodian: Arthur Georges -- Post to 
#' \url{https://groups.google.com/d/forum/dartr}
# @examples 
# \dontrun{
# #Test run with a focal individual from the Macleay River (EmmacMaclGeor) 
# test <- gl.assign.pca(testset.gl, unknown='UC_00146',verbose=3) 
# }

gl.assign.pca <- function(x,
                          unknown,
                          nmin=10,
                          plevel = 0.001,
                          plot.out=TRUE,
                          verbose = NULL) {
  # SET VERBOSITY
    verbose <- gl.check.verbosity(verbose)
    if(verbose==0){
        plot.out <- FALSE
    }

    # FLAG SCRIPT START
    funname <- match.call()[[1]]
    utils.flag.start(func = funname,
                     build = "v2025-05-29",
                     verbose = verbose)
 
    # check if package is installed
    pkg <- "SIBER"
    if (!(requireNamespace(pkg, quietly = TRUE))) {
      cat(error(
        "Package",
        pkg,
        " needed for this function to work. Please install it.\n"
      ))
      return(-1)
    }

    # CHECK DATATYPE
    datatype <- utils.check.datatype(x, verbose = 0)
    
    if (!is(x, "dartR")) {
      class(x) <- "dartR"  
      if (verbose>2) {
        cat(warn("Warning: Standard adegenet genlight object encountered. Converted to compatible dartR genlight object\n"))
        cat(warn("                    Should you wish to convert it back to an adegenet genlight object for later use outside dartR, 
                 please use function dartR2gl\n"))
      }
    }
    
    if (nPop(x) < 2) {
        stop(
            error(
                "Fatal Error: Only one population, including the unknown, no 
                putative source\n"
            )
        )
    }

    # FUNCTION SPECIFIC ERROR CHECKING
    
    if (any(duplicated(indNames(x)))) {
      stop(error("Fatal Error: Duplicate individual names in genlight object.\n"))
    }
    if (any(is.na(indNames(x)))) {
      stop(error("Fatal Error: NA values found in individual names.\n"))
    }

    if (!(unknown %in% indNames(x))) {
        stop(
            error(
                "Fatal Error: Unknown must be listed among the individuals in 
                the genlight object!\n"
            )
        )
    }
    if (plevel > 1 || plevel < 0) {
        cat(warn(
            "  Warning: Value of plevel must be between 0 and 1, set to 0.95\n"
        ))
        plevel <- 0.95
    }
    
    if (is.null(pop(x))) {
      stop(error("Fatal Error: Population assignments are NULL.\n"))
    }
    if (any(is.na(pop(x)))) {
      stop(error("Fatal Error: NA values found in population assignments.\n"))
    }

    if (nLoc(x) < nPop(x)) {
        stop(error(
            "Fatal Error: Number of loci less than number of populations\n"
        ))
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
    
    # DO THE JOB
    
    plevel <- 1 - plevel
    
    # Assign the unknown individual to its own population
    
    tmp <- as.character(pop(x))
    tmp[indNames(x)==unknown] <- "Unknown"
    pop(x) <- as.factor(tmp)
    
    unknown.ind <- gl.keep.pop(x,pop.list="Unknown")
    
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
    
    # Remove all known populations with less than nmin individuals (including the unknown)
    pop.keep <- levels(pop(x))[table(pop(x)) >= nmin]
    pop.keep <- c(pop.keep,"Unknown")
    pop.toss <- levels(pop(x))[table(pop(x)) < nmin]
    pop.toss <- pop.toss[pop.toss!="Unknown"]
    
    if (verbose >= 2) {
      cat("  Discarding",length(pop.toss),"populations with sample size < nmin =",nmin,":\n")
      if (verbose >= 3) {
        cat(paste(pop.toss, collapse = ", "), "\n")
      }
    }
    x <- gl.keep.pop(x, pop.list = c(pop.keep), verbose = 0) 
    
    if (length(pop.keep) == 0) {
      stop(error("Fatal Error: All target populations excluded based on minimum sample size.\n"))
    }
    
    # Remove loci scored as NA for the unknown
    # Fuck, it changed the locus names, replaced hyphens with periods, use check.names=FALSE
    tmp1 <- data.frame(as.matrix(unknown.ind),check.names=FALSE)  
    tmp2 <- names(tmp1)[is.na(tmp1)]
    if (length(tmp2) > 0) {
      x <- gl.drop.loc(x, loc.list = tmp2, verbose = 0)
    }
    
    # Impute remaining missing values
    x <- gl.impute(x, verbose=0)
    
    # Ordinate a reduced space of 2 dimensions
    if (verbose >= 2){
        cat(report("  Calculating a PCA to represent the unknown in the context
                   of putative sources\n"))
    }    
    # Perform the PCA
    pcoa <- gl.pcoa(x, nfactors = 2, verbose = 0)
    
    # Plot
    if (verbose >= 0 && plot.out==TRUE) {
        suppressWarnings(suppressMessages(gl.pcoa.plot(
            pcoa, x, ellipse = TRUE, plevel = plevel, verbose=0
        )))  # Because the unknown pop throws an ellipse error
    }
    
    if (any(!is.finite(pcoa$scores))) {
      stop(error("Fatal Error: PCA ordination failed. Check for missing data or insufficient variation.\n"))
    }
    
    # Combine Pop names and pca scores
    df <- data.frame(pcoa$scores)
    df <- cbind(as.character(pop(x)), df)
    names(df) <- c("pop", "x", "y")
    # Determine if the unknown lies within the confidence ellipses specified by 
    #plevel
    #result <- data.frame()
    count <- 0
    
    # Remove 'Unknown' for training
    
    tmp <- gl.drop.pop(x, pop.list = "Unknown", verbose = 0)
    result <- data.frame(pop = character(), hit = logical(), stringsAsFactors = FALSE)
    
    for (i in popNames(tmp)) {
      # cat("  -> Population", i, "\n")
      # cat("     N = ", sum(df$pop == i), "\n")
      if (sum(df$pop == i) <= 1) {
        result <- rbind(result, data.frame(pop = i, hit = NA))
        next()
      }
      
      A <- pcoa$scores[df$pop == i, ]
      mu <- colMeans(A)
      sigma <- stats::cov(A)
      
      if (det(sigma) == 0 || any(is.na(sigma))) {
        warning(paste("Skipping population", i, ": singular covariance matrix"))
        result <- rbind(result, data.frame(pop = i, hit = NA))
        next()
      }
      
      testset <- rbind(pcoa$scores[df$pop == "Unknown", ], A)
      transform <- SIBER::pointsToEllipsoid(testset, sigma, mu)
      inside.or.out <- SIBER::ellipseInOut(transform, p = plevel)
      
      result <- rbind(result, data.frame(pop = i, hit = inside.or.out[1]))
      
    } #End Loop
    
    if (any(is.na(result$hit))) {
      warning("Some populations could not be evaluated due to singular covariance matrices: ",
              paste(result$pop[is.na(result$hit)], collapse = ", "))
    }
    
    names(result) <- c("pop", "hit")
    nhits <- length(result$pop[result$hit])
    nohits <- length(result$pop[!result$hit])
    if(verbose >= 2){
        cat(report("  Eliminating populations for which the unknown is outside
                   their confidence envelope\n"))
    }
    if (verbose >= 3) {
        if (nhits > 0) {
            cat("  Putative source populations:",
                paste(result[result$hit == TRUE, "pop"], collapse = ", "),
                "\n")
        } else {
            cat("  No putative source populations identified\n")
        }
        if (all(result$hit == TRUE)) {
            cat("  No populations eliminated from consideration\n")
        } else {
            cat(
                "  Populations eliminated from consideration:",
                paste(result[result$hit == FALSE, "pop"], collapse = ", "),
                "\n"
            )
        }
    }
    
    if (all(result$hit == TRUE,na.rm=TRUE)) {
        x2 <- x
    } else {
        x2 <- gl.drop.pop(x, pop.list = result[result$hit == FALSE, "pop"],verbose = 0)
    }
    
    if(verbose >= 2){
        cat(report("  Returning a genlight object with remaining putative source
                   populations plus the unknown\n"))
    }
    if (nPop(x2) <= 1) {
      cat(warn("  Warning: No putative source populations remain after filtering.\n"))
    }
    
    # FLAG SCRIPT END
    
    if (verbose > 0) {
        cat(report("Completed:", funname, "\n"))
    }
    
    return(x2)
}
