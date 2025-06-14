#' @name gl.assign.mahalanobis
#' @title Calculate probabilities of assignment of an individual of unknown provenance 
#' to population based on Mahalanobis Distance
#' @description
#' This script assigns an individual of unknown provenance to one or more target
#' populations based on the unknown individual's proximity to population 
#' centroids; proximity is estimated using Mahalanobis Distance and a z score and
#' probability of assignment is calculated.
#'
#' The following process is followed:
#' \enumerate{
#' \item An ordination is undertaken on the populations to again yield a
#' series of orthogonal (independent) axes.
#' \item A workable subset of dimensions is chosen, that specified as dim.limit, or
#' the number of dimensions with substantive eigenvalues (Kaiser-Guttman criterion),
#' whichever is the smaller.
#' \item The Mahalobalis Distance is calculated for the unknown against each
#' population and probability of membership of each population is calculated.
#' The assignment probabilities are listed in support of a decision.
#' }
#' @details
#' There are three considerations to assignment. First, consider only those
#' populations for which the unknown has no private alleles. Private alleles are
#' an indication that the unknown does not belong to a target population
#' (provided that the sample size is adequate, say >=10). This can be evaluated
#'  with gl.assign.pa().
#'
#' A next step is to consider the PCoA plot for populations remaining after step 1.
#' The position of the unknown in relation to the
#' confidence ellipses is plotted by this script as a basis for narrowing down
#' the list of putative source populations. This can be evaluated with 
#' gl.assign.pca().
#' 
#' The third step (delivered by this script) is to consider the assignment 
#' probabilities based on the squared Generalised Linear Distance 
#' (Mahalanobis distance) of the unknown from the centroid for each population, 
#' then to consider the probability associated with its quantile using the 
#' Chisquare approximation. In effect, this index takes into account position 
#' of the unknown in relation to the confidence envelope in all selected 
#' dimensions of the ordination. The larger the assignment probability, 
#' the greater the confidence in the assignment. 
#' 
#' If dim.limit is set to 2, to correspond with the dimensions used in
#' gl.assign.pa(), then the output provides a ranking of the set
#' of putative source populations selected after the PCoA selection step.
#' 
#' If dim.limit is set to be > 2, then this script provides a basis for
#' further narrowing the set of putative populations.If the unknown individual
#' is an extreme outlier, say at less than 0.001 probability of population 
#' membership (0.999 confidence envelope), then the associated population 
#' can be eliminated from further consideration.
#' 
#' Warning: gl.assign.mahalanobis() treats each specified dimension equally, without
#' regard to the percentage variation explained after ordination. If the 
#' unknown is an outlier in a lower dimension with an explanatory variance of,
#' say, 0.1%, the putative population will be eliminated. This is why the script
#' only uses substantive dimensions from the ordination.
#'
#' Each of these above approaches provides evidence, none are 100% definitive. 
#' They need to be interpreted cautiously.
#' 
#' In deciding the assignment, the script considers an individual to be an
#' outlier with respect to a particular population at alpha = 0.001 as default
#'
#' @param x Name of the input genlight object [required].
#' @param unknown Identity label of the focal individual whose provenance is
#' unknown [required].
#' @param nmin Minimum sample size for a target population to be included in the
#' analysis [default 10].
#' @param dim.limit Maximum number of dimensions to consider for the
#' confidence ellipses [default nPop(x)-1]
#' @param plevel Probability level for bounding ellipses
#' [default 0.001].
#' @param n.best If given a value, dictates the best n=n.best populations to
#' retain for consideration (or more if their are ties). If not
#' specified, then the putative source populations identified as possibilities
#' by the PCA are retained. [default NULL].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default 2 or as specified using gl.set.verbosity].
#'
#' @return A data frame with the results of the assignment analysis. 
#'
#' @importFrom stats dnorm qnorm
#' @export
#'
#' @author Script: Arthur Georges. Custodian: Arthur Georges --
#' Post to \url{https://groups.google.com/d/forum/dartr}
#'
#@examples
# # Test run with a focal individual from the Macleay River (EmmacMaclGeor)
# test <- gl.assign.mahalanobis(testset.gl,unknown='UC_00146',verbose=3)
#
# @examples 
# \dontrun{
# #Test run with a focal individual from the Macleay River (EmmacMaclGeor) 
# test <- gl.assign.mahalanobis(testset.gl,unknown='UC_00146',verbose=3)
# }

gl.assign.mahalanobis <- function(x,
                                  nmin=10,
                                  dim.limit=NULL,
                                  plevel=0.001,
                                  n.best = NULL,
                                  unknown,
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
    
    # FUNCTION SPECIFIC ERROR CHECKING
    
    if(is.null(dim.limit)){
      dim.limit <- nPop(x) - 1
    }
    
    if (any(is.na(indNames(x))) || any(duplicated(indNames(x)))) {
      stop(error("Fatal Error: Individual names must be unique and non-missing.\n"))
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
            "  Warning: Value of plevel must be between 0 and 1, set to 0.001\n"
        ))
        plevel <- 0.001
    }
    
    if (is.null(pop(x)) || any(is.na(pop(x)))) {
      stop(error("Fatal Error: Population assignments are missing or incomplete.\n"))
    }
    
    if (nLoc(x) < nPop(x)) {
        stop(error(
            "Fatal Error: Number of loci less than number of populations"
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
    # 
    # # Hold the unknown
    # if("unknown" %in% popNames(x)){
    #   hold.unknown <- gl.keep.pop(x,pop.list="unknown",verbose=0)
    #   unk.flag <- TRUE
    # } else {
    #   unk.flag <- FALSE
    # }
    
    # DO THE JOB
    
    # Remove all known populations with less than nmin individuals (including the unknown)
    pop.keep <- levels(pop(x))[table(pop(x)) >= nmin]
    pop.toss <- levels(pop(x))[table(pop(x)) < nmin]
    if (verbose >= 2) {
      if(length(pop.toss != 0)){
        cat("  Discarding",length(pop.toss),"populations with sample size < nmin =",nmin,":\n")
        if (verbose >= 3) {
          cat(paste(pop.toss, collapse = ", "), "\n")
        }
      }
    }
    
    # Assign the unknown individual to its own population
    
    tmp <- as.character(pop(x))
    tmp[indNames(x)==unknown] <- "Unknown"
    pop(x) <- as.factor(tmp)
    
    # Keep populations with n >= nmin and the unknown
    
    if(length(pop.keep != 0)){
      x <- gl.keep.pop(x, pop.list = c(pop.keep,"unknown"), verbose = 0) 
    } else{
      x <- gl.keep.pop(x, pop.list = "unknown", verbose = 0) 
      cat(warn(" Warning: No populations have sufficient individuals\n"))
    }
    
    if (length(pop.keep) == 0) {
      stop(error("Fatal Error: All target populations excluded based on minimum sample size.\n"))
    }
    
    # # Add the unknown back in
    # if(unk.flag){
    #   x <- gl.join(x1=x,x2=hold.unknown,method="end2end",verbose=0)
    # }

    plevel <- 1- plevel
    
    if (nPop(x) < 2) {
        if(verbose >= 2){
            cat(warn("  Only one population, including the unknown, no putative 
                     source identified.\n"))
        }
        # Add a row
        df[i,1] <- unknown
        df[i,2] <- NA
        df[i,3] <- NA
        df[i,4] <- NA
        df[i,5] <- NA
        df[i,6] <- NA
        # FLAG SCRIPT END
        if (verbose > 0) {cat(report("Completed:", funname, "\n"))}
        return(df)
        
    } else {
    
    # Assign the unknown to population unknown    
    vec <- as.vector(pop(x))
    vec[indNames(x) == unknown] <- "unknown"
    pop(x) <- as.factor(vec)
    
    # Run the pcoa 
    if (nInd(x) < 2) {
        df <- NULL
        if(verbose >= 2){
          cat(warn("  Only one individual, no putative 
                     source identified.\n"))
        }
        return(df)
    }    
    
    # Tidy up the genlight object
      # Remove loci that are missing in the unknown
        b <- gl.keep.pop(x,pop.list="unknown", verbose=0)
        c <- locNames(b)[is.na(as.matrix(b))]
        # # Fuck, it changed the locus names, replaced hyphens with periods, use check.names=FALSE
        # b <- data.frame(as.matrix(unknown.ind),check.names=FALSE)  
        # c <- names(b)[is.na(b)]
        if (length(c) > 0) {
          x <- gl.drop.loc(x, loc.list = c, verbose = 0)
          #unknown.ind <- gl.drop.loc(unknown.ind, loc.list=c,verbose=0)
        }
        if(verbose >= 3){
          cat(report("  Rendering the data matrix dense by imputation\n"))
        }
        x <- gl.impute(x,verbose=0)
    # This is necessary to avoid difficulties with Mahalanobis Distances. Needs dense matrix.
    
    # Run the PCA
        if(verbose >= 3){
          cat(report("  Undertaking a PCA\n"))
        }
        #pcoa <- gl.pcoa(x,nfactors=dim.limit,correction='cailliez',verbose=0)
        pcoa <- gl.pcoa(x,nfactors=dim.limit,verbose=0)

        #   if(verbose >= 2){
        #     cat(report("  Plotting the unknown against putative populations\n"))
        #   }
        #     suppressWarnings(suppressMessages(gl.pcoa.plot(pcoa,
        #                                                    x,
        #                                                    ellipse=TRUE,
        #                                                    plevel=plevel,
        #                                                    verbose=0)))
        
    # Additional tests, should be OK with dense input matrix
        if (any(!is.finite(pcoa$scores))) {
          stop(error("Fatal Error: PCA ordination returned non-finite values.\n"))
        }
        
        if(any(pcoa$eig <= 0)){
          stop(error("Fatal Error: PCA ordination returned negative eigenvalues.This will likely result in negative Mahalanobis Distances\n"))
          cat("  Run again with a smaller dim.limit informed by the Broken-Stick plot\n")
        }
        
        # Determine the number of dimensions for confidence envelope (the 
        # ordination and dimension reduction) From the eigenvalue
        # distribution
        
        # FUNCTION FOR BROKEN STICK
        bs.statistics <- function(eigenvalues,
                                  gap_threshold = 2) {
          
          ### Step 1: Initialization ###
          # Number of eigenvalues and their indices
          n <- length(eigenvalues)
          index <- 1:n
          
          ### Step 2: Calculate Broken Stick Thresholds ###
          # Broken stick model for expected eigenvalue lengths
          broken.sticks <- sum(eigenvalues) * rev(cumsum(1 / n:1)) / n
          
          # Create a data frame to hold eigenvalues and their thresholds
          df <- data.frame(index = index, 
                           eigenvalues = eigenvalues, 
                           broken.sticks = broken.sticks)
          
          ### Step 3: Classify Eigenvalues as "Structured" or "Noisy" ###
          # Logical vector: TRUE if eigenvalue > broken stick threshold
          idx <- eigenvalues > broken.sticks
          idx2 <- which(idx)  # Indices of eigenvalues greater than the thresholds
          
          # Handle edge case: No eigenvalues greater than thresholds
          if (length(idx2) == 0) {
            # Return all eigenvalues as "noisy" if no structure is detected
            return(list(struct = data.frame(), noise = df))
          }
          
          # Identify gaps in structured indices
          # - A "gap" occurs if the difference between consecutive indices > gap_threshold
          gaps <- which(diff(idx2) > gap_threshold)
          
          # If gaps exist, mark all eigenvalues after the gap as "noisy"
          if (length(gaps) > 0) {
            idx[(gaps[1] + 1):length(idx)] <- FALSE
          }
          
          # Separate structured and noisy eigenvalues based on updated idx
          struc <- df[idx, ]
          noise <- df[!idx, ]
          
          # Add classification labels to each subset
          struc$structure <- 'structured'
          noise$structure <- 'noisy'
          
          return(length(struc$structure))

          }
          
        # Set the number of dimensions for confidence envelope 
        first.est <- bs.statistics(pcoa$eig)
        if(verbose >= 3){
          cat(report("  Dimensions retained:",first.est,"\n"))
        }
        
        # s <- sum(pcoa$eig)
        # e <- round(pcoa$eig * 100 / s, 1)
        # e.sign <- e[e > mean(e,na.rm=TRUE)]
        # first.est <- length(e.sign)
        # # From the number of populations, including the unknown 
        # # sec.est <- nPop(x)
        
        # cat(' Number of populations, including the unknown:',sec.est,'\n')
        if (verbose >= 2) {
            cat(
                report(
                    "  Number of dimensions with substantial eigenvalues (Broken-Stick Criterion):",
                    first.est,
                    ". Hardwired limit",
                    dim.limit,
                    "\n"
                )
            )
            cat(report("    Selecting the smallest of the two\n"))
        }
        dim <- min(first.est, dim.limit)
        if (verbose >= 2) {
            cat(report("    Dimension of confidence envelope set at", dim,
                       "\n"))
        }
        pcoa$scores <- pcoa$scores[, 1:dim]
        
        
        # Add population names to the scores
        c <-
            data.frame(cbind(pcoa$scores, as.character(pop(x))), 
                       stringsAsFactors = FALSE)
        colnames(c)[dim + 1] <- "pop"
        
        # Create a set of data without the unknown
        knowns <- c[c[, "pop"] != "unknown",]
        Unknown <- c[c[, "pop"] == "unknown",]

        # For each population
        p <- as.factor(unique(knowns[, "pop"]))
        # Create a dataframe to hold the results
        df <- data.frame(
          unknown = character(),
          pop = character(),
          MahalD = numeric(),
          pval = numeric(),
          critval = numeric(),
          assign = character(),
          stringsAsFactors = FALSE
        )
        for (i in 1:length(levels(p))) {
            # Pull out population i
           # if(verbose>=2){
           #   cat(" Considering population",levels(p)[i],"\n")
           # }
            m <- knowns[knowns[, "pop"] == levels(p)[i],]
            
            # Discard the population labels
            m <- m[, 1:dim]
            hold <- row.names(m)
            row.names(m) <- NULL
            Unknown <- Unknown[1:dim]
            row.names(Unknown) <- NULL
            
            # Convert to numeric, and reformat as a matrix
            n <- as.matrix(sapply(m, as.numeric))
            Unknown <- matrix(as.numeric(Unknown[1:dim]), nrow = 1)
            
            # Calculate statistics for the data without the unknown
            means <- colMeans(n)
            covariance <- stats::cov(n)
            
            # Add back in the unknown
            all <- rbind(n,Unknown)
            
            # Calculate Mahalanobis Distances
            error.flag <- FALSE
            D <- tryCatch({
              stats::mahalanobis(all, means, covariance, tol = 1e-20)
            }, error = function(e) {
              cat(error("Error: Mahalanobis distance failed. Covariance matrix may be singular.\n"))
              cat(  "Unable to consider population",levels(p)[i],"-- discarded\n")
              error.flag <- TRUE
            })
            if(error.flag){next()}
            
#            wtD <- WMDB::wmahalanobis(all,means,covariance,weight=e)
            names(D) <- c(hold,"unknown")
            
            # Calculate the associated probabilities
            pval <- (pchisq(D, df=length(D)-1, lower.tail=FALSE))
            # Is the result non-significant, then assign=yes
            if (pval["unknown"] >= 1-plevel) {
                assign <- "yes"
            } else {
                assign <- "no"
            }
            
            # df.known <-  data.frame(
            # unknown = unknown,
            # pop = levels(p)[i],
            # MahalD = D["unknown"],
            # pval = pval["unknown"],
            # critval = 1 - plevel,
            # assign = assign,
            # stringsAsFactors = FALSE
            # )
   
            new.entry <- data.frame(
              #unknown = unknown,
              pop = levels(p)[i],
              MahalD = D["unknown"],
              pval = pval["unknown"],
              #critval = plevel,
              assign = assign,
              stringsAsFactors = FALSE
            )
            row.names(new.entry) <- NULL
            
            df <- rbind(df, new.entry)
        }

        # Order the dataframe in descending order on pval
        df <- df[order(df$pval,decreasing=TRUE),]
        
        # Dump unwanted variables
        df <- df[c("pop","MahalD","pval","assign")]
        row.names(df) <- NULL
        
        # Display
        if (verbose >= 3) {
            cat("Assignment of unknown individual:",unknown,"\n")
            cat("Alpha level of significance:",1-plevel,"\n")
            print(df)
        }
        
        # Extract the best match, and report
        best <- as.character(df[df$assign == "yes",][1,1])
       
        if (verbose >= 3) {
            cat("  Best assignment is the population with the largest probability
                of assignment, in this case",
                    best,
                    "\n"
                )
        }

        
        # Retain only those n.best populations, or if n.best not set, the non-significant
        # populations
        
        if(!is.null(n.best)){
          pop.keep <- df$pop[1:n.best]
        } else {
          pop.keep <- df$pop[df$assign == "yes"]
        }
        
        gl.out <- gl.keep.pop(x,pop.list=c(pop.keep,"unknown"),verbose=0)
       # gl.out <- gl.join(gl.out,unknown.ind,method="end2end",verbose=0)
        
        gl.out <- gl.filter.monomorphs(gl.out,verbose=0)
        
        if (nInd(gl.out) == 1) {
          cat(warn("Warning: Final genlight object contains only the unknown individual, no populations assigned.\n"))
        }
        
    if(verbose >= 2){
      cat(report("  Returning a genlight object with the putative source populations and the unknown\n"))
    }
    
    # FLAG SCRIPT END
    
    if (verbose > 0) {
        cat(report("Completed:", funname, "\n"))
    }
  }
    
    return(gl.out)
}
