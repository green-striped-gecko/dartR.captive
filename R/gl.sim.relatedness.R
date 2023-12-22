#' @name gl.sim.relatedness
#' @title Simulate relatedness estimates.
#' @description
#' A simulation based tool to estimate different degrees of relatedness using
#' genlight object to bootstrap the results of kinship estimates. This method 
#' uses EMIBD9 (Wang, J. 2022).
#' 
#'Below is a table modified from Speed & Balding (2015) showing kinship values,
#'and their confidence intervals (CI), for different relationships that could 
#'be used to guide the choosing of the relatedness threshold in the function.
#'
#'|Relationship                               |Kinship  |     95% CI       |
#'
#'|Identical twins/clones/same individual     | 0.5     |        -         |
#'
#'|Sibling/Parent-Offspring                   | 0.25    |    (0.204, 0.296)|
#'
#'|Half-sibling                               | 0.125   |    (0.092, 0.158)|
#'
#'|First cousin                               | 0.062   |    (0.038, 0.089)|
#'
#'|Unrelated                                  | 0       |        -         | 
#' @param x Name of the genlight object containing the SNP data [required].
#' @param rel The degree of relatedness you wish to simulate. One of, 'full.sib', 
#' 'half.sib','first.cousin' [default 'full.sib']. 
#' @param nboots The number of simulation replicates you wish to perform [default 10].
#' @param emibd9.path The location of all necessary files to run EMIBD9 (read more at gl.run.EMIBD9)
#' @param conf The specified threshold for confidence interval calculation from simulated relatedness values [default 0.95]
#' @param nclusters The number of cores used to perform calculations [default 1].
#' @param plot.out A boolean that indicates whether to plot the results [default TRUE].
#' @param plot.dir Directory to save the plot RDS files [default as specified by the global working directory or tempdir()]
#' @param plot.file Name for the RDS binary file to save (base name only, exclude extension) [default NULL]
#' @param iseed An integer used to seed the random number generator [default 42].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#'  [default NULL, unless specified using gl.set.verbosity]
#' @return Summary statistics of chosen relatedness relationship and a histogram of relatedness values showing the mean. 
#' @author Custodian: Sam Amini -- Post to
#' \url{https://groups.google.com/d/forum/dartr}
#' @examples
#' \dontrun{
#' #To run this function needs EMIBD9 installed in your computer
#' }
#'
#' @references
#' \itemize{
#' \item Wang, J. (2022). A joint likelihood estimator of relatedness and allele
#'  frequencies from a small sample of individuals. Methods in Ecology and
#'  Evolution, 13(11), 2443-2462.
#'  
#'  Speed, D., Balding, D. Relatedness in the post-genomic era: is it still useful?. Nat
#'  Rev Genet 16, 33–44 (2015). 
#' }
#' @importFrom stringr str_split
#' @export

gl.sim.relatedness <- function(x,
                               rel = "full.sib", 
                               nboots = 10,
                               emibd9.path = getwd(),
                               conf = 0.95,
                               nclusters = 1,
                               ISeed = 42, 
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
  
  #Check dartR.captive is installed so we can use EMIBD9
  pkg <- "dartR.captive"
  if (!(requireNamespace(pkg, quietly = TRUE))) {
    cat(error(
      "Package",
      pkg,
      " needed for this function to work. Please install it.\n"
    ))
    return(-1)
  }
  
  #Do the job
  
  if (rel == "full.sib") {
    full.sib <- function(x) {
    parents <- sample(1:nInd(x), 2, replace = T)
    
    mother <- x[parents[1],]
    father <- x[parents[2],]
    
    off <- gl.sim.offspring(father, mother, noffpermother = 2, sexratio = 0.5)
    ppoff <- rbind(x,off[1,])
    
    res <- gl.run.EMIBD9(ppoff, 
                         emibd9.path = emibd9.path)
    
    
    relo <- (res$rel[parents[1], nInd(ppoff)] + res$rel[parents[2], nInd(ppoff)])/2
    return(relo)
    }
    } else if (rel == "half.sib") {
      half.sib <- function(x) {
        parents <- sample(1:nInd(x),3, replace = F)
        mother1 <- x[parents[1],]
        mother2 <- x[parents[2],]
        
        father <- x[parents[3],]
        
        off1 <- gl.sim.offspring(father, mother1, noffpermother = 2, sexratio = 1)
        off2 <- gl.sim.offspring(father, mother2, noffpermother = 2, sexratio = 1)

        ppoff <- rbind(x,off1[1,],off2[2,])
        
        res <- gl.run.EMIBD9(ppoff, 
                             emibd9.path = emibd9.path)
        
  
        
        relo <- (res$rel["Po_1","Po_2"])
        return(relo)
      }
      } else if (rel == "first.cousin") {
        first.cousin <- function(x) {
          parents1 <- sample(1:nInd(x),2, replace = F)
          mother <- x[parents1[1],]
          father <- x[parents1[2],]
          
          off1 <- gl.sim.offspring(father, mother, noffpermother = 2, sexratio = 0.5)
        
          parents2 <- sample(1:nInd(x),2, replace = F)
          mother2 <- x[parents2[1],]
          father2 <- x[parents2[2],]
          
          cousin1 <- gl.sim.offspring(off1[1,], mother2, noffpermother = 2, sexratio = 0.5)
          cousin1 <- cousin1[1,]
          indNames(cousin1) <- ("cousin1")
          
          cousin2 <- gl.sim.offspring(father2, off1[2,], noffpermother = 2, sexratio = 0.5)
          cousin2 <- cousin2[1,]
          indNames(cousin2) <- ("cousin2")
          
          ppoff <- rbind(x,off1, cousin1, cousin2)
          
          res <- gl.run.EMIBD9(ppoff, 
                               emibd9.path = emibd9.path)
          
          relo <- (res$rel["cousin1","cousin2"])
          return(relo)
        }
      }
  

#Now we repeat the simulation 
      if (rel == "full.sib") {
      rr <- replicate(n = nboots,
                      expr = full.sib(x))
      } else if (rel == "half.sib") {
          rr <- replicate(n = nboots,
                          half.sib(x))
        } else if (rel == "first.cousin") {
          rr <- replicate(n = nboots, 
                          first.cousin(x))
        }
        

        rr <- data.frame(rr)
        colnames(rr) <- c("Relatedness")
        sum <- summary(rr)
        mean_rel <- mean(rr$Relatedness)
        
#Add 95% CI's for simulated relatedness (does this match the propiosed fs,hs,fc?)
        l.model <- lm(Relatedness ~ 1, rr)
        CI <- confint(l.model, level = (conf))
        l.ci <- CI[,1]
        u.ci <- CI[,2]

#Set the verbosity
        if (verbose>0) {
          cat(
            report(
              "Returning summary statistics for simulated relatedness values:\n")
          )
          }
        
#Print plot 
        
        if (rel == "full.sib") {
          p1 <- ggplot(data = rr, aes(x = Relatedness)) +
            geom_histogram(binwidth = 0.005, colour = "black") +
            ggtitle("Histogram of simulated full sibling relatedness values and CI's") +
            theme(plot.title = element_text(hjust = 0.5)) +
            theme_dartR()+
            geom_vline(aes(xintercept = mean_rel), color = "red", linewidth = 1) + #mean
            geom_vline(aes(xintercept = l.ci), color = "green", linewidth = 1, linetype = 2) + #lower CI
            geom_vline(aes(xintercept = u.ci), color = "green", linewidth = 1, linetype = 2) #upper CI
          labs(y = "Count", x = "Relatedness") 
          if (plot.out) print (p1)
        } else if (rel == "half.sib") {
          p1 <- ggplot(data = rr, aes(x = Relatedness)) +
            geom_histogram(binwidth = 0.005, colour = "black") +
            ggtitle("Histogram of simulated half sibling relatedness values and CI's") +
            theme(plot.title = element_text(hjust = 0.5)) +
            labs(y = "Count", x = "Relatedness") + 
            theme_dartR() +
            geom_vline(aes(xintercept = mean_rel), color = "red", linewidth = 1) +
            geom_vline(aes(xintercept = l.ci), color = "green", linewidth = 1, linetype = 2) + #lower CI
            geom_vline(aes(xintercept = u.ci), color = "green", linewidth = 1, linetype = 2) #upper CI
          if (plot.out) print (p1)
        } else if (rel == "first.cousin") {
          p1 <- ggplot(data = rr, aes(x = Relatedness)) +
            geom_histogram(binwidth = 0.005, colour = "black") +
            ggtitle("Histogram of simulated first cousin relatedness values and CI's") +
            theme(plot.title = element_text(hjust = 0.5)) +
            theme_dartR()+
            geom_vline(aes(xintercept = mean_rel), color = "red", linewidth = 1) + 
            geom_vline(aes(xintercept = l.ci), color = "green", linewidth = 1, linetype = 2) + #lower CI
            geom_vline(aes(xintercept = u.ci), color = "green", linewidth = 1, linetype = 2) #upper CI
          labs(y = "Count", x = "Relatedness") 
          if (plot.out) print (p1)
        }

#Optionally save the plot 

  if(!is.null(plot.file)) {
  tmp <- utils.plot.save(p1,
                         dir=plot.dir,
                         file=plot.file,
                         verbose=verbose)
  }
        
        print(sum)
        print(CI)
}

##Things to add: 
# - parrallelise the function 
# - add more to the outputs i.e. return the object as a list
# - make the code tidier? - mmm nah 

#Simulations to run: 
# - 

















