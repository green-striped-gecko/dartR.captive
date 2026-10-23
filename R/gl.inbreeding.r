#' @name gl.inbreeding
#' @title Estimates individual inbreeding coefficients with PLINK and GCTA
#' @family inbreeding functions
#' @description
#' Estimates the inbreeding coefficient of each individual from SNP data by
#' running PLINK 1.9 (--het) and GCTA (--ibc), and plots the estimates of
#' each population, individuals sorted by value.
#'
#' @param x Name of the genlight object containing the SNP data [required].
#' @param method Programs to run: "plink", "gcta" or both [default
#' c("plink", "gcta")].
#' @param by.pop If TRUE, each population is analysed on its own, so allele
#' frequencies come from that population; if FALSE, from all individuals
#' pooled. See Details [default FALSE].
#' @param plink.path Path to the folder holding the PLINK 1.9 executable
#' (plink, or plink.exe on Windows), or the path to the executable itself.
#' PLINK is needed for both methods [default getwd()].
#' @param gcta.path Path to the folder holding the GCTA executable (gcta64,
#' or gcta64.exe on Windows), or the path to the executable itself. Only
#' used when method includes "gcta" [default getwd()].
#' @param plot.display If TRUE, the plot is displayed in the plot window
#' [default TRUE].
#' @param plot.theme Theme for the plot. See Details for options
#' [default theme_dartR()].
#' @param plot.colors A vector of colours, one per population (named by
#' population, or in the order of the population levels), or a palette
#' function [default gl.select.colors(x)].
#' @param plot.file Name for the RDS binary file to save (base name only,
#' exclude extension) [default NULL].
#' @param plot.dir Directory to save the plot RDS files [default as specified
#' by the global working directory or tempdir()].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default 2, unless specified using gl.set.verbosity].
#'
#' @details
#' Both programs measure the excess of homozygosity of each individual over
#' the homozygosity expected from allele frequencies estimated in the
#' sample. PLINK --het reports F = (O(HOM) - E(HOM)) / (N(NM) - E(HOM)),
#' where O(HOM) is the observed and E(HOM) the expected number of
#' homozygous loci, and N(NM) the number of loci called in the individual.
#' GCTA --ibc reports three estimators (Yang et al. 2011), with x the
#' allele count (0, 1, 2) at a locus, p the allele frequency and
#' h = 2p(1 - p), each averaged over the loci called in the individual:
#' Fhat1 = (x - 2p)^2 / h - 1, the diagonal of the genomic relationship
#' matrix minus 1; Fhat2 = 1 - x(2 - x) / h, based on excess homozygosity;
#' and Fhat3 = (x^2 - (1 + 2p)x + 2p^2) / h, based on the correlation
#' between uniting gametes. Fhat3 gives a large weight to rare alleles: an
#' individual homozygous for an allele at frequency p adds (1 - p) / p at
#' that locus, so a few such loci can push Fhat3 well above 1. Removing
#' rare alleles beforehand (gl.filter.maf) reduces this effect.
#'
#' Loci that are monomorphic or not called in the individuals analysed are
#' removed before the programs run, so PLINK and GCTA estimate on the same
#' loci (GCTA would otherwise count monomorphic loci in its denominator and
#' pull every estimate towards 0). All loci are treated as autosomal and
#' placed on one chromosome; neither estimator uses map positions. Remove
#' sex-linked loci first (for example with dartR.sexlinked
#' gl.drop.sexlinked), because loci on the X chromosome are hemizygous in
#' males and would make males appear more inbred.
#'
#' With by.pop = FALSE and more than one population, the expected
#' homozygosity comes from the pooled allele frequencies, so the
#' differentiation between populations (the Wahlund effect) is counted as
#' inbreeding in every individual. With by.pop = TRUE, each population is
#' run separately: populations of one individual get NA (every homozygous
#' locus is monomorphic, so F would be -1 by construction), and populations
#' of fewer than 10 individuals get a warning that their allele
#' frequencies, and hence the estimates, are imprecise.
#'
#' PLINK 1.9 can be downloaded from \url{https://www.cog-genomics.org/plink/}
#' and GCTA from \url{https://yanglab.westlake.edu.cn/software/gcta/}. The
#' files for both programs are written to a new folder under tempdir(),
#' which is deleted on exit.
#'
#' The plot shows PLINK F and GCTA Fhat3 side by side, one row of panels
#' per population. A colour vector can be obtained with gl.select.colors()
#' and passed with the plot.colors parameter. Themes can be obtained from
#' \url{https://ggplot2.tidyverse.org/reference/ggtheme.html}. If a
#' plot.file is given, the ggplot is saved as an RDS binary file using
#' saveRDS(); it can be reloaded with readRDS(). If a plot directory
#' (plot.dir) is specified, the ggplot binary is saved to that directory;
#' otherwise to the tempdir().
#'
#' @return A data frame, returned invisibly, with one row per individual in
#' the order of x: id, pop, and for PLINK the columns O.HOM, E.HOM, N.NM and
#' F, and for GCTA the columns NOMISS (loci used), Fhat1, Fhat2 and Fhat3.
#' Individuals that could not be estimated have NA.
#'
#' @author Author(s): Luis Mijangos. Custodian: Luis Mijangos -- Post to
#' \url{https://groups.google.com/d/forum/dartr}
#'
#' @examples
#' \dontrun{
#' # needs PLINK 1.9 and GCTA installed; point to their folders
#' res <- gl.inbreeding(platypus.gl, plink.path = "~/programs",
#'                      gcta.path = "~/programs")
#' # allele frequencies from each population
#' res.pop <- gl.inbreeding(platypus.gl, by.pop = TRUE,
#'                          plink.path = "~/programs",
#'                          gcta.path = "~/programs")
#' # PLINK only
#' res.plink <- gl.inbreeding(platypus.gl, method = "plink",
#'                            plink.path = "~/programs")
#' }
#'
#' @references
#' \itemize{
#' \item Chang, C. C., Chow, C. C., Tellier, L. C., Vattikuti, S., Purcell,
#' S. M., & Lee, J. J. (2015). Second-generation PLINK: rising to the
#' challenge of larger and richer datasets. GigaScience, 4(1), 7.
#' \item Purcell, S., Neale, B., Todd-Brown, K., Thomas, L., Ferreira, M. A.,
#' Bender, D., ... & Sham, P. C. (2007). PLINK: a tool set for whole-genome
#' association and population-based linkage analyses. The American Journal
#' of Human Genetics, 81(3), 559-575.
#' \item Yang, J., Lee, S. H., Goddard, M. E., & Visscher, P. M. (2011). GCTA:
#' a tool for genome-wide complex trait analysis. The American Journal of
#' Human Genetics, 88(1), 76-82.
#' }
#'
#' @seealso \code{\link{gl.report.kinship}}, \code{\link{gl.run.EMIBD9}}
#' @export

gl.inbreeding <- function(x,
                          method = c("plink", "gcta"),
                          by.pop = FALSE,
                          plink.path = getwd(),
                          gcta.path = getwd(),
                          plot.display = TRUE,
                          plot.theme = theme_dartR(),
                          plot.colors = NULL,
                          plot.file = NULL,
                          plot.dir = NULL,
                          verbose = NULL) {
  # PRELIMINARIES -- checking ----------------
  # SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)
  if (verbose == 0) {
    plot.display <- FALSE
  }

  # SET WORKING DIRECTORY
  plot.dir <- gl.check.wd(plot.dir, verbose = 0)

  # FLAG SCRIPT START
  funname <- match.call()[[1]]
  utils.flag.start(func = funname,
                   build = "v.2026.1",
                   verbose = verbose)

  # CHECK DATATYPE
  datatype <- utils.check.datatype(x, accept = "SNP", verbose = verbose)

  # FUNCTION SPECIFIC ERROR CHECKING
  method <- match.arg(method, several.ok = TRUE)
  if (!is.logical(by.pop) || length(by.pop) != 1 || is.na(by.pop)) {
    stop(error("  by.pop must be TRUE or FALSE\n"))
  }

  # a folder holding the executable, or the executable itself (GCTA's macOS
  # builds are not always named gcta64)
  find.exe <- function(path, name, prog) {
    if (.Platform$OS.type == "windows") {
      name <- paste0(name, ".exe")
    }
    exe <- if (utils::file_test("-f", path)) path else file.path(path, name)
    if (!utils::file_test("-f", exe)) {
      stop(error("  Cannot find", prog, "at", exe, "\n"))
    }
    normalizePath(exe)
  }
  # PLINK also builds the binary files GCTA reads, so it is always needed
  plink.exe <- find.exe(plink.path, "plink", "PLINK (plink.path)")
  if ("gcta" %in% method) {
    gcta.exe <- find.exe(gcta.path, "gcta64", "GCTA (gcta.path)")
  }

  # Assign a population if none is specified
  if (is.null(pop(x)) | is.na(length(pop(x))) | length(pop(x)) <= 0) {
    if (verbose >= 2) {
      cat(important("  Population assignments not detected, individuals",
                    "assigned to a single population labelled 'pop1'\n"))
    }
    pop(x) <- factor(rep("pop1", nInd(x)))
  }
  pop(x) <- droplevels(pop(x))

  if (!by.pop && nPop(x) > 1 && verbose >= 1) {
    cat(warn("  Warning: allele frequencies are pooled across", nPop(x),
             "populations, so F includes the differentiation between",
             "populations; use by.pop = TRUE for F relative to each",
             "population\n"))
  }

  # DO THE JOB ----------------------
  plink.cols <- c("O.HOM", "E.HOM", "N.NM", "F")
  gcta.cols <- c("NOMISS", "Fhat1", "Fhat2", "Fhat3")
  res <- data.frame(id = indNames(x),
                    pop = as.character(pop(x)),
                    stringsAsFactors = FALSE)
  if ("plink" %in% method) {
    res[plink.cols] <- NA_real_
  }
  if ("gcta" %in% method) {
    res[gcta.cols] <- NA_real_
  }

  if (by.pop) {
    units <- split(seq_len(nInd(x)), pop(x))
  } else {
    units <- list(all = seq_len(nInd(x)))
  }

  # a new folder per call, so no output of an earlier run can be read back
  rundir <- tempfile("gl.inbreeding_")
  dir.create(rundir)
  on.exit(unlink(rundir, recursive = TRUE), add = TRUE)
  # the programs' own logs are shown only in a full report
  verbose.prog <- if (verbose >= 5) 3 else 0

  run.gcta <- function(dir.unit) {
    args <- c("--bfile", shQuote(file.path(dir.unit, "gl_plink")),
              "--autosome", "--ibc",
              "--out", shQuote(file.path(dir.unit, "gcta_ibc")))
    gcta.out <- tryCatch(
      suppressWarnings(system2(gcta.exe, args, stdout = TRUE, stderr = TRUE)),
      error = function(e) structure(conditionMessage(e), status = 127L)
    )
    status <- attr(gcta.out, "status")
    if (verbose >= 5) {
      cat(report(paste(gcta.out, collapse = "\n"), "\n"))
    }
    ibc.file <- file.path(dir.unit, "gcta_ibc.ibc")
    if ((!is.null(status) && status != 0) || !file.exists(ibc.file)) {
      stop(error(
        "  GCTA failed running", gcta.exe, "\n  Last lines of the GCTA",
        "output:\n", paste(utils::tail(gcta.out, 5), collapse = "\n"), "\n"
      ))
    }
    utils::read.table(ibc.file, header = TRUE,
                      colClasses = c(FID = "character", IID = "character"))
  }

  n.loci <- rep(NA_integer_, length(units))
  names(n.loci) <- names(units)
  too.small <- NULL
  imprecise <- NULL
  no.loci <- NULL

  for (u in seq_along(units)) {
    idx <- units[[u]]
    label <- names(units)[u]
    if (length(idx) < 2) {
      too.small <- c(too.small, label)
      next
    }
    if (length(idx) < 10) {
      imprecise <- c(imprecise, label)
    }
    sub <- x[idx, ]
    # loci monomorphic or uncalled in these individuals carry no information,
    # and GCTA would count them in its denominator
    p <- glMean(sub)
    poly <- which(!is.na(p) & p > 0 & p < 1)
    if (length(poly) == 0) {
      no.loci <- c(no.loci, label)
      next
    }
    sub <- sub[, poly]
    n.loci[u] <- nLoc(sub)
    # one autosome for all loci: gl2plink codes loci at position 0 as
    # chromosome 0, which GCTA --autosome drops; neither estimator uses
    # positions
    sub$chromosome <- factor(rep("1", nLoc(sub)))
    sub$position <- seq_len(nLoc(sub))

    if (verbose >= 2) {
      cat(report(paste0("  Estimating inbreeding",
                        if (by.pop) paste0(" in population ", label),
                        " (", length(idx), " individuals, ", nLoc(sub),
                        " loci)\n")))
    }

    dir.unit <- file.path(rundir, paste0("unit", u))
    dir.create(dir.unit)
    gl2plink(sub, bed.files = FALSE, outfile = "gl_plink",
             outpath = dir.unit, verbose = 0)
    utils.plink.run(dir.in = dir.unit,
                    plink.cmd = basename(plink.exe),
                    plink.path = dirname(plink.exe),
                    out = "gl_plink",
                    syntax = "--file gl_plink --allow-no-sex --make-bed",
                    verbose = verbose.prog)
    # PLINK and GCTA keep the order of the .fam file, so rows are matched
    # by position; the IIDs confirm it
    iid <- gsub(" ", "_", indNames(sub))

    if ("plink" %in% method) {
      utils.plink.run(dir.in = dir.unit,
                      plink.cmd = basename(plink.exe),
                      plink.path = dirname(plink.exe),
                      out = "plink_het",
                      syntax = "--bfile gl_plink --allow-no-sex --het",
                      verbose = verbose.prog)
      het <- utils::read.table(file.path(dir.unit, "plink_het.het"),
                               header = TRUE,
                               colClasses = c(FID = "character",
                                              IID = "character"))
      if (!identical(het$IID, iid)) {
        stop(error("  PLINK returned individuals in an unexpected order\n"))
      }
      res[idx, plink.cols] <- het[, c("O.HOM.", "E.HOM.", "N.NM.", "F")]
    }

    if ("gcta" %in% method) {
      ibc <- run.gcta(dir.unit)
      if (!identical(ibc$IID, iid)) {
        stop(error("  GCTA returned individuals in an unexpected order\n"))
      }
      res[idx, gcta.cols] <- ibc[, gcta.cols]
    }
  }

  # these warnings change which estimates exist or how far to trust them,
  # so they print at verbose >= 1
  if (length(too.small) > 0 && verbose >= 1) {
    cat(warn("  Warning: populations with one individual have no allele",
             "frequencies of their own; estimates set to NA:",
             paste(too.small, collapse = ", "), "\n"))
  }
  if (length(no.loci) > 0 && verbose >= 1) {
    cat(warn("  Warning: no polymorphic loci; estimates set to NA:",
             paste(no.loci, collapse = ", "), "\n"))
  }
  if (length(imprecise) > 0 && verbose >= 1) {
    cat(warn("  Warning: fewer than 10 individuals, so allele frequencies",
             "and estimates are imprecise:",
             paste(imprecise, collapse = ", "), "\n"))
  }

  # PRINTING OUTPUTS -----------
  if (verbose >= 3) {
    est.cols <- intersect(c("F", "Fhat3"), names(res))
    summ <- stats::aggregate(res[est.cols],
                             by = list(pop = factor(res$pop,
                                                    levels = levels(pop(x)))),
                             FUN = mean, na.rm = TRUE)
    summ <- cbind(summ[1], n = as.vector(table(pop(x))), summ[-1])
    if (by.pop) {
      summ$n.loci <- n.loci[as.character(summ$pop)]
    } else {
      cat(report("  Loci used:", n.loci[1], "\n"))
    }
    cat(report("  Mean inbreeding per population:\n"))
    print(summ, row.names = FALSE)
  }

  # colours follow the population names when they cover every population,
  # otherwise the order of the population levels
  if (is.null(plot.colors)) {
    plot.colors <- gl.select.colors(x, plot.display = FALSE, verbose = 0)
  }
  if (is.function(plot.colors)) {
    plot.colors <- plot.colors(nPop(x))
  }
  if (!all(levels(pop(x)) %in% names(plot.colors))) {
    plot.colors <- rep_len(plot.colors, nPop(x))
    names(plot.colors) <- levels(pop(x))
  }

  value <- NULL  # silence R CMD check note for ggplot NSE
  plot.panel <- function(col, title, strips) {
    df <- res[!is.na(res[[col]]), c("id", "pop", col)]
    names(df)[3] <- "value"
    df$pop <- factor(df$pop, levels = levels(pop(x)))
    # individuals sorted by value within each population
    df$order <- stats::ave(-df$value, df$pop,
                           FUN = function(v) rank(v, ties.method = "first"))
    p <- ggplot(df, aes(x = value, y = order, colour = pop)) +
      geom_point(size = 3, alpha = 1 / 2) +
      facet_grid(rows = vars(pop), scales = "free_y") +
      scale_colour_manual(values = plot.colors) +
      ggtitle(title) +
      plot.theme +
      theme(legend.position = "none",
            axis.title = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank())
    if (strips) {
      # horizontal labels, so long population names stay readable
      p <- p + theme(strip.text.y = element_text(angle = 0))
    } else {
      p <- p + theme(strip.background = element_blank(),
                     strip.text.y = element_blank())
    }
    p
  }

  p1 <- NULL
  has.plink <- "plink" %in% method && any(!is.na(res$F))
  has.gcta <- "gcta" %in% method && any(!is.na(res$Fhat3))
  if (has.plink && has.gcta) {
    # population labels on the right-hand panel only
    p1 <- plot.panel("F", "PLINK F", strips = FALSE) +
      plot.panel("Fhat3", "GCTA Fhat3", strips = TRUE)
  } else if (has.plink) {
    p1 <- plot.panel("F", "PLINK F", strips = TRUE)
  } else if (has.gcta) {
    p1 <- plot.panel("Fhat3", "GCTA Fhat3", strips = TRUE)
  } else if (verbose >= 1) {
    cat(warn("  Warning: no estimates to plot\n"))
  }

  if (plot.display && !is.null(p1)) {
    print(p1)
  }

  # Optionally save the plot ---------------------
  if (!is.null(plot.file) && !is.null(p1)) {
    tmp <- utils.plot.save(p1,
                           dir = plot.dir,
                           file = plot.file,
                           verbose = verbose)
  }

  # FLAG SCRIPT END ---------------
  if (verbose >= 1) {
    cat(report("Completed:", funname, "\n"))
  }
  # ----------------------

  # RETURN
  invisible(res)
}
