coanct_clean <- function(input, coanctTests = NULL){
  
  # Coancestry estimators from gl.relatedness (relatedness scale, symmetric
  # matrices with diagonal NA); pairs sharing no called locus are NA
  res <- gl.relatedness(input, estimators = coanctTests, plot.out = FALSE,
                        verbose = 0)
  
  new_x <- NULL
  
  for(i in coanctTests){
    # same order as GRM_clean, which is bound to this table by position
    mat_coan <- res[[i]]
    attr(mat_coan, "scale") <- NULL
    ord <- colnames(mat_coan)[order(colnames(mat_coan))]
    coan_col <- mat_coan[ord, ord]
    coan_col[upper.tri(coan_col)] <- NA
    coan_col <- as.data.frame(as.table(coan_col))
    # relatedness to kinship
    coan_col$Freq <- coan_col$Freq/2
    names(coan_col)[names(coan_col) == 'Freq'] <- i
    
    if(is.null(new_x)){
      new_x <- coan_col
    }else{
      new_x[i] <- coan_col[i]
    }
    
  }
  
  return(new_x)
}

# Clean up gl.grm output 
GRM_clean <- function(input){
  GRM <- gl.grm(input, plotheatmap = FALSE, verbose = 0)
  order_grm <- colnames(GRM)[order(colnames(GRM))]
  GRM <- GRM[order_grm, order_grm]
  
  GRM_col <- GRM
  GRM_col[upper.tri(GRM_col)] <- NA
  GRM_col <- as.data.frame(as.table(as.matrix(GRM_col)))
  GRM_col$Freq <- GRM_col$Freq/2
  
  return(GRM_col$Freq)
}


# Cleanup relatedness output 
cleanup_rel <- function(input_val, testSelect=NULL){
  rel_cal <- cbind(coanct_clean(input_val, coanctTests=testSelect), GRM_clean(input_val))
  rel_cal <- rel_cal[complete.cases(rel_cal[testSelect]),]
  colnames(rel_cal) <- c("ind1","ind2", testSelect, "rrBLUP")
  rel_plot_2 <- reshape2::melt(rel_cal, id.vars = c("ind1","ind2"))
  
  #rel_plot_2 <- rel_cal %>%
  #  pivot_longer(cols = c("rrBLUP",which_tests),
  #               values_to = "value") %>%
  #  {colnames(.) <- c("ind1", "ind2", "variable", "value");.}
  
  return(rel_plot_2)
  
}

# Returns relatedness plot 
plot_rel <- function(cleanup_out){
  value <- NULL
  variable <- NULL
  yintercept<- NULL
  rel_plotgg <- ggplot(cleanup_out,aes(y=value,x=variable,color = variable,fill = variable)) +
    geom_point(position=position_jitterdodge(dodge.width=1),show.legend = F) +
    geom_violin(alpha=0.1,scale = "count")+
    geom_boxplot(alpha=0.35)+
    theme_bw(base_size = 16) +
    theme( legend.position = "bottom",
           axis.ticks.x=element_blank() ,
           axis.text.x = element_blank(),
           axis.title.x=element_blank(),
           legend.title=element_blank(),
           legend.text=element_text(size = 18),
           legend.key.height=unit(1, "cm"))+
    ylab("Kinship") 
  
  return(rel_plotgg)
  
}


# Extracts parents from iteration output 
# Inbreeding of each population of x (1 - sum Ho / sum uHe, as dartR.sim's
# real_inbreeding) on loci called in at least `threshold` of individuals,
# where lost heterozygote calls inflate F less. When fewer than min.loci
# loci pass, the threshold is lowered in steps of 0.01 until they do.
# Returns F per population (levels(pop(x)) order, negative values kept)
# with the threshold and number of loci used as attributes
inbreedingHighCallrate <- function(x, threshold = 0.99, min.loci = 100) {
  called <- colMeans(!is.na(as.matrix(x)))
  while (sum(called >= threshold) < min(min.loci, nLoc(x)) &&
         threshold > 0) {
    threshold <- round(threshold - 0.01, 2)
  }
  xs <- x[, called >= threshold]
  f <- vapply(seppop(xs), function(p) {
    g <- as.matrix(p)
    n <- colSums(!is.na(g))
    q <- colMeans(g, na.rm = TRUE) / 2
    ho <- colMeans(g == 1, na.rm = TRUE)
    he <- 2 * q * (1 - q) * 2 * n / (2 * n - 1)
    keep <- n > 0 & is.finite(he) & is.finite(ho)
    if (sum(he[keep]) == 0) return(NA_real_)
    1 - sum(ho[keep]) / sum(he[keep])
  }, numeric(1))
  f <- f[levels(pop(x))]
  attr(f, "threshold") <- threshold
  attr(f, "n.loci") <- nLoc(xs)
  f
}

# Sizes of the full-sib families (2 or more individuals) of each population
# of x, largest first, as a list in the order of levels(pop(x)); families
# come from an ind.metrics column or from COLONY ("colony")
resolveFamilies <- function(x, families, colony.path, verbose) {
  if (identical(families, "colony")) {
    if (is.null(colony.path)) {
      stop(error("  families = 'colony' needs colony.path, the folder with",
                 "the COLONY executable\n"))
    }
    col <- gl.run.colony(x, colony.path = colony.path, verbose = 0)
    bc <- col$best.config
    fam <- paste(bc$FatherID, bc$MotherID)[match(indNames(x),
                                                  bc$OffspringID)]
  } else {
    if (!is.character(families) || length(families) != 1 ||
        !families %in% colnames(x@other$ind.metrics)) {
      stop(error("  families must be 'colony' or the name of a column of",
                 "x@other$ind.metrics\n"))
    }
    fam <- as.character(x@other$ind.metrics[[families]])
  }
  fam[fam %in% ""] <- NA
  pops <- levels(pop(x))
  sizes <- lapply(pops, function(p) {
    tb <- table(fam[as.character(pop(x)) == p])
    sort(as.integer(tb[tb >= 2]), decreasing = TRUE)
  })
  names(sizes) <- pops
  if (verbose >= 2) {
    cat(report("  Full-sib family sizes of x:",
               paste0(pops, " ", vapply(sizes, function(v)
                 if (length(v)) paste(v, collapse = "/") else "none",
                 character(1)), collapse = "; "), "\n"))
  }
  sizes
}

# An LDNe estimate is unreliable when its jackknife upper limit is infinite
# or it is more than 10 times the sample size
neUnreliable <- function(est, ci.high, n) {
  is.finite(est) & (!is.finite(ci.high) | est > 10 * n)
}

# Effective population size of each population of x (in the order of
# levels(pop(x))): Ne given by the user, recycled to one per population, or
# estimated with gl.LDNe (critical allele frequency 0.05)
resolveNe <- function(x, Ne, neest.path, verbose) {
  pops <- levels(pop(x))
  if (!is.null(Ne)) {
    if (!length(Ne) %in% c(1, length(pops))) {
      stop(error("  Ne must be one value or one per population (",
                 length(pops), ")\n"))
    }
    return(rep_len(as.numeric(Ne), length(pops)))
  }
  if (is.null(neest.path)) {
    stop(error(
      "  The simulation needs the effective population size of each",
      "population of x: give Ne, or neest.path (folder with the NeEstimator",
      "binary) to estimate it with gl.LDNe\n"))
  }
  pkg <- "dartR.popgen"
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(error("  Estimating Ne needs the package dartR.popgen; install it",
               "or give Ne\n"))
  }
  ldne <- getExportedValue(pkg, "gl.LDNe")
  res <- tryCatch(
    ldne(x, neest.path = neest.path, critical = 0.05,
         singleton.rm = TRUE, plot.out = FALSE, verbose = 0),
    error = function(e) {
      stop(error("  gl.LDNe failed (", conditionMessage(e), "); give Ne",
                 "instead\n"), call. = FALSE)
    })
  est <- vapply(pops, function(p) {
    d <- res[[p]]
    if (is.null(d)) return(NA_real_)
    suppressWarnings(as.numeric(d[d$Statistic == "Estimated Ne^", 2]))
  }, numeric(1))
  # the jackknife upper limit is infinite, or Ne is far above the sample
  # size, when the sample holds too little LD information for a usable Ne
  ci.high <- vapply(pops, function(p) {
    d <- res[[p]]
    if (is.null(d)) return(NA_real_)
    suppressWarnings(as.numeric(d[d$Statistic == "CI high JackKnife", 2]))
  }, numeric(1))
  n.pop <- as.vector(table(pop(x))[pops])
  unreliable <- pops[neUnreliable(est, ci.high, n.pop)]
  if (length(unreliable) > 0 && verbose >= 1) {
    cat(warn(
      "  Warning: the Ne estimated with gl.LDNe is unreliable for",
      paste0(unreliable, " (Ne ", signif(est[unreliable], 3), ", n ",
             n.pop[match(unreliable, pops)], ", jackknife upper limit ",
             signif(ci.high[unreliable], 3), ")", collapse = "; "),
      "; consider giving Ne\n"))
  }
  bad <- pops[!is.finite(est) | est <= 0]
  if (length(bad) > 0) {
    stop(error("  gl.LDNe could not estimate a finite Ne for",
               paste(bad, collapse = ", "), "; give Ne instead\n"))
  }
  if (verbose >= 2) {
    cat(report("  Ne estimated with gl.LDNe:",
               paste0(pops, " ", signif(est, 3), collapse = "; "), "\n"))
  }
  unname(est)
}

# Rows of a pedigree (id, dad, mom, ...) for the individuals in ids and all
# their ancestors; the kinship of the ids needs no one else
ancestorPedigree <- function(ped, ids) {
  keep <- ids
  current <- ids
  while (length(current) > 0) {
    rows <- ped[ped$id %in% current, , drop = FALSE]
    parents <- setdiff(unique(c(rows$dad, rows$mom)), c(keep, NA))
    keep <- c(keep, parents)
    current <- parents
  }
  ped[ped$id %in% keep, , drop = FALSE]
}

# Binds ind.metrics tables whose columns differ (dartR.sim adds F_founder
# only to generation_0), filling missing columns with NA
bindIndMetrics <- function(tables) {
  cols <- unique(unlist(lapply(tables, colnames)))
  do.call(rbind, lapply(unname(tables), function(im) {
    for (cn in setdiff(cols, colnames(im))) im[[cn]] <- rep(NA, nrow(im))
    im[, cols, drop = FALSE]
  }))
}

ExtractParents <- function(inputClass, iteration=1){
  
  indDf <- bindIndMetrics(lapply(inputClass[[iteration]], function(g) {
    im <- g@other$ind.metrics
    rownames(im) <- indNames(g)
    im
  }))
  
  parental.df <- indDf %>%
    {. <- .[,c(3,4)]; .} %>%
    {colnames(.) <- c("dad", "mom"); .} %>%
    {.["id"] <- rownames(.); .} 
  # realised inbreeding of founders stored by dartR.sim (store_founders)
  if (!is.null(indDf$F_founder)) {
    parental.df$F <- as.numeric(indDf$F_founder)
  }
  
  
  
  return(parental.df)
  
}

################################################################################
#  Kinship-Classifier for Wright-Fisher Simulations (non-overlapping generations)
#
#  Purpose
#  -------
#  •  Build a pedigree from the three most recent generations of a Wright-Fisher
#     simulation and label every *detectable* pair of individuals as:
#       "parent_offspring", "full_sibs", "half_sibs",
#       "full_first_cousins", "second_cousins", or "unrelated".
#  •  Designed for very large simulated populations: avoids constructing the
#     full N × N Cartesian product by generating pairs *only* inside shared-
#     ancestor sets and by keying intermediate tables for O(log M) look-ups.
#
#  Key Features
#  ------------
#  ▸ **Linear climbs, local combinations** – three ancestor-climbing joins
#    (depth 1–3) are O(N); sibling/cousin pairs are formed inside each small
#    family group, not across the whole population.
#  ▸ **data.table back-end** – lightning-fast joins, grouping and keyed binary
#    searches; multi-threaded if `data.table` was compiled with OpenMP.
#  ▸ **Memory-safe** – keeps only essential columns, drops helpers on the fly,
#    and returns a single tidy table `related` plus an on-demand query function
#    `get_relationship(a, b)`.
#  ▸ **Configurable depth** – increase the `step_up()` loop to classify third
#    cousins, etc., without touching the rest of the logic.
#  ▸ **Pedigree-agnostic I/O** – expects three data.tables (`gen1`, `gen2`,
#    `gen3`) each with columns `id`, `dad`, `mom`; IDs must be unique across
#    generations; missing parents are `NA`.
#
#  Output Objects
#  --------------
#  ♦ `related`          – data.table [id1, id2, relationship] (deduplicated)
#  ♦ `get_relationship` – helper: fast O(log M) look-up for any two IDs
#
#  Dependencies
#  ------------
#    data.table ≥ 1.14.0   (base R ≥ 4.0 recommended for `\(x)` lambda syntax)
#
#  Usage
#  -----
#    source("kinship_classifier.R")            # after defining gen1 … gen3
#    related                                   # view all labelled pairs
#    get_relationship(1042, 2198)              # quick ad-hoc query
################################################################################


CleanupExtractParents <- function(parentalTable){
  
  child1 <- NULL
  child2 <- NULL
  dad <- NULL
  id <- NULL
  id1 <- NULL
  id2 <- NULL
  ind1 <- NULL
  ind2 <- NULL
  mom <- NULL
  child <- parent <- grandparent <- ggparent <- s1 <- s2 <- NULL
  
  ped <- as.data.table(parentalTable)
  
  # Ensure required columns exist
  stopifnot(all(c("id","dad","mom") %in% names(ped)))
  
  # Estimate number of generations
  gen_depth <- max(rle(!is.na(ped$dad) | !is.na(ped$mom))$lengths)
  
  # -----------------------------
  # 1. Parent-Offspring
  # -----------------------------
  po_dad <- ped[!is.na(dad), .(id1 = dad, id2 = id, relationship = "parent_offspring", r = 0.25)]
  po_mom <- ped[!is.na(mom), .(id1 = mom, id2 = id, relationship = "parent_offspring", r = 0.25)]
  parent_offspring <- unique(rbind(po_dad, po_mom, fill=TRUE))
  
  # -----------------------------
  # 2. Siblings
  # -----------------------------
  sibs <- merge(
    ped[, .(id, dad, mom)], 
    ped[, .(id2 = id, dad, mom)], 
    by = c("dad","mom"), allow.cartesian = TRUE
  )[id < id2]
  
  full_sibs <- sibs[!is.na(dad) & !is.na(mom),
                    .(id1 = id, id2, relationship = "full_sibs", r = 0.25)]
  
  # missing parents are unknown, not shared: data.table joins match NA to
  # NA, so they are excluded before joining
  half_sibs <- rbind(
    merge(ped[!is.na(dad),.(id,dad)], ped[!is.na(dad),.(id2=id,dad)], by="dad", allow.cartesian=TRUE)[id<id2,
                                                                                .(id1=id,id2,relationship="half_sibs", r=0.125)],
    merge(ped[!is.na(mom),.(id,mom)], ped[!is.na(mom),.(id2=id,mom)], by="mom", allow.cartesian=TRUE)[id<id2,
                                                                                .(id1=id,id2,relationship="half_sibs", r=0.125)],
    fill=TRUE
  )
  half_sibs <- fsetdiff(half_sibs, full_sibs[,.(id1,id2,relationship="half_sibs", r=0.125)])
  
  # -----------------------------
  # 3. First cousins
  # -----------------------------
  full_first_cousins <- data.table()
  half_first_cousins <- data.table()
  second_cousins <- data.table()
  
  if (gen_depth >= 3) {
    get_children <- rbind(
      ped[!is.na(dad),.(child=id, parent=dad)],
      ped[!is.na(mom),.(child=id, parent=mom)],
      fill=TRUE
    )
    
    # Full first cousins
    cousin_pairs_full <- merge(get_children, full_sibs[,.(p1=id1,p2=id2)], 
                               by.x="parent", by.y="p1", allow.cartesian=TRUE)
    cousin_pairs_full <- merge(cousin_pairs_full, get_children, 
                               by.x="p2", by.y="parent", suffixes=c("1","2"), allow.cartesian=TRUE)
    full_first_cousins <- unique(
      cousin_pairs_full[child1<child2, .(id1=child1, id2=child2, relationship="full_first_cousins", r=0.0625)]
    )
    
    # Half first cousins
    cousin_pairs_half <- merge(get_children, half_sibs[,.(p1=id1,p2=id2)], 
                               by.x="parent", by.y="p1", allow.cartesian=TRUE)
    cousin_pairs_half <- merge(cousin_pairs_half, get_children, 
                               by.x="p2", by.y="parent", suffixes=c("1","2"), allow.cartesian=TRUE)
    half_first_cousins <- unique(
      cousin_pairs_half[child1<child2, .(id1=child1, id2=child2, relationship="half_first_cousins", r=0.03125)]
    )
  }
  
  # -----------------------------
  # 4. Second cousins
  # -----------------------------
  if (gen_depth >= 4 && (nrow(full_first_cousins) > 0 | nrow(half_first_cousins) > 0)) {
    fc_all <- rbind(full_first_cousins[,.(id1,id2)], half_first_cousins[,.(id1,id2)], fill=TRUE)
    
    fc_parents <- unique(rbind(
      ped[!is.na(dad),.(child=id, parent=dad)],
      ped[!is.na(mom),.(child=id, parent=mom)],
      fill=TRUE
    ))
    
    sc <- merge(fc_parents, fc_all[,.(p1=id1,p2=id2)], by.x="parent", by.y="p1", allow.cartesian=TRUE)
    sc <- merge(sc, fc_parents, by.x="p2", by.y="parent", suffixes=c("1","2"), allow.cartesian=TRUE)
    
    second_cousins <- unique(
      sc[child1<child2, .(id1=child1, id2=child2, relationship="second_cousins", r=0.015625)]
    )
  }
  
  # -----------------------------
  # 5. Relatives across generations
  # -----------------------------
  # without these classes grandparents and aunts or uncles of an individual
  # would be labelled unrelated
  pc <- unique(rbind(
    ped[!is.na(dad), .(child = id, parent = dad)],
    ped[!is.na(mom), .(child = id, parent = mom)]
  ))
  gp <- unique(merge(pc, pc[, .(parent = child, grandparent = parent)],
                     by = "parent", allow.cartesian = TRUE)[
                       , .(child, grandparent)])
  grandparents <- gp[, .(id1 = grandparent, id2 = child,
                         relationship = "grandparent_grandchild", r = 0.125)]
  great_grandparents <- unique(
    merge(gp, pc[, .(grandparent = child, ggparent = parent)],
          by = "grandparent", allow.cartesian = TRUE)[
            , .(id1 = ggparent, id2 = child,
                relationship = "great_grandparent_grandchild", r = 0.0625)])
  # an aunt or uncle is a sibling of a parent
  uncles <- function(sibs, label, r) {
    both <- rbind(sibs[, .(s1 = id1, s2 = id2)], sibs[, .(s1 = id2, s2 = id1)])
    unique(merge(pc, both, by.x = "parent", by.y = "s1",
                 allow.cartesian = TRUE)[
                   , .(id1 = s2, id2 = child, relationship = label, r = r)])
  }
  avuncular <- uncles(full_sibs, "avuncular", 0.125)
  half_avuncular <- uncles(half_sibs, "half_avuncular", 0.0625)

  # -----------------------------
  # Combine All
  # -----------------------------
  all_rel <- rbindlist(list(
    parent_offspring,
    full_sibs,
    half_sibs,
    full_first_cousins,
    half_first_cousins,
    second_cousins,
    grandparents,
    great_grandparents,
    avuncular,
    half_avuncular
  ), fill=TRUE)
  
  return(all_rel[])
  
}


printCorVals <- function(corValues, whichTests){
  
  cat("Number iterations:", length(corValues), "\n")
  for(i in 1:length(corValues)){
    cat("Iteration:", i, "\n")
    cat("Correlation between `related` and:","\n")
    colSelect <- colnames(corValues[[i]])[-which(colnames(corValues[[i]])=="rel")]
    for(j in 1:length(colSelect)){
      colSelect <- colnames(corValues[[i]])[-which(colnames(corValues[[i]])=="rel")]
      cat(colSelect[j],": ", corValues[[i]]["rel", colSelect[j]], "\n", sep = "")
    }
  }
}

# Bias (estimate minus pedigree kinship) and RMSE of each estimator by
# relationship class, as tiles coloured by bias and labelled
# "bias (RMSE)"; classes absent from the data are left out
biasTilePlot <- function(relatedDf, which_tests) {
  estimator <- RelDegree <- bias <- label <- NULL
  b <- calcBias(list(relatedDf), which_tests)[[1]]
  r <- calcRMSE(list(relatedDf), which_tests)[[1]]
  df <- data.frame(
    estimator = factor(rep(rownames(b), times = ncol(b)),
                       levels = which_tests),
    RelDegree = factor(rep(colnames(b), each = nrow(b)),
                       levels = rev(relationshipClasses)),
    bias = unlist(b), rmse = unlist(r))
  df <- df[!is.na(df$bias), ]
  df$label <- sprintf("%.3f\n(%.3f)", df$bias, df$rmse)
  lim <- max(abs(df$bias), 0.01)
  ggplot(df, aes(x = estimator, y = RelDegree, fill = bias)) +
    geom_tile(color = "white") +
    geom_text(aes(label = label), size = 2.8) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                         midpoint = 0, limits = c(-lim, lim)) +
    theme_bw() +
    labs(x = "Estimator", y = "Relationship class",
         fill = "Bias",
         title = "Bias (RMSE) against the pedigree kinship")
}

relatedLevelPlots <- function(relatedDf, which_tests, pedSim=F){
  
  value <- NULL
  variable <- NULL
  yintercept<- NULL
  RelDegree <- NULL
  
  if(pedSim){
    
    df1 <- reshape2::melt(relatedDf, id.vars = "RelDegree", measure.vars = which_tests) %>%
      {.$RelDegree <- factor(.$RelDegree, 
                             levels = relationshipClasses); .}
    
    # other_relatives has no single expected kinship
    lines_df <- data.frame(
      RelDegree = relationshipClasses,
      yintercept = unname(relationshipKinship)
    )
    lines_df <- lines_df[!is.na(lines_df$yintercept), ]
    
    lines_df$RelDegree <- factor(lines_df$RelDegree,
                                 levels = relationshipClasses)
    
    outputBoxPlot <- ggplot(df1, aes(x=variable,y=value,color=variable,
                                     fill=variable))+
      geom_boxplot(alpha=0.5,show.legend = F)+
      geom_hline(
        data = lines_df, 
        aes(yintercept = yintercept), 
        color = "red",
        linetype = "dashed", 
        linewidth = 1             
      ) +
      facet_wrap(~ RelDegree,   scales = "free", ncol =2) + 
      theme_bw() + 
      labs(
        x = "Estimator",
        y = "Relatedness Value"
      )
    
    outputDensityPlot <- ggplot(df1, aes(x = value, color= RelDegree,
                                         fill = RelDegree)) + 
      geom_density(alpha=0.5) + 
      geom_vline(
        data = lines_df, 
        aes(xintercept = yintercept), 
        color = "red",
        linetype = "dashed", 
        linewidth = 1             
      ) +
      facet_wrap(~ variable,scales ="fixed") + 
      theme_bw() + 
      labs(
        x = "Relatedness Value",
        y = "Count"
      )
    
    asf <- NULL
    asf[[1]] <- outputBoxPlot
    asf[[2]] <- outputDensityPlot
    asf[[3]] <- biasTilePlot(relatedDf, which_tests)
    
    return(asf)
  }else{
    
    df1 <- relatedDf
    
    outputBoxPlot <- ggplot(df1, aes(x=variable,y=value,color=variable,
                                     fill=variable))+
      geom_boxplot(alpha=0.5,show.legend = F) + 
      labs(
        x="Estimator", 
        y="Relatedness Value"
      )
    
    return(outputBoxPlot)
  }
}

runE9 <- function(inputObj, e9Path, numCores, e9parallel=e9parallel, E9Inbreed=F){
  e9Name <- "E9"
  if(E9Inbreed){
    e9Name <- "E9_Inbred"
  }
  e9runObj <- gl.run.EMIBD9(inputObj,
                            emibd9.path = e9Path, 
                            parallel=e9parallel, 
                            ncores=numCores,
                            Inbreed = E9Inbreed, 
                            plot.out = F) %>%
    {. <- .$rel}%>%
    {.[upper.tri(.)] <- NA; .} %>% 
    as.matrix() %>%
    as.table() %>%
    as.data.frame() %>%
    {colnames(.) <- c("ind1", "ind2", e9Name); .} %>%
    {. <- na.omit(.);.}
  
  return(e9runObj)
}


# Exact kinship coefficients from a pedigree (id, dad, mom; missing parents
# NA or 0; optional column F, the inbreeding of founders). Founders are
# unrelated to each other; a founder with F gets K[i, i] = (1 + F) / 2, the
# others (and parents that are not listed, which are added as founders) are
# not inbred. Individuals are ordered so that parents come before their
# offspring, then the standard recursion is applied:
#   K[i, i] = (1 + K[dad, mom]) / 2
#   K[i, j] = (K[dad, j] + K[mom, j]) / 2   (j before i)
pedigreeKinship <- function(ped) {
  founder.F <- if (!is.null(ped$F)) as.numeric(ped$F) else NA_real_
  ped <- data.frame(id = as.character(ped$id),
                    dad = as.character(ped$dad),
                    mom = as.character(ped$mom),
                    F = founder.F,
                    stringsAsFactors = FALSE)
  ped$dad[ped$dad %in% c("0", "")] <- NA
  ped$mom[ped$mom %in% c("0", "")] <- NA
  if (anyDuplicated(ped$id)) {
    stop(error("  The pedigree lists some individuals more than once\n"))
  }
  founders <- setdiff(unique(c(ped$dad, ped$mom)), c(ped$id, NA))
  if (length(founders) > 0) {
    ped <- rbind(data.frame(id = founders, dad = NA_character_,
                            mom = NA_character_, F = NA_real_,
                            stringsAsFactors = FALSE),
                 ped)
  }
  ordered.ids <- character(0)
  remaining <- ped
  while (nrow(remaining) > 0) {
    ready <- (is.na(remaining$dad) | remaining$dad %in% ordered.ids) &
      (is.na(remaining$mom) | remaining$mom %in% ordered.ids)
    if (!any(ready)) {
      stop(error("  The pedigree contains a loop (an individual is its own",
                 "ancestor)\n"))
    }
    ordered.ids <- c(ordered.ids, remaining$id[ready])
    remaining <- remaining[!ready, , drop = FALSE]
  }
  ped <- ped[match(ordered.ids, ped$id), ]
  n <- nrow(ped)
  d <- match(ped$dad, ped$id)
  m <- match(ped$mom, ped$id)
  K <- matrix(0, n, n, dimnames = list(ped$id, ped$id))
  for (i in seq_len(n)) {
    if (i > 1) {
      earlier <- seq_len(i - 1)
      k.dad <- if (!is.na(d[i])) K[d[i], earlier] else 0
      k.mom <- if (!is.na(m[i])) K[m[i], earlier] else 0
      K[i, earlier] <- (k.dad + k.mom) / 2
      K[earlier, i] <- K[i, earlier]
    }
    K[i, i] <- (1 + if (!is.na(d[i]) && !is.na(m[i])) K[d[i], m[i]] else
      if (is.na(d[i]) && is.na(m[i]) && !is.na(ped$F[i])) ped$F[i] else 0) / 2
  }
  K
}

# Relationship classes, from the closest to the most distant, with their
# kinship when the pedigree has no inbreeding; a pair that the classifier
# places in several classes keeps the closest one, and a pair in none of
# them is "other_relatives" when its pedigree kinship is above 0 and
# "unrelated" otherwise
relationshipKinship <- c(parent_offspring = 0.25, full_sibs = 0.25,
                         half_sibs = 0.125, grandparent_grandchild = 0.125,
                         avuncular = 0.125, full_first_cousins = 0.0625,
                         great_grandparent_grandchild = 0.0625,
                         half_avuncular = 0.0625,
                         half_first_cousins = 0.03125,
                         second_cousins = 0.015625, other_relatives = NA,
                         unrelated = 0)
relationshipClasses <- names(relationshipKinship)

# Value of a variable in a dartR.sim variable file, as written there
simVariableValue <- function(file, variable) {
  v <- utils::read.csv(file, stringsAsFactors = FALSE)
  trimws(v$value[v$variable == variable])
}

# Copy of a dartR.sim variable file with some values changed, written to a
# temporary file whose path is returned
simVariableFile <- function(file, changes) {
  v <- utils::read.csv(file, stringsAsFactors = FALSE)
  for (k in names(changes)) {
    v$value[v$variable == k] <- changes[[k]]
  }
  out <- tempfile(fileext = ".csv")
  utils::write.csv(v, out, row.names = FALSE)
  out
}

# Adds the pedigree truth to the relatedness estimates: one row per pair
# with ind1, ind2, RelDegree (closest relationship class), rel (exact
# pedigree kinship) and one column per estimator
mergePedigreeTruth <- function(relatedDf, ped) {

  ind1 <- ind2 <- ID1 <- ID2 <- id1 <- id2 <- relationship <- NULL

  estimates <- as.data.frame(stats::na.omit(relatedDf))
  estimates$ind1 <- as.character(estimates$ind1)
  estimates$ind2 <- as.character(estimates$ind2)
  estimates$ID1 <- pmin(estimates$ind1, estimates$ind2)
  estimates$ID2 <- pmax(estimates$ind1, estimates$ind2)
  estimates <- as.data.frame(tidyr::pivot_wider(
    estimates[, c("ID1", "ID2", "variable", "value")],
    names_from = "variable", values_from = "value"))

  # missing parents may be coded 0 or ""; the classifier needs NA
  ped <- data.frame(id = as.character(ped$id), dad = as.character(ped$dad),
                    mom = as.character(ped$mom),
                    F = if (!is.null(ped$F)) as.numeric(ped$F) else NA_real_,
                    stringsAsFactors = FALSE)
  ped$dad[ped$dad %in% c("0", "")] <- NA
  ped$mom[ped$mom %in% c("0", "")] <- NA

  classes <- as.data.frame(CleanupExtractParents(ped[, c("id", "dad", "mom")]))
  classes$ID1 <- pmin(as.character(classes$id1), as.character(classes$id2))
  classes$ID2 <- pmax(as.character(classes$id1), as.character(classes$id2))
  classes$rank <- match(classes$relationship, relationshipClasses)
  classes <- classes[order(classes$rank), ]
  classes <- classes[!duplicated(classes[, c("ID1", "ID2")]),
                     c("ID1", "ID2", "relationship")]

  out <- merge(estimates, classes, by = c("ID1", "ID2"), all.x = TRUE)

  K <- pedigreeKinship(ped)
  in.ped <- out$ID1 %in% rownames(K) & out$ID2 %in% rownames(K)
  out$rel <- NA_real_
  out$rel[in.ped] <- K[cbind(out$ID1[in.ped], out$ID2[in.ped])]

  none <- is.na(out$relationship)
  out$relationship[none] <- ifelse(!is.na(out$rel[none]) & out$rel[none] > 0,
                                   "other_relatives", "unrelated")

  est.cols <- setdiff(colnames(out), c("ID1", "ID2", "relationship", "rel"))
  out <- data.frame(ind1 = out$ID1, ind2 = out$ID2,
                    RelDegree = out$relationship, rel = out$rel,
                    out[, est.cols, drop = FALSE],
                    stringsAsFactors = FALSE, check.names = FALSE)
  out
}

mergeE9Related <- function(relatedDf, RecodeDf,test_select){
  
  ind1 <- NULL
  ind2 <- NULL
  ID1 <- NULL
  ID2 <- NULL
  
  relatedTransform <- relatedDf %>%
    rbind() %>%
    na.omit() 
  relatedTransform$ind1 <- as.character(relatedTransform$ind1)
  relatedTransform$ind2 <- as.character(relatedTransform$ind2)
  
  recodeBound <- RecodeDf %>%
    as.data.frame() %>%
    na.omit()
  recodeBound$ind1 <- as.character(recodeBound$ind1)
  recodeBound$ind2 <- as.character(recodeBound$ind2)
  
  
  setDT(recodeBound); setDT(relatedTransform)
  # 1. create canonical ID pairs in-place
  recodeBound[ , `:=`(
    ID1 = pmin(ind1, ind2),
    ID2 = pmax(ind1, ind2)
  )]
  relatedTransform[ , `:=`(
    ID1 = pmin(ind1, ind2),
    ID2 = pmax(ind1, ind2)
  )]
  
  # 2. set the join keys (very fast lookup)
  setkey(recodeBound, ID1, ID2)
  setkey(relatedTransform, ID1, ID2)
  
  # 3. merge (inner join; drop unmatched)
  #    this will bring relatedTransform’s stats alongside recodeBound’s
  merged <- recodeBound[relatedTransform, nomatch=0]
  # 4. clean up helper columns if you like
  merged[ , c("ind1","ind2") := NULL]  # drop originals
  # or rename ID1/ID2 back to sample1/sample2
  setnames(merged, c("ID1","ID2"), c("ind1","ind2"))
  
  mergedWider <- as.data.frame(merged) %>%
    pivot_wider(names_from = "variable", values_from = c("value")) %>%
    {. <- .[,c("ind1", "ind2", test_select)];.} %>%
    {. <- reshape2::melt(., id.vars=c("ind1", "ind2")); .}
  
  return(mergedWider)
  
}

# Root mean square error of each estimator against the exact pedigree
# kinship (column rel), by relationship class
calcRMSE <- function(inputDf, which_tests){
  lapply(inputDf, function(df) {
    out <- matrix(NA_real_, nrow = length(which_tests),
                  ncol = length(relationshipClasses),
                  dimnames = list(which_tests, relationshipClasses))
    for (j in which_tests) {
      for (k in relationshipClasses) {
        rows <- df$RelDegree == k & !is.na(df[[j]]) & !is.na(df$rel)
        if (any(rows)) {
          out[j, k] <- sqrt(mean((df[[j]][rows] - df$rel[rows])^2))
        }
      }
    }
    as.data.frame(out)
  })
}

# Bias (mean of estimate minus pedigree kinship) of each estimator, by
# relationship class
calcBias <- function(inputDf, which_tests){
  lapply(inputDf, function(df) {
    out <- matrix(NA_real_, nrow = length(which_tests),
                  ncol = length(relationshipClasses),
                  dimnames = list(which_tests, relationshipClasses))
    for (j in which_tests) {
      for (k in relationshipClasses) {
        rows <- df$RelDegree == k & !is.na(df[[j]]) & !is.na(df$rel)
        if (any(rows)) {
          out[j, k] <- mean(df[[j]][rows] - df$rel[rows])
        }
      }
    }
    as.data.frame(out)
  })
}

# Copies the missing-data pattern of x onto simulated genotypes: each
# simulated individual takes the missing loci of an individual of x drawn at
# random from the same population (from all of x when the population is not
# in x). Needs the same loci in the same order.
copyMissing <- function(sim, x) {
  miss <- is.na(as.matrix(x))
  if (!any(miss)) return(sim)
  pop.x <- as.character(pop(x))
  pop.sim <- as.character(pop(sim))
  donor <- vapply(seq_len(nInd(sim)), function(i) {
    pool <- which(pop.x == pop.sim[i])
    if (length(pool) == 0) pool <- seq_len(nInd(x))
    pool[sample.int(length(pool), 1)]
  }, integer(1))
  g <- as.matrix(sim)
  g[miss[donor, , drop = FALSE]] <- NA
  out <- new("genlight", g, ind.names = indNames(sim),
             loc.names = locNames(sim), pop = pop(sim), ploidy = ploidy(sim))
  out@loc.all <- sim@loc.all
  out@position <- sim@position
  out@chromosome <- sim@chromosome
  out@other <- sim@other
  if (is(sim, "dartR")) out <- methods::as(out, "dartR")
  out
}

# Variance of each estimator, by relationship class
calcVar <- function(inputDf, which_tests){
  lapply(inputDf, function(df) {
    out <- matrix(NA_real_, nrow = length(which_tests),
                  ncol = length(relationshipClasses),
                  dimnames = list(which_tests, relationshipClasses))
    for (j in which_tests) {
      for (k in relationshipClasses) {
        v <- df[[j]][df$RelDegree == k & !is.na(df[[j]])]
        if (length(v) > 1) {
          out[j, k] <- stats::var(v)
        }
      }
    }
    as.data.frame(out)
  })
}

tableColor <- function(dfIn){
  val_to_col <- function(x,
                         col.palette = gl.colors("con")) {
    # scale to 0–1
    rng <- range(x, na.rm = TRUE)
    scaled <- (x - rng[1]) / (rng[2] - rng[1])
    # map to colors (blue low, red high)
    cols <- col.palette(100)
    cols[as.integer(scaled * 99) + 1]
  }
  
  # Apply coloring to the data matrix (not headers)
  cell_fill <- matrix(val_to_col(as.matrix(dfIn)),
                      nrow = nrow(dfIn),
                      ncol = ncol(dfIn))
  
  # Add header row color (e.g., grey)
  header_fill <- rep("grey80", ncol(dfIn))
  
  # Combine header + body fills
  fills <- rbind(header_fill, cell_fill)
  
  # Create table grob with custom fills
  table_grob <- tableGrob(
    dfIn,
    theme = ttheme_default(
      core = list(bg_params = list(fill = cell_fill, col = "black")),
      colhead = list(bg_params = list(fill = header_fill))
    )
  )
  
  # Plot with ggplot background
  plot <- ggplot() +
    annotation_custom(table_grob)
  return(plot)
}

tableOut <- function(valIn){
  listOut <- NULL
  for(i in 1:length(valIn)){
    listOut[[i]] <- tableColor(valIn[[i]])
  }
  return(listOut)
}

# Pedigree attached to the genlight object (ind.metrics columns id, dad,
# mom; missing parents 0 or NA)
attachedPedigree <- function(baseInput) {
  im <- baseInput@other$ind.metrics
  data.frame(id = as.character(im$id), dad = as.character(im$dad),
             mom = as.character(im$mom), stringsAsFactors = FALSE)
}
