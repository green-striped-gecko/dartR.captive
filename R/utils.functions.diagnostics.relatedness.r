coanct_clean <- function(input, coanctTests = NULL){
  
  # Gets gl2realted input or coancestry 
  x_gl2 <- gl2related(input, save = FALSE, verbose = 0)
  
  # Crappy code to avoid having 7 if loops to set tests based 
  # on function input 
  tests <- c("trioml", "wang", "lynchli", "lynchrd", "ritland",
             "quellergt", "dyadml") 
  test_select <- rep(0, 7)
  
  test_select[which(tests %in% coanctTests)] <- 1
  
  # 'related' is not on CRAN so it cannot be a declared dependency; its
  # availability is checked in gl.diagnostics.relatedness, and the namespace
  # lookup here works whether or not the package is attached
  coancestry <- utils::getFromNamespace("coancestry", "related")

  x_coancest <- coancestry(x_gl2, trioml = test_select[1], wang = test_select[2],
                           lynchli = test_select[3], lynchrd = test_select[4], ritland = test_select[5],
                           quellergt = test_select[6], dyadml = test_select[7])
  
  test_col <- which(colnames(x_coancest$relatedness) %in% coanctTests)
  
  x_coancest2 <- x_coancest$relatedness[,c(2,3,test_col)]
  x_coancest3 <- cbind(x_coancest2[,2],
                       x_coancest2[,1],
                       x_coancest2[,coanctTests])
  colnames(x_coancest3) <- colnames(x_coancest2)
  x_coancest4 <- rbind(x_coancest2,x_coancest3)
  x_coancest4[,3] <- as.numeric(x_coancest4[,3])
  new_x <- NULL
  
  for(i in coanctTests){
    mat_coan <- as.matrix(acast(x_coancest4, ind1.id~ind2.id, value.var=i))
    
    mat_coan <- apply(mat_coan, 2, as.numeric)
    mat_coan
    rownames(mat_coan) <- colnames(mat_coan)
    mat_coan
    coan_col <- mat_coan
    coan_col[upper.tri(coan_col)] <- NA
    coan_col
    coan_col <- as.data.frame(as.table(as.matrix(coan_col)))
    coan_col$Freq <- coan_col$Freq/2
    names(coan_col)[names(coan_col) == 'Freq'] <- i
    
    
    if(is.null(new_x) == T){
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
ExtractParents <- function(inputClass, iteration=1){
  
  indDf <- NULL
  
  for(i in 1:length(inputClass[[iteration]])){
    rownames(inputClass[[iteration]][[i]]@other$ind.metrics) <- 
      indNames(inputClass[[iteration]][[i]])
    indDf <- rbind(indDf, inputClass[[iteration]][[i]]@other$ind.metrics)
  }
  
  parental.df <- indDf %>%
    {. <- .[,c(3,4)]; .} %>%
    {colnames(.) <- c("dad", "mom"); .} %>%
    {.["id"] <- rownames(.); .} 
  
  
  
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
  # Combine All
  # -----------------------------
  all_rel <- rbindlist(list(
    parent_offspring,
    full_sibs,
    half_sibs,
    full_first_cousins,
    half_first_cousins,
    second_cousins
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

relatedLevelPlots <- function(relatedDf, which_tests, pedSim=F){
  
  value <- NULL
  variable <- NULL
  yintercept<- NULL
  RelDegree <- NULL
  
  if(pedSim){
    
    df1 <- reshape2::melt(relatedDf, id.vars = "RelDegree", measure.vars = which_tests) %>%
      {.$RelDegree <- factor(.$RelDegree, 
                             levels = relationshipClasses); .}
    
    lines_df <- data.frame(
      RelDegree = relationshipClasses,
      yintercept = c(0.25, 0.25, 0.125, 0.0625, 0.03125, 0.015625, 0)
    )
    
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
# NA or 0). Parents that are not listed as individuals are added as unrelated,
# non-inbred founders. Individuals are ordered so that parents come before
# their offspring, then the standard recursion is applied:
#   K[i, i] = (1 + K[dad, mom]) / 2
#   K[i, j] = (K[dad, j] + K[mom, j]) / 2   (j before i)
pedigreeKinship <- function(ped) {
  ped <- data.frame(id = as.character(ped$id),
                    dad = as.character(ped$dad),
                    mom = as.character(ped$mom),
                    stringsAsFactors = FALSE)
  ped$dad[ped$dad %in% c("0", "")] <- NA
  ped$mom[ped$mom %in% c("0", "")] <- NA
  if (anyDuplicated(ped$id)) {
    stop(error("  The pedigree lists some individuals more than once\n"))
  }
  founders <- setdiff(unique(c(ped$dad, ped$mom)), c(ped$id, NA))
  if (length(founders) > 0) {
    ped <- rbind(data.frame(id = founders, dad = NA_character_,
                            mom = NA_character_, stringsAsFactors = FALSE),
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
    K[i, i] <- (1 + if (!is.na(d[i]) && !is.na(m[i])) K[d[i], m[i]] else 0) / 2
  }
  K
}

# Relationship classes, from the closest to the most distant; a pair that
# the classifier places in several classes keeps the closest one, and a pair
# in none of them is "unrelated"
relationshipClasses <- c("parent_offspring", "full_sibs", "half_sibs",
                         "full_first_cousins", "half_first_cousins",
                         "second_cousins", "unrelated")

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
                    mom = as.character(ped$mom), stringsAsFactors = FALSE)
  ped$dad[ped$dad %in% c("0", "")] <- NA
  ped$mom[ped$mom %in% c("0", "")] <- NA

  classes <- as.data.frame(CleanupExtractParents(ped))
  classes$ID1 <- pmin(as.character(classes$id1), as.character(classes$id2))
  classes$ID2 <- pmax(as.character(classes$id1), as.character(classes$id2))
  classes$rank <- match(classes$relationship, relationshipClasses)
  classes <- classes[order(classes$rank), ]
  classes <- classes[!duplicated(classes[, c("ID1", "ID2")]),
                     c("ID1", "ID2", "relationship")]

  out <- merge(estimates, classes, by = c("ID1", "ID2"), all.x = TRUE)
  out$relationship[is.na(out$relationship)] <- "unrelated"

  K <- pedigreeKinship(ped)
  in.ped <- out$ID1 %in% rownames(K) & out$ID2 %in% rownames(K)
  out$rel <- NA_real_
  out$rel[in.ped] <- K[cbind(out$ID1[in.ped], out$ID2[in.ped])]

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
