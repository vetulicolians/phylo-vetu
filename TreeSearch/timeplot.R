minEdge <- 0.7

# Print current working directory, which should contain the scripts,
# matrix and trees.
getwd()

# If getwd() does not contain the relevant files, set wd to working directory
wd <- "./"
outgroup <- c("") # Specify taxa on which to root tree

source(paste0(wd, "/common.R"))
source(paste0(wd, "/plot.R"))

ages <- data.frame(readxl::read_xlsx("../Ages.xlsx", sheet = "taxa"),
                   row.names = "morphoname")


latest <- LatestMatrix(wd)
dat <- ReadAsPhyDat(latest)
if (outgroup == "") {
  outgroup <- names(dat)[1]
}
treeFiles <- list.files(
  path = wd,
  pattern = paste0(".+_", sub("^.*/", "", latest), ".trees"),
  full.names = TRUE
)
treeFile <- treeFiles[1]

for (treeFile in treeFiles) {
  trees <- read.nexus(treeFile)
  
  # Ignore outgroup taxa that aren't in tree
  og <- intersect(outgroup, TipLabels(trees)[[1]])
  if (length(og)) {
    # Root trees on outgroup
    trees <- RootTree(trees, og)
  }
  # rogues <- Rogue::QuickRogue(trees, p = 1)
  rogues <- Rogue::RogueTaxa(trees, threshold = 100)
  cons <- SortTree(ConsensusWithout(trees, rogues[-1, "taxon"]))
  
  pdf(gsub(".trees", ".timed.pdf", treeFile, fixed = TRUE), 
      width = 8, height = 10)
  
  plotted <- ColPlot(cons, ec = "black", minEdge = minEdge)
  pp <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  outgroupEdge <- plotted$tree$edge[, 2] == 
    match(ages[outgroup, "taxon"], 
          plotted$tree$tip.label)
  outgroupTime <- ages[outgroup, "time"] 
  rootAge <- outgroupTime + 
    plotted$tr$edge.length[outgroupEdge]
  # axisPhylo(line = -2, cex = 0.4, root.time = rootAge)
  
  outgroupX <- pp$xx[pp$edge[outgroupEdge]]
  TimeToX <- function(time) {
    (rootAge - time) * 
    (rootAge - outgroupTime) / (outgroupX[2] - outgroupX[1])
  }
  
  ts <- timescales$ICS2015
  plottedPeriods <- c("Cambrian", "Ordovician")
  epochs <- ts[ts$Type %in% "Epoch" & ts$Part_of %in% plottedPeriods, ]
  subepochs <- ts[ts$Type %in% "Age" & ts$Part_of %in% plottedPeriods, ]
  
  rect(
    xleft = TimeToX(subepochs$Start),
    ybottom = 0,
    xright = TimeToX(subepochs$End),
    ytop = 0.5,
    col = rgb(subepochs$Col_R, subepochs$Col_G, subepochs$Col_B, 255, NULL, 255)
  )
  text(TimeToX(subepochs$Midpoint), 0.25, subepochs$Abbrev, cex = 0.5)
  rect(
    xleft = TimeToX(epochs$Start),
    ybottom = -0.5,
    xright = TimeToX(epochs$End),
    ytop = 0,
    col = rgb(epochs$Col_R, epochs$Col_G, epochs$Col_B, 255, NULL, 255)
  )
  text(TimeToX(epochs$Midpoint), -0.25, epochs$Abbrev, cex = 0.5)
  
  periods <- ts[ts$Type %in% "Period" & ts$Name %in% plottedPeriods, ]
  abline(
    v = TimeToX(periods$End),
    lty = "dotted",
    col = rgb(periods$Col_R, periods$Col_G, periods$Col_B, 180, NULL, 255)
  )
  
  gps <- !duplicated(ages$group)
  legend("bottomright", ages[gps, "group"][-1],
         bty = "n", cex = 0.8,
         pch = 15, col = TipCol(rownames(ages)[gps][-1]))
  if (nrow(rogues) > 1) {
    legend("topleft", ages[rogues[-1, "taxon"], "taxon"], bty = "n", lty = 2,
           text.font = 3, text.col = TipCol(rogues[-1, "taxon"]))
  }
  k <- KValue(treeFile)
  legend(
    "topright",
    c(
      sub("^(?:.*/)*([^/_]+)_.+", "\\1", treeFile, perl = TRUE),
      paste("Score:", signif(TreeLength(trees[1], dat, concavity = k)))
    ),
    bty = "n" # No bounding box
  )
  
  
  distances <- TreeDist::ClusteringInfoDistance(trees)
  whenHit <- gsub("_\\d+$", "", names(trees), perl = TRUE)
  firstHit <- table(whenHit)
  searchStages <- length(firstHit)
  map <- cmdscale(distances, k = 3)
  cols <- hcl.colors(searchStages, alpha = 0.8)
  presOrder <- c("seed", "start", paste0("ratch", 1:10000), "final")
  presOrder <- c(presOrder, setdiff(names(firstHit), presOrder))
  treeCols <- cols[match(whenHit, intersect(presOrder, whenHit))]
  
  # Prepare plotting area
  par(mar = rep(0, 4))
  plot(map, type = "n", axes = FALSE, xlab = "", ylab = "", asp = 1)
  
  # Add minimum spanning tree
  TreeTools::MSTEdges(distances, plot = TRUE, map[, 1], map[, 2],
                      col = '#00000030', lty = 2)
  
  # Connect trees by order found
  lines(map[, 1], map[, 2], col = "#ffccaa", lty = 1)
  
  # Add points
  TreeDist::Plot3(map,
                  col = treeCols,
                  pch = 16, cex = 2,
                  add = TRUE)
  
  # Add legends
  legend("topright", 
         intersect(presOrder, names(firstHit)),
         col = cols, pch = 16, bty = "n")
  legend("topleft", 
         c("Minimum spanning tree (mapping distortion)",
           "Order in which trees found"),
         lty = c(2, 1),
         col = c("#00000030", "#ffccaaaa"),
         bty = "n")
  
  dev.off()
}
