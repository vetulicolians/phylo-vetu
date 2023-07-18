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

influence <- tryCatch(as.matrix(read.table(InfluenceFile(latest))),
                      error = function(e) NULL)

treeFiles <- list.files(
  path = wd,
  pattern = paste0(".+_", sub("^.*/", "", latest), ".trees"),
  full.names = TRUE
)


for (treeFile in treeFiles) {
  trees <- read.nexus(treeFile)
  prefix <- strsplit(basename(treeFile), "_")[[1]][1]
  
  # Ignore outgroup taxa that aren't in tree
  og <- intersect(outgroup, TipLabels(trees)[[1]])
  if (length(og)) {
    # Root trees on outgroup
    trees <- RootTree(trees, og)
  }
  rogues <- Rogue::QuickRogue(trees, p = 1)
  cons <- SortTree(ConsensusWithout(trees, rogues[-1, "taxon"]))
  
  pdf(gsub(".trees", ".pdf", treeFile, fixed = TRUE), 
      width = 8, height = 10)
  
  stabCol <- Rogue::ColByStability(trees)
  Plot(cons, tip.col = stabCol)
  k <- KValue(treeFile)
  legend(
    "topright",
    c(
      sub("^(?:.*/)*([^/_]+)_.+", "\\1", treeFile, perl = TRUE),
      paste("Score:", signif(TreeLength(trees[1], dat, concavity = k)))
    ),
    bty = "n" # No bounding box
  )
  PlotTools::SpectrumLegend(
    "bottomright",
    palette = hcl.colors(131, "inferno")[1:101],
    title = paste0("Tip stability",
                   if(nrow(rogues) > 1) "\n(before rogue removal)"),
    legend = c("Least stable", "-", "Most stable"),
    bty = "n"
  )
  if (nrow(rogues) > 1) {
    rogueTaxa <- rogues[-1, "taxon"]
    legend("topleft", rogueTaxa, bty = "n", lty = 2, col = stabCol[rogueTaxa])
  }
  
  influenceAvailable <- !is.null(influence) &&
    paste0(prefix, "_max") %in% rownames(influence)
  
  tip.col <- if (influenceAvailable) {
    TipCol(cons$tip.label)
  } else {
    maxPossible <- ClusteringEntropy(PectinateTree(NTip(dat) - 1)) * 2
    upperBound <- max(influence[paste0(prefix, "_max"), ])
    nBin <- 128
    bin <- cut(
      influence[paste0(prefix, "_dwMean"), ],
      breaks = seq(0, upperBound, length.out = nBin),
      include.lowest = TRUE
    )
    palette <- hcl.colors(nBin, "inferno")
    palette[bin]
  }
  Plot(cons, tip.col = tip.col)
  
  if (nrow(rogues) > 1) {
    rogueTaxa <- rogues[-1, "taxon"]
    legend("topleft", rogueTaxa, bty = "n", lty = 2, col = tip.col[rogueTaxa])
  }

  if (influenceAvailable) {
    PlotTools::SpectrumLegend(
      "bottomright",
      palette = palette,
      title = paste("Tip influence\n Max:", signif(maxPossible, 3), "bits"),
      legend = signif(seq(upperBound, 0, length.out = 4), 3),
      bty = "n"
    )
    
    influencers <- colnames(influence)[
      order(influence[paste0(prefix, "_dwMean"), ], decreasing = TRUE)[1:3]]
    for (influencer in influencers) {
      cache <- paste0(InfluenceCache(latest), "/",
                      prefix, "_", fs::path_sanitize(influencer),
                      ".nex")
      if (file.exists(cache)) {
        infTrees <- read.nexus(cache)
        # Ignore outgroup taxa that aren't in tree
        infOG <- intersect(outgroup, TipLabels(infTrees)[[1]])
        if (length(infOG)) {
          # Root trees on outgroup
          infTrees <- RootTree(infTrees, infOG)
        }
        infRogues <- Rogue::QuickRogue(infTrees, p = 1)
        infCons <- SortTree(ConsensusWithout(infTrees, infRogues[-1, "taxon"]))
        
        infStabCol <- Rogue::ColByStability(infTrees)
        Plot(infCons, tip.col = infStabCol)
        
        legend(
          "topright",
          paste("Search without", influencer),
          bty = "n" # No bounding box
        )
        
        if (nrow(infRogues) > 1) {
          rogueTaxa <- infRogues[-1, "taxon"]
          legend("topleft", rogueTaxa, bty = "n", lty = 2,
                 col = stabCol[rogueTaxa])
        }
        PlotTools::SpectrumLegend(
          "bottomright",
          palette = hcl.colors(131, "inferno")[1:101],
          title = paste0("Tip stability",
                         if(nrow(infRogues) > 1) "\n(before rogue removal)"),
          legend = c("Least stable", "--", "Most stable"),
          bty = "n"
        )
      }
    }
    
  }
  
  
  
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
