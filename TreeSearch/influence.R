devtools::load_all("../../TreeSearch")
source("common.R")

kValues <- c(10, 40, 3, 20, 6) # Concavity constants for implied weighting
timeout <- 60 # Minutes after which to terminate each search
ratchets <- 8 # Ratchet iterations
hits <- 100 # Maximum times to hit best tree

kValues <- c(10) # Concavity constants for implied weighting
timeout <- 3 # Minutes after which to terminate each search
ratchets <- 8 # Ratchet iterations
hits <- 60 # Maximum times to hit best tree

# Load data from locally downloaded matrix
latest <- LatestMatrix()
message("* Reading ", latest)
dat <- ReadAsPhyDat(latest)
#dat <- dat[names(dat)[substr(names(dat), 1, 11) != "Yanjiahella"]]
dat <- dat[!names(dat) %in% c(
  "Yanjiahella_biscarpa_cs1_aulacophore",
  "Mitrate_Jaekelocarpus_oklahomensis"
)]

resultsFile <- ResultsFile(latest, "ew")
infFile <- InfluenceFile(latest)
infDir <- InfluenceCache(latest)

startTree <- LatestTree(dat, "ew")
if (is.null(startTree)) {
  startTree <- AdditionTree(dat)
}

ew <- TaxonInfluence(
  dataset = dat,
  tree = startTree,
  maxHits = hits,
  ratchIter = ratchets,
  maxTime = timeout,
  savePath = paste0(infDir, "/ew_"),
  useCache = TRUE
)
results <- ew
rownames(results) <- paste0("ew_", rownames(ew))
write.table(results, file = infFile)

# ew <- as.matrix(read.table(infFile))

par(mar = rep(0, 4), cex = 0.8)
maxPossible <- ClusteringEntropy(PectinateTree(NTip(dat) - 1)) * 2
upperBound <- max(ew)
nBin <- 128
bin <- cut(
  ew["dwMean", ],
  breaks = seq(0, upperBound, length.out = nBin),
  include.lowest = TRUE
)
palette <- hcl.colors(nBin, "inferno")

plot(startTree, tip.color = palette[bin])
PlotTools::SpectrumLegend(
  "bottomleft",
  palette = palette,
  title = paste("Tip influence\n Max:", signif(maxPossible, 3), "bits"),
  legend = signif(seq(upperBound, 0, length.out = 4), 3),
  bty = "n"
)

for (k in kValues) {
  startTree <- LatestTree(dat, paste0("iw", k))
  if (is.null(startTree)) {
    startTree <- AdditionTree(dat, concavity = k)
  }
  kRes <- cbind(TaxonInfluence(
    dataset = dat,
    tree = startTree,
    concavity = k,
    maxHits = hits,
    ratchIter = ratchets,
    savePath = paste0(infDir, "/iw", k, "_"),
    useCache = TRUE,
    maxTime = timeout
  ))
  rownames(kRes) <- paste0("iw", k, "_", rownames(kRes))
  results <- rbind(results, kRes)
  write.table(results, file = infFile)
}

if (interactive()) {
  par(mar = rep(0, 4), cex = 0.8)
  maxPossible <- ClusteringEntropy(PectinateTree(NTip(dat) - 1)) * 2
  upperBound <- max(results["iw10_max", ])
  nBin <- 128
  bin <- cut(
    results["iw10_dwMean", ],
    breaks = seq(0, upperBound, length.out = nBin),
    include.lowest = TRUE
  )
  palette <- hcl.colors(nBin, "inferno")
  
  plot(startTree, tip.color = palette[bin])
  PlotTools::SpectrumLegend(
    "bottomleft",
    palette = palette,
    title = paste("Tip influence\n Max:", signif(maxPossible, 3), "bits"),
    legend = signif(seq(upperBound, 0, length.out = 4), 3),
    bty = "n"
  )
}
  
