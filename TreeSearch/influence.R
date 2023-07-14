devtools::load_all("../../TreeSearch")
source("common.R")

kValues <- c(10, 40, 3, 20, 6) # Concavity constants for implied weighting
timeout <- 60 # Minutes after which to terminate each search
ratchets <- 8 # Ratchet iterations
hits <- 100 # Maximum times to hit best tree

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

startTree <- LatestTree(dat, "ew")
if (is.null(startTree)) {
  startTree <- AdditionTree(dat)
}

ew <- TaxonInfluence(
  dataset = dat,
  tree = startTree,
  maxHits = hits,
  ratchIter = ratchets,
  maxTime = timeout
)
results <- cbind(ew = ew)
write.csv(results, file = "influence.csv")


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
    maxTime = timeout
  ))
  setNames(kRes, paste0("k", k))
  results <- cbind(results, kRes)
  write.nexus(best, file = resultsFile)
}
