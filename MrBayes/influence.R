wd <- if (basename(getwd()) == "TreeSearch") "../MrBayes/" else 
  if (basename(getwd()) == "MrBayes") "./" else "/MrBayes/"
exclusions <- list.files(wd, "^.*-no-.*\\.nex$")
base <- substr(exclusions[1], 0, regexpr("-no-", exclusions[1]) - 1)
nTrees <- 20
Distance <- TreeDist::ClusteringInfoDistance

#' @importFrom ape read.nexus
.ReadMrBayesTrees <- function(treeFile, burninFrac) {
  trees <- read.nexus(treeFile)
  trees[-seq_len(burninFrac * length(trees))]
}

ReadMrBayes <- function(filename, n = NULL, burninFrac = NULL) {
  lines <- readLines(filename, warn = FALSE)
  if (is.null(burninFrac)) {
    burninPattern <- ".*\\bmcmcp?\\b.*burninf?r?a?c?\\s*=\\s*([\\d\\.]+)\\b.*"
    burninFrac <- rev(as.numeric(
      gsub(burninPattern, "\\1", lines[grep(burninPattern, lines, TRUE, TRUE)],
           perl = TRUE, ignore.case = TRUE)
    ))[1]
  }
  treeFiles <- list.files(dirname(filename),
                          paste0(basename(filename), "\\.run\\d+\\.t"),
                          full.names = TRUE)
  trees <- tryCatch(do.call(c, lapply(treeFiles, .ReadMrBayesTrees,
                                      burninFrac = burninFrac)),
                    error = function(e) NULL)
  if (is.null(trees)) {
    warning("Could not read trees from ", filename)
    return(NULL)
  }
  if (!is.null(n)) {
    trees <- trees[seq(1, length(trees), length.out = n)]
  }
}

baseTrees <- ReadMrBayes(paste0(wd, base, ".nex"), n = nTrees)

res <- vapply(exclusions, function(file) {
  trees <- ReadMrBayes(paste0(wd, file), n = nTrees)
  if (is.null(trees)) {
    c(NA_real_, NA, NA)
  } else {
    thinnedTrees <- KeepTip(baseTrees, TipLabels(trees[[1]]))
    self <- Distance(baseTrees)
    d <- Distance(baseTrees, trees)
    c(mean(d) - mean(self), mean(d), sd(d) - sd(self))
  }
}, c(meanIncrease = 0, originalMean = 0, sdIncrease = 0))

colnames(res) <- substr(exclusions, nchar(base) + 5, nchar(exclusions) - 4)
write.table(res, file = "influence.txt")