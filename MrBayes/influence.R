library("TreeTools") # for median.multiPhylo
library("TreeDist") # for median.multiPhylo

wd <- if (basename(getwd()) == "TreeSearch") "../MrBayes/" else 
  if (basename(getwd()) == "MrBayes") "./" else "/MrBayes/"
exclusions <- list.files(wd, "^.*-no-.*\\.nex$")
base <- substr(exclusions[1], 0, regexpr("-no-", exclusions[1]) - 1)
nTrees <- 20
nTrees <- 100
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

RogueCons <- function(trees) {
  rogues <- Rogue::RogueTaxa(trees)[-1, "taxon"]
  
  cons <- ConsensusWithout(trees, rogues, p = 0.5)
  stabCol <- Rogue::ColByStability(trees)
  plot(cons, tip.col = stabCol[cons$tip.label])
  PlotTools::SpectrumLegend("bottomright", bty = "n",
                            title = "Leaf stability",
                            legend = c("Unstable", "", "", "Stable"),
                            palette = hcl.colors(131, "inferno")[1:101])
  
  splitP <- SplitFrequency(cons, trees) / length(trees)
  LabelSplits(cons, frame = "none", pos = 3, signif(splitP * 100, 2),
              unit = "%", col = SupportColor(splitP), cex = 0.8)
  if (length(rogues)) {
    legend("topright", lty = "dotted", gsub("_", " ", rogues),
           text.col = stabCol[rogues],
           text.font = 3, bty = "n", cex = 0.8)
  }
}

res <- vapply(exclusions, function(file) {
  trees <- ReadMrBayes(paste0(wd, file), n = nTrees)
  if (is.null(trees)) {
    rep(NA_real_, 7)
  } else {
    taxon <- gsub("_", " ", fixed = TRUE,
                  substr(file, nchar(base) + 5, nchar(file) - 4))
    cli::cli_h2(taxon)
    
    stopifnot(length(baseTrees) == length(trees))
    
    thinnedTrees <- KeepTip(baseTrees, TipLabels(trees[[1]]))
    d <- as.matrix(Distance(c(trees, thinnedTrees)))
    
    without <- d[seq_len(nTrees), seq_len(nTrees)]
    with <- d[-seq_len(nTrees), -seq_len(nTrees)]
    
    mdnWithout <- median(trees, index = TRUE)
    mdnWith <- median(thinnedTrees, index = TRUE)
    
    dIn <- median(vapply(seq_len(nTrees), function (i) {
      theseTrees <- seq_len(nTrees)[-i]
      median(d[i, theseTrees]) / median(d[i, nTrees + seq_len(nTrees)])
    }, double(1)))
    
    pdf(paste0(wd, sub(".nex", ".pdf", fixed = TRUE, file)), 7, 8.8)
    par(mar = c(0, 0, 1, 0), cex = 0.9)
    
    RogueCons(trees)
    title(paste0(gsub("_", " ", taxon), " removed a priori"),
          font.main = 3)
    RogueCons(thinnedTrees)
    title(paste0(gsub("_", " ", taxon), " removed a posteriori"),
          font.main = 3)
    
    
    map <- cmdscale(d, k = 2)
    plot(map, asp = 1, ann = FALSE, axes = FALSE,
         pch = c(rep(16, length(trees)), rep(1, length(thinnedTrees))),
         col = c(rep(3, length(trees)), rep(8, length(thinnedTrees))))
    points(map[c(mdnWithout, length(trees) + mdnWith), ],
           col = c(3, 8), pch = c(3, 4), cex = 2.5)
    legend("topleft", bty = "n", text.font = 3, taxon)
    legend("topright",
           c("a priori", "a posteriori", "median"),
           pch = c(16, 1, 3), col = c(3, 8, 1),
           title = "Taxon removed:", bty = "n")
    qual <- TreeDist::MappingQuality(d, dist(map))[1:2]
    legend("bottom", bty = "n", paste(names(qual), "=", signif(qual, 3)),
           title = "Tree space mapping quality")
    dev.off()
    
    # Return:
    c(mean(d[seq_len(nTrees), nTrees + seq_len(nTrees)]), # 10.1002/ece3.3941
      median(without), median(with),
      Distance(thinnedTrees[[mdnWith]], trees[[mdnWithout]]),
      dIn,
      sd(without), sd(with))
  }
}, c(tii = 0, mdnPre = 0, mdnPost = 0, distMdns = 0, withOwn = 0,
     sdPre = 0, sdPost = 0))
colnames(res) <- substr(exclusions, nchar(base) + 5, nchar(exclusions) - 4)

write.table(res, file = paste0(wd, base, "-influence.txt"))

ShowMe <- function(what, title = "", sub = "", diverging = FALSE) {
  nBin <- 128
  palRange <- if (diverging) {
    extreme <- max(abs(range(what, na.rm = TRUE)))
    c(-extreme, extreme)
  } else {
    range(what, na.rm = TRUE)
  }
  palette <- if (diverging) {
    bin <- cut(what, seq(palRange[1], palRange[2], length.out = nBin),
               include.lowest = TRUE)
    hcl.colors(nBin, "Green-Orange", rev = TRUE)
  } else {
    bin <- cut(what, breaks = nBin, include.lowest = TRUE)
    hcl.colors(nBin * 1.1, "plasma")
  }
  
  cons <- Consensus(baseTrees, p = 0.5)
  tipIndex <- TreeDist::LAPJV(adist(TipLabels(cons), colnames(res)))$matching
  
  oPar <- par(mar = rep(0, 4), cex = 0.8)
  on.exit(par(oPar))
  plot(cons, tip.color = palette[bin][tipIndex])
  PlotTools::SpectrumLegend(
    "bottomright",
    palette = palette[seq_len(nBin)],
    legend = signif(seq(palRange[2], palRange[1], length.out = 5), 3),
    bty = "n"
  )
  title(title, font.main = 1, line = -1, sub = sub)
}

# load from file:
# res <- as.matrix(read.table(paste0(wd, base, "-influence.txt")))

pdf(paste0(wd, base, "-influence.pdf"), 7, 9, title = "Taxon influence")
par(mar = c(0, 0, 1, 0), cex = 0.9)
# Bright: Including this taxon changes the central tendency of trees
ShowMe(res["distMdns", ], "Distance between median trees",
       sub = "High value: including taxon changes typical tree")

# Bright: Including this taxon changes the central tendency of trees
ShowMe(-log(res["withOwn", ]), diverging = TRUE,
       "log(Similarity to self: similarity to other)",
       sub = "High value: Distinct results if taxon removed a priori")

# Bright: Including this taxon results in different trees
ShowMe(res["tii", ], "Taxon influence index")

# Green: Including this taxon INCREASES uncertainty
ShowMe(res["mdnPost", ] - res["mdnPre", ], diverging = TRUE,
       "Increase in median distance when taxon removed a priori",
       sub = "Positive value: Including the taxon in analysis decreases uncertainty")
# ShowMe(res["sdPre", ] - res["sdPost", ], diverging = TRUE)
dev.off()

extra <- c("vetu", "vetulicolans")
res["mdnPost", extra] - res["mdnPre", extra]
-log(res["withOwn", extra])

cli::cli_h1("Complete")