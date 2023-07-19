wd <- if (basename(getwd()) == "TreeSearch") "../MrBayes/" else 
  if (basename(getwd()) == "MrBayes") "./" else "/MrBayes/"
exclusions <- list.files(wd, "^.*-no-.*\\.nex$")
base <- substr(exclusions[1], 0, regexpr("-no-", exclusions[1]) - 1)
nTrees <- 20


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
  trees <- do.call(c, lapply(treeFiles, .ReadMrBayesTrees,
                             burninFrac = burninFrac))
  if (!is.null(n)) {
    trees <- trees[seq(1, length(trees), length.out = n)]
  }
}

baseTrees <- ReadMrBayes(paste0(wd, base, ".nex"), n = nTrees)

