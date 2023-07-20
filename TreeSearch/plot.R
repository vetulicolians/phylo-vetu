bold <- character(0)

TipCol <- function (tip.label) {
  c("black" = "#222222", palette.colors(pal = "Tableau 10"))[
    match(ages[tip.label, "group"], unique(ages$group))]
}

# REQUIRE tr, a phylo object.
Plot <- function (tr, pdf = FALSE, direction = "rightwards", font = 3,
                  plotw = 3, ploth = 2.5, minEdge = 2,
                  pts = 10, ec = "black", bi = FALSE, annot = FALSE,
                  bi.nudge = 0,
                  col.factor = 1, brightest = 0.9, filename = "plot/plot.pdf",
                  tip.col, fig = FALSE) {
  if (pdf) {
    pdf(file = paste(filename, "pdf", sep = "."), width = plotw,
        height = ploth, pointsize = pts, colormodel = 'rgb')
    on.exit(dev.off())
  }
  tip.id <- tr$tip.label
  tip.label <- ages[tip.id, "taxon"]
  if (any(is.na(tip.label))) {
    warning("Keys not found in ages.xlsx: ",
            paste(tip.id[is.na(tip.label)], collapse = ", "))
  }
  nTip <- length(tip.id)
  nNode <- tr$Nnode
  
  if (missing(tip.col)) {
    if (fig == FALSE) {
      tip.col <- TipCol(tip.id)
    } else {
      tip.col <- "black";
    }
  }
  
  if (bi) label.offset <- 2 * min(tr$edge.length)
  oPar <- par(cex = 0.8)  # Character expansion
  on.exit(par(oPar), add = TRUE)
  if (direction == 'rightwards') {
    align <- 0
  } else {
    align <- 0.5
  }
  tr$tip.label <- tip.label
  tr <- TipTimedTree(tr, ages[tip.id, "time"], minEdge)
  
  plotted <- plot(
    tr,
    edge.color = ec,
    edge.width = 2,
    font = font,
    #cex = 1,
    tip.col = tip.col,
    adj = align,
    label.offset = 0.5,
    #use.edge.length = bi,
    direction = direction,
    no.margin = TRUE,
    root.edge = TRUE,
    underscore = TRUE
  )
  
  if (annot) {
    labex <- regexpr("([0-9]+)", tr$node.label);
    lablen <- attr(labex, 'match.length')
    lab <- ifelse(lablen > 0 & lablen < 3,
                  substr(tr$node.label, labex, labex + lablen - 1),
                  " ")
    nodelabels(lab, adj = c(nudgel + bi.nudge, -0.5), frame = 'none', cex = 0.8)
  }
  
  list(
    plotted = plotted,
    tree = tr
  )
}

ColPlot <- function (tr, taxnames = "", direction = "rightwards",
                     ec = 0, ...) {
  tip.id <- tr$tip.label
  tip.label <- ages[tip.id, "taxon"]
  nTip <- length(tip.id)
  nNode <- tr$Nnode
  
  tip.col <- TipCol(tip.id)
  
  for (tax in names(taxnames)) {
    taxa <- taxnames[[tax]]
    tr <- drop.tip(tr, which(tr$tip.label%in%taxnames[[tax]]), subtree = TRUE)
    new.clade.name <- paste(tax, " (", length(taxnames[[tax]]), ")", sep = "")
    tr$tip.label[length(tr$tip.label)] <- new.clade.name
  }
  
  nTip <- length(tr$tip.label)
  nNode <- tr$Nnode
  font <- 3 + (tr$tip.label %in% bold)
  Plot(tr, direction = direction, font = font, ec = ec, tip.col = tip.col, ...)
}

RoguePlot <- function(trees, outgroup, p = 1) {
  # Ignore outgroup taxa that aren't in tree
  outgroup <- intersect(outgroup, TipLabels(c(trees)[[1]]))
  if (length(outgroup)) {
    # Root trees on outgroup
    trees <- RootTree(trees, outgroup)
  }
  rogues <- Rogue::QuickRogue(trees, p = p)
  cons <- SortTree(ConsensusWithout(trees, rogues[-1, "taxon"], p = p))
  
  ColPlot(cons, ec = "black")
  if (nrow(rogues) > 1) {
    legend("topleft", rogues[-1, "taxon"], bty = "n", lty = 2)
  }
  invisible(cons)
}
