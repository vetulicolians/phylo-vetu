source(paste0(wd, "/common.R"))
source(paste0(wd, "/plot.R"))


ages <- data.frame(readxl::read_xlsx("../Ages.xlsx", sheet = "taxa"),
                   row.names = "morphoname")

latest <- LatestMatrix(wd)

mbt <- MrBayesTrees(paste0("../MrBayes/", basename(latest)), n = 100)

cons <- ConsensusWithout(mbt, p = 0.5,
                         c(
                           "Herpetogaster_collinsi",
                           # "Yanjiahella_biscarpa",
                           # "Amiskwia_sagittiformis",
                           "Rotadiscus"
                           )
                         ); Plot(cons)
sort(cons$tip.label)

