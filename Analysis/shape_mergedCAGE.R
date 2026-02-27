packages<-c("rtracklayer","tidyverse","BSgenome.Hsapiens.UCSC.hg38","ggridges", "patchwork")
invisible(lapply(packages, library, character.only = TRUE))
lib_gr=import.bed("Library_data/res/library.bed")
epd=rtracklayer::import.bed("External_data/EPD/human38_epdnew.bed") %>% 
  promoters(downstream=16, upstream = 236)

theme_set(theme_bw())

filenames=paste0("External_data/FANTOM5/",
                 unique(cell_ont$`Facet ontology term`),".bed")
filenames=filenames[filenames%in%list.files("External_data/FANTOM5/", full.names = T)]
sample_name=str_remove(filenames, "External_data/FANTOM5/")%>%
  str_remove(".bed")
cage_data=map(filenames, ~import.bed(.x) %>% plyranges::join_overlap_inner_within(lib_gr))
names(cage_data)=sample_name

summ_cage=cage_data %>%
  imap(~as_tibble(.x) %>% mutate(sample=.y)) %>% 
  list_rbind() %>% group_by(seqnames, start,strand, name) %>% 
  summarise(score=sum(score))

summ_cage_gr=summ_cage %>% mutate(end=start) %>% GRanges()
export.bed(summ_cage_gr, "Analysis/Tables/merged_cage.bed")

summ_cager=CAGEr::CAGEexp( genomeName     = "BSgenome.Hsapiens.UCSC.hg38"
                           , inputFiles     = "Analysis/Tables/merged_cage.bed"
                           , inputFilesType = "bedScore"
                           , sampleLabels   = "merged")
cager=CAGEr::getCTSS(summ_cager)
lib_gr$gene_name=lib_gr$name
cager2=CAGEr::annotateCTSS(cager, lib_gr, upstream = 0, downstream = 0)
CAGEr::plotReverseCumulatives(cager2, fitInRange = c(100, 1E6))
cager2 <- CAGEr::normalizeTagCount(cager2, method = "powerLaw", fitInRange = c(100, 1E6), alpha = 1.04, T = 1E9)
cager2 <- CAGEr::filterLowExpCTSS(cager2, thresholdIsTpm = F, nrPassThreshold = 1, threshold = 10)
cager2 <- CAGEr::distclu(cager2, maxDist = 5, keepSingletonsAbove = 10)
cager2 <- CAGEr::cumulativeCTSSdistribution(cager2, clusters = "tagClusters", useMulticore = T)
cager2 <- CAGEr::quantilePositions(cager2, clusters = "tagClusters", qLow = 0.05, qUp = 0.95)
cage_df=CAGEr::tagClustersGR(cager2,1,qLow = 0.05,qUp = 0.95) %>%
  plyranges::join_overlap_inner_directed(., epd) %>%
  as_tibble()
write_tsv(cage_df, "Analysis/Tables/shape_merged.tsv")