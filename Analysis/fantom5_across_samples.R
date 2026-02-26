packages<-c("rtracklayer","tidyverse","BSgenome.Hsapiens.UCSC.hg38","ggridges", "patchwork")
invisible(lapply(packages, library, character.only = TRUE))
lib_gr<- import.bed("Library_data/res/library.bed")
epd<- rtracklayer::import.bed("External_data/EPD/human38_epdnew.bed") %>% promoters(downstream=16, upstream = 236)

theme_set(theme_bw())

load_cage<- function(df){
  ls=map(df$filename, ~import.bed(.x)%>% #load bed files
           plyranges::join_overlap_left_directed(lib_gr,.)%>% #overlap with the library
           values()%>%as_tibble()%>%group_by(name.x)%>%
           dplyr::summarise(counts=sum(score.y)%>%replace_na(0)))%>% #sum the counts overlapping each promoter
    list_rbind()%>%group_by(name.x)%>%dplyr::summarise(counts=sum(counts)) #join all the samples
  return(ls)
}

set<- "primary_cell"

cell_ont=readxl::read_xls(paste0("External_data/FANTOM5/",set,"_ontology_FANTOM5.xlsx")) %>% 
  dplyr::rename("sample_id"=`Sample ID`)
inputFiles<- list.files(paste0("External_data/FANTOM5/hg38_",set),
                        full.names = T,pattern=".nobarcode.ctss.bed.gz$")
sample_id<- inputFiles%>%str_extract("CNhs.{5}")
cell_ont<- tibble(filename=inputFiles, sample_id=sample_id)%>%inner_join(cell_ont, by=c("sample_id"))
cell_ont_list<- cell_ont%>%group_split(`Facet ontology term`)
names(cell_ont_list)<-group_keys(cell_ont%>%group_by(`Facet ontology term`))%>%as_vector()

libsize<- read.delim("External_data/FANTOM5/library_size_pc.txt", sep = " ", header = F)%>%mutate(sample_id=str_extract(V2, "CNhs.{5}"))
names(libsize)<- c("size", "filename", "sample_id")

cell_ont=cell_ont%>%left_join(libsize, by="sample_id")

libsize_tissues<- cell_ont%>%group_by(`Facet ontology term`)%>%summarise(libsize=sum(size, na.rm=T))

ls<- cell_ont_list%>%imap(~load_cage(.x)%>%mutate(sample=.y))%>%list_rbind()%>%dplyr::rename(name=name.x)
ls %>% write_tsv(paste0("Analysis/Tables/",set,"_CAGE_activity.tsv"))

long_ls<- ls%>%mutate(TPM=1E6*counts/libsize)

long_ls_filt<- long_ls%>%mutate(min_signal=TPM>1)%>%group_by(name)%>%mutate(sum_counts=sum(counts), n_min_signal=sum(min_signal), sum_TPM=sum(TPM), median_TPM=median(TPM),mean_TPM=mean(TPM))