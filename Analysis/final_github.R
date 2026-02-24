library(rtracklayer)
library(tidyverse)
library(patchwork)
library(ggpubr)
library(rcartocolor)

theme_set(theme_bw())

# From counts per gate to mean activity and noise
# sampling effort
perc_cells=tibble(sample=factor(1:6), perc=c(14.59,8.87,6.04,3.09, 1.41, 1.03))
perc_cells$perc=perc_cells$perc/sum(perc_cells$perc)
perc_cells$rel_factor=perc_cells$perc/min(perc_cells$perc)
perc_cells%>%ggplot(aes(sample, perc))+
  geom_col(fill="#FFB400")+xlab("Gate")+
  ylab("Sample effort relativization factor")+
  ggpubr::theme_pubclean()
ggsave("Plots/Fig1/sample_effort.pdf")

# data loading, cleaning and normalization
files<-list.files("Data_processing/Count_table","counts_ampliconseq2023_.-..tsv", 
                 full.names = T)
data<-map(files, ~read_tsv(.x, col_names = c("seq_id", "counts")))
names(data)<-files %>% str_remove("^.+2023_")%>%str_remove(".tsv$")
data<-imap(data, ~mutate(.x, counts=counts/2, sample_rep=.y))
data_long<-data %>% list_rbind()%>%
  mutate(sample=factor(as.numeric(str_remove(sample_rep, "-.$"))-1), 
         rep=factor(str_replace(sample_rep, "^.-", "Rep ") %>%
                      str_replace("3","2")),
         name=str_remove(seq_id, "^FP.{6}_"))

counts_tot<- data_long %>%
  group_by(sample,rep) %>%
  summarise(sumCounts=sum(counts)) %>% ungroup() %>% 
  group_by(rep) %>% mutate(mean_size_rep=mean(sumCounts)) %>% ungroup()
data_long<- data_long %>% 
  left_join(counts_tot) %>% 
  mutate(counts_norm_libsize=counts/(sumCounts/mean_size_rep)) %>%
  select(-mean_size_rep)
data_long<- data_long %>% 
  left_join(perc_cells, by="sample") %>%
  mutate(counts_norm=counts_norm_libsize*rel_factor) %>%
  select(-c(sumCounts,perc, sample_rep, rel_factor))


# promoter statistics
get_stats<-function(df,col_counts){
  if(!all(df[[col_counts]]==0)){
    k<-sum(df[[col_counts]]!=0)
    df$sample<-as.numeric(df$sample)
    a<-sum(df[[col_counts]])
    vctr<-map(1:length(unique(df$sample)), ~rep(df$sample[.x],times=df[[col_counts]][.x]))%>%list_c()
    mean=mean(vctr)
    median=median(vctr)
    mode<-df$sample[which.max(df[[col_counts]])]
    var<-var(vctr)
    return(data.frame(seq_id=unique(df$seq_id),name=unique(df$name),positive_gates=k,mean,median, mode,var))
  }}

stats_rep1=data_long%>%filter(rep=="Rep 1")%>%group_split(seq_id)%>%map(~get_stats(.x, "counts_norm"))%>%list_rbind() %>% mutate(rep="Rep 1")
stats_rep2=data_long%>%filter(rep=="Rep 2")%>%group_split(seq_id)%>%map(~get_stats(.x, "counts_norm"))%>%list_rbind()%>% mutate(rep="Rep 2")
stats=bind_rows(stats_rep1, stats_rep2)
stats=data_long %>% group_by(seq_id, rep) %>%
  summarise(sum_counts=sum(counts), 
            sum_norm_counts_lib=sum(counts_norm_libsize), 
            sum_norm_counts=sum(counts_norm))%>%
  left_join(stats,.,by=c("seq_id","rep"))
write_tsv(stats, "Analysis/Tables/stats.tsv")

# spike-in
spike_id=c("spikeIn_SV40", "spikeIn_CMVe", "spikeIn_CMVeMut")
spikes=data_long%>%filter(seq_id%in%spike_id) %>% 
  mutate(spike=ifelse(seq_id=="spikeIn_CMVeMut", "Spike-In 1",
                      ifelse(seq_id=="spikeIn_SV40", "Spike-In 2", "Spike-In 3")))

spikes%>%filter(spike%in%c("Spike-In 2"))%>%
  pivot_wider(names_from = "spike", values_from = "counts") %>%
  select(`Spike-In 2`, sample, rep)%>%
  left_join(spikes,.)%>%mutate(rel_spike=counts/`Spike-In 2`) %>% 
  ggplot(aes(as_factor(as.numeric(sample)), rel_spike, col=spike, group=spike)) +
  geom_point(size=2)+
  geom_smooth(method="lm", se=F, linewidth=0.5)+
  scale_y_log10(labels=scales::label_number(), breaks=10^(-2:2))+
  scale_color_manual(values=clrs[1:3])+facet_wrap(~rep)+
  ylab("Relative counts")+xlab("Sample")+labs(col="Spike") +
  theme(text=element_text(size=15))
ggsave("Plots/Fig1/spike_in.pdf")

# validation
files=list.files(path="Experimental_data/Citometry/Stable_validation/Sample Group - 1/Tables/", pattern=".cells", full.names = T, recursive=T)
prom_name=files %>% str_remove("^.+Tables/(lvs-)?") %>% str_remove("_Data .+$") #otra opcion con str_match(files, "^.+Tables/(?:lvs-)?(.*?)(_Data .+$)")[,2]
df=map2(files, prom_name, ~read_csv(.x, ) %>%
          mutate(sample_name=.y)) %>% list_rbind()
names(df)=str_replace_all(names(df), "-","_")
df$prom_name=str_remove(df$sample_name, " -.+$")

name_df=tibble(short_name=unique(df$prom_name), name=c("BTG1_1","ETS1_1",  "KIAA0753_1", "LSM1_1", "METAP2_1", "PPP1R14B_3", "TMEM87A_1","ZKSCAN2_1", "Control", "Strong" ))
df=inner_join( name_df, df, by=c("short_name"="prom_name"))
max_EGFP=max(df$`Comp_FL2_A`)
act_mean<- stats %>% group_by(name) %>%  
  summarise(mean_hist=mean(mean), 
            sd_mean_hist=sd(mean) %>% replace_na(0),
            var_hist=mean(var), 
            sd_var_hist=sd(var) %>% replace_na(0), 
            cv_hist=mean(sqrt(var)/mean), 
            sd_cv_hist=sd(sqrt(var)/mean) %>% replace_na(0), 
            fano_hist=mean(var/mean)%>% replace_na(0), 
            sd_fano_hist=sd(var/mean)%>% replace_na(0))
density_mean<-df %>% group_by(sample_name, name) %>% 
  summarise(mean_EGFP=mean(log10(Comp_FL2_A), na.rm=T)) %>% 
  ungroup() %>%group_by(name) %>% 
  summarise(mean_density=mean(mean_EGFP), 
            sd_density=sd(mean_EGFP))

density_mean %>% left_join(act_mean) %>% ggplot(aes(mean_hist, mean_density, col=name)) + 
  geom_point(size=3.5)+
  geom_errorbar(aes(ymin=mean_density-sd_density, ymax=mean_density+sd_density), width=.125)+
  geom_errorbarh(aes( xmin=mean_hist-sd_mean_hist, xmax=mean_hist+sd_mean_hist), height=.05)+
  labs(x="Mean (Histogram)", y="Mean (Density)")+
  theme_bw()+
  theme(text=element_text(size=20))+
  labs(col="Promoter")+
  scale_color_manual(values = c(KIAA0753_1="#AD343E", TMEM87A_1="#FFB400", BTG1_1="#7FB800", METAP2_1="#D6741F",PPP1R14B_3="#1B8C8E" , ZKSCAN2_1="#9EDBCB", LSM1_1="#0D2C54", ETS1_1="#25afe9"))+
  ylim(c(1.2,1.8))
ggsave("Plots/Fig1/stable_validation_mean.pdf", width = 9, height = 6.75, units="in")

df %>% filter(name%in%c( "KIAA0753_1","TMEM87A_1")) %>% 
  ggplot(aes(`Comp_FL2_A`, fill=name, group=name))+
  geom_density( alpha=0.6)+
  scale_x_log10(limits=c(1,1000))+
  ggpubr::theme_pubr()+
  labs(x="Gate", y="Relative Counts", fill="Promoter", alpha=NULL)+
  scale_fill_manual(values=c(KIAA0753_1="#AD343E", TMEM87A_1="#FFB400"))
ggsave("Plots/Fig1/stable_validation_density.pdf", width = 9, height = 6.75, units="in")

#qPCR
qpcr=readxl::read_xlsx("Experimental_data/qPCR_Gates.xlsx")
qpcr=qpcr%>%filter(Task=="Unknown", Sample!="NTC") %>%
  mutate(Cq=as.numeric(Cq), 
         Gate=str_remove(Sample,"_rep_.") %>% str_remove("Gate_") %>% as.numeric() %>% as_factor(), 
         Rep=str_remove(Sample, "Gate_._"))
qpcr$tech=rep(c(1,2,3), times=7*3*2)
qpcr_wide=qpcr%>%select(Cq,Target,Sample, Rep,tech,Gate) %>%
  pivot_wider(values_from=Cq, names_from=Target) %>%
  mutate(delta=EGFP-GAPDH)

df<- qpcr_wide %>% filter(Rep!="rep_1", Gate!=1) %>%
  group_by(Rep,Gate) %>%
  mutate(mean_GAPDH=mean(GAPDH), 
         delta=2^-(EGFP-mean_GAPDH))
mean_g2<- df %>% filter(Gate==2) %>%
  select(Rep,delta)%>%group_by(Rep)%>%
  summarise(mean_g2=mean(delta))  
df<- df %>% left_join(mean_g2) %>% 
  mutate(FC=delta/mean_g2) %>%
  mutate(Gate=as.factor(as.numeric(Gate)-1))
p2<- df%>%filter(Rep=="rep_2")%>%
  ggplot(aes(Gate, FC, fill=Gate, col=Gate))+
  geom_line(linewidth=2)+geom_point(size=5)+
  geom_point(shape = 1,size=5,colour = "black")+
  ggtitle("Rep 1")+scale_y_log10()+geom_point()+
  ylab("EGFP RNA (Fold change vs Gate 1)")+
  rcartocolor::scale_color_carto_d(palette="Emrld")+
  ggpubr::theme_pubr()+theme(legend.position = "none")
p3=df %>% filter(Rep=="rep_3") %>%
  ggplot(aes(Gate, FC, fill=Gate, col=Gate))+
  geom_line(linewidth=2)+geom_point(size=5)+
  geom_point(shape = 1,size=5,colour = "black")+
  ggtitle("Rep 2")+
  scale_y_log10()+
  ylab("")+
  rcartocolor::scale_color_carto_d(palette="Emrld")+
  ggpubr::theme_pubr()+theme(legend.position = "none")
ggsave("Plots/Fig1/qpcr.pdf", plot=cowplot::plot_grid(p2,p3))


# Filters 
stats<- stats%>%filter(sum_counts>=100)
inconsistent_median<- stats  %>% select(seq_id, median, rep) %>% 
  pivot_wider(names_from = rep, values_from = median) %>%
  mutate(median_dif=abs(as.numeric(`Rep 1`)-as.numeric(`Rep 2`))) %>% 
  filter(median_dif>1) %>% select(seq_id) %>% as_vector()
stats<- stats %>% filter(!seq_id%in%inconsistent_median)

max_var=function(k){
  if(k%%2==1){
    return((sum(seq(0, (k-1)/2)^2)*2)/k)
  }
  if(k%%2==0){
    return((sum((seq(1, k/2)-0.5)^2)*2)/k)
  }
}
stats$max_var=sapply(stats$positive_gates,max_var)
stats %>% mutate(bimodal=var>max_var*1.2) %>%
  arrange(bimodal)%>% ggplot(aes(mean, var, col=bimodal))+
  geom_point(size=0.5)+facet_wrap(~rep)+ylab("Variance")+
  xlab("Mean")+ggpubr::theme_pubr()+
  scale_color_manual(values=c( "FALSE"="#0D2C54", "TRUE"="#1B8C8E"))+
  labs(col="Bimodal")+
  guides(colour = guide_legend(override.aes = list(size=4)))+
  theme(text = element_text(size=20))
ggsave("Plots/Fig4/mean_var_bimodal.pdf", width = 9, height = 6.75, units="in")
stats<- stats%>%filter(var<=max_var*1.2)

# Mean correlation
p<- stats %>% select(seq_id, mean, rep) %>% 
  pivot_wider(values_from = "mean", names_from = "rep") %>% 
  ggplot(aes(`Rep 1`, `Rep 2`))+
  geom_point(size=0.2, alpha=0.5, col="#0D2C54")+
  labs(x="Mean activity (Rep 1)", y="Mean activity (Rep 2)")
ggExtra::ggMarginal(p, type="density", margins = "both", col="#0D2C54")
ggsave("Plots/Fig1/mean_replicates.pdf")
stats %>% select(seq_id, mean, rep) %>% 
  pivot_wider(values_from = "mean", names_from = "rep") %>% 
  summarize(corre=cor(`Rep 1`, `Rep 2`, use="complete.obs", method="pearson"))

# activity analysis
prom_df=read_tsv("Analysis/Tables/prom_df.tsv") %>% filter(type=="promoter")
data=inner_join(stats, prom_df) 
excess=data %>% group_by(rep) %>% summarise(excess=n()%%100)
data=data %>% ungroup() %>% 
  group_split(rep) %>% map2(excess$excess, ~.x %>% 
                              mutate(mean_rank_sw=row_number(mean)) %>% 
                              filter(mean_rank_sw>.y) %>% 
                              mutate(mean_rank_sw=row_number(mean), 
                                     mean_sw=factor(ceiling(mean_rank_sw/100))) %>% 
                              group_by(rep, mean_sw) %>% mutate(var_rank_sw=row_number(var)) %>% 
                              ungroup()) %>% bind_rows()

#TATA activity
data %>% group_by(mean_sw,rep) %>% 
  summarise(TATA_EPD=sum(TATA_EPD)/n()) %>% 
  ggplot(aes(mean_sw %>% as.numeric(),TATA_EPD))+
  geom_point(col="#14AFB2")+geom_smooth(col="#216869")+
  theme(text=element_text(size=20))+facet_wrap(~rep)+ylim(0,1)+
  labs(x="Mean activity (binned by rank)",y="Proportion of promoters")+
  ggtitle("TATA-box")+theme_bw()
ggsave("Plots/Fig2/TATA_scatter.pdf", width = 9, height = 6.75, units="in")

#CGI activity
data %>% group_by(mean_sw,rep) %>% 
  summarise(CGI=sum(CGI)/n()) %>% 
  ggplot(aes(mean_sw %>% as.numeric(),CGI))+
  geom_point(col="#14AFB2")+geom_smooth(col="#216869")+
  theme(text=element_text(size=20))+facet_wrap(~rep)+ylim(0,1)+
  labs(x="Mean activity (binned by rank)",y="Proportion of promoters")+
  ggtitle("CpG island overlap")
ggsave("Plots/Fig2/CGI_scatter.pdf", width = 9, height = 6.75, units="in")

p=ggviolin(data %>% filter(rep=="Rep 1"), 
           x="CGI", y="mean", fill="CGI", draw_quantiles = 0.5, 
           add="median_q1q3", palette=c("#14AFB2", "#216869"), alpha=0.7)
p+stat_compare_means()+
  theme_pubclean()+labs(fill="Cpg island", x=NULL, y="Transcriptional mean activity")+
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(), 
        text=element_text(size=20))+
  ylim(c(0.5,6.2))
ggsave("Plots/Fig2/CGI_violin_rep1.pdf", width = 9, height = 6.75, units="in")
p<- ggviolin(data %>% filter(rep=="Rep 2"), 
             x="CGI", y="mean", fill="CGI", draw_quantiles = 0.5, 
             add="median_q1q3", palette=c("#14AFB2", "#216869"), alpha=0.7)
p + stat_compare_means()+
  theme_pubclean()+labs(fill="CpG isalnd", x=NULL, y="Transcriptional mean activity")+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), 
        text=element_text(size=20))+
  ylim(c(0.5,6.2))
ggsave("Plots/Fig2/CGI_violin_rep2.pdf", width = 9, height = 6.75, units="in")

#tissue specificity
data %>% filter(!is.na(sample_specificity_gini), type=="promoter") %>%
  group_by(rep) %>% arrange(rep, mean) %>%  
  mutate(mean_sw=row_number(mean) %>% cut_width(100) %>% as.numeric(), 
         high_sp=(cut_number(sample_specificity_gini, n=3) %>% as.numeric())==3) %>% 
  group_by(mean_sw, rep) %>% summarise(high_sp=sum(high_sp)/n()) %>%
  ggplot(aes(mean_sw,high_sp))+geom_point(col="#14AFB2")+
  geom_smooth(col="#216869")+theme(text=element_text(size=20))+
  facet_wrap(~rep)+ylim(0,1)+
  labs(x="Mean activity (binned by rank)",y="Proportion of promoters")+
  ggtitle("High tissue specificity")
ggsave("Plots/Fig2/high_tissue_sp_scatter.pdf", width = 9, height = 6.75, units="in")

df=data %>% filter(!is.na(sample_specificity_gini), type=="promoter") %>%group_by(rep) %>%  mutate(high_sp=(cut_number(sample_specificity_gini, n=3) %>% as.numeric())==3)

p=ggviolin(df %>% filter(rep=="Rep 1"), 
           x="high_sp", y="mean", fill="high_sp", draw_quantiles = 0.5, 
           add="median_q1q3", palette=c("#14AFB2", "#216869"), alpha=0.7)
p+stat_compare_means()+theme_pubclean()+
  labs(fill="High tissue specificity", x=NULL, y="Transcriptional mean activity")+
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),text=element_text(size=20))+ylim(c(0.5,6.2))
ggsave("Plots/Fig2/high_tissue_sp_violin_rep1.pdf", width = 9, height = 6.75, units="in")
p=ggviolin(df %>% filter(rep=="Rep 2"), 
           x="high_sp", y="mean", fill="high_sp", draw_quantiles = 0.5, 
           add="median_q1q3", palette=c("#14AFB2", "#216869"), alpha=0.7)
p+stat_compare_means()+theme_pubclean()+
  labs(fill="High tissue specificity", x=NULL, y="Transcriptional mean activity")+
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), text=element_text(size=20))+
  ylim(c(0.5,6.2))
ggsave("Plots/Fig2/high_tissue_sp_violin_rep2.pdf", width = 9, height = 6.75, units="in")

# activity summary
exclude=c("mean", "mean_sw", "var_rank_sw", "rep")

tidy_data=data %>% fastDummies::dummy_cols("TE_superclass", ignore_na = T, omit_colname_prefix=T, remove_selected_columns = F) %>%
  fastDummies::dummy_cols("sample_specificity_class", ignore_na = T, omit_colname_prefix=T, remove_selected_columns = T) %>% 
  group_by(rep) %>% 
  mutate(`High G+C content`=(cut_number(g_c, n=3) %>% as.numeric())==3,
         across(c(DNA, SINE, LINE,LTR), ~replace_na(.x, 0)),
         `No activity in mouse`=turnover%in%c("expression-turnover", "mouse-dimished"),
         `Deleted in mouse`=turnover=="mouse-deleted",
         `Human inserted`=turnover=="human-inserted",
         `High TSS predictability`=(cut_number(`Sum of total effect`, n=3) %>% as.numeric())==3,
         `Transposable elements`=!is.na(TE_superclass),
         `Low Complexity Repeats`=LCR_overlap>0,
         `High conservation (16 to -50pb)`=(cut_number(phylop100_close, n=3) %>% as.numeric())==3,
         `High conservation (-50 to -150)`=(cut_number(phylop100_intermediate, n=3) %>% as.numeric())==3,,
         `High conservation (-150 to -235)`=(cut_number(phylop100_far, n=3) %>% as.numeric())==3,
         `High tissue specificity`=(cut_number(sample_specificity_gini, n=3) %>% as.numeric())==3,
         `Low tissue specificity`=(cut_number(sample_specificity_gini, n=3) %>% as.numeric())==1,
         shape_class= cut_number(interquantile_width, 3) %>% as_factor() %>% as.numeric(),
         `Narrow promoters`=shape_class==1,
         `Broad promoters`=shape_class==3,
         `High activity in HEK293`=hek_tpm>median(data$hek_tpm[data$hek_tpm>0]),
         `None Cis Regulatory Module`=N_TF_CRM==0,
         `High chromatin accesibility (DNase-seq)`=(cut_number(mean_dnase, n=3) %>% as.numeric())==3,
         across(starts_with("enh"), ~.x>0)
  ) %>%
  dplyr::rename(
    `TATA-box`=TATA_EPD,
    `CCAAT` = CCAAT_EPD,
    `GC-box` = GCbox_EPD,
    `CpG islands`=CGI,
    Retrotransposon=DNA,
    TCT=TCT_TSS,
    `CG at TSS`=CG_TSS,
    `TA at TSS`=TA_TSS,
    `TG at TSS`=TG_TSS,
    `CA at TSS`=CA_TSS,
    `Non canonical TSS`=other_TSS,
    `Strong TSS`=INR_strong_TSS,
    `Ubiquitous (FANTOM5)`=all_samples_detected,
    `Sample-specific (FANTOM5)` = group_enrichment,
    `Sample-enhanced`=sample_enhanced,
    `Non detected activity (FANTOM5)`= non_detected, 
    `Enhancers at 10kb window`=enh10kb,
    `Enhancers at 50kb window`=enh50kb,
    `Enhancers at 100kb window`=enh100kb
  ) %>% 
  select(mean, rep,
         LINE,
         SINE,
         LTR,
         `High G+C content`,
         `No activity in mouse`,
         `Human inserted`,
         `Transposable elements`,
         `Low Complexity Repeats`,
         `High conservation (16 to -50pb)`,
         `High conservation (-50 to -150)`,
         `High conservation (-150 to -235)`,
         `TATA-box`,
         `CCAAT`,
         `GC-box`,
         `CpG islands`,
         Retrotransposon,
         TCT,
         `CG at TSS`,
         `TA at TSS`,
         `TG at TSS`,
         `CA at TSS`,
         `Non canonical TSS`,
         `Strong TSS`,
         `Non detected activity (FANTOM5)`,
         `High tissue specificity`,
         `Low tissue specificity`,
         `Narrow promoters`,
         `Broad promoters`,
         `High activity in HEK293`,
         `High chromatin accesibility (DNase-seq)`,
         `Enhancers at 50kb window`,
         `None Cis Regulatory Module`)%>%
  
  mutate(across(c(where(~ is.numeric(.x)), -contains("mean"),-contains("var")), as.logical), rep=as.factor(rep), across(where(is.logical), ~.x%>% factor(level=c("TRUE","FALSE")))) %>% ungroup()
tidy_data %>% write_tsv("Analysis/Tables/tidy_data.tsv")
tidy_names=read_tsv( "Analysis/Tables/tidy_names.tsv")
n_features=tidy_data %>% select(-mean) %>% group_by(rep) %>% summarise(across(everything(), function(x)as.logical(x) %>% sum(na.rm=T) %>% magrittr::divide_by(n()))) %>% pivot_longer(-rep, names_to = "feature", values_to = "prop")
vbles_split=map(c("seq","endo"), ~tidy_names$feature[(tidy_names$group==.x %>% na.exclude())&(tidy_names$feature%in%names(tidy_data))])

wilcox=map(tidy_data %>% select(-rep,-contains("mean"),-contains("var")), ~coin::wilcox_test(formula=mean~.x|rep, data=tidy_data) %>% pvalue())
wilcox=tibble(feature=names(wilcox), pval=list_c(wilcox), pval_corr=p.adjust(pval, "BH", length(wilcox)))
wilcox_rep1=map(tidy_data %>% filter(rep=="Rep 1")%>%  select(-mean,-rep), ~coin::wilcox_test(formula=mean~.x, data=tidy_data%>% filter(rep=="Rep 1"), conf.int=T))
wilcox_rep2=map(tidy_data %>% filter(rep=="Rep 2")%>%  select(-mean,-rep), ~coin::wilcox_test(formula=mean~.x, data=tidy_data%>% filter(rep=="Rep 2"), conf.int=T))

wilcox_df2<- map2(list(wilcox_rep1, wilcox_rep2),c("Rep 1","Rep 2"),
                ~tibble(feature=names(.x),estimate=map(.x,~confint(.x)$estimate) %>% list_c(), 
                        P2.5=map(.x,~confint(.x)$conf.int[1]) %>% list_c(),P97.5=map(.x,~confint(.x)$conf.int[2]) %>% list_c() ,
                        rep=.y)%>% 
                  pivot_longer(c(starts_with("estimate"),starts_with("P2.5"), starts_with("P97.5")),names_to = "val", values_to="estimate") %>% arrange(desc(estimate)) %>% mutate(feature=fct_inorder(feature))) %>% list_rbind() %>% left_join(wilcox)
output_file<-c("Plots/Fig2/summary_seq.pdf", "Plots/Fig2/summary_endo.pdf")
for(i in 1:2){
  wilcox_df2 %>% filter(feature%in% vbles_split[[i]], pval_corr<0.05) %>% 
    mutate(act=ifelse(sign(estimate)==1,"High activity", "Low activity")) %>% 
    group_by(feature) %>%
    mutate(n=length(unique(act))) %>% filter(n==1) %>%
    arrange(estimate, desc(rep)) %>% 
    pivot_wider(names_from = val, values_from = estimate) %>% 
    filter(pval_corr<0.05) %>%
    ggplot(aes(x=estimate, y=fct_inorder(feature), group=fct_inorder(rep)))+
    geom_col(orientation="y", position="dodge", aes( fill=act, alpha=rep))+
    scale_alpha_manual(values=c("Rep 1"=1, "Rep 2"=0.7))+
    theme(text=element_text(size=20), legend.position = "top")+
    labs(fill="Effect", x="Activity Effect", y="Features", alpha="")+
    geom_errorbarh(aes(xmax=P2.5 , xmin=P97.5), position=position_dodge(1), height=0.05, col="#777777",  linewidth=1.5)+
    xlim(c(-0.3, 0.55))+
    scale_fill_manual(values=c("High activity"="#D6741F", "Low activity"="#7FB800"))+
    guides(alpha="none")+theme_pubr()
  #ggsave(output_file[i])
  }

#hek
data %>% mutate(as_factor(mean_sw)) %>% group_by(mean_sw, rep) %>% 
  summarise(hek_tpm=mean(hek_tpm)) %>% 
  ggplot(aes(mean_sw %>% as.numeric(), hek_tpm))+
  geom_point(col="#AD343E")+
  geom_smooth(method="lm", col="#216869")+
  facet_wrap(~rep)+labs(x="Mean activity (binned by rank)", 
                        y="HEK293 endogenous activity (mean TPM)")+
  theme_bw()+theme(text=element_text(size=20))+scale_y_log10()
ggsave("Plots/Fig3/HEK_mean.pdf", width = 9, height = 6.75, units="in")

data %>% group_split(rep) %>% map(~cor(.x$mean,.x$hek_tpm, method="spearman" ))
data %>% group_split(rep) %>% map(~cor.test(.x$mean,.x$hek_tpm, method="spearman" ))
data %>% filter(!is.na(gini_tissue)) %>% mutate(tissue_sp=cut_number(gini_tissue, n=3)) %>% 
  group_split(rep, tissue_sp) %>% map(~cor.test(.x$mean,.x$hek_tpm, method="spearman" ))

data %>% mutate(as_factor(mean_sw)) %>% 
  mutate(mean_sw=row_number(mean) %>% cut_number(20) %>% as.numeric()) %>%
  ggplot(aes(factor(mean_sw), hek_tpm+.1))+
  geom_boxplot(fill="#AD343E", alpha=0.5)+
  #geom_smooth(method="lm", col="#216869")+
  facet_wrap(~rep)+labs(x="Mean activity (binned by rank)", 
                        y="HEK293 endogenous activity (TPM)")+
  theme_bw()+theme(text=element_text(size=20))+scale_y_log10()
ggsave("Plots/Fig3/HEK_mean_boxplot.pdf", width = 9, height = 6.75, units="in")

data %>% filter(!is.na(gini_tissue)) %>% 
  mutate(sample_specificity_gini=cut_number(gini_tissue, n=3) %>% as.numeric() %>% as_factor(),
         sample_specificity_gini=ifelse(sample_specificity_gini==1, "Housekeeping", 
                                        ifelse(sample_specificity_gini==2, "Intermediate", "Tissue-specific"))) %>% 
  filter(sample_specificity_gini!="Intermediate") %>% 
  mutate(hek_tpm=replace_na(hek_tpm,0)) %>% 
  group_by(rep, sample_specificity_gini) %>% 
  mutate(mean_sw=row_number(mean) %>% cut_number(15) %>% as.numeric() %>% as_factor()) %>%
  ggplot(aes(mean_sw, hek_tpm+.1, fill=sample_specificity_gini))+
  geom_boxplot( alpha=0.5, outliers = T)+
  facet_wrap(~rep)+
  labs(x="Mean reporter activity (binned)", 
       y="Endogenous activity (TPM)", fill="")+
  theme_bw()+
  theme(text=element_text(size=20))+
  scale_y_log10()+
  ggtitle("HEK293")+
  scale_fill_manual(values=c("Housekeeping"="#1B8C8E", "Tissue-specific"="#0D2C54"))#+theme(legend.position = "none")
ggsave("Plots/Fig3/HEK_mean_boxplot_by_tissue_sp.pdf", width = 9, height = 6.75, units="in")


#PUFFIN
data %>% group_by(rep) %>% filter(!is.na(selectivity_dudnik)) %>%
  mutate(mean_sw=row_number(mean) %>% cut_number(100) %>% as.numeric()) %>%
  group_by(rep, mean_sw)%>%  summarise(puffin=median(selectivity_dudnik, na.rm=T)) %>% 
  ggplot(aes(mean_sw,puffin))+geom_point(col="#14AFB2")+
  geom_smooth(col="#216869", method="lm")+
  facet_wrap(~rep)+
  labs(x="Mean activity (binned by rank)",y="Selectivity index")+
  ggtitle("PUFFIN selectivity")+theme_bw()+
  theme(text=element_text(size=20))
ggsave("Plots/Fig3/PUFFIN_selectivity_scatter.pdf", width = 9, height = 6.75, units="in")



data %>% group_by(rep) %>% filter(!is.na(selectivity_dudnik)) %>%
  mutate(mean_sw=row_number(mean) %>% cut_number(20) %>% as.numeric()) %>%
  ggplot(aes(mean_sw,selectivity_dudnik, group=mean_sw))+
  geom_boxplot(alpha=0.6, fill="#216869", outliers = F)+
  facet_wrap(~rep, scale="free_x")+
  labs(x="Mean activity (binned by rank)",y="Selectivity index")+
  ggtitle("PUFFIN selectivity")+theme_bw()+
  theme(text=element_text(size=20))
ggsave("Plots/Fig3/PUFFIN_selectivity_boxplot.pdf", width = 9, height = 6.75, units="in")


data %>% group_by(mean_sw,rep) %>% summarise(puffin=median(teoTSS_score)) %>%
  ggplot(aes(mean_sw %>% as.numeric(),puffin))+
  geom_point(col="#14AFB2")+geom_smooth(col="#216869")+
  theme(text=element_text(size=20))+facet_wrap(~rep)+
  labs(x="Mean activity (binned by rank)",y="Median prediction score")+
  ggtitle("PUFFIN prediction over TSS position")+theme_bw()+
  theme(text=element_text(size=20))
ggsave("Plots/Fig3/PUFFIN_predscore_scatter.pdf", width = 9, height = 6.75, units="in")

data %>% filter(!is.na(teoTSS_score)) %>% 
  mutate(mean_sw=row_number(mean) %>% cut_number(20) %>% as.numeric()) %>%
  ggplot(aes(mean_sw,teoTSS_score, group=mean_sw))+
  geom_boxplot(alpha=0.6, fill="#216869", outliers = F)+
  facet_wrap(~rep, scale="free_x")+
  labs(x="Mean activity (binned by rank)",y="Prediction score")+
  ggtitle("PUFFIN prediction over TSS position")+theme_bw()+
  theme(text=element_text(size=20))
ggsave("Plots/Fig3/PUFFIN_predscore_boxplot.pdf", width = 9, height = 6.75, units="in")


## Noise
library(pROC)
#CGI
data %>% ggplot(aes(as.numeric(mean_sw), var_rank_sw, col=CGI))+
  geom_point(size=0.2, alpha=0.5)+facet_wrap(~rep)+
  xlab("Mean rank group")+ylab("Variance rank (SW)")+
  labs(col="CpG island")+
  scale_color_manual(values = c("grey", "#1B8C8E"),
                     guide = guide_legend(override.aes = list(size = 3)))+
  ggnewscale::new_scale_color() +
  geom_smooth(se=F, aes(color=CGI), show.legend = F)+
  scale_color_manual(values = c("#525252", "#136869ff"))+
  theme_pubr(base_size = 20)+labs(col="")
ggsave("Plots/Fig4/noise_scatter_CGI.pdf")

data %>% ungroup() %>% group_split(rep) %>% 
  map(~roc_curve(.x, "CGI") %>% 
        mutate(rep=unique(.x$rep), AUC=auc(.))) %>% 
  list_rbind()%>% 
  ggplot(aes(FPR,TPR, col=rep))+
  geom_line()+
  geom_abline(linetype="dashed")+
  geom_text(data=.%>% select(AUC, rep) %>% 
              mutate(AUC=paste("AUC:", round(AUC,3))) %>% 
              unique() %>% bind_cols(tibble(TPR=c(0.9,0.8), FPR=0.02)), 
            aes(label=AUC), hjust="left")+ggtitle("CGI ROC curve")+
  scale_color_manual(values=c("#0D2C54", "#AD343E"))+
  theme_pubr(base_size = 20)+labs(linetype="", col="")
ggsave("Plots/Fig4/ROC_CGI.pdf")

# TATA
data %>%filter(type=="promoter") %>% 
  ggplot(aes(as.numeric(mean_sw), var_rank_sw, col=TATA_EPD))+
  geom_point(size=0.2, alpha=0.5)+facet_wrap(~rep)+
  xlab("Mean rank group")+ylab("Variance rank (SW)")+
  scale_color_manual(values = c("grey", "#1B8C8E"),
                     guide = guide_legend(override.aes = list(size = 3)))+
  ggnewscale::new_scale_color() +
  geom_smooth(se=F, aes(color=TATA_EPD), show.legend = F)+
  scale_color_manual(values = c("#525252", "#136869ff"))+
  theme_pubr(base_size = 20)
ggsave("Plots/Fig4/noise_scatter_TATA.pdf")

data  %>%  ungroup() %>% group_split(rep) %>% 
  map(~roc_curve(.x, "TATA_EPD") %>% 
        mutate(rep=unique(.x$rep), AUC=auc(.))) %>% 
  list_rbind()%>% 
  ggplot(aes(FPR,TPR, col=rep))+
  geom_line()+geom_abline(linetype="dashed")+
  geom_text(data=.%>% select(AUC, rep) %>% 
              mutate(AUC=paste("AUC:", round(AUC,3)))%>%
              unique() %>% bind_cols(tibble(TPR=c(0.9,0.8), FPR=0.02)), 
            aes(label=AUC), hjust="left")+ggtitle("TATA-box")+
  scale_color_manual(values=c("#0D2C54", "#AD343E"))+
  theme_pubr(base_size = 20)+labs(linetype="", col="")
ggsave("Plots/Fig4/ROC_TATA.pdf", create.dir = T)

data  %>%  mutate(expr_level=ifelse(as.numeric(mean_sw)<(0.75*max(as.numeric(mean_sw))), "Low expression", "High expression")) %>% 
  ungroup() %>% group_split(rep, expr_level) %>% 
  map(~roc_curve(.x, "TATA_EPD") %>% 
        mutate(rep=unique(.x$rep), expr_level=unique(.x$expr_level),AUC=auc(.))) %>% 
  list_rbind()%>% 
  ggplot(aes(FPR,TPR, col=rep, linetype=expr_level))+
  geom_line()+geom_abline(linetype="dashed")+
  geom_text(data=.%>% select(AUC, rep,expr_level) %>% 
              mutate(AUC=paste("AUC:", round(AUC,3)))%>% unique() %>%
              bind_cols(tibble(TPR=c(0.9,0.8, 0.7,0.6), FPR=0.02)), 
            aes(label=AUC), hjust="left")+ggtitle("TATA-box") +
  scale_color_manual(values=c("#14AFB2", "#216869")) + theme_pubr() +
  labs(linetype="", col="")
ggsave("Plots/Fig4/ROC_TATA_expr.pdf", create.dir = T, height=5, width=5, units="in")

# summary (try to reuse the one for activity!)
tidy_data=data %>% 
  fastDummies::dummy_cols("TE_superclass", ignore_na = T, omit_colname_prefix=T, remove_selected_columns = F) %>% 
  fastDummies::dummy_cols("sample_specificity_class", ignore_na = T, omit_colname_prefix=T, remove_selected_columns = T) %>% 
  group_by(rep) %>% 
  mutate(`High G+C content`=(cut_number(g_c, n=3) %>% as.numeric())==3,
         across(c(DNA, SINE, LINE,LTR), ~replace_na(.x, 0)),
         `No activity in mouse`=turnover%in%c("expression-turnover", "mouse-dimished"),
         `Deleted in mouse`=turnover=="mouse-deleted",
         `Human inserted`=turnover=="human-inserted",
         `High TSS predictability`=(cut_number(`Sum of total effect`, n=3) %>% as.numeric())==3,
         `Transposable elements`=!is.na(TE_superclass),
         `Low Complexity Repeats`=LCR_overlap>0,
         `High conservation (16 to -50pb)`=(cut_number(phylop100_close, n=3) %>% as.numeric())==3,
         `High conservation (-50 to -150)`=(cut_number(phylop100_intermediate, n=3) %>% as.numeric())==3,,
         `High conservation (-150 to -235)`=(cut_number(phylop100_far, n=3) %>% as.numeric())==3,
         `High tissue specificity`=(cut_number(sample_specificity_gini, n=3) %>% as.numeric())==3,
         `Low tissue specificity`=(cut_number(sample_specificity_gini, n=3) %>% as.numeric())==1,
         `Narrow promoters`=shape_class=="Narrow",
         `Broad promoters`=shape_class=="Broad",
         `High activity in HEK293`=hek_tpm>median(data$hek_tpm[data$hek_tpm>0]),
         `None Cis Regulatory Module`=N_TF_CRM==0,
         `High chromatin accesibility (DNase-seq)`=(cut_number(mean_dnase, n=3) %>% as.numeric())==3,
         across(starts_with("enh"), ~.x>0), 
  ) %>%
  dplyr::rename(
    `TATA-box`=TATA_EPD,
    `CCAAT` = CCAAT_EPD,
    `GC-box` = GCbox_EPD,
    `CpG islands`=CGI,
    Retrotransposon=DNA,
    TCT=TCT_TSS,
    `CG at TSS`=CG_TSS,
    `TA at TSS`=TA_TSS,
    `TG at TSS`=TG_TSS,
    `CA at TSS`=CA_TSS,
    `Non canonical TSS`=other_TSS,
    `Strong TSS`=INR_strong_TSS,
    #`Ubiquitous (FANTOM5)`=all_samples_detected,
    `Sample-specific (FANTOM5)` = group_enrichment,
    `Sample-enhanced`=sample_enhanced,
    `Non detected activity (FANTOM5)`= non_detected, 
    `Enhancers at 10kb window`=enh10kb,
    `Enhancers at 50kb window`=enh50kb,
    `Enhancers at 100kb window`=enh100kb
  ) %>% 
  select(mean, rep,mean_sw, var_rank_sw,
         LINE,
         SINE,
         LTR,
         `High G+C content`,
         `No activity in mouse`,
         `Human inserted`,
         `Transposable elements`,
         `Low Complexity Repeats`,
         `High conservation (16 to -50pb)`,
         `High conservation (-50 to -150)`,
         `High conservation (-150 to -235)`,
         `TATA-box`,
         `CCAAT`,
         `GC-box`,
         `CpG islands`,
         Retrotransposon,
         TCT,
         `CG at TSS`,
         `TA at TSS`,
         `TG at TSS`,
         `CA at TSS`,
         `Non canonical TSS`,
         `Strong TSS`,
         `Non detected activity (FANTOM5)`,
         `High tissue specificity`,
         `Low tissue specificity`,
         `Narrow promoters`,
         `Broad promoters`,
         `High activity in HEK293`,
         `None Cis Regulatory Module`,
         `High chromatin accesibility (DNase-seq)`,
         `Enhancers at 50kb window`
  )%>%
  mutate(across(c(where(~ is.numeric(.x)), 
                  -c(mean, mean_sw, var_rank_sw)), as.logical), 
         rep=as.factor(rep), 
         across(where(is.logical), ~.x%>% factor(level=c("TRUE","FALSE")))) %>% 
  ungroup()

tidy_names=read_tsv( "Analysis/Tables/tidy_names.tsv")
n_features=tidy_data %>% select(-mean) %>% group_by(rep) %>% summarise(across(everything(), function(x)as.logical(x) %>% sum(na.rm=T) %>% magrittr::divide_by(n()))) %>% pivot_longer(-rep, names_to = "feature", values_to = "prop")
vbles_split=map(c("seq","endo"), ~tidy_names$feature[(tidy_names$group==.x %>% na.exclude())&(tidy_names$feature%in%names(tidy_data))])
vbles=tidy_data  %>% select( -mean_sw, -var_rank_sw, -rep, -mean) %>% names()
vbles_split=map(c("seq","endo"), ~tidy_names$new_name[tidy_names$choose&(tidy_names$group==.x)]) 

auc_l=list()
for(j in unique(tidy_data$rep)){
  df=tidy_data %>% filter(rep==j)
  for(i in vbles){
    auc_l[[j]][[i]]=roc(df[[i]],df[["var_rank_sw"]], ci=T, boot.n=2000, direction=">")
  }
}
auc_df=list()
for (i in c(1,2)){
  auc_df[[i]]=map(auc_l[[i]], ~tibble(rep=unique(tidy_data$rep)[i], ci2.5=(.x$ci %>% as.numeric)[1], ci97.5=(.x$ci %>% as.numeric)[3], AUC=(.x$ci %>% as.numeric)[2])) %>% 
    list_rbind() %>% mutate(feature=names(auc_l[[1]]))
}
auc_df=list_rbind(auc_df)


plot_auc2=function(df){
  df %>%  arrange(desc(rep), AUC) %>% mutate(noise=ifelse(AUC>0.5, "High noise", "Low noise")) %>% filter((ci2.5>0.5&ci97.5>0.5)|(ci2.5<0.5&ci97.5<0.5)) %>% add_count(feature) %>% filter(n==2) %>% 
    ggplot(aes(x=AUC-0.5, y=fct_inorder(feature), group=fct_inorder(rep))) +
    geom_col(orientation="y", position="dodge", aes( fill=noise, alpha=rep))+
    scale_alpha_manual(values=c("Rep 1"=1, "Rep 2"=0.7))+
    theme(text=element_text(size=20), legend.position = "top")+labs(fill="Effect", x="AUC (Noise effect)", y="Features", alpha="")+
    scale_x_continuous(labels = function(x) x + 0.5)+
    geom_errorbarh(aes(xmax=abs(ci2.5)-0.5 , xmin=abs(ci97.5)-0.5), position=position_dodge(1), height=0.05, col="#777777",  linewidth=1.5)+
    scale_fill_manual(values=c("High noise"="#D6741F", "Low noise"="#7FB800"))+scale_alpha_manual(values = c("Rep 1"=1, "Rep 2"=0.5))
}+guides(alpha=F)
plot_auc2(auc_df)
ggsave("Plots/Fig4/summary_noise.pdf", width = 13.5, height = 6.75, units="in")

plot_auc2(auc_df %>% filter(feature%in%vbles_split[[1]]))
ggsave("Plots/Fig4/summary_noise_seq.pdf", width = 13.5, height = 6.75, units="in")

plot_auc2(auc_df %>% filter(feature%in%vbles_split[[2]]))
ggsave("Plots/Fig4/summary_noise_endo.pdf", width = 13.5, height = 6.75, units="in")

# Tissue specificity:
# Max CAGE expression along samples vs. exogenous activity (to better consider tissue specific promoters)
tissue=read_tsv("Analysis/Tables/tissue_CAGE_activity.tsv") %>% mutate(group="tissue")
cell=read_tsv("Analysis/Tables/primary_cell_CAGE_activity.tsv")%>% mutate(group="primary_cell")
fantom_max_expresion=rbind(tissue ,cell) %>% mutate(tpm=1E6*counts/libsize) %>% 
  arrange(desc(tpm)) %>% group_by(name) %>%
  mutate(order=row_number(desc(tpm))) %>%
  filter(order==1) %>% 
  left_join(data,., by=c("seq_id"="name")) %>%
  mutate(tpm=replace_na(tpm, 0)) 

fantom_max_expresion %>% group_split(rep) %>%
  map(~cor(.x$mean,.x$tpm, method="pearson" ))

#check how to be the x-axis numeric not factor, while grouping for fill color
fantom_max_expresion %>% 
  filter(!is.na(gini_tissue)) %>% 
  group_by(rep) %>% 
  mutate(sample_specificity_gini=cut_number(gini_tissue, n=3) %>% as.numeric()) %>% 
  group_by(rep, sample_specificity_gini) %>% 
  mutate(mean=row_number(mean) %>% cut_number(15) %>% as.numeric()) %>%
  filter( sample_specificity_gini!=2) %>%
  mutate(`Expression pattern`=ifelse(sample_specificity_gini==1, "Housekeeping", "Tissue-specific")) %>% 
  ggplot(aes(factor(mean), tpm, fill=`Expression pattern`))+
  geom_boxplot( outliers = F, alpha=0.6)+
  facet_wrap(~rep, scale="free_x")+
  scale_y_log10()+
  labs(x="Ranked groups of mean exogenous activity (bin size=200)",
       y="TPM (max value along CAGE FANTOM5)")+
  scale_fill_manual(values=c("#1B8C8E", "#0D2C54"))+
  theme_bw()+
  theme(text=element_text(size=15), legend.position = "top")
ggsave("Plots/Fig3/tissue_specificty_boxplot.pdf", width = 9, height = 6.75, units="in")

lm_l<- fantom_max_expresion %>%  filter(!is.na(gini_tissue)) %>% 
  mutate(sample_specificity_gini=cut_number(gini_tissue, n=3) %>% as.numeric()) %>% 
  group_by(rep, sample_specificity_gini) %>% 
  mutate(mean=dense_rank(mean) , 
         tpm=log10(tpm+0.1),
         group=paste(rep, sample_specificity_gini, sep="_")) %>%  
  split(.$group) %>%  
  map(~lm(tpm~mean, data=.x)) 

lm_l %>% imap(~tibble(residuals=residuals(.x), group=.y) %>% separate(group, into=c("rep", "Expression pattern"), sep="_")) %>% 
  list_rbind()%>% filter(!is.na(`Expression pattern`), `Expression pattern`!=2) %>% mutate(`Expression pattern`=ifelse(`Expression pattern`==1, "Housekeeping", "Tissue-specific")) %>% 
  ggplot(aes(residuals, fill=`Expression pattern`))+geom_density(alpha=0.5)+facet_wrap(~rep)+
  scale_fill_manual(values=c("#1B8C8E", "#0D2C54"))+
  ggpubr::theme_pubclean()#theme_bw()+theme(text=element_text(size=20))
ggsave("Plots/Fig3/tissue_specificity_residuals.pdf", width = 9, height = 6.75, units="in")

fantom_max_expresion_obs=fantom_max_expresion  %>% filter(!is.na(gini_tissue)) %>% mutate(sample_specificity_gini=cut_number(gini_tissue, n=3) %>% as.numeric()) %>% 
  filter(sample_specificity_gini!=2) %>% 
  group_by(rep) %>% 
  mutate(rank_reporter=row_number(mean),
         rank_endo=row_number(tpm)) 

# fantom_max_expresion_obs%>% 
#   group_split(rep) %>% map(
#     ~tibble(Control=abs(sample(nrow(.x))-sample(nrow(.x))),
#             Observed=abs(.x$rank_endo-.x$rank_reporter),
#             rep=.x$rep, sample_specificity_gini=.x$sample_specificity_gini)) %>% list_rbind() %>% 
#   pivot_longer(-c(rep, sample_specificity_gini)) %>% 
#   mutate(sample_specificity_gini=ifelse(name=="Control", "All",
#                                         ifelse(sample_specificity_gini==1, "Housekeeping", 
#                                                ifelse(sample_specificity_gini==2, "Intermediate", "Tissue-specific")))) %>% 
#              ggplot(aes(value, linetype=name, col=factor(sample_specificity_gini), 
#              fill=factor(sample_specificity_gini),
#              alpha=factor(sample_specificity_gini)))+
#   geom_density()+
#   facet_wrap(~rep)+xlab("Reporter vs Endogenous rank difference")+
#   labs(col="", linetype="", fill="")+
#   scale_fill_manual(values=c("Housekeeping"="#1B8C8E", "Tissue-specific"="#0D2C54", "All"="darkred"))+
#   scale_color_manual(values=c("Housekeeping"="#1B8C8E", "Tissue-specific"="#0D2C54", "All"="darkred"))+
#   scale_alpha_manual(values=c("Housekeeping"=0.5, "Tissue-specific"=0.5, "All"=0))+
#   theme_bw()+
#   guides(col="none", alpha="none")+
#   #ggpubr::theme_pubr(legend = "right")+
#   scale_linetype_manual(values=c("Control"=2, "Observed"=1)) 
# ggsave("Plots/Fig3/reporter_endo_rank_dif.pdf", width = 12, height = 6.75, units="in")

fantom_max_expresion_ctl<-fantom_max_expresion_obs %>% 
  group_split(rep) %>% map(
    ~tibble(Control=abs(sample(1:nrow(.x), size=nrow(.x)*1000, replace = T)-sample(1:nrow(.x), size=nrow(.x)*1000, replace=T)),
            rep=unique(.x$rep))) %>% list_rbind()
  
fantom_max_expresion_obs%>% 
  mutate(Observed=abs(rank_endo-rank_reporter)) %>% 
  mutate(sample_specificity_gini=ifelse(sample_specificity_gini==1, "Housekeeping", 
                                               ifelse(sample_specificity_gini==2, "Intermediate", "Tissue-specific"))) %>% 
  ggplot(aes(Observed))+
  geom_density(alpha=0.5, aes(col=factor(sample_specificity_gini), 
                              fill=factor(sample_specificity_gini)))+
  geom_density(data=fantom_max_expresion_ctl, aes(Control), col="darkred", linetype="dashed")+
  facet_wrap(~rep)+xlab("Reporter vs Endogenous rank difference")+
  labs(col="", fill="")+
  scale_fill_manual(values=c("Housekeeping"="#1B8C8E", "Tissue-specific"="#0D2C54", "All"="darkred"))+
  scale_color_manual(values=c("Housekeeping"="#1B8C8E", "Tissue-specific"="#0D2C54", "All"="darkred"))+
  theme_bw()+
  guides(col="none")
ggsave("Plots/Fig3/reporter_endo_rank_dif.pdf", width = 12, height = 6.75, units="in")



fantom_mean_expresion<- rbind(tissue ,cell) %>%
  mutate(tpm=1E6*counts/libsize) %>% 
  arrange(desc(tpm)) %>% group_by(name) %>%
  summarise(tpm=mean(tpm)) %>% 
  inner_join(data %>% select(seq_id, gene_sym, mean,mean_sw, rep), by=c("name"="seq_id"))
fantom_mean_expresion %>% 
  group_by(rep, mean_sw) %>% 
  summarise(median_endo=median(tpm)) %>% 
  ggplot(aes(mean_sw %>% as.numeric(), median_endo))+
  geom_point(col="#AD343E")+
  geom_smooth(col="#216869", method="lm")+
  facet_wrap(~rep)+labs(x="Mean activity (binned by rank)", 
                        y="Median of mean endogenous activity (TPM)")+
  theme_bw()+theme(text=element_text(size=20))+scale_y_log10()
ggsave("Plots/Fig3/reporter_endo_mean_fantom.pdf", width = 9, height = 6.75, units="in")

fantom_mean_expresion %>% 
  group_by(rep) %>% 
  mutate(mean_sw=row_number(mean) %>% cut_number(20) %>% as.numeric() %>% as_factor()) %>%
  ggplot(aes(mean_sw, tpm))+
  geom_boxplot(fill="#AD343E", alpha=0.5, outliers = F)+
  facet_wrap(~rep)+labs(x="Mean reporter activity (binned)", 
                        y="Mean Endogenous activity across tissues (TPM)")+
  theme_bw()+
  theme(text=element_text(size=20))+scale_y_log10()
ggsave("Plots/Fig3/reporter_endo_mean_fantom_boxplot.pdf", width = 9, height = 6.75, units="in")

fantom_max_expresion %>% 
  group_by(rep) %>% 
  mutate(mean_sw=row_number(mean) %>% cut_number(20) %>% as.numeric() %>% as_factor()) %>%
  ggplot(aes(mean_sw, tpm))+
  geom_boxplot(fill="#AD343E", alpha=0.5, outliers = F)+
  facet_wrap(~rep)+labs(x="Mean reporter activity (binned)", 
                        y="TPM (max value along CAGE FANTOM5)")+
  theme_bw()+
  theme(text=element_text(size=20))+scale_y_log10()
ggsave("Plots/Fig3/reporter_endo_max_fantom_boxplot.pdf", width = 9, height = 6.75, units="in")

### ReMap dataset

library(clusterProfiler)
library(ReactomePA)
library(org.Hs.eg.db)
hs <- org.Hs.eg.db
# WARNING: HIGH TIME AND MEMORY DEMAND:
# nr=vroom::vroom("Analysis/Tables/library_remap_NR.bed", 
#                 col_names = c("seqnames", "start", "end", "seq_id", "score", "strand", "seqnames_peak", 
#                               "start_peak",  "end_peak", "sample", "N_tissue", "strand_peak", "thickstart_peak", "thickend_peak", "rgb"))
# nr=nr %>% separate(sample, into=c("TF", "sample"), sep=":")
# nr_N=nr %>% group_split(TF) %>% 
#   map(~tibble(TF=unique(.x$TF),N=nrow(.x))) %>% list_rbind()
# 
# TF_data=data
# for(TFs in unique(nr$feature)){
#   subdata=nr %>% filter(feature==TFs)
#   if(any(subdata$seq_id%in%TF_data$seq_id)){
#     TF_data[[paste0("TF_",TFs)]]=TF_data$seq_id%in%(subdata$seq_id)
#   }}
# rm(nr)
# gc()
# 
# write_tsv(TF_data, "Analysis/Tables/data_TF_remap.tsv")
TF_data=read_tsv("Analysis/Tables/data_TF_remap.tsv")
nTF=TF_data %>% select(-TF_knownmotif) %>% group_by(rep) %>% summarise(across(starts_with("TF_"),sum))%>% 
  pivot_longer(starts_with("TF_"), names_to = "feature",values_to = "counts") %>% 
  filter(counts>100) %>% add_count(feature) %>% filter(n==2)

tidy_TF_data=TF_data %>% dplyr::select(mean, mean_sw, var_rank_sw, rep, any_of(unique(nTF$feature))) %>% 
  mutate(rep=as.factor(rep),across(where(is.logical), ~.x%>% factor(level=c("TRUE","FALSE"))))
repname=unique(tidy_TF_data$rep)
# TF Binding site activity association
TF_data %>% group_by(rep, mean_sw) %>% summarise(NFYA=sum(TF_NFYA)/n()) %>% 
  ggplot(aes(mean_sw %>% as.numeric(),NFYA))+geom_point(col="#14AFB2")+theme_bw()+
  geom_smooth(col="#216869")+theme(text=element_text(size=20))+
  facet_wrap(~rep)+ylim(0,1)+labs(x="Mean activity (binned by rank)",y="Proportion of promoters")+ggtitle("NFYA binding")

ggsave("Plots/Fig2/NFYA_scatter.pdf", width = 9, height = 6.75, units="in")

TF_data %>% group_by(rep, mean_sw) %>% summarise(SP1=sum(TF_SP1)/n()) %>% 
  ggplot(aes(mean_sw %>% as.numeric(),SP1))+geom_point(col="#14AFB2")+
  geom_smooth(col="#216869")+theme_bw()+theme(text=element_text(size=20))+
  facet_wrap(~rep)+ylim(0,1)+labs(x="Mean activity (binned by rank)",y="Proportion of promoters")+ggtitle("SP1 binding")

ggsave("Plots/Fig2/SP1_scatter.pdf", width = 9, height = 6.75, units="in")

TF_data %>% group_by(rep, mean_sw) %>% summarise(SP2=sum(TF_SP2)/n()) %>% 
  ggplot(aes(mean_sw %>% as.numeric(),SP2))+geom_point(col="#14AFB2")+
  geom_smooth(col="#216869")+theme_bw()+theme(text=element_text(size=20))+
  facet_wrap(~rep)+ylim(0,1)+labs(x="Mean activity (binned by rank)",y="Proportion of promoters")+ggtitle("SP2 binding")

ggsave("Plots/Fig2/SP2_scatter.pdf", width = 9, height = 6.75, units="in")

# Whole Remap - activity analysis
wilcox=map(tidy_TF_data %>% select(-mean, -mean_sw,-var_rank_sw, -rep),
           ~coin::wilcox_test(formula=mean~.x|rep, data=tidy_TF_data) %>% pvalue())
wilcox=tibble(feature=names(wilcox), 
              pval=list_c(wilcox), 
              pval_corr=p.adjust(pval, "BH", length(wilcox)))
wilcox_rep1=map(tidy_TF_data %>% 
                  filter(rep==repname[1])%>%  
                  select(-mean, -mean_sw,-var_rank_sw, -rep), 
                ~coin::wilcox_test(formula=mean~.x, 
                                   data=tidy_TF_data%>% filter(rep==repname[1]), conf.int=T))
wilcox_rep2=map(tidy_TF_data %>% filter(rep==repname[2])%>%  select(-mean, -mean_sw,-var_rank_sw, -rep), ~coin::wilcox_test(formula=mean~.x, data=tidy_TF_data%>% filter(rep==repname[2]), conf.int=T))

wilcox_df2=map2(list(wilcox_rep1, wilcox_rep2),repname,
                ~tibble(feature=names(.x),estimate=map(.x,~confint(.x)$estimate) %>% list_c(), 
                        P2.5=map(.x,~confint(.x)$conf.int[1]) %>% list_c(),
                        P97.5=map(.x,~confint(.x)$conf.int[2]) %>% list_c() ,
                        rep=.y)%>% 
                  pivot_longer(c(starts_with("estimate"),starts_with("P2.5"), starts_with("P97.5")),names_to = "val", values_to="estimate") %>% 
                  arrange(desc(estimate)) %>% 
                  mutate(feature=fct_inorder(feature))) %>% 
  list_rbind() %>% left_join(wilcox)

wilcox_df2 %>% filter(pval_corr<0.05)%>% 
  mutate(act=ifelse(sign(estimate)==1,"High activity", "Low activity")) %>% group_by(feature) %>%
  mutate(n=length(unique(act))) %>% filter(n==1) %>%
  arrange(desc(rep),estimate) %>% 
  pivot_wider(names_from = val, values_from = estimate) %>% 
  filter(pval_corr<0.05) %>% 
  ggplot(aes(x=estimate, y=fct_inorder(feature), group=fct_inorder(rep)))+
  geom_col(orientation="y", position="dodge", aes( fill=act, alpha=rep))+
  scale_alpha_manual(values=c("Rep 1"=1, "Rep 2"=0.5))+
  ggpubr::theme_pubr()+
  theme(text=element_text(size=20), legend.position = "top", axis.text.y = element_blank(), axis.ticks.y = element_blank())+labs(fill="Effect", x="Activity Effect", y="Features", alpha="")+
  scale_fill_manual(values=c("High activity"="#D6741F", "Low activity"="#7FB800"))
ggsave("Plots/Fig2/remap_act.pdf", width = 9, height = 6.75, units="in")

wilcox_df2 %>% filter(pval_corr<0.05)%>% 
  mutate(act=ifelse(sign(estimate)==1,"High activity", "Low activity")) %>% group_by(feature) %>%
  mutate(n=length(unique(act))) %>% filter(n==1, act=="Low activity") %>%arrange(desc(rep),estimate) %>% 
  pivot_wider(names_from = val, values_from = estimate) %>% 
  filter(pval_corr<0.05) %>%  ggplot(aes(x=estimate, y=fct_inorder(feature), group=fct_inorder(rep)))+
  geom_col(orientation="y", position="dodge", aes( fill=act, alpha=rep))+
  scale_alpha_manual(values=c("Rep 1"=1, "Rep 2"=0.5))+
  ggpubr::theme_pubr()+
  theme(text=element_text(size=20), legend.position = "top")+labs(fill="Effect", x="Activity Effect", y="Features", alpha="")+
  geom_errorbarh(aes(xmax=P2.5 , xmin=P97.5), position=position_dodge(1), height=0.05, col="#777777",  linewidth=1.5)+
  scale_fill_manual(values=c("High activity"="#D6741F", "Low activity"="#7FB800"))
ggsave("Plots/Fig2/remap_act_low.pdf", width = 9, height = 6.75, units="in")

# GSEA Remap activity
pregsea_act=function(df, rep){
  ids<-AnnotationDbi::select(hs, 
                             keys = (df %>% filter(rep==rep) %>% mutate(feature=str_remove(feature,"TF_")))$feature,
                             columns = c("ENTREZID", "SYMBOL"),
                             keytype = "SYMBOL") %>% unique()
  
  pregsea=df%>% mutate(feature=str_remove(feature,"TF_")) %>% 
    filter(val=="estimate") %>% arrange(desc(estimate))%>% left_join(ids, by=c("feature"="SYMBOL"))%>% 
    filter(!is.na(ENTREZID)) %>% group_split(rep)
  
  ordered_genes=pregsea[[rep]]$estimate
  names(ordered_genes)=pregsea[[rep]]$ENTREZID
  return(ordered_genes)
}

ordered_genes=pregsea_act(wilcox_df2, repname[1])
GSE_GO_rep1<-gseGO(ordered_genes, ont="all", OrgDb = "org.Hs.eg.db")
GSE_GO_rep1_df<-as.data.frame(GSE_GO_rep1) %>% arrange(desc(NES))

ordered_genes=pregsea_act(wilcox_df2, repname[2])
GSE_GO_rep2<-gseGO(ordered_genes, ont="all", OrgDb = "org.Hs.eg.db")
GSE_GO_rep2_df<-as.data.frame(GSE_GO_rep2) %>% arrange(desc(NES))
writexl::write_xlsx(list(Rep1=GSE_GO_rep1_df,Rep2=GSE_GO_rep2_df), "Analysis/Tables/GSEA_activity_remap.xlsx")

pol2_entrez=GSE_GO_rep1@geneSets[[GSE_GO_rep1$ID[1]]] 
pol2_tf=pregsea[[1]]$feature[pregsea[[1]]$ENTREZID%in%pol2_entrez]
p1=enrichplot::gseaplot2(GSE_GO_rep1, GSE_GO_rep1_df$ID[1], title=GSE_GO_rep1_df$Description[1], subplots = 1, color = "darkred")

p2=wilcox_df2 %>% filter(pval_corr<0.05, rep=="Rep 1")%>% 
  mutate(act=ifelse(feature%in%paste0("TF_", pol2_tf), "set",
                    ifelse(sign(estimate)==1,"High activity", "Low activity"))) %>% 
  filter(val=="estimate") %>% 
  group_by(feature) %>%arrange(desc(estimate)) %>% 
  pivot_wider(names_from = val, values_from = estimate) %>% 
  ggplot(aes(x=estimate, y=fct_inorder(feature)))+
  geom_col(orientation="y", position="dodge", aes( fill=act, alpha=act))+
  ggpubr::theme_pubr()+
  theme(text=element_text(size=10), legend.position = "none", 
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank())+
  labs(fill="Effect", x="Activity Effect", y="Transcription factors")+
  scale_fill_manual(values=c("High activity"="#D6741F", "Low activity"="#7FB800", "set"="darkred"))+
  scale_alpha_manual(values=c("High activity"=0.4, "Low activity"=0.4, "set"=1))+
  coord_flip()
aplot::gglist(list(p1,p2), ncol = 1, heights = c(.75, 0.5))
ggsave("Plots/Fig2/GSEA_POL2.pdf", width = 9, height = 6.75, units="in" )

## Remap noise analysis

vbles=tidy_TF_data  %>% select( -mean_sw, -var_rank_sw, -rep, -mean) %>% names()

auc_l=list()
for(j in unique(tidy_TF_data$rep)){
  df=tidy_TF_data %>% filter(rep==j)
  for(i in vbles){
    auc_l[[j]][[i]]=pROC::roc(df[[i]],df[["var_rank_sw"]],ci=T, n.boots=1000, verbose=F, direction=">")
  }
}
auc_df=list()
for (i in c(1,2)){
  auc_df[[i]]=map(auc_l[[i]], ~tibble(rep=unique(tidy_TF_data$rep)[i], 
                                      ci2.5=(.x$ci %>% as.numeric)[1], 
                                      ci97.5=(.x$ci %>% as.numeric)[3],
                                      AUC=(.x$ci %>% as.numeric)[2])) %>% 
    list_rbind() %>% mutate(feature=names(auc_l[[1]]))
}
auc_df=list_rbind(auc_df)

plot_auc2(auc_df)
ggsave("Plots/Fig4/Remap_noise.pdf", width = 9, height = 6.75, units="in")

ids<-AnnotationDbi::select(hs, 
                           keys = (auc_df %>% filter(rep==rep) %>% mutate(feature=str_remove(feature,"TF_")))$feature,
                           columns = c("ENTREZID", "SYMBOL"),
                           keytype = "SYMBOL") %>% unique()

pregsea_noise= auc_df %>%  
  mutate(feature=str_remove(feature,"TF_")) %>% arrange(desc(AUC))%>% 
  left_join(ids, by=c("feature"="SYMBOL"))%>% filter(!is.na(ENTREZID)) %>% 
  group_split(rep)
ordered_noise=pregsea_noise[[1]]$AUC
names(ordered_noise)=pregsea_noise[[1]]$ENTREZID
GSE_GO_noise_rep1<-gseGO(ordered_noise, ont="CC", OrgDb = "org.Hs.eg.db")
GSE_GO_noise_rep1_df<-as.data.frame(GSE_GO_noise_rep1) %>% arrange(desc(NES))

ordered_noise=pregsea_noise[[2]]$AUC
names(ordered_noise)=pregsea_noise[[2]]$ENTREZID
GSE_GO_noise_rep2<-gseGO(ordered_noise, ont="CC", OrgDb = "org.Hs.eg.db")
GSE_GO_noise_rep2_df<-as.data.frame(GSE_GO_noise_rep2) %>% arrange(desc(NES))

writexl::write_xlsx(list(Remap_activity_rep1=GSE_GO_rep1_df, Remap_activity_rep2=GSE_GO_rep2_df, Remap_noise_rep1=GSE_GO_noise_rep1_df, Remap_noise_rep2=GSE_GO_noise_rep2_df),"Analysis/Tables/GSEA_remap.xlsx")

MLL1_entrez=GSE_GO_noise_rep1@geneSets[["GO:0071339"]] 
MLL1_tf=pregsea_noise[[1]]$feature[pregsea_noise[[1]]$ENTREZID%in%MLL1_entrez]

p1=enrichplot::gseaplot2(GSE_GO_noise_rep1, "GO:0071339", title="MLL1 complex", subplots = 1, color = "darkred")
p2=auc_df %>% filter(rep=="Rep 1", ((AUC<0.5 & ci97.5<0.5)|(AUC>0.5 &ci2.5>0.5))) %>% 
  mutate(act=ifelse(feature%in%paste0("TF_", MLL1_tf), "set",
                                                                ifelse(AUC>0.5,"High noise", "Low noise"))) %>% 
  arrange(desc(AUC)) %>% 
  ggplot(aes(x=AUC, y=fct_inorder(feature)))+
  geom_col(orientation="y", position="dodge", aes( fill=act, alpha=act))+
  ggpubr::theme_pubr()+
  theme(text=element_text(size=10), legend.position = "none", 
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank())+
  labs(fill="Effect", x="Noise Effect", y="Transcription factors")+
  scale_fill_manual(values=c("High noise"="#D6741F", "Low noise"="#7FB800", "set"="darkred"))+
  scale_alpha_manual(values=c("High noise"=0.4, "Low noise"=0.4, "set"=1))+
  coord_flip()
aplot::gglist(list(p1,p2), ncol = 1, heights = c(.75, 0.5))
#ggsave("Plots/Fig4/GSEA_MLL1.pdf", width = 9, height = 6.75, units="in" )

# ChIP-Atlas 
chipatlas=read_tsv("External_data/allPeaks_light.hg38.50_lib_unique.tsv", col_names = c("seq_id", "feature"))
explist=read_tsv("External_data/experimentList_chipatlas.tab", col_names = F) %>% select(X3, X4) %>% unique()
explist=explist[duplicated(explist$X4),]
chipatlas=chipatlas %>% count(seq_id, feature, name = "n_samples") %>% 
  left_join(explist, by=c("feature"="X4")) %>% rename(group=X3) 
chipatlas %>% write_tsv("Analysis/Tables/allPeaks_light.hg38.50_lib_final.tsv")

histones=chipatlas %>% filter(group=="Histone")
hist_data=data
for(TFs in unique(histones$feature)){
  subdata=histones %>% filter(feature==TFs)
  if(any(subdata$seq_id%in%data$seq_id)){
    hist_data[[paste0("Hist_",TFs)]]=data$seq_id%in%(subdata$seq_id)
  }}
write_tsv(hist_data, "Analysis/Tables/data_Hist_chipatlas.tsv")

nTF=hist_data  %>% group_by(rep) %>% 
  summarise(across(starts_with("Hist_"),sum)) %>% 
  pivot_longer(starts_with("Hist_"), names_to = "feature",values_to = "counts") %>% 
  filter(counts>100) %>% add_count(feature) %>% filter(n==2)

tidy_hist_data=hist_data %>% 
  select(mean, mean_sw, var_rank_sw, rep, any_of(unique(nTF$feature))) %>%
  mutate(rep=as.factor(rep),across(where(is.logical), 
                                   ~.x%>% factor(level=c("TRUE","FALSE"))))
repname=unique(tidy_hist_data$rep)

vbles=tidy_hist_data  %>% select( -mean_sw, -var_rank_sw, -rep, -mean) %>% names()
auc_l=list()
for(j in unique(tidy_hist_data$rep)){
  df=tidy_hist_data %>% filter(rep==j)
  for(i in vbles){
    auc_l[[j]][[i]]=roc(df[[i]],df[["var_rank_sw"]], ci=T, boot.n=2000, direction=">")
  }
}
auc_df=list()
for (i in c(1,2)){
  auc_df[[i]]=map(auc_l[[i]], ~tibble(rep=unique(tidy_hist_data$rep)[i], ci2.5=(.x$ci %>% as.numeric)[1], ci97.5=(.x$ci %>% as.numeric)[3], AUC=(.x$ci %>% as.numeric)[2])) %>% 
    list_rbind() %>% mutate(feature=names(auc_l[[1]]))
}
auc_df=list_rbind(auc_df)

plot_auc2(auc_df)
ggsave("Plots/Fig2/chipatlas_noise.pdf", width = 9, height = 6.75, units="in")

### Alternative promoters (poner en otro script como armo df sample_activity)

epd<- rtracklayer::import.bed("External_data/EPD/human38_epdnew.bed") %>% 
  IRanges::promoters(downstream=16, upstream = 236) #full EPD database
lib<- rtracklayer::import.bed("Library_data/res/library.bed")
lib$name<- str_remove(lib$name, "^FP.{6}_")
full_lib<- c(epd, BiocGenerics::subset(lib, !lib$name%in%epd$name)) %>% unique() #add the rest of EPD to library
full_lib$score=NULL

sample_activity<- read_tsv("Analysis/Tables/sample_CAGE_activity.tsv") # if it was not generated
sample_activity<- mutate(sample_activity, tpm=counts*1E6/libsize)
sample_activity_prom<- sample_activity %>% filter(name%in%epd$name) #solo los promotores
sample_activity_prom$gene_sym<- str_remove(sample_activity_prom$name, "_.{1,2}$") 

# primmary classification
non_detected<- unique(sample_activity_prom$name)[!unique(sample_activity_prom$name)%in%unique(sample_activity_prom$name[sample_activity_prom$tpm>1])]

tot_counts=sample_activity_prom %>% 
  group_by(name, gene_sym) %>% 
  summarise(tpm=sum(tpm)) %>% 
  group_by(gene_sym) %>% 
  mutate(max_tpm=tpm==max(tpm), 
         tpm_gene=sum(tpm), 
         perc_tpm_gene=tpm/tpm_gene) %>%
  ungroup()

unique_prom<- tot_counts %>% 
  filter(perc_tpm_gene>0.999, !name%in%non_detected) %>% 
  select(name) %>% unique() %>% as_vector()

epd_unique<- sample_activity_prom %>%
  select(name, gene_sym) %>% unique() %>% 
  count(gene_sym) %>% filter(n==1) %>% 
  select(gene_sym) %>% as_vector()

sample_activity_alt<- sample_activity_prom %>% 
  filter(!name%in%c(non_detected, unique_prom)) %>% 
  group_by(gene_sym, sample) %>% 
  mutate(tpm_gene_sample=sum(tpm), 
         perc_tpm_gene_sample=tpm/tpm_gene_sample, 
         is_highest_in_sample=tpm==max(tpm)) %>% 
  ungroup()
sample_activity_alt<- sample_activity_alt %>% 
  filter(tpm_gene_sample>1) %>% 
  group_by(name) %>%
  summarise(highest_in_N_samples= sum(is_highest_in_sample), total_tpm_prom=sum(tpm)) %>% 
  inner_join(sample_activity_alt,.) %>% ungroup() %>% 
  group_by(gene_sym,sample) %>% 
  mutate(main_N_samples= highest_in_N_samples/sum(highest_in_N_samples)) %>%
  ungroup()
main_prom<- sample_activity_alt %>% 
  left_join(tot_counts %>% select(-tpm), by=c("name", "gene_sym")) %>% 
  group_by(gene_sym) %>% 
  filter((perc_tpm_gene==max(perc_tpm_gene))&(main_N_samples==max(main_N_samples)))

unclassified<-sample_activity_alt %>% 
  left_join(tot_counts %>% select(-tpm), by=c("name", "gene_sym")) %>% 
  group_by(gene_sym) %>% 
  filter(xor((perc_tpm_gene==max(perc_tpm_gene)),(main_N_samples==max(main_N_samples)))) 

#secondary classification

wide_df<- sample_activity_alt %>% 
  select(name, gene_sym, tpm, sample) %>% 
  inner_join(main_prom %>% 
               select(name, gene_sym, tpm, sample) %>% 
               unique() %>% 
               rename(main_name=name, main_tpm=tpm), by=c("sample","gene_sym")) %>% 
  filter(name!=main_name)

cor_data_pearson<- wide_df %>% ungroup() %>% 
  filter(tpm+main_tpm>0) %>% group_split(name) %>% 
  map(~summarise(.x, cor=cor(tpm, main_tpm, method = "pearson"), 
                 n_samples=n(), name=unique(name), 
                 main_name=unique(main_name))) %>% 
  list_rbind()

cor_data_spearman<- wide_df %>% ungroup() %>% 
  filter(tpm+main_tpm>0) %>% group_split(name) %>% 
  map(~summarise(.x, cor=cor(tpm, main_tpm, method = "spearman"), 
                 n_samples=n(), name=unique(name), 
                 main_name=unique(main_name))) %>% 
  list_rbind()

cor_data<- left_join(cor_data_pearson, cor_data_spearman, 
                     by=c("n_samples", "name", "main_name"), 
                     suffix = c("_pearson", "_spearman"))

wide_df<- wide_df %>% left_join(cor_data)
wide_df<- wide_df %>% mutate(ratio_alt_main=(log10(tpm)/log10(main_tpm)), 
                             switch=(ratio_alt_main>1.5)&(tpm>5), 
                             switch=ifelse(tpm>5&main_tpm<1, T,switch)) %>% 
  group_by(name) %>% mutate(N_switch=sum(switch)) %>% 
  arrange(desc(N_switch))

short_df<- wide_df %>% 
  select(name, main_name, cor_pearson,
         n_samples, cor_spearman, prom_alt, N_switch) %>%
  unique()

short_df %>% mutate(with_switch_cases=ifelse(N_switch>0, "With switch cases", "Without switch cases"),
                    prom_alt=factor(prom_alt, levels = c("independent", "correlated", "switch"))) %>% 
  ggplot(aes(cor_pearson, fill=prom_alt))+
  facet_wrap(~with_switch_cases, nrow=2)+
  geom_density(alpha=0.75)+
  ggpubr::theme_pubclean()+
  labs(x="Pearson correlation between the activity \nof the Main-Secondary pair", 
       fill="Secondary \npromoter class")+
  scale_fill_manual(values=c("#48d6d9","#14AFB2", "#216869"))
ggsave("Plots/Fig Promalt/secondary_classification.pdf", width = 9, height = 6.75, units="in")

p=ggviolin(data %>% 
             filter(!is.na(prom_alt), !prom_alt%in%c("unclassified", "non_detected")) %>% 
             mutate(prom_alt2=ifelse(prom_alt%in%c("independent", "correlated", "switch"),
                                     "Secondary promoter", prom_alt)) %>% 
             mutate(prom_alt2=ifelse(prom_alt2=="unique", 
                                     "No alternative promoters", prom_alt2) %>% 
                      factor(levels = c("No alternative promoters", "Main promoter", "Secondary promoter"))), 
           x="prom_alt2", y="mean", fill="prom_alt2", draw_quantiles = 0.5,
           add="median_q1q3",  alpha=0.85, 
           palette = c(`Main promoter`= "#35B166", `Secondary promoter`="#5D9D1B", `No alternative promoters`="#B2D4FB"))
p=facet(p, facet.by = "rep", ncol=2)
p+stat_compare_means(comparisons = list(c("Main promoter", "Secondary promoter"),
                                        c("Main promoter", "No alternative promoters"), 
                                        c("No alternative promoters", "Secondary promoter")), 
                     label = "p.signif", method="wilcox.test", na.rm = T)+
  theme_pubclean()+
  labs(fill="Alternative promoters", x=NULL, 
       y="Transcriptional mean activity")+
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(), 
        text=element_text(size=15))+ylim(0,7.5)
ggsave("Plots/Fig Promalt/violin_promalt_main_sec_unique.pdf",width = 12, height = 6.75, units="in")


df_split_main<- wide_df %>% select(main_name, prom_alt) %>% unique() %>% 
  rename(name=main_name) %>%
  mutate(prom_alt=paste(prom_alt, "main")) %>% 
  rbind(wide_df %>% select(name, prom_alt) %>% unique()) %>% 
  inner_join(data %>% select(-prom_alt))

tidy_palt<- df_split_main %>% 
  select(prom_alt, rep, mean) %>% 
  fastDummies::dummy_cols("prom_alt", ignore_na = T, omit_colname_prefix=T, remove_selected_columns = T) %>%
  mutate(across(c(where(~ is.numeric(.x)), -contains("mean")), as.logical),
         #`Secondary promoter`=independent|correlated|switch, 
         across(where(is.logical), ~.x%>% factor(level=c("TRUE","FALSE"))))

wilcox=map(tidy_palt %>% select(-rep,-contains("mean"),-contains("var")), ~coin::wilcox_test(formula=mean~.x|rep, data=tidy_palt %>% mutate(rep=as.factor(rep))) %>% pvalue())
wilcox=tibble(feature=names(wilcox), pval=list_c(wilcox), pval_corr=p.adjust(pval, "BH", length(wilcox)))


wilcox_rep1=map(tidy_palt %>% filter(rep=="Rep 1") %>%  
                  select(-mean,-rep), 
                ~coin::wilcox_test(formula=mean~.x, data=tidy_palt%>% filter(rep=="Rep 1"), conf.int=T))
wilcox_rep2=map(tidy_palt %>% filter(rep=="Rep 2")%>%  
                  select(-mean,-rep), 
                ~coin::wilcox_test(formula=mean~.x, data=tidy_palt%>% filter(rep=="Rep 2"), conf.int=T))
wilcox_df2=map2(list(wilcox_rep1, wilcox_rep2),c("Rep 1","Rep 2"),
                ~tibble(feature=names(.x),estimate=map(.x,~confint(.x)$estimate) %>% list_c(), 
                        P2.5=map(.x,~confint(.x)$conf.int[1]) %>% list_c(),
                        P97.5=map(.x,~confint(.x)$conf.int[2]) %>% list_c() ,
                        rep=.y)%>% 
                  pivot_longer(c(starts_with("estimate"),starts_with("P2.5"), starts_with("P97.5")),
                               names_to = "val", values_to="estimate") %>% 
                  arrange(desc(estimate)) %>% 
                  mutate(feature=fct_reorder())) %>% list_rbind()

wilcox_df2 %>% 
  pivot_wider(names_from = val, values_from = estimate) %>% 
  mutate(act=ifelse(sign(estimate)==1,"High activity", "Low activity"), 
         cat=str_remove(feature, " main"), cat2=ifelse(str_detect(feature, "main"), "main", "sec"
         )) %>% 
  arrange(desc(cat), desc(rep), desc(cat2)) %>% 
  ggplot(aes(x=estimate, y=fct_inorder(feature), group=fct_inorder(rep))) +
  geom_col(orientation="y", position="dodge", aes( fill=act, alpha=rep)) +
  scale_alpha_manual(values=c("Rep 1"=1, "Rep 2"=0.7)) +
  theme(text=element_text(size=20), legend.position = "top") + 
  labs(fill="Effect", x="Activity Effect", y="Features", alpha="") +
  geom_errorbarh(aes(xmax=P2.5 , xmin= P97.5), 
                 position=position_dodge(1), height=0.05, col="#777777",  
                 linewidth=1.5)+
  xlim(c(-0.3, 0.3)) +
  scale_fill_manual(values=c("High activity"="#D6741F", "Low activity"="#7FB800")) +
  guides(alpha="none") + theme_pubr()
ggsave("Plots/Fig Promalt/summary_act_promalt_with_mains_v2.pdf", width = 13.5, height = 6.75, units="in")


data %>% 
  mutate(prom_alt2=ifelse(prom_alt%in%c("independent", "correlated", "switch"), "Secondary promoter", prom_alt)) %>% 
  mutate(prom_alt2=ifelse(prom_alt2=="unique", "No alternative promoters", prom_alt2)) %>% 
  filter(!is.na(prom_alt2), !prom_alt2%in%c("non_detected", "unclassified")) %>%
  arrange(rep, mean) %>% 
  group_by(rep) %>% 
  mutate(mean_rank=row_number(mean),
         mean_sw= mean_rank %>% cut_width(100)) %>% ungroup() %>% 
  group_by(rep,mean_sw) %>% 
  mutate(var_rank_sw=row_number(var) %>% as.numeric()) %>% ungroup() %>% 
  mutate(mean_sw=as.numeric(mean_sw)) %>% 
  ggplot(aes(mean_sw,var_rank_sw, col=prom_alt2)) +
  geom_point(size=0.2, alpha=0.8)+facet_wrap(~rep) +
  labs(x="Mean rank group", y="Variance rank (SW)", col="") +
  scale_color_manual(values = c(`Main promoter`= "#16703A", `Secondary promoter`="#5D9D1B", `No alternative promoters`="#B2D4FB"),
                     guide = guide_legend(override.aes = list(size = 3))) +
  geom_smooth(se=F, aes(color=prom_alt2), show.legend = F) +
  theme_pubr(base_size = 20)
ggsave("Plots/Fig Promalt/scatter_noise_promalt.pdf", width = 9, height = 6.75, units="in")

df_split_main<-wide_df %>% select(main_name, prom_alt) %>% unique() %>% rename(name=main_name) %>% mutate(prom_alt=paste(prom_alt, "main")) %>% rbind(wide_df %>% select(name, prom_alt) %>% unique()) %>% inner_join(data %>% select(-prom_alt))

tidy_palt<- df_split_main %>% 
  select(prom_alt, rep, var_rank_sw) %>% 
  fastDummies::dummy_cols("prom_alt", ignore_na = T, omit_colname_prefix=T, remove_selected_columns = T) %>%
  mutate(across(c(where(~ is.numeric(.x)), -contains("var_rank_sw")), as.logical),
         across(where(is.logical), ~.x%>% factor(level=c("TRUE","FALSE"))))

library(pROC)
auc_l=list()
for(j in unique(tidy_palt$rep)){
  df=tidy_palt %>% filter(rep==j)
  for(i in names(tidy_palt %>% select(-rep,-var_rank_sw))){
    auc_l[[j]][[i]]=roc(df[[i]],df[["var_rank_sw"]], ci=T, boot.n=2000, direction=">")
  }
}
auc_df=list()
for (i in c(1,2)){
  auc_df[[i]]=map(auc_l[[i]], ~tibble(rep=unique(tidy_palt$rep)[i], ci2.5=(.x$ci %>% as.numeric)[1], ci97.5=(.x$ci %>% as.numeric)[3], AUC=(.x$ci %>% as.numeric)[2])) %>% 
    list_rbind() %>% mutate(feature=names(auc_l[[1]]))
}
auc_df=list_rbind(auc_df)



plot_auc3=function(df){
  df %>% mutate(noise=ifelse(AUC>0.5, "High noise", "Low noise"), cat=str_remove(feature, " main"),
  ) %>% 
    mutate(a=ifelse(str_detect(feature, "independent"), 6, ifelse(str_detect(feature, "correlated"), 4, 2)), 
           b=ifelse(str_detect(feature, "main"), 1, 0), c=a+b) %>% 
    arrange(c, desc(rep)) %>% 
    ggplot(aes(x=AUC-0.5, y=fct_inorder(feature), group=fct_inorder(rep))) +
    geom_col(orientation="y", position="dodge", aes( fill=noise, alpha=rep))+
    scale_alpha_manual(values=c("Rep 1"=1, "Rep 2"=0.7))+
    theme(text=element_text(size=20), legend.position = "top")+
    labs(fill="Effect", x="AUC (Noise effect)", y="Features", alpha="")+
    scale_x_continuous(labels = function(x) x + 0.5)+
    geom_errorbarh(aes(xmax=abs(ci2.5)-0.5 , xmin=abs(ci97.5)-0.5), 
                   position=position_dodge(1), height=0.05, col="#777777",  linewidth=1.5)+
    scale_fill_manual(values=c("High noise"="#D6741F", "Low noise"="#7FB800"))+ 
    theme_pubr()+guides(alpha=F)
}
plot_auc3(auc_df)
ggsave("Plots/Fig Promalt/summary_noise_promalt_with_mains.pdf", width = 13.5, height = 6.75, units="in")

# same-gene promoters:
pair_proms=wide_df %>% select(name, main_name, starts_with("cor"),, N_switch, prom_alt) %>% 
  unique() %>% 
  inner_join(data %>% select(name, rep, mean_rank_sw, mean, var, var_rank_sw)) %>% 
  inner_join(data %>% select(name, rep,  mean_rank_sw, mean, var, var_rank_sw) %>%
               rename(main_mean_rank=mean_rank_sw,
                      main_mean=mean,
                      main_var=var, 
                      main_var_rank=var_rank_sw), 
             by=c("main_name"="name", "rep")) %>% 
  mutate(true_dif_mean_rank=abs(mean_rank_sw-main_mean_rank),
         true_dif_mean=abs(mean-main_mean),
         true_dif_var=abs(var-main_var),
         true_dif_var_rank=abs(var_rank_sw-main_var_rank)) %>% 
  group_split(rep)

permut_pairs<-function(list,main_var, alt_var, n=10){
  permut<- vector(mode="list", length = 2)
  for(rep in c(1,2)){
    permut[[rep]]=vector(mode="list", length = n)
    for(i in seq(n)){
      permut[[rep]][[i]]=abs(pair_proms[[rep]][[main_var]]-sample(pair_proms[[rep]][[alt_var]]))
    }
    permut[[rep]]<- do.call(cbind, permut[[rep]]) %>% as_tibble() %>% pivot_longer(everything())%>%  mutate(dif="Permuted")
  }
  names(permut)<-map(list, ~unique(.x$rep))
  return(permut)
}


permut_main_mean=permut_pairs(pair_proms, "main_mean", "mean")
pairwise_mean_dif_df<-map2(pair_proms, permut_main_mean, ~ .x%>% ungroup() %>% 
                             select(true_dif_mean) %>% set_names("value") %>% 
                             mutate(dif="Observed")  %>% 
                             rbind(.y %>% select(-name)) %>% 
                             mutate(rep=unique(.x$rep))) %>% 
  list_rbind() 
pairwise_mean_dif_df%>% 
  ggplot(aes(value, fill=dif, linetype=rep))+geom_density(alpha=0.5)+
  labs(x="Pairwise Mean Diference", fill="")+
  scale_fill_manual(values=c("#3D518C", "#14AFB2"),
                    labels=c("true_dif_mean"="Observed", "random_dif"="Permuted"))+
  ggpubr::theme_pubclean()+
  ggtitle("Alternative promoters activity similarity")
ggsave("Plots/Fig Promalt/pairwise_dif_promalt.pdf", width = 10, height = 10, units="in")

wilcox_test(value~factor(dif),data=pairwise_mean_dif_df %>% filter(rep=="Rep 1"))
wilcox_test(value~factor(dif),data=pairwise_mean_dif_df %>% filter(rep=="Rep 2"))
wilcox_test(value~factor(dif)|factor(rep),data=pairwise_mean_dif_df)
effsize::cliff.delta(value~factor(dif),data=pairwise_mean_dif_df %>% filter(rep=="Rep 1"))
effsize::cliff.delta(value~factor(dif),data=pairwise_mean_dif_df %>% filter(rep=="Rep 2"))

permut_main_var_rank=permut_pairs(pair_proms, "main_var_rank", "var_rank_sw")

pairwise_varrank_dif_df<-map2(pair_proms, permut_main_var_rank, ~ .x%>% ungroup() %>% 
                                select(true_dif_var_rank) %>% set_names("value") %>% 
                                mutate(dif="Observed")  %>% 
                                rbind(.y %>% select(-name)) %>% 
                                mutate(rep=unique(.x$rep))) %>% 
  list_rbind() 
pairwise_varrank_dif_df %>% 
  ggplot(aes(value, fill=dif, linetype=rep))+geom_density(alpha=0.5)+
  labs(x="Pairwise Var Rank Diference", fill="")+
  scale_fill_manual(values=c("#3D518C", "#14AFB2"),
                    labels=c("true_dif_var_rank"="Observed", "random_dif"="Permuted"))+
  ggpubr::theme_pubclean()+
  ggtitle("Alternative promoters Noise similarity")
ggsave("Plots/Fig Promalt/pairwise_dif_promalt_noise.pdf",, width = 10, height = 10, units="in")

wilcox_test(value~factor(dif),data=pairwise_varrank_dif_df %>% filter(rep=="Rep 1"))
wilcox_test(value~factor(dif),data=pairwise_varrank_dif_df %>% filter(rep=="Rep 2"))
effsize::cliff.delta(value~factor(dif),data=pairwise_varrank_dif_df %>% filter(rep=="Rep 1"))
effsize::cliff.delta(value~factor(dif),data=pairwise_varrank_dif_df %>% filter(rep=="Rep 2"))

# upstream vs downstream

ordered_pair_proms<-map(pair_proms, ~.x %>% 
                          left_join(as_tibble(lib) %>% select(name, strand, start)) %>% 
                          left_join(as_tibble(lib) %>% select(name, start), 
                                    by=c("main_name"="name"), 
                                    suffix = c("", "_main")) %>% 
                          mutate(order_from_main=ifelse(strand=="+",
                                                        ifelse(start<start_main, "upstream", "downstream"),
                                                        ifelse(start>start_main, "upstream", "downstream")),
                                 rep=unique(.x$rep)))


p1<-ggviolin(ordered_pair_proms %>% list_rbind(), x="order_from_main", y="mean", fill="order_from_main", draw_quantiles = 0.5, add="median_q1q3",  alpha=0.7, palette = c("#3D518C", "#14AFB2"))
p1<-facet(p1,facet.by="rep", ncol=2)
p1+stat_compare_means(comparisons = list(c("upstream", "downstream")), 
                      label = "p.signif", method="wilcox.test", na.rm = T)+
  theme_pubclean()+
  labs( y="Transcriptional mean activity", x="Relative position to main promoter", fill=NULL)+
  theme( 
    text=element_text(size=15))+ylim(0,7.5)
ggsave("Plots/Fig Promalt/promalt_position_mean.pdf", , width = 12, height = 6.75, units="in")

coin::wilcox_test(formula=mean~factor(order_from_main), data=ordered_pair_proms[[1]])
coin::wilcox_test(formula=mean~factor(order_from_main), data=ordered_pair_proms[[2]])
effsize::cliff.delta(formula=mean~factor(order_from_main), data=ordered_pair_proms[[1]])
effsize::cliff.delta(formula=mean~factor(order_from_main), data=ordered_pair_proms[[2]])

p1<-ggviolin(ordered_pair_proms %>% list_rbind(), x="order_from_main", y="var_rank_sw", fill="order_from_main", draw_quantiles = 0.5, add="median_q1q3",  alpha=0.7, palette = c("#3D518C", "#14AFB2"))
p1<-facet(p1,facet.by="rep", ncol=2)
p1+stat_compare_means(comparisons = list(c("upstream", "downstream")), 
                      label = "p.signif", method="wilcox.test", na.rm = T)+
  theme_pubclean()+
  labs( y="Transcriptional noise", x="Relative position to main promoter", fill=NULL)+
  theme( 
    text=element_text(size=15))
ggsave("Plots/Fig Promalt/promalt_position_noise.pdf", , width = 12, height = 6.75, units="in")
coin::wilcox_test(formula=var_rank_sw~factor(order_from_main), data=ordered_pair_proms[[1]])
coin::wilcox_test(formula=var_rank_sw~factor(order_from_main), data=ordered_pair_proms[[2]])
effsize::cliff.delta(formula=var_rank_sw~factor(order_from_main), data=ordered_pair_proms[[1]])
effsize::cliff.delta(formula=var_rank_sw~factor(order_from_main), data=ordered_pair_proms[[2]])

coin::wilcox_test(formula=var_rank_sw~factor(order_from_main)|factor(rep), data=list_rbind(ordered_pair_proms))

#falta corer y guardar:
rank_dif_order<-as_tibble(lib) %>% add_count(gene_sym) %>% filter(n==2) %>% mutate(start=ifelse(strand=="-", -start,start)) %>% 
  arrange(gene_sym, start) %>% 
  group_by(gene_sym) %>% 
  mutate(order=row_number(start) %>% as_factor()) %>% 
  ungroup() %>% 
  inner_join(tot_counts %>% select(name, tpm)) %>% 
  inner_join(data %>% select(name, gene_sym, mean, var, rep)) %>% 
  group_by(rep) %>% 
  mutate(rank_reporter=row_number(mean), 
         rank_endo=row_number(tpm))

map(group_split(rank_dif_order,rep), 
    ~tibble(Control=abs(sample(nrow(.x))-sample(nrow(.x))), 
            Observed=abs(.x$rank_endo-.x$rank_reporter), 
            rep=.x$rep, order=.x$order)) %>% 
  list_rbind() %>% pivot_longer(-c(rep, order)) %>%
  mutate(order=ifelse(name=="Control", "All", 
                      ifelse(order==1, "Upstream", "Downstream"))) %>% 
  ggplot(aes(value, linetype=name, col=factor(order), 
           fill=factor(order),
           alpha=factor(order)))+
  geom_density()+
  facet_wrap(~rep)+xlab("Reporter vs Endogenous rank difference")+
  labs(col="", linetype="", fill="")+
  scale_fill_manual(values=c("Upstream"="#1B8C8E", "Downstream"="#0D2C54", "All"="darkred"))+
  scale_color_manual(values=c("Upstream"="#1B8C8E", "Downstream"="#0D2C54", "All"="darkred"))+
  scale_alpha_manual(values=c("Upstream"=0.5, "Downstream"=0.5, "All"=0))+
  theme_bw()+
  guides(col="none", alpha="none")+
  #ggpubr::theme_pubr(legend = "right")+
  scale_linetype_manual(values=c("Control"=2, "Observed"=1)) 

ggsave("Plots/Fig Promalt/promoter_order_rank_dif.pdf", width = 12, height = 6.75, units="in")
