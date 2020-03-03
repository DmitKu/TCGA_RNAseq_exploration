#loading libraries
devtools::install_github("mariodeng/FirebrowseR")
library(FirebrowseR)
library(reshape2)
library(Hmisc)
library(plyr)
#for scter
library(ggpubr)
library(ggpmisc)
library(ggplot2)
library(ds4psy)
library(ggrepel)

setwd("XXXXXXXXX")#set working directory
#==================================================================================== 

###load TCGA mRNAseq data for PDCD1 and CD274 mRNA
mRNA_Exp = Samples.mRNASeq(format = "csv",
                           gene = c("CD279","PDCD1", "CD274"),
                           protocol = "RSEM",
                           page = "1", page_size = "100000", sort_by = "cohort")


###making log2 values numeric============================== 
mRNA_Exp$expression_log2<-as.numeric(mRNA_Exp$expression_log2)


#### investigating data structure============================== 
summary(mRNA_Exp) ### checking data types in df
unique(mRNA_Exp$cohort) ###list of cohorts
unique(mRNA_Exp$sample_type) ### list of sample type
unique(table(mRNA_Exp$tcga_participant_barcode)) ### number of entries for a tcga_participant_barcode
table(mRNA_Exp$cohort)


####Preparing wide data format============================== 
mRNA_Exp_wide <- dcast(mRNA_Exp, tcga_participant_barcode  + cohort + 
                         protocol + sample_type~ gene, value.var="expression_log2")


####ploting scater plots, density of points is visualized
cols <-  colorRampPalette(c("#000099", "#00FEFF", "#45FE4F",
                            "#FCFF00", "#FF9400", "#FF3100"))(256)

ggplot(subset(mRNA_Exp_wide), aes(x=CD274, y=PDCD1)) +
  facet_wrap(~cohort, scales="fixed")+
  stat_density_2d(aes(fill = stat(level)), geom = 'polygon')+
  #geom_point(shape = '.')+
  scale_fill_gradientn(colours=cols)+theme_ds4psy()+ coord_fixed()+
  stat_smooth(method = "lm", color="black") +
  stat_cor(label.y = 10.5,size = 4)

ggsave("20200210_Pearson_PD-1_PD-L1_corr_in_each_cohort.pdf",
       units = "mm", width = 250, height = 200)

#Calculating pearson correlation for all cohorts
R_count_function<-function(df=mRNA_Exp_wide,r_method="pearson"){
  final_data = data.frame()
  for (cohort in unique(df$cohort)){
    df_cohort<-df[df$cohort==cohort,]
    r_cohort<-cor.test(df_cohort[,"PDCD1"],df_cohort[,"CD274"],method =r_method)
    data<-as.data.frame(cohort)
    data$corr<-r_cohort$estimate
    data$p.value<-r_cohort[3]
    data$method<-r_method
    final_data<-rbind(final_data,data)
  }
  names(final_data) <- c("cohort","Corr","P.value","Corr_Method")
  return (final_data)
}


####calculating  R
R_cor_df=R_count_function(df=mRNA_Exp_wide,r_method="pearson")
names(mRNA_Exp_wide)
unique(mRNA_Exp_wide$cohort)


####calculating median
length(table(mRNA_Exp_wide$tcga_participant_barcode))# number of all cases
table_median<-ddply(subset(mRNA_Exp_wide,!sample_type=="NT"),.variables=.(cohort),
                  summarise,Median_PDCD1=median(PDCD1,na.rm =T),
                  Median_CD274=median(CD274,na.rm =T),
                  SEM_PDCD1=sd(PDCD1,na.rm =T)/sqrt(length(PDCD1)),
                  SEM_CD274=sd(CD274,na.rm =T)/sqrt(length(CD274)),
                  .progress="text")

table_median_NT<-ddply(subset(mRNA_Exp_wide,sample_type=="NT"),.variables=.(cohort),
                    summarise,Median_PDCD1=median(PDCD1,na.rm =T),
                    Median_CD274=median(CD274,na.rm =T),
                    SEM_PDCD1=sd(PDCD1,na.rm =T)/sqrt(length(PDCD1)),
                    SEM_CD274=sd(CD274,na.rm =T)/sqrt(length(CD274)),
                    .progress="text")

table_median$Tumor<-"mixed"
table_median_NT$Tumor<-"no"
table_median_all<-rbind(table_median,
                        table_median_NT)

summary(table_median)
ggplot(subset(table_median_all,Tumor=="mixed"), aes(x=Median_CD274, y=Median_PDCD1, fill=Tumor))+
  stat_smooth(method = "lm", color="black") +
  stat_cor(label.y = 10,size = 7)+
  geom_point(shape = 21,size=3)+
  geom_errorbar(aes(ymin = Median_PDCD1+SEM_PDCD1,ymax = Median_PDCD1-SEM_PDCD1))+ 
  geom_errorbarh(aes(xmin = Median_CD274+SEM_CD274,xmax = Median_CD274-SEM_CD274))+
  theme_ds4psy()+ #coord_fixed()+
  geom_text_repel(aes(label=cohort))+ expand_limits(x = c(0, 11), y = c(0, 11))+ coord_fixed()
ggsave("20200210_PD-1_PD-L1_corr_cohort_only_tumors.pdf",
       units = "mm", width = 170, height = 130)


#Combine median data and corr data  
data_r_med<-merge(table_median,R_cor_df,by="cohort")

table_median$Combined_expession<-table_median$Median_PDCD1+table_median$Median_CD274
tp<-read.csv("tumor_type.csv")
table_median<-merge(table_median,tp,by="cohort")
table_median <- table_median[order(-table_median$Combined_expession),]
write.csv(table_median[,c(1:3,7,8)], "20200210_results_PD1_PD1L_kuchenov.CSV")
