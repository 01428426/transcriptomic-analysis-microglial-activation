library(dplyr)
library(tidyr)
require(data.table)
library(DescTools)
library(ggplot2)


files = list.files(path = "/rds/general/project/hda_students_data/live/Group7/General/Eleonore/thesis/new_aim11/samples_umap1_500", 
                   pattern = "clara_res", full.names = TRUE)
files=files[order(as.numeric(gsub("[^0-9]", "",files )),decreasing = FALSE)]
dat_list = lapply(files, function (x) data.table(readRDS(x)))

# optimising best k ranging 2-15
clara=list()
for(num in seq(1, 14, by=1)){
  clara[[num]]=as.numeric(dat_list[[num]][["V1"]][[9]][["widths"]][,3])
}

# Median and 95% CI
CI=list()
for(list in seq(1, 14, by=1)){
  CI[[list]]=MedianCI(clara[[list]])
}

df=as.data.frame((CI))
colnames(df)=as.character(seq(2, 15, by=1))
df=t(df)
df=data.frame(df)
colnames(df)=c('median','lwr.ci','upr.ci') # median, low CI, high CI
df$names=rownames(df)

# Find thresholds
threshold <- max(df['lwr.ci']) #max lower confidence interval
choice=as.character(max(as.numeric(rownames(df[which(df$median>=threshold),]))))

clara_df=plyr::ldply(clara, rbind)
clara_df=data.frame(t(clara_df))
colnames(clara_df)=as.character(seq(2, 15, by=1))
data=clara_df%>% pivot_longer(everything()) 

# Plotting avg silhouette widths 
avg_sil=list()
for(num in seq(1, 14, by=1)){
  avg_sil[[num]]=as.numeric(dat_list[[num]][["V1"]][[9]][["clus.avg.widths"]])
}

avg_sil_df=plyr::ldply(avg_sil, rbind)
avg_sil_df=data.frame(t(avg_sil_df))
colnames(avg_sil_df)=as.character(seq(2, 15, by=1))

data2=avg_sil_df%>% pivot_longer(everything()) 


# Plotting Silhouette distribution

ggplot(df, aes(reorder(factor(names)), median)) +
  geom_crossbar(
    aes(ymin = lwr.ci, ymax = upr.ci),
    fill = "grey",
    size = 0.25
  )+
  xlab('k')+
  ylab('Silhouette Score')+
  geom_hline(aes(yintercept = threshold), colour = "blue")+
  geom_vline(aes(xintercept = choice), colour = "red")+
  scale_y_continuous(
    "Silhouette Score",
    expand = c(0, 0),
    limits = c(0, 1),
    breaks = seq(-1, 1, 0.25),
    oob = scales::squish
  ) +
  cowplot::theme_minimal_hgrid() +
  theme(
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 7),
    axis.line.x = element_line(colour = "black"),
    axis.line.y = element_line(colour = "black"),
    axis.ticks = element_line(colour = "black"),
  ) +  geom_jitter(
    data = data2,
    aes(x=reorder(factor(name)),y=value),
    size = 0.05,
    width = 0.1
  ) 
  


ggsave(
  filename ="/rds/general/project/hda_students_data/live/Group7/General/Eleonore/thesis/new_aim11/samples_umap1_500/silhouette_distribution_plot.png",
  dpi = 300,
  height = 3.5,
  width = 3.5,
  units = "in"
)

