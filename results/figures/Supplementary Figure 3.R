#####################
# Import libraries
#####################
library(dplyr) 
library(ggplot2)
library("readxl")
library(patchwork)
library(rstatix)
library(scales)
library(stringr)
library(ggpubr)

#####################
# Import date
#####################
timesorg <- read_excel("/Users/Andreas/Desktop/Project/compare_profiles/Times/Times.xlsx", sheet="Merged")
profiles <- read_excel("/Users/Andreas/Desktop/ProfileStats.xlsx", sheet="ProfileSize")
#create new column to join on
timesorg <- transform(timesorg, Prof=Profile)
timesorg <- timesorg %>% mutate(across('Prof', str_replace, '_HeuristicsOff', '')) %>% mutate(across('Prof', str_replace, '_HeuristicsOn', ''))
timesorg <- transform(timesorg, ID=paste(Protein,"_",Prof,sep=''))
profiles <- transform(profiles, ID=paste(str_trim(Protein),"_",Profile,sep=''))

df <- merge(timesorg, profiles, by="ID", all.x=TRUE)
df2 <- df[,-which(names(df) %in% c("Organism","Start", "End", "Profile.y", "Protein.y", "Prof"))]

df3 <- df2 %>% group_by(Profile.x, Protein.x) %>% mutate(med=median(Time))
write.csv(df3,"/Users/Andreas/Desktop/df.csv",row.names=FALSE)

df4 <- profiles <- read_excel("/Users/Andreas/Desktop/Median.xlsx", sheet="Sheet1")



pssm<-ggplot(df4[df4$Profile.x == "PSSM",]) +
          geom_point(aes(x=QueryLength,y=med)) +
          geom_segment(aes(x=QueryLength, xend=QueryLength, y=med, yend=med)) +
  
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
  
  )+
  scale_y_continuous(breaks = seq(0, 3.5, by=0.5), limits=c(0,3.5))+
scale_x_continuous(breaks = seq(0, 2800, by=200), limits=c(0,2800)) +
  geom_smooth(method=lm, aes(x=QueryLength, y=med), se=TRUE, fullrange=FALSE, level=0.95) +
  ggtitle("PSSM")+ 
  labs(y = "Median of Time (s)",x="Query Length") +
  stat_regline_equation(aes(x=QueryLength, y=med, label=paste(..eq.label.., ..rr.label.., sep = "~~~~")))
  


NUPoff<-ggplot(df4[df4$Profile.x == "NUP-HMM_HeuristicsOff",]) +
  geom_point(aes(x=QueryLength,y=med)) +
  geom_segment(aes(x=QueryLength, xend=QueryLength, y=med, yend=med)) +
  
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
  )+
  scale_y_continuous(breaks = seq(0, 135, by=10), limits=c(0,135))+
  scale_x_continuous(breaks = seq(0, 2800, by=200), limits=c(0,2800)) +
  geom_smooth(method=lm, aes(x=QueryLength, y=med), se=TRUE, fullrange=FALSE, level=0.95) +
  ggtitle("NUP-HMM (Heuristics off)")+ 
  labs(y = "Median of Time (s)",x="Query Length") +
  stat_regline_equation(aes(x=QueryLength, y=med, label=paste(..eq.label.., ..rr.label.., sep = "~~~~")))

NUPon<-ggplot(df4[df4$Profile.x == "NUP-HMM_HeuristicsOn",]) +
  geom_point(aes(x=QueryLength,y=med)) +
  geom_segment(aes(x=QueryLength, xend=QueryLength, y=med, yend=med)) +
  
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
  )+
  scale_y_continuous(breaks = seq(0, 20, by=5), limits=c(0,20))+
  scale_x_continuous(breaks = seq(0, 2800, by=200), limits=c(0,2800)) +
  geom_smooth(method=lm, aes(x=QueryLength, y=med), se=TRUE, fullrange=FALSE, level=0.95) +
  ggtitle("NUP-HMM (Heuristics on)")+ 
  labs(y = "Median of Time (s)",x="Query Length") +
  stat_regline_equation(aes(x=QueryLength, y=med, label=paste(..eq.label.., ..rr.label.., sep = "~~~~")))

PFAMon<-ggplot(df4[df4$Profile.x == "PFAM-HMM_HeuristicsOn",]) +
  geom_point(aes(x=QueryLength,y=med)) +
  geom_segment(aes(x=QueryLength, xend=QueryLength, y=med, yend=med)) +
  
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
  )+
  scale_y_continuous(breaks = seq(0, 2.5, by=0.5), limits=c(0,2.5))+
  scale_x_continuous(breaks = seq(0, 1800, by=200), limits=c(0,1800)) +
  geom_smooth(method=lm, aes(x=QueryLength, y=med), se=TRUE, fullrange=FALSE, level=0.95) +
  ggtitle("PFAM-HMM (Heuristics on)")+
  stat_regline_equation(aes(x=QueryLength, y=med, label=paste(..eq.label.., ..rr.label.., sep = "~~~~")))

PFAMoff<-ggplot(df4[df4$Profile.x == "PFAM-HMM_HeuristicsOff",]) +
  geom_point(aes(x=QueryLength,y=med)) +
  geom_segment(aes(x=QueryLength, xend=QueryLength, y=med, yend=med)) +
  
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),

  )+
  scale_y_continuous(breaks = seq(0, 60, by=10), limits=c(0,60))+
  scale_x_continuous(breaks = seq(0, 1800, by=200), limits=c(0,1800)) +
  geom_smooth(method=lm, aes(x=QueryLength, y=med), se=TRUE, fullrange=FALSE, level=0.95) +
  stat_regline_equation(aes(x=QueryLength, y=med, label=paste(..eq.label.., ..rr.label.., sep = "~~~~"))) +
  ggtitle("PFAM-HMM (Heuristics off)")+ 
  labs(y = "Median of Time (s)",x="Query Length")





layout <- "
##BBCC
AABBCC
AADDEE
##DDEE
"
patch <- pssm + PFAMoff + PFAMon + NUPoff + NUPon + plot_layout(design=layout) 
tiff("Supplementary Figure 3.tiff", units="in", width=10, height=8, res=300)
patch
dev.off()
