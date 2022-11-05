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

                                          #####################
                                          # Per nup family
                                          #####################
#####################
# Import data
#####################
times <- read_excel("/Users/Andreas/Desktop/Project/compare_profiles/Times/Execution_times.xlsx", sheet="Merged")
df = subset(times, select = -c(Start,End) )

#####################
# Remove outliers
#####################
df2 <- df %>%
  group_by(Protein, Profile) %>%
  #mutate(Time = remove_outliers(Time))
  identify_outliers("Time") 
df_wo_outliers <- 
  df %>% 
  anti_join(df2, by = c("Run","Organism","Protein","Profile")) 
df3 <-df_wo_outliers
#####################
# Create the plots
#####################
e1 <- ggplot(df3[df3$Profile == "PSSM",], aes(x = Protein, y = Time)) +
  geom_boxplot(aes(fill = Profile), outlier.shape = NA,
  position = position_dodge(0.9) ) +
  scale_fill_manual(values = c("#999999", "#E69F00")) + 
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5),
        #axis.title.x = element_blank()
        ) +
  scale_y_continuous(breaks = seq(0, 4, by=0.5), limits=c(0,4)) +
  ggtitle("PSSM")+ 
  labs(y = "Time (s)")  + 
  coord_flip()

e2 <- ggplot(df3 %>% filter(str_detect(Profile,"PFAM-HMM_HeuristicsOff")) , aes(x = Protein, y = Time)) +
  geom_boxplot(aes(fill = Profile), outlier.shape = NA,
               position = position_dodge(0.9)) +
  scale_fill_manual(values = c("#999999", "#E69F00")) + 
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5),
        axis.text.y=element_blank(),  
        axis.ticks.y=element_blank(),
        axis.title.y = element_blank(),
        #axis.title.x = element_blank()
        )+
  scale_y_continuous(breaks = seq(0, 90, by=10), limits=c(0,90)) +
  ggtitle("PFAM-HMM (Heuristics off)")+ 
  labs(y = "Time (s)")  + 
  coord_flip()

e3 <- ggplot(df3 %>% filter(str_detect(Profile,"PFAM-HMM_HeuristicsOn")) , aes(x = Protein, y = Time)) +
  geom_boxplot(aes(fill = Profile), outlier.shape = NA,
               position = position_dodge(0.9)) +
  scale_fill_manual(values = c("#999999", "#E69F00")) + 
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5),
        axis.text.y=element_blank(),  
        axis.ticks.y=element_blank(),
        axis.title.y = element_blank(),
        #axis.title.x = element_blank()
  )+
  scale_y_continuous(breaks = seq(0, 4.5, by=0.5), limits=c(0,4.5)) +
  ggtitle("PFAM-HMM (Heuristics on)")+ 
  labs(y = "Time (s)")  + 
  coord_flip()

e4 <- ggplot(df3 %>% filter(str_detect(Profile,"NUP-HMM_HeuristicsOff")) , aes(x = Protein, y = Time)) +
  geom_boxplot(aes(fill = Profile), outlier.shape = NA,
               position = position_dodge(0.9)) +
  scale_fill_manual(values = c("#999999", "#E69F00")) + 
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5),
        #axis.text.y=element_blank(),  
        #axis.ticks.y=element_blank(),
        axis.title.y = element_blank(),
        #axis.title.x = element_blank()
  )+
  scale_y_continuous(breaks = seq(0, 90, by=10), limits=c(0,90)) +
  ggtitle("NUP-HMM (Heuristics off)")+ 
  labs(y = "Time (s)")  + 
  coord_flip()

e5 <- ggplot(df3 %>% filter(str_detect(Profile,"NUP-HMM_HeuristicsOn")) , aes(x = Protein, y = Time)) +
  geom_boxplot(aes(fill = Profile), outlier.shape = NA,
               position = position_dodge(0.9)) +
  scale_fill_manual(values = c("#999999", "#E69F00")) + 
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5),
        axis.text.y=element_blank(),  
        axis.ticks.y=element_blank(),
        axis.title.y = element_blank(),
        #axis.title.x = element_blank()
  )+
  scale_y_continuous(breaks = seq(0, 4.5, by=0.5), limits=c(0,4.5)) +
  ggtitle("NUP-HMM (Heuristics on)")+ 
  labs(y = "Time (s)")  + 
  coord_flip()

layout <- "
##BBCC
AABBCC
AADDEE
##DDEE
"

#####################
# Merge plots and export PSSM/PFAM
#####################
patch <- e1 + e2 + e3
tiff("Figure 3.tiff", units="in", width=10, height=5, res=300)
patch
dev.off()

#####################
# Merge plots and export PSSM/PFAM/NUPHMM
#####################

# Do not hide Y axis for this plot
e2 <- ggplot(df3 %>% filter(str_detect(Profile,"PFAM-HMM_HeuristicsOff")) , aes(x = Protein, y = Time)) +
  geom_boxplot(aes(fill = Profile), outlier.shape = NA,
               position = position_dodge(0.9)) +
  scale_fill_manual(values = c("#999999", "#E69F00")) + 
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5),
        #axis.text.y=element_blank(),  
        #axis.ticks.y=element_blank(),
        axis.title.y = element_blank(),
        #axis.title.x = element_blank()
  )+
  scale_y_continuous(breaks = seq(0, 90, by=10), limits=c(0,90)) +
  ggtitle("PFAM-HMM (Heuristics off)")+ 
  labs(y = "Time (s)")  + 
  coord_flip()

patch <- e1 + e2 + e3 + e4 + e5 + plot_layout(design=layout) 
tiff("Supplementary Figure 1.tiff", units="in", width=10, height=8, res=300)
patch
dev.off()
###########################################################################

                                            #####################
                                            # Per organism
                                            #####################
#####################
# Import data
#####################
timesorg <- read_excel("/Users/Andreas/Desktop/Project/compare_profiles/Times/Execution_times.xlsx", sheet="Organisms")
df3 <- timesorg
#####################
# Create the plots
#####################
p1 <- ggplot(df3[df3$Profile == "PSSM",], aes(x = Organism, y = Time)) +
  geom_boxplot(aes(fill = Profile), outlier.shape = NA,
               position = position_dodge(0.9) ) +
  scale_fill_manual(values = c("#999999", "#E69F00")) + 
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)
  ) +
  scale_y_continuous(breaks = seq(0, 35, by=5), limits=c(0,35)) +
  ggtitle("PSSM")+ 
  labs(y = "Time (s)")  + 
  coord_flip()

p2 <- ggplot(df3 %>% filter(str_detect(Profile,"PFAM-HMM_HeuristicsOff")), aes(x = Organism, y = Time)) +
  geom_boxplot(aes(fill = Profile), outlier.shape = NA,
               position = position_dodge(0.9)) +
  scale_fill_manual(values = c("#999999", "#E69F00")) + 
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5),
        axis.text.y=element_blank(),  
        axis.ticks.y=element_blank(),
        axis.title.y = element_blank()
  )+
  scale_y_continuous(breaks = seq(0, 1500, by=300), limits=c(0,1500)) +
  ggtitle("PFAM-HMM  (Heuristics off)")+ 
  labs(y = "Time (s)") + 
  coord_flip()

p3 <- ggplot(df3 %>% filter(str_detect(Profile,"PFAM-HMM_HeuristicsOn")), aes(x = Organism, y = Time)) +
  geom_boxplot(aes(fill = Profile), outlier.shape = NA,
               position = position_dodge(0.9)) +
  scale_fill_manual(values = c("#999999", "#E69F00")) + 
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5),
        axis.text.y=element_blank(),  
        axis.ticks.y=element_blank(),
        axis.title.y = element_blank()
  )+
  scale_y_continuous(breaks = seq(0, 15, by=5), limits=c(0,15)) +
  ggtitle("PFAM-HMM (Heuristics on)")+ 
  labs(y = "Time (s)") + 
  coord_flip()


p4 <- ggplot(df3 %>% filter(str_detect(Profile,"NUP-HMM_HeuristicsOff")), aes(x = Organism, y = Time)) +
  geom_boxplot(aes(fill = Profile), outlier.shape = NA,
               position = position_dodge(0.9)) +
  scale_fill_manual(values = c("#999999", "#E69F00")) + 
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5),
        #axis.text.y=element_blank(),  
        #axis.ticks.y=element_blank(),
        axis.title.y = element_blank()
  )+
  scale_y_continuous(breaks = seq(0, 2500, by=300), limits=c(0,2500)) +
  ggtitle("NUP-HMM  (Heuristics off)")+ 
  labs(y = "Time (s)") + 
  coord_flip()

p5 <- ggplot(df3 %>% filter(str_detect(Profile,"NUP-HMM_HeuristicsOn")), aes(x = Organism, y = Time)) +
  geom_boxplot(aes(fill = Profile), outlier.shape = NA,
               position = position_dodge(0.9)) +
  scale_fill_manual(values = c("#999999", "#E69F00")) + 
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5),
        axis.text.y=element_blank(),  
        axis.ticks.y=element_blank(),
        axis.title.y = element_blank()
  )+
  scale_y_continuous(breaks = seq(0, 500, by=100), limits=c(0,500)) +
  ggtitle("NUP-HMM (Heuristics on)")+ 
  labs(y = "Time (s)") + 
  coord_flip()


layout2 <- "
##BBCC
AABBCC
AADDEE
##DDEE
"

#####################
# Merge plots and export PSSM/PFAM
#####################
patch2 <- p1 + p2 + p3
tiff("Figure 4.tiff", units="in", width=10, height=5,  res=300)
patch2
dev.off()

#####################
# Merge plots and export PSSM/PFAM/NUPHMM
#####################

# Do not hide Y axis for this plot
p2 <- ggplot(df3 %>% filter(str_detect(Profile,"PFAM-HMM_HeuristicsOff")), aes(x = Organism, y = Time)) +
  geom_boxplot(aes(fill = Profile), outlier.shape = NA,
               position = position_dodge(0.9)) +
  scale_fill_manual(values = c("#999999", "#E69F00")) + 
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5),
        #axis.text.y=element_blank(),  
        #axis.ticks.y=element_blank(),
        axis.title.y = element_blank()
  )+
  scale_y_continuous(breaks = seq(0, 1500, by=300), limits=c(0,1500)) +
  ggtitle("PFAM-HMM  (Heuristics off)")+ 
  labs(y = "Time (s)") + 
  coord_flip()

patch2 <- p1 + p2 + p3 + p4 + p5 + plot_layout(design=layout2) 
tiff("Supplementary Figure 2.tiff", units="in", width=10, height=5,  res=300)
patch2
dev.off()
#############################################
                  