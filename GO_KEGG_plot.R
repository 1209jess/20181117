rm(list=ls())

# Libraries
library(ggplot2)
library(dplyr)
library(stringr)
library(tibble)
library(ggpubr)
library(patchwork)
library(ggrepel)


#read data
df <- read.delim("ALA_24hr_cluster(-489)_GOTERM_BP.txt")
colnames(df)[10]  <- "Fold"
df$Term <- str_sub(df$Term, start = 12)

# Add new column from existing column
log_pvalue <- -log(df$PValue)
df2 <- df %>%
  add_column("log_pvalue" = log_pvalue)
df3 <- df2[order(df2$log_pvalue,decreasing=TRUE),]

#select the column needed
df4 <- df3[,c(2,3,10,14)]

#plot chart
p1 <- ggplot(df4, aes(x = Fold, y = reorder(Term, log_pvalue), size = Count,colour = log_pvalue)) +
  geom_point()+
  theme_bw()+
  guides(size = guide_legend(order = 1),log_pvalue = guide_legend(order = 2))+
  labs(title="ALA_24hr_cluster(-489)_GOTERM_BP", x="Fold Enrichment", y="Term", size="# genes", colour="-log(Pvalue)")+
  scale_colour_gradient(guide=guide_colourbar(reverse = FALSE),
                        low="blue", high="red")+
  scale_y_discrete(labels = function(y) str_wrap(y, width = 45))+
  scale_x_continuous(expand = c(0.1, 0))+
  theme(aspect.ratio = 2, plot.title = element_text(size=12,face = "bold",hjust = 0.5),
        axis.text.y = element_text(colour = "black", size = 11, vjust = 0.5),
        axis.text.x = element_text(colour = "black", size = 11, hjust = 0,
                                   margin = margin(r = 5)),
        legend.key.size = unit(0.5, 'cm'))
geom_point()
p1


#read data
df_2 <- read.delim("ALA_24hr_cluster(-489)_GOTERM_CC.txt")
colnames(df_2)[10]  <- "Fold"
df_2$Term <- str_sub(df_2$Term, start = 12)

# Add new column from existing column
log_pvalue_2 <- -log(df_2$PValue)
df2_2 <- df_2 %>%
  add_column("log_pvalue_2" = log_pvalue_2)
df3_2 <- df2_2[order(df2_2$log_pvalue_2,decreasing=TRUE),]

#select the column needed
df4_2 <- df3_2[,c(2,3,10,14)]

#plot chart
p2 <- ggplot(df4_2, aes(x = Fold, y = reorder(Term, log_pvalue_2), size = Count,colour = log_pvalue_2)) +
  geom_point()+
  theme_bw()+
  guides(size = guide_legend(order = 1),log_pvalue_2 = guide_legend(order = 2))+
  labs(title="ALA_24hr_cluster(-489)_GOTERM_CC", x="Fold Enrichment", y="Term", size="# genes", colour="-log(Pvalue)")+
  scale_colour_gradient(guide=guide_colourbar(reverse = FALSE),
                        low="blue", high="red")+
  scale_y_discrete(labels = function(y) str_wrap(y, width = 45))+
  scale_x_continuous(expand = c(0.1, 0))+
  theme(aspect.ratio = 2, plot.title = element_text(size=12,face = "bold",hjust = 0.5),
        axis.text.y = element_text(colour = "black", size = 11, vjust = 0.5),
        axis.text.x = element_text(colour = "black", size = 11, hjust = 0,
                                   margin = margin(r = 5)),
        legend.key.size = unit(0.5, 'cm'))
geom_point()
p2


#read data
df_3 <- read.delim("ALA_24hr_cluster(-489)_GOTERM_MF.txt")
colnames(df_3)[10]  <- "Fold"
df_3$Term <- str_sub(df_3$Term, start = 12)

# Add new column from existing column
log_pvalue_3 <- -log(df_3$PValue)
df2_3 <- df_3 %>%
  add_column("log_pvalue_3" = log_pvalue_3)
df3_3 <- df2_3[order(df2_3$log_pvalue_3,decreasing=TRUE),]

#select the column needed
df4_3 <- df3_3[,c(2,3,10,14)]

#plot chart
p3 <- ggplot(df4_3, aes(x = Fold, y = reorder(Term, log_pvalue_3), size = Count,colour = log_pvalue_3)) +
  geom_point()+
  theme_bw()+
  guides(size = guide_legend(order = 1),log_pvalue_3 = guide_legend(order = 2))+
  labs(title="ALA_24hr_cluster(-489)_GOTERM_MF", x="Fold Enrichment", y="Term", size="# genes", colour="-log(Pvalue)")+
  scale_colour_gradient(guide=guide_colourbar(reverse = FALSE),
                        low="blue", high="red")+
  scale_y_discrete(labels = function(y) str_wrap(y, width = 45))+
  scale_x_continuous(expand = c(0.1, 0))+
  theme(aspect.ratio = 2, plot.title = element_text(size=12,face = "bold",hjust = 0.5),
        axis.text.y = element_text(colour = "black", size = 11, vjust = 0.5),
        axis.text.x = element_text(colour = "black", size = 11, hjust = 0,
                                   margin = margin(r = 5)),
        legend.key.size = unit(0.5, 'cm'))
geom_point()
p3


#read data
df_4 <- read.delim("ALA_24hr_cluster(-489)_KEGG.txt")
colnames(df_4)[10]  <- "Fold"
df_4$Term <- str_sub(df_4$Term, start = 10)

# Add new column from existing column
log_pvalue_4 <- -log(df_4$PValue)
df2_4 <- df_4 %>%
  add_column("log_pvalue_4" = log_pvalue_4)
df3_4 <- df2_4[order(df2_4$log_pvalue_4,decreasing=TRUE),]

#select the column needed
df4_4 <- df3_4[,c(2,3,10,14)]

#plot chart
p4 <- ggplot(df4_4, aes(x = Fold, y = reorder(Term, log_pvalue_4), size = Count, colour = log_pvalue_4)) +
  geom_point()+
  theme_bw()+
  guides(size = guide_legend(order = 1),log_pvalue_4 = guide_legend(order = 2))+
  labs(title="ALA_24hr_cluster(-489)_KEGG", x="Fold Enrichment", y="Term", size="# genes", colour="-log(Pvalue)")+
  scale_colour_gradient(guide=guide_colourbar(reverse = FALSE),
                        low="blue", high="red")+
  scale_y_discrete(labels = function(y) str_wrap(y, width = 45))+
  scale_x_continuous(expand = c(0.1, 0))+
  theme(aspect.ratio = 2, plot.title = element_text(size=12,face = "bold",hjust = 0.5),
        axis.text.y = element_text(colour = "black", size = 11, vjust = 0.5),
        axis.text.x = element_text(colour = "black", size = 11, hjust = 0,
                                   margin = margin(r = 5)),
        legend.key.size = unit(0.5, 'cm'))
geom_point()
p4



#combine the graph
p1+p2+p3+p4
p1
p2
p3
p4
