rm(list=ls())

## Libraries

library(tidyverse)
library(ggplot2)
library(ggpubr)
library(corrplot)
library(ggcorrplot)
library(ggstatsplot)


#### FUNCTIONS

## Boxplot

box.plot = function(df, axis.x, axis.y, title.x, title.y, m.y, n){
  ggboxplot(df, x = axis.x, y = axis.y,
            color = axis.x, palette =c("skyblue", "tomato", "#E7B801"),
            add = "jitter", shape = axis.x) + 
    stat_compare_means(method = "wilcox.test", label.y = m.y + n) + # Add significance levels
    stat_summary(fun=mean, geom="point", shape=20, size=2, color="black", fill="black") +
    scale_y_continuous(name = title.y) +
    scale_x_discrete(name = title.x) +
    theme (axis.text.x=element_text(colour="black", face = "bold", size = 9),
           axis.text.y=element_text(colour="black", size = 9)) + 
    theme(axis.line = element_line(size = 0.5, linetype = "solid"),
          axis.title = element_text(face = "bold", size = 10)) 
}


## Density plot

density.plot = function(df, axis.x, title.x){
  ggdensity(df, x = axis.x,
            color = "msi", fill = "msi", 
            palette =c("skyblue", "tomato"),
            add = "median") + 
    scale_y_continuous(name = "") +
    scale_x_discrete(name = title.x) +
    theme (axis.text.y=element_text(colour="black", size = 9)) + 
    theme(axis.line = element_line(size = 0.5, linetype = "solid"),
          axis.title = element_text(face = "bold", size = 12)) 
}


#### READING FILES

ucec.table = read.table("./data/ucec.final.no.pol.txt",
                        header = T, sep = "\t", quote = "", fill = T)


## Selecting only molecular data

ucec.molecular = select(ucec.table, c(1, 13:27))
ucec.molecular = select(ucec.molecular, -frameshift.fraction)

## Transforming data

ucec.molecular$msi = factor(ucec.molecular$msi,
                            levels = c(0,1),
                            labels = c("MSS", "MSI"))

#### PLOT PARAMETERS

max.mut.load = round(max(ucec.molecular$mutation.load)) + 15
max.nonsynonymous = round(max(ucec.molecular$nonsynonymous)) + 15
max.frameshift = round(max(ucec.molecular$frameshift)) + 15
max.snp = round(max(ucec.molecular$snp)) + 15
max.indel = round(max(ucec.molecular$indel)) + 15
max.indel.fraction = round(max(ucec.molecular$indel.snp.ratio))

max.tmb = round(max(ucec.molecular$tmb)) + 2
max.c.a = round(max(ucec.molecular$C.A)) + 5
max.c.g = round(max(ucec.molecular$C.G)) + 5
max.c.t = round(max(ucec.molecular$C.T)) + 5
max.t.c = round(max(ucec.molecular$T.C)) + 5
max.t.a = round(max(ucec.molecular$T.A)) + 5
max.t.g = round(max(ucec.molecular$T.G)) + 5

#### BOXPLOTS

# Mutation burden

boxplot.mutation.load = box.plot(df = ucec.molecular,
                                 axis.x = "msi", axis.y = "mutation.load", 
                                 title.x = "", title.y = "Mutation Load",
                                 m.y = max.mut.load,
                                 n = 1000)
# Nonsynonymous


boxplot.nonsyn = box.plot(df = ucec.molecular,
                                 axis.x = "msi", axis.y = "nonsynonymous", 
                                 title.x = "", title.y = "Nonsynonymous Mutations",
                                 m.y = max.nonsynonymous, n = 1000)

# Frameshifts


boxplot.frameshift = box.plot(df = ucec.molecular,
                          axis.x = "msi", axis.y = "frameshift", 
                          title.x = "", title.y = "Frameshift Mutations",
                          m.y = max.frameshift, n = 30)

# Snps


boxplot.snp = box.plot(df = ucec.molecular,
                              axis.x = "msi", axis.y = "snp", 
                              title.x = "", title.y = "SNP Mutations",
                              m.y = max.snp, n = 1500)

# Indels


boxplot.indel = box.plot(df = ucec.molecular,
                       axis.x = "msi", axis.y = "indel", 
                       title.x = "", title.y = "InDel Mutations",
                       m.y = max.indel, n = 200)


# Indel/SNP

boxplot.indel.ratio = box.plot(df = ucec.molecular,
                         axis.x = "msi", axis.y = "indel.snp.ratio", 
                         title.x = "", title.y = "Indel/SNP Ratio",
                         m.y = max.indel.fraction, n = 0.001)

# TMB

boxplot.tmb = box.plot(df = ucec.molecular,
                               axis.x = "msi", axis.y = "tmb", 
                               title.x = "", title.y = "TMB",
                               m.y = max.tmb, n = 20)



## Merge graphs

ggarrange(boxplot.tmb,
           boxplot.mutation.load,
           boxplot.nonsyn,
           boxplot.frameshift,
           boxplot.indel.ratio,
           boxplot.snp,
           boxplot.indel,
           labels = c("A", "B", "C", "D", "E", "F", "G"),
           ncol = 4,
           nrow = 2,
           common.legend = TRUE, legend = "right",
           font.label = list(size = 16, color = "black", face = "bold"))


#### DENSITY PLOTS

substitutions = select(ucec.molecular, c(1, 9:15))

density.c.a = density.plot(df = substitutions,
                           axis.x = "C.A",
                           title.x = "C>A")

density.c.t = density.plot(df = substitutions,
                           axis.x = "C.T",
                           title.x = "C>T")

density.c.g = density.plot(df = substitutions,
                           axis.x = "C.G",
                           title.x = "C>G")

density.t.a = density.plot(df = substitutions,
                           axis.x = "T.A",
                           title.x = "T>A")

density.t.c = density.plot(df = substitutions,
                           axis.x = "T.C",
                           title.x = "T>C")

density.t.g = density.plot(df = substitutions,
                           axis.x = "T.G",
                           title.x = "T>G")


ggarrange(density.c.a,
          density.c.g,
          density.c.t,
          density.t.a,
          density.t.c,
          density.t.g,
          labels = c("A", "B", "C", "D", "E", "F"),
          ncol = 3,
          nrow = 2,
          common.legend = TRUE, legend = "bottom",
          font.label = list(size = 16, color = "black", face = "bold"))


#### CORRELATION MATRIX

cor.matrix <- cor(ucec.molecular[, 2:14])
sig.matrix <- cor_pmat(cor.matrix)

ggcorrplot(cor.matrix,
           legend.title = "Correlation",
           outline.color = "white",
           lab = TRUE,
           lab_size = 3,
           p.mat = sig.matrix,
           sig.level = 0.05,
           pch = 4,
           pch.col = "darkslategray4",
           pch.cex = 4)