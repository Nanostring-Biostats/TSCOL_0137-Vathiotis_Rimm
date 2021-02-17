##########################################
# Collaboration: TSCOL_0137
#   Senior PI: David Rimm
#   Lead Investigator: Ioannis Vathiotis
# Jason Reeves

# Copyright (C) 2020, NanoString Technologies, Inc.
#
# This program is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public License along with
# this program.  If not, see https://www.gnu.org/licenses/.

# Inputs needed: graphing_adjusted_log2_expression_data_WITH DSP.xlsx
##########################################

##########################################
# Libraries
##########################################
library(readxl)
library(ggplot2)
library(ggrepel)
library(pheatmap)

##########################
#### Custom Functions ####

# mean_geo: compute the geometric mean of linear count data
mean_geo <- function(x,
                     na.rm=TRUE) {
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

######################
#### Data Loading ####

#### 0.1: Load data ####
annot <- as.data.frame(read_xlsx(path = 'data/graphing_adjusted_log2_expression_data_WITH DSP.xlsx',
                                 sheet = 5))
targetannot <- as.data.frame(read_xlsx(path = 'data/graphing_adjusted_log2_expression_data_WITH DSP.xlsx',
                                       sheet = 6))
targetannot <- targetannot[,1:7]
data <- as.data.frame(read_xlsx(path = 'data/graphing_adjusted_log2_expression_data_WITH DSP.xlsx',
                                sheet = 2))
rownames(data) <- data$Target
data <- data[,-1]

#### 0.2: transform DSP data & update annotations ####
data[targetannot$Overall == 'DSP',] <- log2(data[targetannot$Overall == 'DSP',])
annot$Response.NS <- factor(annot$Response.NS, levels = c('NR','R'))
targetannot$run <- ifelse(grepl('2nd_run_', rownames(data)),'Run 2',
                          ifelse(targetannot$Overall == 'IO360',
                                 'Bulk','Run 1'))

#### 0.3: mean center expression from DSP runs ####
data_orig <- data
targetannot_orig <- targetannot
dsp_genes <- sort(targetannot$Gene[targetannot$run == 'Run 1' & targetannot$Source == 'CD45'])
for(comp in c('CD45','CD68','Melanocyte')) {
  mean_df <- matrix(data = NA, nrow = length(dsp_genes), ncol = ncol(data),
                    dimnames = list(row = dsp_genes, col = colnames(data)))
  for(g in dsp_genes) {
    r1 <- data[targetannot$Gene == g & targetannot$run == 'Run 1' & targetannot$Source == comp, ]
    r2 <- data[targetannot$Gene == g & targetannot$run == 'Run 2' & targetannot$Source == comp, ]
    rep_mean <- colMeans(rbind(r1, r2), na.rm = T)
    mean_df[g, ] <- rep_mean
  }
  rownames(mean_df) <- paste0('Mean_',dsp_genes,'_',comp)
  mean_ann <- targetannot[targetannot$Source == comp & targetannot$run == 'Run 1', ]
  mean_ann$run <- 'Mean DSP'
  mean_ann$Gene <- dsp_genes
  mean_ann$Target <- paste0('Mean_',dsp_genes,'_',comp)
  targetannot <- rbind(targetannot, mean_ann)
  data <- rbind(data, mean_df)
}

rm(list = c('mean_df','mean_ann','r1','r2','rep_mean'))

##############################################
#### Section 1: DE of IO 360 and DSP data ####

#### 1.1 Response Analysis (Univariate) ####
# build data tables and run univariate analysis
out_R <- data.frame(target = targetannot$Target)
out_R$pval <- unlist(apply(data, 1, function(x) {
  mod <- lm(formula = unlist(x) ~ annot$Response.NS)
  coef(summary(mod))[2,4]
}))
out_R$FC <- unlist(apply(data, 1, function(x) {
  mod <- lm(formula = unlist(x) ~ annot$Response.NS)
  coef(summary(mod))[2,1]
}))

# merge with annotations
out_R <- merge(out_R, targetannot, by.x = 1, by.y = 1)
out_R$Source <- factor(out_R$Source, levels = c('IO360','CD45','CD68','Melanocyte'))
# remove technical controls
out_R <- subset(out_R, !grepl('IgG', Gene))

#### 1.2 Graph Objective Response DE Results ####
OR_VP <- ggplot(subset(out_R, run %in% c('Bulk','Mean DSP')),
                aes(x = FC, y = -log10(pval),
                    color = Source, size = -log10(pval),
                    label = Gene)) +
  geom_point(alpha = 0.25) +
  theme_bw(base_size = 17) +
  labs(title = '', y = 'Significance, -log10(P)', x = 'Fold Change, log2(FC)') +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  geom_text_repel(data = subset(out_R, run %in% c('Bulk','Mean DSP') &
                                  (FC < -0.5 | FC > 1 | pval < 0.01)),
                  color = 'black', size = 5, fontface = 'bold',
                  segment.size = 1, box.padding = .45, min.segment.length = .1,
                  point.padding = .3) +
  scale_color_manual(values = c(CD45 = '#DD0000',
                                CD68 = '#FF00FF',
                                Melanocyte = '#009000',
                                IO360 = '#000080')) +
  guides(color = guide_legend(title = 'Analyte', override.aes = list(size = 5)), size = FALSE) +
  geom_hline(yintercept = -log10(0.05), lty = 'dashed', color = 'black') +
  annotate(geom = 'text', x = 1.6, y = 1.4, label = 'P = 0.05', size = 5) +
  labs(tag = "a") +
  theme(legend.position = c(0.15,0.85), aspect.ratio = 0.9,
        legend.background = element_rect(color = 'darkgray', fill = 'white'), 
        plot.tag = element_text(size = 12, face = "bold")) 

write.csv(out_R[which(out_R$FC<=0), ], file = 'output/negative_predictors.csv', row.names = FALSE)

#### 1.3 Save DE Results ####
dir.create("figs/")
dir.create("figs/jpg")
dir.create("figs/eps")

ggsave(OR_VP, filename = 'figs/jpg/Figure2a.jpg', width = 8, height = 8, dpi = 300)
ggsave(OR_VP, filename = 'figs/eps/Figure2a.eps', width = 8, height = 8, dpi = 300, device = cairo_ps)

#############################
#### 2.0 Correlation Map ####

#### 2.1: correlation heatmap of genes * genes ####
gene_cor <- cor(t(data[targetannot$run %in% c('Bulk','Mean DSP'), ]), method = 'spearman', use = 'complete')
gene_ann <- data.frame(Source = as.character(targetannot$Source)[targetannot$run %in% c('Bulk','Mean DSP')],
                       D = as.character(targetannot$run)[targetannot$run %in% c('Bulk','Mean DSP')])
gene_ann$Data <- ifelse(gene_ann$D == 'Bulk', 'Bulk RNA', 'GeoMx DSP')
rownames(gene_ann) <- rownames(gene_cor)
cor_heat <- pheatmap(gene_cor,
         annotation_col = gene_ann[,c('Data','Source')],
         show_rownames = F, show_colnames = F,
         clustering_method = 'average',
         annotation_colors = list(Data = c(`Bulk RNA` = '#000080',`GeoMx DSP` = '#F5C242'),
                                  Source = c(CD45 = '#DD0000', CD68 = '#FF00FF',
                                             IO360 = '#000080', Melanocyte = '#009000')))

#### 2.2: Save Files ####
jpeg(filename = 'figs/jpg/Figure1a.jpg', width = 8, height = 6, units = 'in', res = 300)
grid::grid.newpage()
grid::grid.draw(cor_heat$gtable)
grid::grid.text(expression(bold("a")), x = 0.01, y = 0.99)
dev.off()

setEPS()
postscript('figs/eps/Figure1a.eps', width = 8, height = 6, bg = "white")
grid::grid.newpage()
grid::grid.draw(cor_heat$gtable)
grid::grid.text(expression(bold("a")), x = 0.01, y = 0.99)
dev.off()

########################################
#### 3: Save Workspace Environments ####

out_Surv <- out_R
rm(list = c('comp','cor_heat','dsp_genes','g','gene_ann','gene_cor','mean_geo',
            'OR_VP','out_R','data_orig','targetannot_orig'))
save.image(file = 'data/TSCOL_0137_Initialized.Rdata')
