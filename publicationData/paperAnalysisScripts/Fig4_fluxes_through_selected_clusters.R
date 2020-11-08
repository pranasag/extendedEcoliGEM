require(ggplot2)
require(reshape2)

resultsDir = '' #Path for the results directory
# Define correct paths for these two annotation files
met_clusters <- read.csv('MetabolicClusters_ExtendedModel.csv') 
cond_names <- read.csv('ConditionNames_Schmidt2016.csv')

conds = list.files(resultsDir, pattern = 'fvadata')
setwd(resultsDir)
for(i in seq(1, length(conds))) {
  fva_s <- read.csv(conds[i])
  fva_s$condition <- cond_names$name[grep(paste0("^",conds[i]), cond_names$csv)]
  
  if (i == 1) (compiled = fva_s) else (compiled = rbind(compiled, fva_s))}

compiled <- merge.data.frame(compiled, met_clusters, by.x = 'Jid', by.y = 'Reaction.ID')
compiled$min <- ifelse(is.na(compiled$min), compiled$optval, compiled$min)
compiled$max <- ifelse(is.na(compiled$max), compiled$optval, compiled$max)

clusters <- unique(compiled$Subsystem)
growth_conds <- unique(compiled$condition)
for(i in seq(1, length(clusters))) {
  rxn <- subset(compiled, compiled$Subsystem == clusters[i])
  rxn_clust <- data.frame(condition=growth_conds, cluster=clusters[i])
  for(j in seq(1, length(growth_conds))) {
    rxn_clust$opt[j] <- sum(abs(subset(rxn, rxn$condition == growth_conds[j])$optval))
    rxn_clust$min[j] <- sum(abs(subset(rxn, rxn$condition == growth_conds[j])$min))
    rxn_clust$max[j] <- sum(abs(subset(rxn, rxn$condition == growth_conds[j])$max))
  }

  if (i == 1) (glc_normalized <- rxn_clust) else (glc_normalized <- rbind(glc_normalized, rxn_clust))
}

glc_normalized <- unique(glc_normalized)
glc_normalized <- subset(glc_normalized, glc_normalized$condition %in% 
                           c('Ac', 'Fru', 'Fum', 'Gal', 'Glc', 'GlcN', 'Glyc', 'Glyc + AA', 'Man', 'Pyr', 'Succ', 'Xyl'))
glc_normalized$condition <- factor(glc_normalized$condition, 
                                   levels=c('Fru','Gal','Glc','Man','Ac','Fum','Pyr','Succ','GlcN','Glyc','Glyc + AA','Xyl', ordered=TRUE))
glc_normalized$cluster <- gsub('Citric Acid Cycle', 'Tricarboxylic Acid Cycle', glc_normalized$cluster)
preferred_cluster = c('Tricarboxylic Acid Cycle', 'Val, Leu, and Ile Metabolism', 'Pentose Phosphate Pathway', 
                      'Glycolysis/Gluconeogenesis', 'Glyoxylate Metabolism', 'Methylglyoxal Metabolism')
glc_cluster <- subset(glc_normalized, (glc_normalized$cluster %in% preferred_cluster))

my_breaks <- function(x) { if (max(x) > 1000) seq(0, 1600, 200) 
  else (if (max(x) < 20) (seq(0, 10, 2)) else (seq(0, 200, 50)))}
my_limits <- function(x) { if (max(x) > 1000) c(-80,1600) else (if (max(x) < 20) (c(-0.5, 10)) else (c(-10, 200)))}

g <- ggplot(glc_cluster, aes(ymin = min, ymax = max, y=min, x = condition, fill = condition)) + 
  geom_crossbar(position = position_dodge(1), fatten = 0) + 
  scale_y_continuous(limits = my_limits, breaks = my_breaks, expand = c(0,0)) +
  scale_x_discrete(breaks = NULL) + 
  facet_wrap(cluster ~ ., ncol = 6, scales = 'free_y') +
  labs(y = 'Sum of fluxes through\n metabolic cluster', 
       x = 'Growth medium',
       fill = 'Growth medium:') +
  theme(aspect.ratio = 1,
        plot.title = element_text(hjust = 0.5, size = 24), 
        axis.title.y = element_text(size = 24),
        strip.text.x = element_text(size = 18),
        strip.background = element_blank(),
        axis.text.y  = element_text(size = 20, margin = margin(r=5), hjust = 1),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 20, margin = margin(t=5), hjust = 0.5),
        axis.ticks = element_blank(),
        legend.title = element_text(size = 24, hjust = 0, vjust = 0.5, face='bold'),
        legend.text = element_text(size = 24, hjust = 0, vjust = 0.5),
        legend.key = element_rect(fill = NA),
        legend.position="bottom",
        legend.box = "horizontal",
        panel.background = element_rect(fill = NA),
        panel.border = element_rect(colour = "black", fill=NA, size=2),
        panel.grid.major = element_line(colour = "grey50"),
        plot.margin = margin(20, 20, 20, 20)) +
  guides(fill=guide_legend(nrow=4,byrow=FALSE))
plot(g)

