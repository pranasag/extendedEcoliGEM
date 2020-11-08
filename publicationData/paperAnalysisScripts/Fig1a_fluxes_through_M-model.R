require(ggplot2)
require(reshape2)

setwd('resultsDir')
# Define correct paths for these two annotation files
met_clusters <- read.csv('MetabolicClusters_MModel.csv') 
cond_names <- read.csv('ConditionNames_Li2014.csv')

conds = list.files(pattern = 'fvadata')
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
    rxn_clust$opt[j] <- sum(subset(rxn, rxn$condition == growth_conds[j])$optval)
    rxn_clust$min[j] <- sum(subset(rxn, rxn$condition == growth_conds[j])$min)
    rxn_clust$max[j] <- sum(subset(rxn, rxn$condition == growth_conds[j])$max)
  }

  if (i == 1) (glc_normalized <- rxn_clust) else (glc_normalized <- rbind(glc_normalized, rxn_clust))
}


glc_normalized$cluster <- factor(glc_normalized$cluster, levels = sort(levels(glc_normalized$cluster)))

glc_normalized <- subset(glc_normalized, (glc_normalized$cluster %in%
                                            c('Thr and Lys Metabolism', 'Val, Leu, and Ile Metabolism',
                                              'Pyruvate Metabolism', 'Transport, Inner Membrane', 'Folate Metabolism', 'Pentose Phosphate Pathway')))
glc_normalized$cluster <- gsub('Citric Acid Cycle', 'Tricarboxylic Acid Cycle', glc_normalized$cluster)

g <- ggplot(glc_normalized, aes(ymin = min, ymax = max, y=min, x = condition, fill = condition)) + 
  geom_crossbar(position = position_dodge(1), fatten = 0, size = 1) + 
  scale_fill_manual(values= c("#0084A9", "#007A60", "#32A669", "#33ADAC", "#0089CF", "#404040")) +
  scale_y_continuous(trans = 'pseudo_log', limits = c(-2e4, +2e4), breaks = c(-2000,-500,-100,-10,0,10,100,500,2000), expand = c(0,0)) +
  scale_x_discrete(breaks = NULL) + 
  facet_wrap(cluster ~ ., ncol = 6, labeller = label_wrap_gen(multi_line = TRUE)) +
  labs(y = 'Sum of fluxes through\n metabolic cluster', 
       x = 'Growth medium',
       fill = 'Growth medium:', 
       title = 'iML1515') +
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
        legend.position="none",
        legend.box = "horizontal",
        panel.background = element_rect(fill = NA),
        panel.border = element_rect(colour = "black", fill=NA, size=2),
        panel.grid.major = element_line(colour = "grey50"),
        plot.margin = margin(20, 20, 20, 20)) +
  guides(fill=guide_legend(nrow=1,byrow=FALSE))
plot(g)

