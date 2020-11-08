require(ggplot2)
require(reshape2)

resultsDir = '' #Path for the results directory
# Define correct paths for these two annotation files
met_clusters <- read.csv('MetabolicClusters_ExtendedModel.csv') 
cond_names <- read.csv('ConditionNames_Li2014.csv')

conds = list.files(resultsDir, pattern = 'fvadata')
setwd(resultsDir)
for(i in seq(1, length(conds))) {
  fva_s <- read.csv(conds[i])
  fva$glucose <- strsplit(strsplit(conds[i], '_')[[1]][3], '.fvadata.csv')[[1]][1]
  fva$growth <- fva$optval[grep("R_BIOMASS_Ec_iML1515_core_75p37M", fva$Jid, fixed = 1)]
  fva_s$condition <- cond_names$name[grep(paste0(strsplit(conds[i], '_')[[1]][1], '.fvadata.csv'), cond_names$csv)]
  
  if (i == 1) (compiled = fva_s) else (compiled = rbind(compiled, fva_s))}

compiled <- merge.data.frame(compiled, met_clusters, by.x = 'Jid', by.y = 'Reaction.ID')
compiled$min <- ifelse(is.na(compiled$min), compiled$optval, compiled$min)
compiled$max <- ifelse(is.na(compiled$max), compiled$optval, compiled$max)

compiled$glucose <- compiled$glucose/(-1)
compiled_backup <- compiled
compiled = subset(compiled, (compiled$Subsystem == 'Citric Acid Cycle') | (compiled$Subsystem == 'Pyruvate Metabolism'))

## For the growth rates
glcs <- sort(unique(compiled$glucose))
growth_conds <- unique(compiled$condition)
clusters <- unique(compiled$Subsystem)
for(k in seq(1, length(clusters))) {
  for(i in seq(1, length(glcs))) {
    rxn <- subset(compiled, compiled$glucose == glcs[i] & compiled$Subsystem == clusters[k])
    rxn_clust <- data.frame(condition=growth_conds, cluster = clusters[k], glucose=glcs[i])
      for(j in seq(1, length(growth_conds))) {
        rxn_clust$opt[j] <- sum(abs(subset(rxn, rxn$condition == growth_conds[j])$optval))
        rxn_clust$min[j] <- sum(abs(subset(rxn, rxn$condition == growth_conds[j])$min))
        rxn_clust$max[j] <- sum(abs(subset(rxn, rxn$condition == growth_conds[j])$max))
      }
      
      if ((k == 1)&(i == 1)) (glc_normalized <- rxn_clust) else (glc_normalized <- rbind(glc_normalized, rxn_clust))
    }
}

g <- ggplot(glc_normalized, aes(ymin = min, ymax = max, y=min, x = condition, fill = condition, col = cluster)) + 
  geom_crossbar(fatten = 0, size = 1) + 
  scale_y_continuous(breaks = c(0,50,100,150,200), limits = c(0,200), expand = c(0,0)) +
  scale_x_discrete(breaks = NULL) + 
  facet_wrap(glucose ~ ., ncol = 5, labeller = label_wrap_gen(), strip.position = 'bottom') +
  scale_color_manual(values = c('black', "grey50")) +
  scale_fill_manual(values= c("#0084A9", "#007A60", "#32A669", "#33ADAC", "#0089CF", "#404040")) +
  labs(y = 'Sum of fluxes through\n metabolic cluster', 
       fill = 'Growth medium:', 
       col = 'Metabolic cluster:',
       x = expression(Growth~rate~(h^{-1}))) +
  theme(aspect.ratio = 1,
        plot.title = element_text(hjust = 0.5, size = 24), 
        axis.title.x = element_text(size = 24),
        axis.title.y = element_text(size = 24),
        strip.text.x = element_text(size = 18),
        strip.background = element_blank(),
        axis.text.y  = element_text(size = 20, margin = margin(r=5), hjust = 1),
        axis.text.x = element_text(size = 20, margin = margin(t=5), hjust = 0.5),
        axis.ticks = element_blank(),
        legend.title = element_text(size = 24, hjust = 0, vjust = 0.5, face='bold'),
        legend.text = element_text(size = 24, hjust = 0, vjust = 0.5),
        legend.key = element_rect(fill = NA),
        legend.position="bottom",
        legend.box = "vertical",
        panel.background = element_rect(fill = NA),
        panel.border = element_rect(colour = "black", fill=NA, size=2),
        panel.grid.major = element_line(colour = "grey50"),
        plot.margin = margin(20, 20, 20, 20)) +
  guides(fill=guide_legend(nrow=1,byrow=FALSE), col=guide_legend(nrow=1,byrow=FALSE), linetype = 'none')
plot(g)
