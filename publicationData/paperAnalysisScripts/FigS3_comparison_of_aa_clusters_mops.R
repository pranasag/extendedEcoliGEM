require(ggplot2)
require(reshape2)

resultsDir = '' #Path for the results directory
# Define correct paths for these two annotation files
met_clusters <- read.csv('MetabolicClusters_ExtendedModel.csv') 
cond_names <- read.csv('ConditionNames_Li2014.csv')

conds = list.files(paste0(resultsDir, './withoutProteomics_sameAsFig1B/'), pattern = 'fvadata')
setwd(resultsDir)
for(i in seq(1, length(conds))) {
  read.csv(paste0('./withProteomics/', fileNames[j]))
  fva_s$condition <- cond_names$name[grep(paste0("^",conds[i]), cond_names$csv)]
  if (i == 1) (compiled = fva_s) else (compiled = rbind(compiled, fva_s))}

for(i in seq(1, length(conds))) {
  fva_wo <- read.csv(paste0('./withoutProteomics_sameAsFig1B/', fileNames[j]))
  fva_wo$condition <- cond_names$name[grep(paste0("^",conds[i]), cond_names$csv)]
  if (i == 1) (compiled_wo = fva_wo) else (compiled_wo = rbind(compiled_wo, fva_wo))}

compiled <- merge.data.frame(compiled, met_clusters, by.x = 'Jid', by.y = 'Reaction.ID')
compiled$min <- ifelse(is.na(compiled$min), compiled$optval, compiled$min)
compiled$max <- ifelse(is.na(compiled$max), compiled$optval, compiled$max)

compiled_wo <- merge.data.frame(compiled_wo, met_clusters, by.x = 'Jid', by.y = 'Reaction.ID')
compiled_wo$min <- ifelse(is.na(compiled_wo$min), compiled_wo$optval, compiled_wo$min)
compiled_wo$max <- ifelse(is.na(compiled_wo$max), compiled_wo$optval, compiled_wo$max)

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

clusters <- unique(compiled_wo$Subsystem)
growth_conds <- unique(compiled_wo$condition)
for(i in seq(1, length(clusters))) {
  rxn <- subset(compiled_wo, compiled_wo$Subsystem == clusters[i])
  rxn_clust <- data.frame(condition=growth_conds, cluster=clusters[i]) 
  for(j in seq(1, length(growth_conds))) {
    rxn_clust$opt[j] <- sum(abs(subset(rxn, rxn$condition == growth_conds[j])$optval))
    rxn_clust$min[j] <- sum(abs(subset(rxn, rxn$condition == growth_conds[j])$min))
    rxn_clust$max[j] <- sum(abs(subset(rxn, rxn$condition == growth_conds[j])$max))
  }
  
  if (i == 1) (glc_normalized_wo <- rxn_clust) else (glc_normalized_wo <- rbind(glc_normalized_wo, rxn_clust))
}

glc_normalized <- unique(glc_normalized)
glc_normalized$Proteomics <- 'With'
glc_normalized_wo <- unique(glc_normalized_wo)
glc_normalized_wo$Proteomics <- 'Without'

glc_normalized <- rbind(glc_normalized, glc_normalized_wo)
glc_normalized <- subset(glc_normalized, (glc_normalized$cluster %in%
                                            c('Thr and Lys Metabolism', 'Arg and Pro Metabolism', 'Met Metabolism',
                                              'Ala and Asp Metabolism', 'Tyr, Trp, and Phe Metabolism', 'Val, Leu, and Ile Metabolism',
                                              'Gly and Ser Metabolism', 'Cys Metabolism', 'His Metabolism', 
                                              'Glu Metabolism')))
glc_normalized$condition <- factor(glc_normalized$condition, levels=c('MOPS -AA', 'MOPS -Met', 'MOPS +AA', ordered=TRUE))
glc_normalized <- glc_normalized[order(glc_normalized$cluster),]
glc_normalized$cluster <- factor(glc_normalized$cluster, levels = sort(levels(glc_normalized$cluster)))
glc_normalized$Proteomics <- factor(glc_normalized$Proteomics, levels=c('Without', 'With', ordered=TRUE))

my_breaks <- function(x) { if (max(x) < 15) seq(0, 10, 2) else (seq(0, 60, 10))}
my_limits <- function(x) { if (max(x) < 15) c(-0.05,10) else (c(-3, 60))}

g <- ggplot(glc_normalized, aes(ymin = min, ymax = max, y=min, x = condition, col = condition, fill = Proteomics)) + 
  geom_crossbar(position = position_dodge(1), fatten = 0, size = 1.5) + 
  scale_y_continuous(limits = my_limits, breaks = my_breaks, expand = c(0,0)) +
  scale_x_discrete(breaks = unique(glc_normalized$Condition)) + 
  scale_color_manual(values= c("#0084A9", "#007A60", "#32A669", "#33ADAC", "#0089CF", "#404040")) +
  facet_wrap(cluster ~ ., scales = 'free_y', ncol = 5) +
  scale_fill_manual(values = c('grey90', "grey60")) +
  labs(y = 'Sum of fluxes through \n metabolic cluster', 
       col = 'Growth medium:',
       fill = 'Proteomics input:') +
  theme(aspect.ratio = 1,
        plot.title = element_text(hjust = 0.5, size = 24), 
        #axis.title.x = element_text(size = 24),
        axis.title.y = element_text(size = 24),
        strip.text.x = element_text(size = 18),
        strip.background = element_blank(),
        axis.text.y  = element_text(size = 20, margin = margin(r=5), hjust = 1),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 20, margin = margin(t=5), hjust = 0.5),
        axis.ticks = element_blank(),
        legend.title = element_text(size = 24, hjust = 0.5, vjust = 0.5, face='bold'),
        legend.text = element_text(size = 24, hjust = 0.5, vjust = 0.5),
        legend.key = element_rect(fill = NA),
        legend.position="bottom",
        legend.box = "horizontal",
        panel.background = element_rect(fill = NA),
        panel.border = element_rect(colour = "black", fill=NA, size=2),
        panel.grid.major = element_line(colour = "grey50"),
        plot.margin = margin(20, 20, 20, 20)) +
  guides(col=guide_legend(nrow=1,byrow=TRUE), fill=guide_legend(nrow=1,byrow=TRUE))
plot(g)
