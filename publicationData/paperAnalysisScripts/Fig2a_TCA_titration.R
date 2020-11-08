require(ggplot2)
require(reshape2)

resultsDir = '' #Path for the results directory
# Define correct paths for these two annotation files
met_clusters <- read.csv('MetabolicClusters_ExtendedModel.csv') 
cond_names <- read.csv('ConditionNames_Li2014.csv')

experiments = read.csv('acetateExcretionExperimental.csv')

conds = list.files(resultsDir, pattern = 'fvadata')
setwd(resultsDir)
for(i in seq(1, length(conds))) {
  fva <- read.csv(conds[i])
  fva_s <- subset(fva, fva$Jid == 'R_EX_ac_e')
  fva_s$glucose <- fva$optval[grep("R_EX_glc__D_e", fva$Jid, fixed = 1)]
  fva_s$mu <- fva$optval[grep("R_BIOMASS_Ec_iML1515_core_75p37M", fva$Jid, fixed = 1)]
  fva_s$condition <- cond_names$name[grep(paste0(strsplit(conds[i], '_')[[1]][1], '.fvadata.csv'), cond_names$csv)]
  
  if (i == 1) (compiled = fva_s) else (compiled = rbind(compiled, fva_s))}

compiled$min <- ifelse(is.na(compiled$min), compiled$optval, compiled$min)
compiled$max <- ifelse(is.na(compiled$max), compiled$optval, compiled$max)

compiled$glucose <- compiled$glucose/(-1)
compiled$mu_2 <- compiled$mu * 50

g <- ggplot(compiled, aes(x = mu, y = optval)) + 
  geom_line(aes(col = condition), size=1.75) + 
  geom_point(data = experiments, aes(x=mu,y=Acetate, shape=Dataset), size = 7, col = 'black') +
  geom_smooth(data = subset(experiments, experiments$mu>0.38), aes(x=mu,y=Acetate, linetype=Dataset), 
              method = lm, size = 1, col = 'black', se=FALSE) +
  scale_x_continuous(limits = c(0,0.6), breaks=seq(0,0.60,0.10), expand = c(0,0)) +
  scale_y_continuous(limits = c(-1,20), breaks = seq(0,20,5), expand = c(0,0)) +
  scale_color_manual(values= c("#0084A9", "#007A60", "#32A669", "#33ADAC", "#0089CF", "#404040")) +
  labs(y = expression(Acetate~exchange~flux~(mmol~gDW^{-1}~h^{-1})),
       x = expression(Growth~rate~(h^{-1})),
       col = 'Growth medium:', 
       shape = 'Experimental dataset:') +
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
        legend.position=c(0.3,0.6),
        legend.box = "vertical",
        panel.background = element_rect(fill = NA),
        panel.border = element_rect(colour = "black", fill=NA, size=2),
        panel.grid.major = element_line(colour = "grey50"),
        plot.margin = margin(20, 20, 20, 20)) +
  guides(col=guide_legend(nrow=3,byrow=FALSE), shape=guide_legend(nrow=3,byrow=FALSE), linetype = 'none')
plot(g)
