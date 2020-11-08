require(ggplot2)
require(reshape2)

resultsDir = '' #Path for the results directory
# Define correct paths for the annotation file 
cond_names <- read.csv('ConditionNames_Li2014.csv')

conds = list.files(resultsDir, pattern = 'fvadata')
setwd(resultsDir)
for(i in seq(1, length(conds))) {
  fva_s <- read.csv(conds[i])
  fva_s <- subset(fva_s, (grepl("_s$", fva_s$Jid)))
  fva_s$condition <- conds[i]
  
  fva_s$condition<- cond_names$name[grep(paste0("^",conds[i]), cond_names$csv)]
  if (i == 1) (compiled = fva_s) else (compiled = rbind(compiled, fva_s))}

for(i in seq(1,nrow(compiled))) {
   if (is.nan(compiled$span[i])) (compiled$span[i] <- 0)
 }

g <- ggplot(compiled, aes(condition, span)) + 
  geom_boxplot(aes(col=condition), outlier.size = 5) + 
  scale_y_continuous(limits = c(-1e-5,5e-4), breaks = seq(0,5e-4,1e-4), labels= seq(0,5,1), expand = c(0,0)) +
  scale_color_manual(values= c("#0084A9", "#007A60", "#32A669", "#33ADAC", "#0089CF", "#404040")) +
  labs(y = expression(atop('Absolute flux interval of the protein', 'synthesis reaction ('~'\U00D7'~10^-4~')')), 
       x = 'Growth medium') +
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
        legend.position= 'none',
        legend.box = "vertical",
        panel.background = element_rect(fill = NA),
        panel.border = element_rect(colour = "black", fill=NA, size=2),
        panel.grid.major = element_line(colour = "grey50"),
        plot.margin = margin(20, 20, 20, 20)) +
  guides(fill=guide_legend(nrow=3,byrow=FALSE), col=guide_legend(nrow=1,byrow=FALSE), linetype = 'none')
plot(g)
