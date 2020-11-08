require(reshape2)
require(ggplot2)

resultsDir = '' #Path for the results directory
expressionDir = '' #Path for the expression data
# Define correct path for the annotation file
cond_names <- read.csv('ConditionNames_Li2014.csv')

conds = list.files(resultsDir, pattern = 'fvadata')
setwd(resultsDir)
for(i in seq(1, length(conds))) {
  fva <- read.csv(conds[i])
  fva$optval <- fva$optval/fva$optval[grep('R_BIOMASS_Ec_iML1515_core_75p37M', fva$Jid, fixed = 1)]
  fva_s <- subset(fva, (grepl("_s$", fva$Jid)))
  fva_s$condition <- cond_names$name[grep(paste0("^",conds[i]), cond_names$csv)]
  
  if (i == 1) (compiled = fva_s) else (compiled = rbind(compiled, fva_s))}

compiled$Jid <- gsub("R_","",compiled$Jid)
compiled$Jid <- gsub("_s","",compiled$Jid)
compiled <- subset(compiled, compiled$optval != 0.0)
compiled$condition <- factor(compiled$condition, levels=c('MOPS -AA', 'MOPS -Met', 'MOPS +AA', ordered=TRUE))

expr_data = list.files(expressionDir, pattern = '.csv')
setwd(expressionDir)
for(i in seq(1, length(expr_data))) {
  fva <- read.csv(expr_data[i])
  fva$condition <- cond_names$name[grep(paste0("^",expr_data[i],'.fvadata.csv'), cond_names$csv)]
  
  if (i == 1) (experimental = fva) else (experimental = rbind(experimental, fva))}

experimental <- subset(experimental, experimental$conc_mmolgDW != '#VALUE!')
experimental$conc_mmolgDW <- gsub("E", "e", experimental$conc_mmolgDW)
experimental$conc_mmolgDW <- as.numeric(experimental$conc_mmolgDW)

df = merge(compiled, experimental, by.x = c('Jid', 'condition'), by.y = c('Entry', 'condition'))
df$ratio <- log2(df$optval/df$conc_mmolgDW)
df$abs_ratio <- abs(log2(df$optval/df$conc_mmolgDW))

g <- ggplot(df, aes(x = optval, y = conc_mmolgDW, col = condition)) +
  geom_abline(slope = 1) +
  scale_x_continuous(trans='log10', limits = c(1e-10,1e-02), expand = c(0,0)) +
  scale_y_continuous(trans='log10', limits = c(1e-10,1e-02), expand = c(0,0)) +
  geom_point(size = 2) +
  scale_color_manual(values= c("#0084A9", "#007A60", "#32A669", "#33ADAC", "#0089CF", "#404040")) +
  facet_wrap(condition~.) +
labs(y = expression(atop("Experimentally determined", abundance~(mmol~gDW^-1))), 
     x = expression(Predicted~abundance~(mmol~gDW^-1)),
     col = 'Growth medium: ') +
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
        panel.spacing = unit(5, 'lines'),
        plot.margin = margin(20, 20, 20, 20)) +
  guides(fill=guide_legend(nrow=1,byrow=FALSE), col=guide_legend(nrow=1,byrow=FALSE), linetype = 'none')
plot(g)
