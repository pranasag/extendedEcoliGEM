require(reshape2)
require(ggplot2)

resultsDir = '' #Path for the results directory
# Define correct paths for the annotation file
cond_names <- read.csv('ConditionNames_Li2014.csv')

fileNames = list.files(paste0(resultsDir, './withoutProteomics_sameAsFig1B/'), pattern = 'fvadata')
setwd(resultsDir)

for(j in seq(1, length(fileNames))) {
wo_prot <- read.csv(paste0('./withoutProteomics_sameAsFig1B/', fileNames[j]))
w_prot <- read.csv(paste0('./withProteomics/', fileNames[j]))

wo_prot <- subset(wo_prot, !(grepl("_s$", wo_prot$Jid) | grepl("_f$", wo_prot$Jid) | grepl("_inh_d$", wo_prot$Jid) | grepl("_inhibition$", wo_prot$Jid)))
wo_prot$span[is.na(wo_prot$span)] <- 0

w_prot <- subset(w_prot, !(grepl("_s$", w_prot$Jid) | grepl("_f$", w_prot$Jid) | grepl("_inh_d$", w_prot$Jid) | grepl("_inhibition$", w_prot$Jid)))
w_prot$span[is.na(w_prot$span)] <- 0

wo_spans <- sort(wo_prot$span)
w_spans <- sort(w_prot$span)

start_val = 0
for(i in seq(1,length(wo_spans))) {
  start_val <- start_val + wo_spans[i]
  wo_spans[i] <- start_val
}

start_val = 0
for(i in seq(1,length(w_spans))) {
  start_val <- start_val + w_spans[i]
  w_spans[i] <- start_val
}
ks.test(wo_spans, w_spans)

df = data.frame(index=seq(1,length(wo_spans)), wo=wo_spans, w=w_spans)
df$condition<- cond_names$name[grep(paste0("^",fileNames[j]), cond_names$csv)]

if (j==1) (results <- df) else (results <- rbind(results, df))
}

df_melt <- melt(results, measure.vars = c('wo', 'w'))
df_melt$condition <- factor(df_melt$condition, levels=c('MOPS -AA', 'MOPS -Met', 'MOPS +AA', ordered=TRUE))
df_melt$variable <- gsub('wo', 'Without', df_melt$variable)
df_melt$variable <- gsub('w', 'With', df_melt$variable)
df_melt <- subset(df_melt, df_melt$index >5000)
df_melt <- df_melt[order(df_melt$variable),]
df_melt$value <- df_melt$value/(max(subset(df_melt$value, df_melt$condition == 'MOPS +AA' & df_melt$variable == 'With')))

g <- ggplot(df_melt, aes(x = index, y = value, col = condition)) +
  geom_line(aes(linetype = variable), size=1.5) +
  scale_y_continuous(limits = c(-0.05,1), breaks = seq(0, 1, 0.2), expand = c(0,0)) +
  scale_x_continuous(limits = c(5e3,6e3), breaks = c(5000,5200,5400,5600,5800,6000), expand = c(0,0)) +
  scale_color_manual(values= c("#0084A9", "#007A60", "#32A669", "#33ADAC", "#0089CF", "#404040")) +
labs(y = 'Cumulative flux variability', 
     x = 'Rank of the metabolic reaction',
     col = 'Growth medium: ', 
     linetype = 'Proteomics input: ') +
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
        legend.position=c(0.3,0.6),
        legend.key = element_rect(fill = NA),
        legend.box = "vertical",
        panel.background = element_rect(fill = NA),
        panel.border = element_rect(colour = "black", fill=NA, size=2),
        panel.grid.major = element_line(colour = "grey50"),
        plot.margin = margin(20, 20, 20, 20)) +
  guides(col=guide_legend(nrow=3,byrow=TRUE), linetype=guide_legend(nrow=2,byrow=TRUE))
plot(g)
