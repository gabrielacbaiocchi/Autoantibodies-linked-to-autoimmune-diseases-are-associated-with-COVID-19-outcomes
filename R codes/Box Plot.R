library(ggpubr)
library(lemon)
library(ggplot2)
library(openxlsx)
library(stringr)

IgA.all<- "input.data"


my_comparisons <- list( c("Control-Young", "Control-Elderly"), c("Mild-Young", "Mild-Elderly"), c ("Moderate-Young", "Moderate-Elderly"), c("Severe-Young", "Severe-Elderly"))#define comparisons between groups

symnum.args <- list(
  cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1.01), #1.01 Because in some rare cases, 1 is shown as NS for no reason whatsoever 
  symbols = c("****", "***", "**", "*", waiver()))

IgA.all$counts <- as.numeric(IgA.all$counts)

IgA.all$Group <- factor(IgA.all$Group, levels= c("Control-Young","Control-Elderly", "Mild-Young", "Mild-Elderly", "Moderate-Young", "Moderate-Elderly", "Severe-Young", "Severe-Elderly")) #define groups
IgA.all$gene <- factor(IgA.all$gene, levels=str_sort(unique(IgA.all$gene)))

png(file = "Fig01-IgA.png", bg = "transparent", width = 4500, height = 2500, units = "px", res = 300)
ggboxplot(IgA.all, x = "Group", y = "counts",
          color = "Group", palette = c("#bcbddc", "#9ecae1", "#807dba", "#6baed6", "#6a51a3", "#08519c", "#3f007d", "#08306b"),
          add = c("jitter", "median_iqr"), #"mean_sd"
          add.params = list(size = 0.5))+
  scale_y_continuous(breaks = pretty(c(0,2.5), n = 5), #Values showed from 0 to 8 with 4 breaks will show 0, 2, 4, 6 and 8
                     limits = c(0,4.5), name="IgA(Units/mL)")+ # y scale + axis' name
  facet_rep_wrap(~gene, strip.position="bottom", ncol=5)+ #if you have more than 2 lines, use ( facet_rep_wrap(~gene, strip.position="bottom") ) instead of facet_Wrap
  stat_compare_means(method= "wilcox.test", label = "p.format", 
                     label.y = c(2.6, 2.9, 3.2, 3.5, 3.8, 4.0), bracket.size = 0.4, tip.length = 0.006, #location and brackets customization
                     hide.ns = F, aes(Group = Group), comparisons = my_comparisons, 
                     vjust = 0.5, size = 5, # * customization: vjust == 0 -> in the middle of the bracket; accepts negative values
                     symnum.args = symnum.args)+
  theme(strip.text.x = element_text(size = 12))+
  theme(axis.line.x=element_line(), axis.ticks.x=element_blank(),axis.title.x=element_blank(), axis.text.x=element_blank(),
        strip.background = element_rect(color="white", fill="white"), strip.placement = "outside",
        panel.border = element_blank(), panel.background = element_blank(), panel.grid = element_blank(), 
        panel.spacing.x = unit(0,"line"), panel.spacing.y = unit(0,"line"))
dev.off()

