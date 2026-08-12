

## Create a dummy volcano plot for methodological figure
saliva_combined_results <- read.csv("S11_2Saliva_Correlation_Results_Top40Features.csv",header = TRUE,stringsAsFactors = FALSE,check.names = FALSE)

## Subset the df
saliva_combined_results22 <- saliva_combined_results[1:25,]


library(ggplot2)
library(ggrepel)

## Distinct colours for functional themes
theme_colors <- c(
  "Phosphate acquisition and regulation"                = "#E41A1C",
  "Ion transport and homeostasis"                       = "#377EB8",
  "Proteolysis and hydrolysis"                           = "#4DAF4A",
  "Sugar transport and exchange"                        = "#984EA3",
  "Stress resistance and methylation"                   = "#FF7F00",
  "DNA repair and replication"                          = "#A65628",
  "Cell envelope and surface-associated functions"      = "#F781BF",
  "Signal transduction and transcriptional regulation" = "#00A6A6",
  "Protein interaction domains"                         = "#FFD700",
  "Mobile genetic elements"                             = "#6A3D9A",
  "Oxidoreductase activity"                             = "#1B9E77",
  "ATPase and ATP-dependent functions"                  = "#D95F02",
  "Translation and tRNA metabolism"                     = "#7570B3",
  "Cell-envelope stress response"                       = "#66A61E",
  "Domains of unknown function"                         = "#E7298A",
  "Unannotated or insufficiently characterized"         = "#666666",
  "Other functional category"                           = "#8C510A"
)

pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S11_FunctionalAnalysis/S11_6Dummy_VolcanoPlot.pdf",width = 12,height = 10)
ggplot(saliva_combined_results2,aes(x = Correlation_Coefficient,y = -log10(P_value),color = Functional_Theme)) +
  geom_point(size = 5, alpha = 0.85) +
  geom_text_repel(
    data = saliva_combined_results2[saliva_combined_results2$P_value <= 0.05, ,drop = FALSE],
    aes(label = Functional_Theme),
    size = 4,
    max.overlaps = Inf,
    box.padding = 0.5,
    point.padding = 0.3,
    min.segment.length = 0,
    show.legend = FALSE
  ) +
  geom_hline(yintercept = -log10(0.05),linewidth = 0.8,linetype = "dashed",color = "black") +
  geom_vline(xintercept = 0,linewidth = 0.8,linetype = "solid",color = "black") +
  scale_color_manual(values = theme_colors,na.value = "grey50") +
  guides(color = guide_legend(title = "Functional theme",nrow = 3,byrow = TRUE,override.aes = list(size = 4, alpha = 1))) +
  labs(x = "Correlation coefficient",y = expression(-log[10]("P-value")),title = "Correlation of Top 40 Functional Features with sHACK Score") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5,face = "bold",size = 12),
    axis.title.x = element_text(color = "#117A65",size = 12),
    axis.title.y = element_text(color = "#117A65",size = 12),
    panel.border = element_rect(color = "black",fill = NA,linewidth = 1),
    legend.position = "top",
    legend.text = element_text(size = 8),
    legend.title = element_text(face = "bold")
  )
dev.off()


save.image("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S11_FunctionalAnalysis/S11_6Dummy_Figure_Methodology_Workspace.RData")