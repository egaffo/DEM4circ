library(data.table)
library(ggplot2)

# #' A function to mimic the ggplot default palette
# gg_color_hue <- function(n) {
#   hues = seq(15, 375, length = n + 1)
#   hcl(h = hues, l = 65, c = 100)[1:n]
# }

method_types <-
  data.table(Method = c("circMeta_Dt_PZT", "circMeta_Lc_PZT", "DESeq2_BP_WaT",
                        "DESeq2_Dt_LRT", "DESeq2_Dt_WaT", "DESeq2_GP_LRT",
                        "DESeq2_Lc_LRT", "DESeq2_Sc_LRT", "DESeq2_Zi_LRT",
                        "DESeq2_ZW_LRT", "DEsingle_Dt_LRT", "edgeR_Dt_LRT",
                        "edgeR_rbst_LRT", "edgeR_rbst_QFT", "edgeR_rbst50df_LRT",
                        "edgeR_rbstEdf_LRT", "edgeR_rbstTMMswp_LRT", "edgeR_ZW_LRT",
                        "limma_St_Vst", "lncDIFF_Dt_LRT", "MAST_CPM_LRT",
                        "metagenomeSeq_Dt_EBT", "monocle_Dt_LRT", "NBID_Dt_LRT",
                        "NBID_Sc_LRT", "NOISeqBIO_TMM_LRT", "PoissonSeq_Dt_TT",
                        "ROTScpm_Dt_TLT", "SAMseq_Dt_TT",
                        # "scCODE_Dt_MT", "scDEA_Dt_MT",
                        # "hytest_TMM_HG",
                        "seurat_Bim_LRT",
                        "seurat_Dt_WLX", "voom_Dt_MFT", "voom_LF_MFT",
                        "voom_Qn_MFT", "voom_Rb_MFT", "voom_Sp_MFT",
                        "voom_ZW_MFT", "wilcoxon_TMM_WLX"),
             Type = c("Bulk RNA-seq", "Bulk RNA-seq", "Bulk RNA-seq",
                      "Bulk RNA-seq", "Bulk RNA-seq", "Single cell / metagenome",
                      "Single cell / metagenome", "Single cell / metagenome", "Single cell / metagenome",
                      "Single cell / metagenome", "Single cell / metagenome", "Bulk RNA-seq",
                      "Bulk RNA-seq", "Bulk RNA-seq", "Bulk RNA-seq",
                      "Bulk RNA-seq", "Single cell / metagenome", "Single cell / metagenome",
                      "Single cell / metagenome", "Bulk RNA-seq", "Single cell / metagenome",
                      "Single cell / metagenome", "Single cell / metagenome", "Single cell / metagenome",
                      "Single cell / metagenome", "Bulk RNA-seq", "Bulk RNA-seq",
                      "Single cell / metagenome", "Bulk RNA-seq",
                      # "Single cell / metagenome", "Single cell / metagenome",
                      # "Bulk RNA-seq",
                      "Single cell / metagenome",
                      "Single cell / metagenome", "Bulk RNA-seq", "Single cell / metagenome",
                      "Bulk RNA-seq", "Bulk RNA-seq", "Bulk RNA-seq",
                      "Single cell / metagenome", "Single cell / metagenome"),
             BaseMethod = c("circMeta", "circMeta", "DESeq2",
                            "DESeq2", "DESeq2", "glmGamPoi",
                            "DESeq2", "DESeq2", "DESeq2",
                            "DESeq2", "DEsingle", "edgeR",
                            "edgeR", "edgeR", "edgeR",
                            "edgeR", "edgeR", "edgeR",
                            "limma", "lncDIFF", "MAST",
                            "metagenomeSeq", "monocle", "NBID",
                            "NBID", "NOISeqBIO", "PoissonSeq",
                            "ROTS", "SAMseq",
                            # "scCODE", "scDEA",
                            # "Hy-test",
                            "Seurat",
                            "Seurat", "Voom", "Voom",
                            "Voom", "Voom", "Voom",
                            "Voom", "Wilcoxon"))

method_types <- method_types[order(toupper(BaseMethod),
                                   toupper(Method))]
method_order_alf <- method_order <- method_types[order(toupper(BaseMethod),
                                                       toupper(Method)),
                                                 Method]

# set colors
source("../utils/color_palettes.R")

baseMethodColor <-  setNames(object = c(custom_old_gdocs),
                             nm = unique(method_types$BaseMethod))

multishade_methods <- merge(method_types[, .N, by = BaseMethod], #[N > 1, ],
                            as.data.table(baseMethodColor, keep.rownames = "BaseMethod"),
                            by = "BaseMethod")

base_color_shades <- mapply(get_color_shades,
                            setNames(multishade_methods$N - 1,
                                     multishade_methods$BaseMethod),
                            multishade_methods$baseMethodColor,
                            USE.NAMES = T)

method_color_shades <-
  rbindlist(sapply(base_color_shades,
                   function(x)data.frame(Color = as.character(x)),
                   simplify = F),
            idcol = "BaseMethod")[order(BaseMethod, Color)]

method_types[order(BaseMethod), Color := method_color_shades$Color]

method_colors <- setNames(object = method_types$Color,
                          nm = method_types$Method)

plot_dt <-
  melt(copy(method_types)[, c("Base", "Parameters",
                              "Test") := tstrsplit(Method, "_")][, .(Method,
                                                                     Type,
                                                                     `Base tool` = BaseMethod,
                                                                     Parameters,
                                                                     Test)],
       id.vars = c("Method", "Type"))

## Fix names
plot_dt[value == "metagenomeSeq", value := "metagenSeq"]
plot_dt[Method == "wilcoxon_TMM_WLX" & variable == "Parameters", value := "RC"]
plot_dt[Method == "edgeR_rbstTMMswp_LRT" & variable == "Parameters", value := "Twsp"]
plot_dt[Method == "edgeR_rbst50df_LRT" & variable == "Parameters", value := "50DF"]
plot_dt[Method == "edgeR_rbstEdf_LRT" & variable == "Parameters", value := "EDF"]
plot_dt[Method == "DESeq2_GP_LRT" & variable == "Parameters", value := "Dt"]
plot_dt[Method == "voom_Rb_MFT" & variable == "Parameters", value := "RBST"]
plot_dt[Method == "limma_St_Vst" & variable == "Parameters", value := "VST"]
plot_dt[Method == "limma_St_Vst" & variable == "Test", value := "MFT"]

plot_dt[Method == "DESeq2_Dt_LRT" & variable == "Parameters", value := "LRT"]
plot_dt[Method == "DESeq2_Dt_WaT" & variable == "Parameters", value := "WAT"]
plot_dt[Method == "edgeR_rbst_QFT" & variable == "Parameters", value := "QFT"]
plot_dt[Method == "seurat_Dt_WLX" & variable == "Parameters", value := "WLX"]

plot_dt[variable == "Parameters", value := toupper(value)]

plot_dt$variable <- factor(plot_dt$variable,
                           levels = rev(sort((levels(plot_dt$variable)))),
                           ordered = T)

plot_base_method_color_text <-
  ggplot(plot_dt[variable != "Test"],
         aes(x = Method,
             y = variable,
             color = Method)) +
  geom_text(aes(label = value),
            angle = 90,
            hjust = 0.3,
            size = 3) +
  labs(x = NULL, y = NULL) +
  coord_cartesian(clip = "off") +
  scale_color_manual(values = method_colors, guide = "none") +
  facet_grid(rows = vars(variable),
             scales = "free",
             space = "free",
             drop = T) +
  scale_x_discrete(limits = method_order_alf) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(angle = 90, hjust = .5, vjust = .5),
        text = element_text(size = 10),
        panel.border = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_rect(fill = "transparent"),
        strip.background = element_rect(fill = NA),
        strip.text = element_blank(),
        plot.background = element_blank(),
        panel.spacing.y = unit(0.05, "lines"),
        panel.spacing.x = unit(.1, "lines"))

# plot_base_method_color_text

plot_base_method_color_text_v <-
  ggplot(plot_dt[variable != "Test"],
         aes(y = Method,
             x = variable,
             color = Method)) +
  geom_text(aes(label = value),
            size = 3) +
  scale_x_discrete(position = "bottom") +
  labs(x = NULL, y = NULL) +
  coord_cartesian(clip = "off") +
  scale_color_manual(values = method_colors, guide = "none") +
  theme_bw() +
  theme(axis.text.y = element_blank(),
        text = element_text(size = 10),
        panel.border = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        strip.background = element_rect(fill = NA),
        panel.spacing.y = unit(0.05, "lines"),
        panel.spacing.x = unit(.1, "lines"))

# plot_base_method_color_text_v

method_labels_table <-
  dcast(plot_dt[, .(Method, variable, value)],
        Method ~ variable)[, .(Method,
                               MethodLabel = paste(`Base tool`,
                                                   Parameters, #Test,
                                                   sep = "\t"))]

library(ggthemes)
plot_dt <- copy(method_types)

plot_base_method_type_tiles_v <-
  ggplot(plot_dt,
         aes(y = Method,
             x = "", fill = Type)) +
  geom_tile(color = "white", width = 0.9, height = 0.9,
            size = .5) +
  labs(x = NULL, y = NULL) +
  scale_fill_manual(values = setNames(tableau_color_pal(palette = "Classic Purple-Gray 6")(2), #"Nuriel Stone"
                                      c("Bulk RNA-seq", "Single cell / metagenome"))) +
  guides(fill = guide_legend(title = "Method purpose",
                             title.position = "left",
                             title.hjust = .5)) +
  theme_dendro() +
  theme(legend.position = "bottom")

# plot_base_method_type_tiles_v
