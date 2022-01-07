


# library -----------------------------------------------------------------

#GeoMX OMD
library(dplyr); library(reshape2); library(tidyr); library(magrittr)
library(ComplexHeatmap); library(RColorBrewer);library(viridis)
library(ggplot2); library(gridExtra)
library(DESeq2); 


#1. load data ---------------------------------------------------------------


# load GeoMX  data
data <- read.csv(file="./DSP_count.csv")

colnames(data)[1] <- "Gene"
colnames(data) <- gsub(".Geometric.Segment", "", colnames(data))
colnames(data) <- gsub(".CD10.", "", colnames(data))

clin.info <- read.csv(file="./ROI_info.csv")
clin.info$Sample_ID<- gsub("Nail 1", "Nail.1", clin.info$Sample_ID)
clin.info$Sample_ID<- gsub("Nail 2", "Nail.2", clin.info$Sample_ID)


library(tibble)
data1 <- data %>% tibble::column_to_rownames("Gene")
data1<-data1[ ,clin.info$Sample_ID]

#color
col_loc <-  c(RColorBrewer::brewer.pal(n=8, "Set2")[2:7],RColorBrewer::brewer.pal(n=10, "Set3")[2:10],RColorBrewer::brewer.pal(n=12, "Set3") ) 
names(col_loc) <- unique(clin.info$Type)
col_loc <- col_loc[1:23]


# 2 Data processing -------------------------------------------------------
#2-1. Pre-processing ----------------------------------------------------------

## 2.1.1  Raw Q3 normalized data 

data2 <-reshape2::melt(data1) 
data2 %<>%  dplyr::inner_join(clin.info, by=c("variable" = "Sample_ID"))
ggplot(data2, aes(x = variable, y=log2(value), colour=Type)) +
  geom_violin()+
  geom_boxplot(width=0.2)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_manual(values = col_loc) +
  facet_wrap(~ Scan_ID+Type, scale="free_x", nrow=1)+
  labs(title ="Rawdata_prefilter")




## 2.1.2 filtering GeneID

library(org.Hs.eg.db)
idfound <- rownames(data1) %in% mappedRkeys(org.Hs.egSYMBOL)
data1 <- data1[idfound,]

data2 <- data1
rownames(data2) <- toTable(org.Hs.egSYMBOL)$gene_id[match(rownames(data1), toTable(org.Hs.egSYMBOL)$symbol)]


# 2-2. TMM normalized ---------------------------------------------------

library(edgeR)
y <- DGEList(counts=data1,
             gene=rownames(data1),
             group=clin.info$Type)
y$genes$ids <- rownames(data2)

keep = filterByExpr(y, group=clin.info$Group,
                    min.count = 5, min.total.count = 15, 
                    large.n = 5, min.prop = 0.7) 


y <- y[keep, ,keep.lib.size=FALSE]
y <- calcNormFactors(y, method = "TMM")



design <- model.matrix(~ 0 + clin.info$Type, data = y$samples) 
y <- estimateDisp(y)


cpm.data <- edgeR::cpm(y, log=FALSE)

cpm.data.symbol <- cpm.data %>% as.data.frame()
cpm.data.symbol$ID <- y$genes$ids
cpm.data.symbol <- cpm.data.symbol[, c("ID", colnames(cpm.data))]




#2-3. CPM heatmap -------------------------------------------------------------

log2.cpm.data <- edgeR::cpm(y, prior.count=2, log=TRUE)

Nail_genes<-c("LGR5", "LGR6","KRT16", "KRT84", 
            "KRT17", "KRT1","KRT10","KRT5", "KRT85", "KRT6A",
            "KRT14","SPINK6", "DSG1", "COL17A1",
            "SBSN", "DSG1", "SPRR1B", "FLG",  "SPINK5","RSPO4", "RSPO3", 
            "CXCL14", "HEY1", "COL1A1", "MSX1", "BMP5",
            "TWIST1", "CRABP1", "MME", "VCAN")

log2.cpm.data <- log2.cpm.data[ , clin.info$Sample_ID ]
log2.cpm.data <- log2.cpm.data[ Nail_genes , ]


library("gplots")

Heatmap(log2.cpm.data, 
        show_row_names = TRUE, cluster_columns = TRUE, cluster_rows = F,
        column_names_gp = gpar(fontsize = 10),
        row_names_gp = gpar(fontsize = 7), 
        column_names_rot = 60, column_split =clin.info$Type,
        top_annotation = HeatmapAnnotation("Location" = clin.info$Type,
                                           col=list(Location = col_loc),
                                           simple_anno_size = unit(0.35, "cm"),
                                           gap = unit(2, "mm")) )




# 3 NichNet ---------------------------------------------------------------

# 3-1 interaction analyses --------------------------------------------------
# 3-2 ligand-receptor database ------------------------------------------------

ligand_target_matrix = readRDS("/Users/manki_lab/Work/Project/_B_onychomatricoma/Polydactyly_SingleCell/ligand_target_matrix.rds")
lr_network = readRDS("/Users/manki_lab/Work/Project/_B_onychomatricoma/Polydactyly_SingleCell/lr_network.rds")

# interactions and their weights in the ligand-receptor + signaling network
weighted_networks = readRDS("/Users/manki_lab/Work/Project/_B_onychomatricoma/Polydactyly_SingleCell/weighted_networks.rds")
weighted_networks_lr = weighted_networks$lr_sig %>% inner_join(lr_network %>% distinct(from,to), by = c("from","to"))



# 3-3 load data DP---------------------------------------------------------------

# load GeoMX  data
data <- read.csv(file="./DSP_count.csv")

colnames(data)[1] <- "Gene"
colnames(data) <- gsub(".Geometric.Segment", "", colnames(data))
colnames(data) <- gsub(".CD10.", "", colnames(data))

clin.info <- read.csv(file="./ROI_info.csv")
clin.info$Sample_ID<- gsub("Nail 1", "Nail.1", clin.info$Sample_ID)
clin.info$Sample_ID<- gsub("Nail 2", "Nail.2", clin.info$Sample_ID)


library(tibble)
data1 <- data %>% tibble::column_to_rownames("Gene")
data1<-data1[ ,clin.info$Sample_ID]


#DP, Hair matrix, Controls
clin.info<- subset(clin.info, Type %in% c("Hair matrix", "Follicular dermal papilla", "Control dermis", "Control epidermis")) 

library(tibble)
data1 <- data %>% tibble::column_to_rownames("Gene")
data1<-data1[ ,clin.info$Sample_ID]

col_loc <-  c(RColorBrewer::brewer.pal(n=12, "Set3"), RColorBrewer::brewer.pal(n=8, "Set2")[1:6],RColorBrewer::brewer.pal(n=10, "Set3")[2:10] ) 
names(col_loc) <- unique(clin.info$Type)
col_loc <- col_loc[1:4]


##"NegProbe-WTX"

library(org.Hs.eg.db)
idfound <- rownames(data1) %in% c(mappedRkeys(org.Hs.egSYMBOL),"NegProbe-WTX")
data1 <- data1[idfound,]

data2 <- data1

#3-4 TMM normalized ---------------------------------------------------

library(edgeR)
y <- DGEList(counts=data1,
             gene=rownames(data1),
             group=clin.info$Type)
y$genes$ids <- rownames(data2)

keep = filterByExpr(y, group=clin.info$Group,
                    min.count = 5, min.total.count = 15, 
                    large.n = 5, min.prop = 0.7) 


y <- y[keep, ,keep.lib.size=FALSE]
y <- calcNormFactors(y, method = "TMM")



design <- model.matrix(~ 0 + clin.info$Type, data = y$samples) 
y <- estimateDisp(y)


cpm.data <- edgeR::cpm(y, log=FALSE)

cpm.data.symbol <- cpm.data %>% as.data.frame()
cpm.data.symbol$ID <- y$genes$ids
cpm.data.symbol <- cpm.data.symbol[, c("ID", colnames(cpm.data))]

log2.cpm.data <- edgeR::cpm(y, prior.count=2, log=TRUE)

#3-5 DEG DP ---------------------------------------------------------

design <- model.matrix(~ 0+ clin.info$Type)
#design <- model.matrix(~0 +clin.info$Group)


fit <- glmQLFit(y,design)
qlf.Contvs1000 <- glmQLFTest(fit,contrast=c(-1, 0,1, 0)) 
#qlf.Contvs1000 <- glmQLFTest(fit,contrast=c(-1,1))

summary(decideTests(qlf.Contvs1000 ))
plotMD(qlf.Contvs1000 )

deg <- topTags(qlf.Contvs1000 )
DEGs <- qlf.Contvs1000$table %>% arrange(PValue)

DEGs.gene <- DEGs$logFC
names(DEGs.gene) <- rownames(DEGs)

DEG.data <- qlf.Contvs1000$table %>% dplyr::filter(PValue <0.05,  abs(logFC) > 1) %>% dplyr::arrange(-logFC) 
DEG.data$fold_change <- 2^DEG.data$logFC
DEG.data$CPM <- 2^DEG.data$logCPM  
DEG.data %<>% dplyr::select(logFC, fold_change, logCPM, CPM, F, PValue)

sig.gene<- DEGs %>% mutate(Gene = rownames(.)) %>% arrange(desc(logFC))  %>% .[ .$PValue < 0.01,]

RP.genes<- grep("^RPL[0-9]|^RPS[0-9]|^RPLP[0-9]|^HNRNP.*", rownames(DEGs), value = T)
RP.genes<-c(RP.genes, "RPSA" )
His.genes<- grep("^H[0-9]", rownames(DEGs), value = T)
KRT.genes<- grep("^KRT", rownames(DEGs), value = T)

sig.gene<- sig.gene[!rownames(sig.gene) %in% c(RP.genes, His.genes, KRT.genes),]
DP.expressed <- sig.gene[sig.gene$logFC >log2(3),]

DP.raw <- sig.gene
DP.all <- DEGs %>% mutate(Gene = rownames(.)) %>% arrange(desc(logFC))


#3-6 Mesenchymal Volcano plot -----------------------------------------------------------------

library(EnhancedVolcano)

EnhancedVolcano(sig.gene,
                lab = rownames(sig.gene),
                x = 'logFC',
                y = 'PValue',
                xlim = c(-5, 5),
                ylim= c(0,19),
                title = 'DP versus Dermis',
                pCutoff = 0.001,
                FCcutoff = 1.5,
                pointSize = 1,
                labSize = 0.2
)



#3-7 DEG HM ---------------------------------------------------------

design <- model.matrix(~ 0+ clin.info$Type)
#design <- model.matrix(~0 +clin.info$Group)

fit <- glmQLFit(y,design)
qlf.Contvs1000 <- glmQLFTest(fit,contrast=c(0,-1,0, 1)) 


summary(decideTests(qlf.Contvs1000 ))
plotMD(qlf.Contvs1000 )

deg <- topTags(qlf.Contvs1000 )
DEGs <- qlf.Contvs1000$table %>% arrange(PValue)

DEGs.gene <- DEGs$logFC
names(DEGs.gene) <- rownames(DEGs)

DEG.data <- qlf.Contvs1000$table %>% dplyr::filter(PValue <0.05,  abs(logFC) > 1) %>% dplyr::arrange(-logFC) 
DEG.data$fold_change <- 2^DEG.data$logFC
DEG.data$CPM <- 2^DEG.data$logCPM  
DEG.data %<>% dplyr::select(logFC, fold_change, logCPM, CPM, F, PValue)

sig.gene<- DEGs %>% mutate(Gene = rownames(.)) %>% arrange(desc(logFC))  %>% .[ .$PValue < 0.01,]

RP.genes<- grep("^RPL[0-9]|^RPS[0-9]|^RPLP[0-9]|^HNRNP.*", rownames(DEGs), value = T)
RP.genes<-c(RP.genes, "RPSA" )
His.genes<- grep("^H[0-9]", rownames(DEGs), value = T)

sig.gene<- sig.gene[!rownames(sig.gene) %in% c(RP.genes, His.genes),]
HM.expressed <- sig.gene[sig.gene$logFC >log2(3),]

HM.raw <- sig.gene
HM.all <- DEGs %>% mutate(Gene = rownames(.)) %>% arrange(desc(logFC))

#3-8 Epithelial Volcano -----------------------------------------------------------------

library(EnhancedVolcano)

#OMD.markers<- FindMarkers(ON.integrated, ident.1 =  "RSPO4 OMD", ident.2 = "fibroblast", min.pct = 0.25, logfc.threshold = log(1.3))
##9,7 landscape 

EnhancedVolcano(sig.gene,
                lab = rownames(sig.gene),
                x = 'logFC',
                y = 'PValue',
                xlim = c(-8.5, 8.5),
                ylim= c(0,20),
                title = 'OT epithelium versus Epidermis',
                pCutoff = 0.001,
                FCcutoff = 1.5,
                pointSize = 1,
                labSize = 0.1
)


# 3-9 ligand-receptor database ------------------------------------------------

#Database from NichNet
ligand_target_matrix = readRDS("./ligand_target_matrix.rds")
lr_network = readRDS("./lr_network.rds")

# interactions and their weights in the ligand-receptor + signaling network
weighted_networks = readRDS("./weighted_networks.rds")
weighted_networks_lr = weighted_networks$lr_sig %>% inner_join(lr_network %>% distinct(from,to), by = c("from","to"))



#3-10 HM-DP Reciever ----------------------------------------------------------------

DE_receiver <-  HM.expressed

colnames(cpm.data) == clin.info$Sample_ID
expression <- t(cpm.data)

#use Negative Probe value instead
expressed_genes_receiver <- expression[clin.info$Type=="Hair matrix",] %>% apply(2,function(x){10*(2**x - 1)}) %>%
  apply(2,function(x){log2(mean(x) + 1)}) %>% 
  .[. >= 5*.["NegProbe-WTX"]] %>% names()

background_expressed_genes = expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]

#3-11 HM-DP Sender ------------------------------------------------------------------

#Define a gene set of interest
DE_sender <-  DP.expressed
geneset_oi <- DP.expressed$Gene
geneset_oi = geneset_oi %>% .[. %in% lr_network$from]

expressed_genes_sender <- expression[clin.info$Type=="Follicular dermal papilla",] %>% apply(2,function(x){10*(2**x - 1)}) %>%
  apply(2,function(x){log2(mean(x) + 1)}) %>% 
  .[. >= 5*.["NegProbe-WTX"]] %>% names()

#if u want to check all set of interaction 
#geneset_oi <- colnames(ligand_target_matrix)

head(geneset_oi)

#Define a set of potential ligands
ligands = lr_network %>% pull(from) %>% unique()
expressed_ligands = intersect(ligands, expressed_genes_sender)

receptors = lr_network %>% pull(to) %>% unique()
expressed_receptors = intersect(receptors,expressed_genes_receiver)

lr_network_expressed = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) 

potential_ligands = lr_network_expressed %>% pull(from) %>% unique()

head(potential_ligands)

#Perform NicheNet ligand activity analysis
ligand_activities = predict_ligand_activities(geneset = geneset_oi, background_expressed_genes = background_expressed_genes, 
                                              ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)


ligand_activities %>% arrange(-pearson)

#best_upstream_ligands = ligand_activities %>% top_n(100, pearson) %>% arrange(-pearson) %>% pull(test_ligand)

best_upstream_ligands = ligand_activities %>% arrange(-pearson) %>% pull(test_ligand)
head(best_upstream_ligands)


#Infer receptors and top-predicted target genes

p_hist_lig_activity = ggplot(ligand_activities, aes(x=pearson)) + 
  geom_histogram(color="black", fill="darkorange")  + 
  # geom_density(alpha=.1, fill="orange") +
  geom_vline(aes(xintercept=min(ligand_activities %>% top_n(20, pearson) %>% pull(pearson))), color="red", linetype="dashed", size=1) + 
  labs(x="ligand activity (PCC)", y = "# ligands") +
  theme_classic()
p_hist_lig_activity


active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 250) %>% bind_rows()

nrow(active_ligand_target_links_df)
head(active_ligand_target_links_df)

active_ligand_target_links_df<- active_ligand_target_links_df[!is.na(active_ligand_target_links_df$weight),] 

active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, ligand_target_matrix = ligand_target_matrix, cutoff = 0.25)


order_ligands = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev()
order_targets = active_ligand_target_links_df$target %>% unique()
vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% t()

p_ligand_target_network = vis_ligand_target %>% 
  make_heatmap_ggplot("Prioritized ligands","Predicted target genes", color = "red",
                      legend_position = "top", x_axis_position = "top",legend_title = "Regulatory potential") + 
  scale_fill_gradient2(low = "whitesmoke",  high = "red", breaks = c(0,0.005,0.01)) + 
  theme(axis.text.x = element_text(face = "italic"))
p_ligand_target_network


#Receptors of top-ranked ligands

lr_network_top = lr_network %>% filter(from %in% best_upstream_ligands & to %in% expressed_receptors) %>% distinct(from,to)
best_upstream_receptors = lr_network_top %>% pull(to) %>% unique()

lr_network_top_df_large = weighted_networks_lr %>% filter(from %in% best_upstream_ligands & to %in% best_upstream_receptors)

lr_network_top_df = lr_network_top_df_large %>% spread("from","weight",fill = 0)
lr_network_top_matrix = lr_network_top_df %>% dplyr::select(-to) %>% as.matrix() %>% magrittr::set_rownames(lr_network_top_df$to)

dist_receptors = dist(lr_network_top_matrix, method = "binary")
hclust_receptors = hclust(dist_receptors, method = "ward.D2")
order_receptors = hclust_receptors$labels[hclust_receptors$order]

dist_ligands = dist(lr_network_top_matrix %>% t(), method = "binary")
hclust_ligands = hclust(dist_ligands, method = "ward.D2")
order_ligands_receptor = hclust_ligands$labels[hclust_ligands$order]

order_receptors = order_receptors %>% intersect(rownames(lr_network_top_matrix))
order_ligands_receptor = order_ligands_receptor %>% intersect(colnames(lr_network_top_matrix))

vis_ligand_receptor_network = lr_network_top_matrix[order_receptors, order_ligands_receptor]
rownames(vis_ligand_receptor_network) = order_receptors %>% make.names()
colnames(vis_ligand_receptor_network) = order_ligands_receptor %>% make.names()

p_ligand_receptor_network = vis_ligand_receptor_network %>% t() %>% make_heatmap_ggplot("Ligands","Receptors", color = "mediumvioletred", x_axis_position = "top",legend_title = "Prior interaction potential")
p_ligand_receptor_network


#Receptors of top-ranked ligands, but after considering only bona fide ligand-receptor interactions documented in literature and publicly available databases

lr_network_strict = lr_network %>% filter(database != "ppi_prediction_go" & database != "ppi_prediction")
ligands_bona_fide = lr_network_strict %>% pull(from) %>% unique()
receptors_bona_fide = lr_network_strict %>% pull(to) %>% unique()

lr_network_top_df_large_strict = lr_network_top_df_large %>% distinct(from,to) %>% inner_join(lr_network_strict, by = c("from","to")) %>% distinct(from,to)
lr_network_top_df_large_strict = lr_network_top_df_large_strict %>% inner_join(lr_network_top_df_large, by = c("from","to"))


DP_network <- lr_network_top_df_large_strict

lr_network_top_df_strict = lr_network_top_df_large_strict %>% spread("from","weight",fill = 0)
lr_network_top_matrix_strict = lr_network_top_df_strict %>% dplyr::select(-to) %>% as.matrix() %>% magrittr::set_rownames(lr_network_top_df_strict$to)

dist_receptors = dist(lr_network_top_matrix_strict, method = "binary")
hclust_receptors = hclust(dist_receptors, method = "ward.D2")
order_receptors = hclust_receptors$labels[hclust_receptors$order]

dist_ligands = dist(lr_network_top_matrix_strict %>% t(), method = "binary")
hclust_ligands = hclust(dist_ligands, method = "ward.D2")
order_ligands_receptor = hclust_ligands$labels[hclust_ligands$order]

order_receptors = order_receptors %>% intersect(rownames(lr_network_top_matrix_strict))
order_ligands_receptor = order_ligands_receptor %>% intersect(colnames(lr_network_top_matrix_strict))

vis_ligand_receptor_network_strict = lr_network_top_matrix_strict[order_receptors, order_ligands_receptor]


rownames(vis_ligand_receptor_network_strict) = order_receptors %>% make.names()
colnames(vis_ligand_receptor_network_strict) = order_ligands_receptor %>% make.names()
p_ligand_receptor_network_strict = vis_ligand_receptor_network_strict %>% t() %>% make_heatmap_ggplot("Ligands","Receptors", color = "mediumvioletred", x_axis_position = "top",legend_title = "Prior interaction potential\n(bona fide)")
p_ligand_receptor_network_strict

p_ligand_receptor_network_strict= vis_ligand_receptor_network_strict %>% t() %>% 
  make_threecolor_heatmap_ggplot("Ligands","Receptors", low_color = "#fff5f0", mid_color = "#ef3b2c", high_color = "#67000d", mid=0.4,
                                 x_axis_position = "top",legend_title = "Prior interaction potential\n(bona fide)")
p_ligand_receptor_network_strict

#4 Make reference matrix ---------------------------------------------------
#scRNAseq of polydactyly. 
#Processed data can be accessed from the NCBI Gene Expression Omnibus database (accession code GSE158970)
library(Seurat)

ON.integrated@meta.data

ref.dt<- ON.integrated@assays$RNA@data
table(colnames(ref.dt) == rownames(ON.integrated@meta.data))
colnames(ref.dt) <- ON.integrated$celltype

ON.integrated$celltype<- gsub("Fibro-[0-9]", "Fibroblast",  ON.integrated$subcelltype)

ON.integrated@meta.data$cell <- rownames(ON.integrated@meta.data)
sel.dt <- ON.integrated@meta.data[, c("cell", "celltype")]

cell_int<- c("Basal KC", "Cornified cells", "EC", "Eccrine duct", "Eccrine gl",  
             "Keratinocyte-1", "Keratinocyte-2", "LEC", "lymphocyte", "Macrophage", "Mast cells",
             "Melanocyte", "Fibroblast", "Myofibroblast", "OMD", 
             "proliferating KC",  "Suprabasal KC")

sel.dt <-sel.dt[sel.dt$celltype %in% cell_int,]
ref.dt <- ref.dt[,colnames(ref.dt) %in% cell_int]

table(colnames(ref.dt) == sel.dt$celltype)
table(names(colnames(ref.dt)) == sel.dt$cell)

Idents(ON.integrated) <-"celltype"
cellclass.markers <- FindAllMarkers(object = ON.integrated, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25)

cellclass.markers2<- cellclass.markers[cellclass.markers$cluster!="unknown",]
class.top100.genes <- cellclass.markers2 %>% group_by(cluster) %>% top_n(n = 75, wt = avg_log2FC)
class.top100.genes<- class.top100.genes$gene%>%unique()

MTgene<- grep("^MT", class.top100.genes,value = T)

class.top100.genes<- class.top100.genes[!class.top100.genes %in% MTgene ]
class.top100.genes <- unique(c(class.top100.genes, "LGR6", "LGR5","LGR4"))

ref.dt<- ref.dt[class.top100.genes,]
ref.dt<- as.data.frame(ref.dt)
ref.dt<-rownames_to_column(ref.dt)
colnames(ref.dt)[1] <- "Gene"
colnames(ref.dt)<- gsub("[.][0-9].*", "", colnames(ref.dt) )

write.table(ref.dt,file = "scRNA_skinappendage.txt", sep = "\t", row.names = F, quote = F)

#5 Spatial decon  -----------------------------------------------------------

library(SpatialDecon)

#load data
norm = as.matrix(Q3data)
raw = as.matrix(data)
annot = clin.info


# use the NegProbe to estimate per-observation background
per.observation.mean.neg = norm["NegProbe-WTX", ]

# and define a background matrix in which each column (observation) is the
# appropriate value of per-observation background:
#bg = sweep(norm * 0, 2, per.observation.mean.neg, "+")

bg<-Q3data

for(i in 1:length(colnames(Q3data))){
  bg[[i]] <- per.observation.mean.neg[[i]]
}

bg
dim(bg)

##Another refSet
##Cell profile matrix derived from single-cell RNA sequencing of polydactyly
refSkin<- read.delim(file = "./CIBERSORTx_Nail_sig2.txt")
refSkin<- refSkin %>% tibble::column_to_rownames("NAME")

refSkin <- as.matrix(refSkin)
dim(refSkin)

heatmap(sweep(refSkin, 1, apply(refSkin, 1, max), "/"),
        labRow = NA, margins = c(10, 5))

res = spatialdecon(norm = norm,
                   bg = bg,
                   X = refSkin,
                   align_genes = TRUE)

heatmap(res$beta, cexCol = 0.5, cexRow = 0.7, margins = c(10,7))

refskin.matches <- list("Melanocytes"=c("Melanocyte"),
                        "onychofibroblast"=c( "OMD"),
                        "LEC" =c("LEC"),
                        "EC" =c("EC"),
                        "NME" =c("Keratinocyte.1"),
                        "NBE" =c("Keratinocyte.2"),
                        "Basal.KC" =c("Basal.KC"),
                        "Suprabasal.KC" =c("Suprabasal.KC"),
                        "Cornified.cells" =c("Cornified.cells"),
                        "Proliferating.KC" =c("proliferating.KC"),
                        "Myofibroblast" =c("Myofibroblast" ),
                        "Fibroblast" =c ("Fibroblast"),
                        "Macrophage" =c("Macrophage"),
                        "Mast.cells" = c("Mast.cells"),
                        "Lymphocyte" = c("lymphocyte"),
                        "Eccrine"=c("Eccrine.duct", "Eccrine.gl"))

unlist(refskin.matches) %in% colnames(refSkin)

restils = spatialdecon(norm = norm,                     # normalized data
                       raw = raw,                       # raw data, used to down-weight low-count observations 
                       bg = bg,                         # expected background counts for every data point in norm
                       X = refSkin,                     # safeTME matrix, used by default
                       cellmerges = refskin.matches,   # safeTME.matches object, used by default
                       cell_counts = annot$nuclei,      # nuclei counts, used to estimate total cells
                       is_pure_tumor = NULL,   # identities of the Tumor segments/observations
                       n_tumor_clusters = 0)            # how many distinct tumor profiles to append to safeTME


heatmap(sweep(restils$X, 1, apply(restils$X, 1, max), "/"),
        labRow = NA, margins = c(10, 5))

data("cellcols")
cellcols
names(cellcols) <-  names(refskin.matches)
cellcols<- cellcols[1:16]

# show just the TME segments, since that's where the immune cells are:
layout(mat = (matrix(c(1, 2), 1)), widths = c(7, 4))

TIL_barplot(restils$cell.counts$cell.counts, draw_legend = TRUE, 
            cex.names = 0.5)

count.all<- restils$cell.counts$cell.counts
colnames(count.all) <-clin.info$AOI.name
TIL_barplot(count.all, draw_legend = TRUE, 
            cex.names = 0.5)

TIL_barplot(count.all, draw_legend = FALSE, 
            cex.names = 0.5)

# or the proportions of cells:
prop.all<-restils$prop_of_all
colnames(prop.all) <- clin.info$AOI.name

TIL_barplot(restils$prop_of_all, 
            draw_legend = TRUE, cex.names = 0.75)

TIL_barplot(prop.all, 
            draw_legend = TRUE, cex.names = 0.75)


prop.all<- prop.all[,order( colnames(prop.all))]
prop.all<- prop.all[,order( colnames(prop.all))]


TIL_barplot(prop.all, 
            draw_legend = TRUE, cex.names = 0.75)


