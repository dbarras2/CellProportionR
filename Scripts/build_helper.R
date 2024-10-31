
# Load metadata
# Filter out list (colnames, names)
# Cell Types and hierarchy

seurat<-readRDS("/Users/dbarras/Documents/CHUV/Bioinformatics_Core/Projects/Denarda/Ovarian_sc_Processing/Results/Seurat/UPENN_Cohort_Final_All.rds")
metadata<-seurat@meta.data


single_cell_data <- metadata
#single_cell_data <- single_cell_data[,c("Cell_Type_1","Cell_Type_CD4vsCD8","Cell_Type_Myeloid","Cell_Type_2","Technology","Patient")]
single_cell_data$Cell_Type_Strat4<-NA
single_cell_data$Cell_Type_Strat1<-single_cell_data$Cell_Type_1
single_cell_data$Cell_Type_Strat1[which(single_cell_data$Cell_Type_Strat1 %in% c("Myeloid_cells","T_cells","NK_cells","B_cells"))]<-"Immune"
single_cell_data$Cell_Type_Strat1[which(single_cell_data$Cell_Type_Strat1 %in% c("Stromal_cells","Malignant"))]<-"Non_Immune"
single_cell_data$Cell_Type_Strat2<-single_cell_data$Cell_Type_1
single_cell_data$Cell_Type_Strat3<-single_cell_data$Cell_Type_Strat2
idx<-which(single_cell_data$Cell_Type_1=="T_cells")
single_cell_data$Cell_Type_Strat3[idx]<-single_cell_data$Cell_Type_CD4vsCD8[idx]
idx<-which(single_cell_data$Cell_Type_1=="Myeloid_cells")
single_cell_data$Cell_Type_Strat3[idx]<-single_cell_data$Cell_Type_Myeloid[idx]
idx<-which(single_cell_data$Cell_Type_1=="Stromal_cells")
single_cell_data$Cell_Type_Strat3[idx]<-single_cell_data$Cell_Type_2[idx]
single_cell_data$Cell_Type_Strat4<-single_cell_data$Cell_Type_2

subset_data <- list("Technology"=c("5prime"),
                    "Patient"=c("1907","1682","1745","1787","1713"))
subset_data <- NULL

sample_colname <- "Patient"
stratification <- "Immunocategory"
group_up<-"purely_inflamed"
group_dn<-"excluded"
statistics = "wilcoxon"


single_cell_data$Cell_Type_Strat4<-factor(single_cell_data$Cell_Type_Strat4,levels=c("Malignant","CD8_Naive-like","CD8_EM-like","CD8_Pex","CD8_Tex","CD8_HSP","CD8_FOXP3","CD8_CX3CR1","CD8_ISG",
                                                                                     "CD8_Low-Quality","CD8_NK-like","CD8_MAIT","NK_prec","NK_cytotoxic","NK_reg","Tgd","DN","CD4_Th1","CD4_CXCL13","CD4_Tregs",
                                                                                     "CD4_ISG","CD4_GZMK","CD4_HSP","CD4_Low-Quality","pDC","DC1","DC2","DC3","DC_CD5","MonoDC","Monocytes_CD16",
                                                                                     "Macro_CXCL9","Macro_ISG","Macro_S100A8","Macro_TREM2","Macro_C1Q","Macro_MMP9","Macro_Low-Quality","Macro_HSP",
                                                                                     "Mast_cells","Neutrophils","B_cells_Naive","B_cells_Memory","B_cells_GC_CD38_MEF2B","B_cells_ISG","B_cells_HSP",
                                                                                     "B_cells_Low-Quality","Plasma_cells","Endothelial_cells","CAFs","Pericytes","Doublet"))


group_up_col<-NULL
group_dn_col<-NULL
include_ratio<-TRUE
pval_range<-NULL
stratification_names<-c("Immune_vs_Not_vs_Doublet","Mal_T_B_Myeloid_Stromal_NK",
                        "Mal_CD8_CD4_DN_Tgd_B_Macro_Mono_DC_Endo_CAFs_Pericyte","Fine_Cell_Types")



single_cell_data <- single_cell_data[,-which(colnames(single_cell_data) %in% c("Cell_Type_Strat1"))]
stratification_names<-c("Fine_Cell_Types")





single_cell_data <- single_cell_data[,c("Patient","Immunocategory","Technology","Cell_Type_Strat1","Cell_Type_Strat2","Cell_Type_Strat3","Cell_Type_Strat4")]
t<-sample(1:nrow(single_cell_data),size = 20000)
single_cell_data <- single_cell_data[t,]
single_cell_data<-single_cell_data[order(single_cell_data$Patient),]
colnames(single_cell_data)[2]<-"Category"
single_cell_data$Category[which(single_cell_data$Category %in% c("desert","excluded"))]<-"group2"
single_cell_data$Category[which(single_cell_data$Category %in% c("purely_inflamed","mixed_inflamed"))]<-"group1"

table(single_cell_data$Patient)
single_cell_data$Patient <- paste0(single_cell_data$Patient ,"231")

patient_names <- unique(single_cell_data$Patient)

# Create a vector of anonymized patient names
anonymized_names <- paste0("Patient", seq_along(patient_names))
names(anonymized_names)<-patient_names
single_cell_data$Patient<-anonymized_names[match(single_cell_data$Patient,names(anonymized_names))]
rownames(single_cell_data)<-paste0("Cell_",seq(1:nrow(single_cell_data)))
save(single_cell_data,file="/Users/dbarras/My_R_packages/CellProportionR/data/example_data.rda")

CellType_Proportion_Heatmap(single_cell_data = single_cell_data,
                            sample_colname = "Patient",
                            subset_data = list("Technology"=c("5prime"),
                                               "Patient"=c("Patient11","Patient27","Patient28","Patient32",
                                                           "Patient31","Patient17","Patient25","Patient30")),
                            stratification = "Category",
                            group_up = c("group1"),
                            group_dn=c("group2"),
                            stratification_names = c("Immune_vs_Not_vs_Doublet",
                                                     "Mal_T_B_Myeloid_Stromal_NK",
                                                     "Mal_CD8_CD4_DN_Tgd_B_Macro_Mono_DC_Endo_CAFs_Pericyte",
                                                     "Fine_Cell_Types"))

table(single_cell_data$Patient[which(single_cell_data$Category=="group2")])

levels(single_cell_data$Cell_Type_Strat4)
single_cell_data$Cell_Type_Strat4 <- droplevels(single_cell_data$Cell_Type_Strat4)
levels(single_cell_data$Cell_Type_Strat4)
table(single_cell_data$Cell_Type_Strat4)

load("/Users/dbarras/My_R_packages/CellProportionR/data/example_data.rda")
check <-Compute_Proportions_Ratios(single_cell_data = single_cell_data,sample_colname = "Patient",stratification = "Category",group_up = "group1",group_dn = "group2")




cell_proportion_object <- Compute_Proportions_Ratios(single_cell_data = single_cell_data,sample_colname = "Patient",stratification = "Category",group_up = "group1",group_dn = "group2")
cell_proportion_object$excluded<-"no"
cell_proportion_object$excluded[which(cell_proportion_object$Patient=="Patient1")]<-"yes"
subset_data = list("excluded"="no")
stratification = "Stratif"
groups = list("up"="up","dn"="dn")
type = "Ratio"
cell_type = "CD8_EM-like"
out_of = NULL
cell_type_ratio = "Neutrophils"
col_var<-c("red","blue")


load("/Users/dbarras/My_R_packages/CellProportionR/data/example_data.rda")
Proportion_object <- Compute_Proportions_Ratios(single_cell_data = single_cell_data,sample_colname = "Patient",stratification = "Category",group_up = "group1",group_dn = "group2")

Proportion_object$excluded<-"no"
Proportion_object$excluded[which(Proportion_object$Patient=="Patient1")]<-"yes"

CellType_Proportion_Boxplot(cell_proportion_object = Proportion_object,
                            subset_data = list("excluded"="no"),
                            stratification="Stratif",
                            groups = list("up"="up","dn"="dn"),
                            type = "Ratio",
                            cell_type = "CD8_EM-like",
                            out_of = NULL,
                            cell_type_ratio = "Neutrophils",
                            col_var = c("red","blue"))

CellType_Proportion_Boxplot(cell_proportion_object = Proportion_object,
                            subset_data = list("excluded"="no"),
                            stratification="Stratif",
                            groups = list("up"="up","dn"="dn"),
                            type = "Proportion",
                            cell_type = "CD8_Tex",
                            out_of = "CD8_T_cells",
                            cell_type_ratio = NULL,
                            col_var = c("red","blue"))

library(devtools)
load_all("/Users/dbarras/My_R_packages/CellProportionR")
document("/Users/dbarras/My_R_packages/CellProportionR")
build("/Users/dbarras/My_R_packages/CellProportionR")




load("/Users/dbarras/Downloads/test.rda")

single_cell_data = object
sample_colname = "Patient"
stratification = "treatment_phase"
group_up="post-NACT"
group_dn = "treatment-naive"
stratification_names = c("Main_Lineage","Fine_Subtypes")
