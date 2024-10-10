#' @export
CellType_Proportion_Heatmap <- function(Cohort_selection = c("ATATIL","NeoTIL","SOLTIL"),
                                        Tumor_Type_selection = NULL,
                                        Timepoint = c("T0"),
                                        Stratification,
                                        group_up,
                                        group_dn,
                                        statistics = "wilcoxon",
                                        excluded_patient = NULL,
                                        breaks = NULL){


  require(ComplexHeatmap)
  require(limma)
  require(reshape2)
  require(crayon)
  require(circlize)
  require(usethis)
  require(scales)

  if(length(which(Cohort_selection  %in% c("ATATIL","NeoTIL","SOLTIL")==F))>0){
    warning("The selected cohorts should be `ATATIL`, `NeoTIL` and/or `SOLTIL`")
    return()
  }
  if(length(which(Timepoint  %in% c("T0","T30")==F))>0){
    warning("The selected time-points should be `T0` (default; baseline), and/or `T30` (post-ACT)")
    return()
  }
  if(is.na(match(statistics,c("wilcoxon","ttest")))){
    warning("The selected statistics should be `wilcoxon` or `ttest`")
    return()
  }



## Load Datasets
#########################################################################

if("NeoTIL" %in% Cohort_selection){
  data("Clin_NeoTIL_all", package = "CellProportionACT")
  data("Clin_NeoTIL_CD45", package = "CellProportionACT")
  data("meta_NeoTIL_all", package = "CellProportionACT")
  data("meta_NeoTIL_CD45", package = "CellProportionACT")
}

if("ATATIL" %in% Cohort_selection){
  data("Clin_ATATIL_all", package = "CellProportionACT")
  data("Clin_ATATIL_CD45", package = "CellProportionACT")
  data("meta_ATATIL_all", package = "CellProportionACT")
  data("meta_ATATIL_CD45", package = "CellProportionACT")
}

if("SOLTIL" %in% Cohort_selection){
  data("Clin_SOLTIL_CD45", package = "CellProportionACT")
  data("meta_SOLTIL_CD45", package = "CellProportionACT")
}

#########################################################################

## Combine Data
#########################################################################

datasets<-c("ATATIL_all","ATATIL_CD45","NeoTIL_all","NeoTIL_CD45","SOLTIL_CD45")
datasets<-datasets[!is.na(match(strsplit2(datasets,split="_")[,1],Cohort_selection))]

if(length(which(is.na(datasets)))>0){datasets<-datasets[-which(is.na(datasets))]}
for(i in 1:length(datasets)){
  if(i==1){
    Clin<-get(paste0("Clin_",datasets[i]))
    meta<-get(paste0("meta_",datasets[i]))
  }
  if(i!=1){
    Clin_tmp<-get(paste0("Clin_",datasets[i]))
    meta_tmp<-get(paste0("meta_",datasets[i]))

    common<-intersect(colnames(Clin),colnames(Clin_tmp))
    Clin<-rbind(Clin[,match(common,colnames(Clin))],
                Clin_tmp[,match(common,colnames(Clin_tmp))])

    common<-intersect(colnames(meta),colnames(meta_tmp))
    meta<-rbind(meta[,match(common,colnames(meta))],
                meta_tmp[,match(common,colnames(meta_tmp))])
    rm(Clin_tmp,meta_tmp)
  }
}
rm(list = ls(pattern="Clin_"))
rm(list = ls(pattern="meta_"))



#########################################################################

## Select wanted cancer type if applicable and exclude potential patients
#########################################################################
# Select Cancer Type
if(!is.null(Tumor_Type_selection)){

  if(length(which(Tumor_Type_selection  %in% Clin$Cancer_Type ==F))>0){
    warning("The selected cancer types should be `Melanoma`, `Lung` and/or `Other_Solid`")
    return()
  }

  patients <- Clin$Patient_ID[which(Clin$Cancer_Type %in% Tumor_Type_selection)]
  Clin<-Clin[which(Clin$Patient_ID %in% patients),]
  meta<-meta[which(meta$Patient %in% patients),]
}
# Select Correct time-points
if(!is.null(Timepoint)){
  meta<-meta[which(meta$Time %in% Timepoint),]
  Clin<-Clin[which(Clin$Patient_ID %in% unique(meta$Patient)),]
}
# remove patient to exclude
if(!is.null(excluded_patient)){
  idx<-as.vector(unlist(sapply(excluded_patient,function(x){grep(x,Clin$Patient_ID)})))
  if(length(idx)==0){
    warning(paste0("In the selected cohort, the patients that are available for exclusion are :",
                   paste(apply(strsplit2(Clin$Patient_ID,split="_")[,c(1,2)],1,function(x){paste(x,collapse="_")}),collapse=", ")))
    return()
  }

  if(length(idx)>0){
    Clin<-Clin[-idx,]
    meta<-meta[which(meta$Patient %in% Clin$Patient_ID),]
  }
}
#########################################################################

## Define categorizing factor and groups
#########################################################################

idx<-match(Stratification,colnames(Clin))
if(is.na(idx)){
  warning(paste0("Stratification factor not found in clinical data; available ones are: ",paste(colnames(Clin)[-which(colnames(Clin) %in% c("Patient_ID","Use_all_cells","Use_CD45_cells"))],collapse=", ")
  ))
  return()
}

if(!is.na(idx)){
  Clin$Stratif<-Clin[,Stratification]

  message(crayon::silver(paste0("Stratification factor found in clinical data")))

  idx_up<-which(Clin$Stratif %in% group_up)
  idx_dn<-which(Clin$Stratif %in% group_dn)
  if((length(idx_up)==0)|(length(idx_dn)==0)){
    warning(paste0("The selected groups within the stratification factor are not present, please chose between these calls: ",
                   paste(as.character(unique(Clin$Stratif)),collapse=", ")))
    return()
  }

  tmp<-Clin[idx_up,]
  l1<-length(unique(apply(strsplit2(tmp$Patient_ID[which(tmp$Use_CD45_cells=="yes")],"_")[,c(1,2)],1,function(x){paste(x,collapse="_")})))
  l2<-length(unique(apply(strsplit2(tmp$Patient_ID[which(tmp$Use_all_cells=="yes")],"_")[,c(1,2)],1,function(x){paste(x,collapse="_")})))
  message(crayon::red(paste0(l1," ",paste(group_up,collapse="/")," will be used for CD45 sorted cells and ",l2," for all-viable cells")))


  tmp<-Clin[idx_dn,]
  l1<-length(unique(apply(strsplit2(tmp$Patient_ID[which(tmp$Use_CD45_cells=="yes")],"_")[,c(1,2)],1,function(x){paste(x,collapse="_")})))
  l2<-length(unique(apply(strsplit2(tmp$Patient_ID[which(tmp$Use_all_cells=="yes")],"_")[,c(1,2)],1,function(x){paste(x,collapse="_")})))
  message(crayon::blue(paste0(l1," ",paste(group_dn,collapse="/")," will be used for CD45 sorted cells and ",l2," for all-viable cells")))
  rm(tmp)


  if(length(intersect(idx_up,idx_dn)) > 0){
    warning("Group1 and Group2 are overlapping - Please provide non-overlapping categories")
    return()
  }

  Clin$Stratif[idx_up]<-"up"
  Clin$Stratif[idx_dn]<-"dn"

  Clin<-Clin[which(Clin$Stratif %in% c("up","dn")),]
  meta<-meta[which(meta$Patient %in% Clin$Patient_ID),]
  Clin$Stratif<-factor(Clin$Stratif,levels = c("up","dn"))
}
#########################################################################

# List and stratification of all cell types
#########################################################################
column_cell_type<-"Cell_Type_Final"
filter_absent_cell_type <- function(celltypes,column_cell_type){
  celltypes<-celltypes[which(celltypes %in% unique(meta[,column_cell_type]))]
  return(celltypes)
}

Malignant<-c("Malignant")
#Malignant<-filter_absent_cell_type(Malignant,column_cell_type)

CD8_T_cells<-c("CD8_Naive-like","CD8_EM-like","CD8_Pex","CD8_Tex","CD8_HSP","CD8_FOXP3","CD8_CX3CR1","CD8_ISG","CD8_Low-Quality","CD8_NK-like","CD8_MAIT")
#CD8_T_cells<-filter_absent_cell_type(CD8_T_cells,column_cell_type)

NK_tgd<-c("NK_cells","Tgd")
#NK_tgd<-filter_absent_cell_type(NK_tgd,column_cell_type)

CD4_T_cells<-c("CD4_Th1","CD4_CXCL13","CD4_Tregs","CD4_ISG","CD4_EM-like","CD4_HSP","CD4_Low-Quality")
#CD4_T_cells<-filter_absent_cell_type(CD4_T_cells,column_cell_type)

DC_cells<-c("pDC","DC1","DC2","DC3","DC_CD5")
#DC_cells<-filter_absent_cell_type(DC_cells,column_cell_type)

Macro_cells<-c("MonoDC","Monocytes_CD16","Macro_CXCL9","Macro_ISG","Macro_S100A8","Macro_TREM2","Macro_C1Q","Macro_MMP9","Macro_Low-Quality","Macro_HSP")
#Macro_cells<-filter_absent_cell_type(Macro_cells,column_cell_type)

Myelo_non_Macro<-c("Mast_cells","Neutrophils")
#Myelo_non_Macro<-filter_absent_cell_type(Myelo_non_Macro,column_cell_type)

B_cells<-c("B_cells_Naive","B_cells_Memory","B_cells_GC_CD38_MEF2B","B_cells_ISG","B_cells_HSP","B_cells_Low-Quality","Plasma_cells")
#B_cells<-filter_absent_cell_type(B_cells,column_cell_type)

Stromal<-c("Endothelial_cells","CAFs")
#Stromal<-filter_absent_cell_type(Stromal,column_cell_type)

#Doublets<-c("Doublets_T_Myelo_cells","Doublets_B_T_cells","Doublets_B_Myelo_cells","Doublets_CD4_CD8_T_cells","Doublet")
#Others<-c("Apoptotic_cells","not_assigned","Low-Quality","Prolif","Blood_cells")


cell_subtypes<-c(Malignant,CD8_T_cells,NK_tgd,CD4_T_cells,DC_cells,Macro_cells,Myelo_non_Macro,B_cells,Stromal)
cell_types<-c("CD8","CD4","T_cells","DC","Macrophages","Myeloid","B_cells","Stromal")
CD45_positive<-c(CD8_T_cells,NK_tgd,CD4_T_cells,DC_cells,Macro_cells,Myelo_non_Macro,B_cells)

# Normalized to
row_norm<-c("All","CD45","T_Myelo_B","CD8_CD4_DC_Macro_B")

collect_stats <- as.data.frame(matrix(nrow=length(c(row_norm,cell_subtypes,cell_types)),ncol=length(c(cell_subtypes,cell_types)),NA))
rownames(collect_stats)<-c(row_norm,cell_subtypes,cell_types)
colnames(collect_stats)<-c(cell_subtypes,cell_types)

meta$Cell_Type_strat_1<-"others"
meta$Cell_Type_strat_1[which(meta[,column_cell_type] %in% c(CD8_T_cells,"Tgd",CD4_T_cells))]<-"T_cells"
meta$Cell_Type_strat_1[which(meta[,column_cell_type] %in% c(DC_cells,Macro_cells,Myelo_non_Macro))]<-"Myeloid"
meta$Cell_Type_strat_1[which(meta[,column_cell_type] %in% c(B_cells))]<-"B_cells"
meta$Cell_Type_strat_1[which(meta[,column_cell_type] %in% c(Stromal))]<-"Stromal"

meta$Cell_Type_strat_2<-"others"
meta$Cell_Type_strat_2[which(meta[,column_cell_type] %in% c(CD8_T_cells))]<-"CD8"
meta$Cell_Type_strat_2[which(meta[,column_cell_type] %in% c(CD4_T_cells))]<-"CD4"
meta$Cell_Type_strat_2[which(meta[,column_cell_type] %in% c(DC_cells))]<-"DC"
meta$Cell_Type_strat_2[which(meta[,column_cell_type] %in% c(Macro_cells))]<-"Macrophages"
meta$Cell_Type_strat_2[which(meta[,column_cell_type] %in% c(B_cells))]<-"B_cells"
meta$Cell_Type_strat_2[which(meta[,column_cell_type] %in% c(Stromal))]<-"Stromal"

#########################################################################


## Split into CD45 and all cells data
#########################################################################
meta_all<-meta[which(meta$Patient %in% Clin$Patient_ID[which(Clin$Use_all_cells == "yes")]),]
meta_CD45<-meta[which(meta$Patient %in% Clin$Patient_ID[which(Clin$Use_CD45_cells == "yes")]),]
rm(meta)
#########################################################################

makeStatisticalTest <- function(Imm1,Imm2,test){
  result<-list()
  if((length(Imm1[!is.na(Imm1)])<1)|(length(Imm2[!is.na(Imm2)])<1)){
    result$statistics<-NA
    result$p.value<-NA
    return(result)
  }
  if(test=="ttest"){
    t<-suppressWarnings(t.test(Imm1,Imm2))
    result$statistics<-t$statistic
    result$p.value<-t$p.value
    return(result)
  }
  if(test=="wilcoxon"){
    t<-suppressWarnings(wilcox.test(Imm1,Imm2))
    result$p.value<-t$p.value
    val<-ifelse((median(Imm1)-median(Imm2))==0,(mean(Imm1)-mean(Imm2)),(median(Imm1)-median(Imm2)))
    result$statistics<-val
    return(result)
  }

}

# Out of CD45 cells
########################################################################

idx_cd45 <- which(meta_CD45[,column_cell_type] %in% CD45_positive)

Proportions<-prop.table(table(meta_CD45[idx_cd45,which(colnames(meta_CD45) %in% c("Patient",column_cell_type))]),1)
Proportions<-melt(Proportions)
Proportions$Stratif<-Clin$Stratif[match(Proportions$Patient,Clin$Patient_ID)]
colnames(Proportions)[2]<-"Cell_Type"

for(i in 1:length(CD45_positive)){
  tmp1<-Proportions[which(Proportions$Cell_Type==CD45_positive[i]),]
  Imm1<-tmp1$value[which(tmp1$Stratif=="up")];names(Imm1)<-as.character(tmp1$Patient[which(tmp1$Stratif=="up")])
  Imm2<-tmp1$value[which(tmp1$Stratif=="dn")];names(Imm2)<-as.character(tmp1$Patient[which(tmp1$Stratif=="dn")])
  t<-makeStatisticalTest(Imm1[!is.na(Imm1)],Imm2[!is.na(Imm2)],statistics)
  collect_stats["CD45",CD45_positive[i]]<-ifelse(t$statistic>0,(-log10(t$p.value)),-(-log10(t$p.value)))

}

strats<-c("Cell_Type_strat_1","Cell_Type_strat_2")
for(strat in strats){
  Proportions<-prop.table(table(meta_CD45[idx_cd45,which(colnames(meta_CD45) %in% c("Patient",strat))]),1)
  Proportions<-melt(Proportions)
  Proportions$Stratif<-Clin$Stratif[match(Proportions$Patient,Clin$Patient_ID)]
  colnames(Proportions)[2]<-"Cell_Type"

  for(ct in as.character(unique(Proportions$Cell_Type))){
    tmp1<-Proportions[which(Proportions$Cell_Type==ct),]
    Imm1<-tmp1$value[which(tmp1$Stratif=="up")];names(Imm1)<-as.character(tmp1$Patient[which(tmp1$Stratif=="up")])
    Imm2<-tmp1$value[which(tmp1$Stratif=="dn")];names(Imm2)<-as.character(tmp1$Patient[which(tmp1$Stratif=="dn")])
    t<-makeStatisticalTest(Imm1[!is.na(Imm1)],Imm2[!is.na(Imm2)],statistics)
    collect_stats["CD45",ct]<-ifelse(t$statistic>0,(-log10(t$p.value)),-(-log10(t$p.value)))
  }

}
########################################################################

# Out of T-Myelo-B
########################################################################

if(1){

  runs<-list(c(CD8_T_cells,"Tgd",CD4_T_cells),
             c(DC_cells,Macro_cells,Myelo_non_Macro),
             c(B_cells))

  for(r in 1:length(runs)){

    idx_cell_type<-which(meta_CD45[,column_cell_type] %in% runs[[r]])

    Proportions<-prop.table(table(meta_CD45[idx_cell_type,which(colnames(meta_CD45) %in% c("Patient",column_cell_type))]),1)
    Proportions<-melt(Proportions)
    Proportions$Stratif<-Clin$Stratif[match(Proportions$Patient,Clin$Patient_ID)]
    colnames(Proportions)[2]<-"Cell_Type"

    for(i in 1:length(runs[[r]])){
      tmp1<-Proportions[which(Proportions$Cell_Type==runs[[r]][i]),]
      Imm1<-tmp1$value[which(tmp1$Stratif=="up")];names(Imm1)<-as.character(tmp1$Patient[which(tmp1$Stratif=="up")])
      Imm2<-tmp1$value[which(tmp1$Stratif=="dn")];names(Imm2)<-as.character(tmp1$Patient[which(tmp1$Stratif=="dn")])
      t<-makeStatisticalTest(Imm1[!is.na(Imm1)],Imm2[!is.na(Imm2)],statistics)
      collect_stats["T_Myelo_B",runs[[r]][i]]<-ifelse(t$statistic>0,(-log10(t$p.value)),-(-log10(t$p.value)))
    }

    strats<-c("Cell_Type_strat_1","Cell_Type_strat_2")
    for(strat in strats){
      Proportions<-prop.table(table(meta_CD45[idx_cell_type,which(colnames(meta_CD45) %in% c("Patient",strat))]),1)
      Proportions<-melt(Proportions)
      Proportions$Stratif<-Clin$Stratif[match(Proportions$Patient,Clin$Patient_ID)]
      colnames(Proportions)[2]<-"Cell_Type"

      sub<-as.character(unique(Proportions$Cell_Type))[which( as.character(unique(Proportions$Cell_Type)) %in% c("CD8","CD4"))]

      for(ct in sub){
        tmp1<-Proportions[which(Proportions$Cell_Type==ct),]
        Imm1<-tmp1$value[which(tmp1$Stratif=="up")];names(Imm1)<-as.character(tmp1$Patient[which(tmp1$Stratif=="up")])
        Imm2<-tmp1$value[which(tmp1$Stratif=="dn")];names(Imm2)<-as.character(tmp1$Patient[which(tmp1$Stratif=="dn")])
        t<-makeStatisticalTest(Imm1[!is.na(Imm1)],Imm2[!is.na(Imm2)],statistics)
        collect_stats["T_Myelo_B",ct]<-ifelse(t$statistic>0,(-log10(t$p.value)),-(-log10(t$p.value)))
      }
    }
  }
}

########################################################################

# Out of CD8-CD4-DC-Myelo-B
########################################################################

runs<-list(list("CD8",CD8_T_cells),
           list("CD4",CD4_T_cells),
           list("DC",DC_cells),
           list("Macrophages",Macro_cells),
           list("B_cells",B_cells))

for(r in 1:length(runs)){

  idx_cells <- which(meta_CD45$Cell_Type_strat_2 %in% runs[[r]][1])

  Proportions<-prop.table(table(meta_CD45[idx_cells,which(colnames(meta_CD45) %in% c("Patient",column_cell_type))]),1)
  Proportions<-melt(Proportions)
  Proportions$Stratif<-Clin$Stratif[match(Proportions$Patient,Clin$Patient_ID)]
  colnames(Proportions)[2]<-"Cell_Type"


  for(i in 1:length(runs[[r]][2][[1]])){
    tmp1<-Proportions[which(Proportions$Cell_Type==runs[[r]][2][[1]][i]),]
    Imm1<-tmp1$value[which(tmp1$Stratif=="up")];names(Imm1)<-as.character(tmp1$Patient[which(tmp1$Stratif=="up")])
    Imm2<-tmp1$value[which(tmp1$Stratif=="dn")];names(Imm2)<-as.character(tmp1$Patient[which(tmp1$Stratif=="dn")])
    t<-makeStatisticalTest(Imm1[!is.na(Imm1)],Imm2[!is.na(Imm2)],statistics)
    collect_stats["CD8_CD4_DC_Macro_B",runs[[r]][2][[1]][i]]<-ifelse(t$statistic>0,(-log10(t$p.value)),-(-log10(t$p.value)))
  }
}

########################################################################


# Ratio of cell types according to CD45
########################################################################

runs<-list(list(Malignant,CD8_T_cells),
           list(Malignant,CD4_T_cells),
           list(Malignant,NK_tgd),
           list(Malignant,DC_cells),
           list(Malignant,Macro_cells),
           list(Malignant,Myelo_non_Macro),
           list(Malignant,B_cells),
           list(Malignant,Stromal),
           list(CD8_T_cells,Malignant),
           list(CD8_T_cells,CD8_T_cells),
           list(CD8_T_cells,CD4_T_cells),
           list(CD8_T_cells,NK_tgd),
           list(CD8_T_cells,DC_cells),
           list(CD8_T_cells,Macro_cells),
           list(CD8_T_cells,Myelo_non_Macro),
           list(CD8_T_cells,B_cells),
           list(CD8_T_cells,Stromal),
           list(NK_tgd,Malignant),
           list(NK_tgd,CD8_T_cells),
           list(NK_tgd,CD4_T_cells),
           list(NK_tgd,NK_tgd),
           list(NK_tgd,DC_cells),
           list(NK_tgd,Macro_cells),
           list(NK_tgd,Myelo_non_Macro),
           list(NK_tgd,B_cells),
           list(NK_tgd,Stromal),
           list(CD4_T_cells,Malignant),
           list(CD4_T_cells,CD8_T_cells),
           list(CD4_T_cells,CD4_T_cells),
           list(CD4_T_cells,NK_tgd),
           list(CD4_T_cells,DC_cells),
           list(CD4_T_cells,Macro_cells),
           list(CD4_T_cells,Myelo_non_Macro),
           list(CD4_T_cells,B_cells),
           list(CD4_T_cells,Stromal),
           list(DC_cells,Malignant),
           list(DC_cells,CD8_T_cells),
           list(DC_cells,CD4_T_cells),
           list(DC_cells,NK_tgd),
           list(DC_cells,DC_cells),
           list(DC_cells,Macro_cells),
           list(DC_cells,Myelo_non_Macro),
           list(DC_cells,B_cells),
           list(DC_cells,Stromal),
           list(Macro_cells,Malignant),
           list(Macro_cells,CD8_T_cells),
           list(Macro_cells,CD4_T_cells),
           list(Macro_cells,NK_tgd),
           list(Macro_cells,DC_cells),
           list(Macro_cells,Macro_cells),
           list(Macro_cells,Myelo_non_Macro),
           list(Macro_cells,B_cells),
           list(Macro_cells,Stromal),
           list(Myelo_non_Macro,Malignant),
           list(Myelo_non_Macro,CD8_T_cells),
           list(Myelo_non_Macro,CD4_T_cells),
           list(Myelo_non_Macro,NK_tgd),
           list(Myelo_non_Macro,DC_cells),
           list(Myelo_non_Macro,Macro_cells),
           list(Myelo_non_Macro,Myelo_non_Macro),
           list(Myelo_non_Macro,B_cells),
           list(Myelo_non_Macro,Stromal),
           list(B_cells,Malignant),
           list(B_cells,CD8_T_cells),
           list(B_cells,CD4_T_cells),
           list(B_cells,NK_tgd),
           list(B_cells,DC_cells),
           list(B_cells,Macro_cells),
           list(B_cells,Myelo_non_Macro),
           list(B_cells,B_cells),
           list(B_cells,Stromal))

idx_cd45 <- which(meta_CD45[,column_cell_type] %in% CD45_positive)

Proportions<-prop.table(table(meta_CD45[idx_cd45,which(colnames(meta_CD45) %in% c("Patient",column_cell_type))]),1)
Proportions<-melt(Proportions)
Proportions$Stratif<-Clin$Stratif[match(Proportions$Patient,Clin$Patient_ID)]
colnames(Proportions)[2]<-"Cell_Type"


for(r in 1:length(runs)){

  ct1<-runs[[r]][[1]]
  ct2<-runs[[r]][[2]]

  for(i in 1:length(ct1)){
    for(j in 1:length(ct2)){

      if((ct1[i] != ct2[j])&(ct1[i] %in% Proportions$Cell_Type)&(ct2[j] %in% Proportions$Cell_Type)){
        tmp1<-Proportions[which(Proportions$Cell_Type==ct1[i]),]
        tmp2<-Proportions[which(Proportions$Cell_Type==ct2[j]),]
        tmp2<-tmp2[match(tmp1$Patient,tmp2$Patient),]

        #tmp1$Ratio<-tmp1$value/tmp2$value
        tmp1$Ratio<-log2(tmp1$value/tmp2$value)

        if(length(which(tmp1$Ratio==(-Inf)))>0){tmp1$Ratio[which(tmp1$Ratio==(-Inf))]<-NA}
        if(length(which(tmp1$Ratio==(Inf)))>0){tmp1$Ratio[which(tmp1$Ratio==(Inf))]<-NA}

        Imm1<-tmp1$Ratio[which(tmp1$Stratif=="up")];names(Imm1)<-as.character(tmp1$Patient[which(tmp1$Stratif=="up")]);Imm1<-Imm1[complete.cases(Imm1)]
        Imm2<-tmp1$Ratio[which(tmp1$Stratif=="dn")];names(Imm2)<-as.character(tmp1$Patient[which(tmp1$Stratif=="dn")]);Imm2<-Imm2[complete.cases(Imm2)]

        if((length(Imm1)>0)&(length(Imm2)>0)){
          t<-makeStatisticalTest(Imm1[!is.na(Imm1)],Imm2[!is.na(Imm2)],statistics)
          collect_stats[ct2[j],ct1[i]] <- ifelse(t$statistic>0,(-log10(t$p.value)),-(-log10(t$p.value)))
        }
      }
    }
  }
}

# Now for the grouped
strats<-c("Cell_Type_strat_1","Cell_Type_strat_2")
for(strat in strats){
  Proportions_grouped<-prop.table(table(meta_CD45[idx_cd45,which(colnames(meta_CD45) %in% c("Patient",strat))]),1)
  Proportions_grouped<-melt(Proportions_grouped)
  Proportions_grouped$Stratif<-Clin$Stratif[match(Proportions_grouped$Patient,Clin$Patient_ID)]
  colnames(Proportions_grouped)[2]<-"Cell_Type"

  runs<-list(list("CD8",CD4_T_cells),
             list("CD8",NK_tgd),
             list("CD8",DC_cells),
             list("CD8",Macro_cells),
             list("CD8",B_cells),
             list("CD4",CD8_T_cells),
             list("CD4",NK_tgd),
             list("CD4",DC_cells),
             list("CD4",Macro_cells),
             list("CD4",B_cells),
             list("DC",CD8_T_cells),
             list("DC",NK_tgd),
             list("DC",CD4_T_cells),
             list("DC",Macro_cells),
             list("DC",B_cells),
             list("Macrophages",CD8_T_cells),
             list("Macrophages",NK_tgd),
             list("Macrophages",CD4_T_cells),
             list("Macrophages",DC_cells),
             list("Macrophages",B_cells),
             list("T_cells",DC_cells),
             list("T_cells",Macro_cells),
             list("T_cells",B_cells),
             list("Myeloid",CD8_T_cells),
             list("Myeloid",NK_tgd),
             list("Myeloid",CD4_T_cells),
             list("Myeloid",B_cells),
             list("B_cells",CD8_T_cells),
             list("B_cells",NK_tgd),
             list("B_cells",CD4_T_cells),
             list("B_cells",DC_cells),
             list("B_cells",Macro_cells))

  keep<-c()
  for(k in 1:length(runs)){
    keep<-c(keep,ifelse(runs[[k]][[1]] %in% unique(Proportions_grouped$Cell_Type),k,NA))
  }
  keep<-keep[complete.cases(keep)]
  runs<-runs[keep]

  for(r in 1:length(runs)){
    ct1<-runs[[r]][[1]]
    ct2<-runs[[r]][[2]]
    for(j in 1:length(ct2)){
      tmp1<-Proportions_grouped[which(Proportions_grouped$Cell_Type==ct1),]
      tmp2<-Proportions[which(Proportions$Cell_Type==ct2[j]),]
      tmp2<-tmp2[match(tmp1$Patient,tmp2$Patient),]

      tmp1$Ratio<-log2(tmp1$value/tmp2$value)

      if(length(which(tmp1$Ratio==Inf))>0){tmp1$Ratio[which(tmp1$Ratio==Inf)]<-NA}
      if(length(which(tmp1$Ratio==(-Inf)))>0){tmp1$Ratio[which(tmp1$Ratio==(-Inf))]<-NA}

      Imm1<-tmp1$Ratio[which(tmp1$Stratif=="up")];names(Imm1)<-as.character(tmp1$Patient[which(tmp1$Stratif=="up")])
      Imm2<-tmp1$Ratio[which(tmp1$Stratif=="dn")];names(Imm2)<-as.character(tmp1$Patient[which(tmp1$Stratif=="dn")])

      if((length(which(!is.na(Imm1)==TRUE))>1)&(length(which(!is.na(Imm2)==TRUE))>1)){
        t<-makeStatisticalTest(Imm1[!is.na(Imm1)],Imm2[!is.na(Imm2)],statistics)
        collect_stats[ct2[j],ct1] <- ifelse(t$statistic>0,(-log10(t$p.value)),-(-log10(t$p.value)))
      }

    }
    for(j in 1:length(ct2)){
      tmp1<-Proportions_grouped[which(Proportions_grouped$Cell_Type==ct1),]
      tmp2<-Proportions[which(Proportions$Cell_Type==ct2[j]),]
      tmp2<-tmp2[match(tmp1$Patient,tmp2$Patient),]

      #tmp1$Ratio<-tmp2$value/tmp1$value
      tmp1$Ratio<-log2(tmp2$value/tmp1$value)

      if(length(which(tmp1$Ratio==Inf))>0){tmp1$Ratio[which(tmp1$Ratio==Inf)]<-NA}
      if(length(which(tmp1$Ratio==(-Inf)))>0){tmp1$Ratio[which(tmp1$Ratio==(-Inf))]<-NA}

      Imm1<-tmp1$Ratio[which(tmp1$Stratif=="up")];names(Imm1)<-as.character(tmp1$Patient[which(tmp1$Stratif=="up")])
      Imm2<-tmp1$Ratio[which(tmp1$Stratif=="dn")];names(Imm2)<-as.character(tmp1$Patient[which(tmp1$Stratif=="dn")])

      if((length(which(!is.na(Imm1)==TRUE))>1)&(length(which(!is.na(Imm2)==TRUE))>1)){
        t<-makeStatisticalTest(Imm1[!is.na(Imm1)],Imm2[!is.na(Imm2)],statistics)
        collect_stats[ct1,ct2[j]] <- ifelse(t$statistic>0,(-log10(t$p.value)),-(-log10(t$p.value)))
      }


    }
  }



}


# Grouped between them
strat<-"Cell_Type_strat_1"
Proportions_gr1<-prop.table(table(meta_CD45[,which(colnames(meta_CD45) %in% c("Patient",strat))]),1)
Proportions_gr1<-melt(Proportions_gr1)
Proportions_gr1$Stratif<-Clin$Stratif[match(Proportions_gr1$Patient,Clin$Patient_ID)]
colnames(Proportions_gr1)[2]<-"Cell_Type"

strat<-"Cell_Type_strat_2"
Proportions_gr2<-prop.table(table(meta_CD45[,which(colnames(meta_CD45) %in% c("Patient",strat))]),1)
Proportions_gr2<-melt(Proportions_gr2)
Proportions_gr2$Stratif<-Clin$Stratif[match(Proportions_gr2$Patient,Clin$Patient_ID)]
colnames(Proportions_gr2)[2]<-"Cell_Type"

Proportions<-rbind(Proportions_gr1,Proportions_gr2)

for(i in 1:length(cell_types)){
  for(j in 1:length(cell_types)){
    if(cell_types[i] != cell_types[j]){

      tmp1<-Proportions[which(Proportions$Cell_Type==cell_types[i]),]
      tmp2<-Proportions[which(Proportions$Cell_Type==cell_types[j]),]
      tmp2<-tmp2[match(tmp1$Patient,tmp2$Patient),]

      #tmp1$Ratio<-tmp1$value/tmp2$value
      tmp1$Ratio<-log2(tmp1$value/tmp2$value)

      if(length(which(tmp1$Ratio==Inf))>0){tmp1$Ratio[which(tmp1$Ratio==Inf)]<-NA}
      if(length(which(tmp1$Ratio==(-Inf)))>0){tmp1$Ratio[which(tmp1$Ratio==(-Inf))]<-NA}

      Imm1<-tmp1$Ratio[which(tmp1$Stratif=="up")];names(Imm1)<-as.character(tmp1$Patient[which(tmp1$Stratif=="up")])
      Imm2<-tmp1$Ratio[which(tmp1$Stratif=="dn")];names(Imm2)<-as.character(tmp1$Patient[which(tmp1$Stratif=="dn")])
      t<-makeStatisticalTest(Imm1[!is.na(Imm1)],Imm2[!is.na(Imm2)],statistics)
      collect_stats[cell_types[j],cell_types[i]] <- ifelse(t$statistic>0,(-log10(t$p.value)),-(-log10(t$p.value)))

      #tmp1$Ratio<-tmp2$value/tmp1$value
      tmp1$Ratio<-log2(tmp2$value/tmp1$value)

      if(length(which(tmp1$Ratio==Inf))>0){tmp1$Ratio[which(tmp1$Ratio==Inf)]<-NA}
      if(length(which(tmp1$Ratio==(-Inf)))>0){tmp1$Ratio[which(tmp1$Ratio==(-Inf))]<-NA}

      Imm1<-tmp1$Ratio[which(tmp1$Stratif=="up")];names(Imm1)<-as.character(tmp1$Patient[which(tmp1$Stratif=="up")])
      Imm2<-tmp1$Ratio[which(tmp1$Stratif=="dn")];names(Imm2)<-as.character(tmp1$Patient[which(tmp1$Stratif=="dn")])
      t<-makeStatisticalTest(Imm1[!is.na(Imm1)],Imm2[!is.na(Imm2)],statistics)
      collect_stats[cell_types[i],cell_types[j]] <- ifelse(t$statistic>0,(-log10(t$p.value)),-(-log10(t$p.value)))

    }

  }
}

########################################################################

# Out of All cells - 3P
########################################################################

Proportions<-prop.table(table(meta_all[,which(colnames(meta_all) %in% c("Patient",column_cell_type))]),1)
Proportions<-melt(Proportions)
Proportions$Stratif<-Clin$Stratif[match(Proportions$Patient,Clin$Patient_ID)]
colnames(Proportions)[2]<-"Cell_Type"

for(i in 1:length(cell_subtypes)){
  tmp1<-Proportions[which(Proportions$Cell_Type==cell_subtypes[i]),]
  if(nrow(tmp1)>0){
    Imm1<-tmp1$value[which(tmp1$Stratif=="up")];names(Imm1)<-as.character(tmp1$Patient[which(tmp1$Stratif=="up")])
    Imm2<-tmp1$value[which(tmp1$Stratif=="dn")];names(Imm2)<-as.character(tmp1$Patient[which(tmp1$Stratif=="dn")])
    t<-makeStatisticalTest(Imm1[!is.na(Imm1)],Imm2[!is.na(Imm2)],statistics)
    collect_stats["All",cell_subtypes[i]]<-ifelse(t$statistic>0,(-log10(t$p.value)),-(-log10(t$p.value)))
  }
}

strats<-c("Cell_Type_strat_1","Cell_Type_strat_2")
for(strat in strats){
  Proportions<-prop.table(table(meta_all[,which(colnames(meta_all) %in% c("Patient",strat))]),1)
  Proportions<-melt(Proportions)
  Proportions$Stratif<-Clin$Stratif[match(Proportions$Patient,Clin$Patient_ID)]
  colnames(Proportions)[2]<-"Cell_Type"

  for(ct in as.character(unique(Proportions$Cell_Type))){
    tmp1<-Proportions[which(Proportions$Cell_Type==ct),]
    Imm1<-tmp1$value[which(tmp1$Stratif=="up")];names(Imm1)<-as.character(tmp1$Patient[which(tmp1$Stratif=="up")])
    Imm2<-tmp1$value[which(tmp1$Stratif=="dn")];names(Imm2)<-as.character(tmp1$Patient[which(tmp1$Stratif=="dn")])
    t<-makeStatisticalTest(Imm1[!is.na(Imm1)],Imm2[!is.na(Imm2)],statistics)
    collect_stats["All",ct]<-ifelse(t$statistic>0,(-log10(t$p.value)),-(-log10(t$p.value)))
  }

}
########################################################################

# Ratio of cell types involving non-CD45 cell types
########################################################################

runs<-list(list(Malignant,CD8_T_cells),
           list(Malignant,CD4_T_cells),
           list(Malignant,NK_tgd),
           list(Malignant,DC_cells),
           list(Malignant,Macro_cells),
           list(Malignant,Myelo_non_Macro),
           list(Malignant,B_cells),
           list(Malignant,Stromal),
           list(CD8_T_cells,Malignant),
           list(CD8_T_cells,Stromal),
           list(NK_tgd,Malignant),
           list(NK_tgd,Stromal),
           list(CD4_T_cells,Malignant),
           list(CD4_T_cells,Stromal),
           list(DC_cells,Malignant),
           list(DC_cells,Stromal),
           list(Macro_cells,Malignant),
           list(Macro_cells,Stromal),
           list(B_cells,Malignant),
           list(B_cells,Stromal),
           list(Myelo_non_Macro,Malignant),
           list(Myelo_non_Macro,CD8_T_cells),
           list(Myelo_non_Macro,CD4_T_cells),
           list(Myelo_non_Macro,NK_tgd),
           list(Myelo_non_Macro,DC_cells),
           list(Myelo_non_Macro,Macro_cells),
           list(Myelo_non_Macro,B_cells),
           list(Myelo_non_Macro,Stromal),
           list(Stromal,Malignant),
           list(Stromal,CD8_T_cells),
           list(Stromal,CD4_T_cells),
           list(Stromal,NK_tgd),
           list(Stromal,DC_cells),
           list(Stromal,Macro_cells),
           list(Stromal,B_cells),
           list(Stromal,Myelo_non_Macro))

Proportions<-prop.table(table(meta_all[,which(colnames(meta_all) %in% c("Patient",column_cell_type))]),1)
Proportions<-melt(Proportions)
Proportions$Stratif<-Clin$Stratif[match(Proportions$Patient,Clin$Patient_ID)]
colnames(Proportions)[2]<-"Cell_Type"


for(r in 1:length(runs)){

  ct1<-runs[[r]][[1]]
  ct2<-runs[[r]][[2]]

  for(i in 1:length(ct1)){
    for(j in 1:length(ct2)){

      if((ct1[i] != ct2[j])&(ct1[i] %in% Proportions$Cell_Type)&(ct2[j] %in% Proportions$Cell_Type)){
        tmp1<-Proportions[which(Proportions$Cell_Type==ct1[i]),]
        #tmp1<-tmp1[complete.cases(tmp1),]

        tmp2<-Proportions[which(Proportions$Cell_Type==ct2[j]),]
        tmp2<-tmp2[match(tmp1$Patient,tmp2$Patient),]
        #tmp2<-tmp2[complete.cases(tmp2),]

        #tmp1$Ratio<-tmp1$value/tmp2$value
        tmp1$Ratio<-log2(tmp1$value/tmp2$value)

        if(length(which(tmp1$Ratio==(-Inf)))>0){tmp1$Ratio[which(tmp1$Ratio==(-Inf))]<-NA}
        if(length(which(tmp1$Ratio==(Inf)))>0){tmp1$Ratio[which(tmp1$Ratio==(Inf))]<-NA}

        Imm1<-tmp1$Ratio[which(tmp1$Stratif=="up")];names(Imm1)<-as.character(tmp1$Patient[which(tmp1$Stratif=="up")]);Imm1<-Imm1[complete.cases(Imm1)]
        Imm2<-tmp1$Ratio[which(tmp1$Stratif=="dn")];names(Imm2)<-as.character(tmp1$Patient[which(tmp1$Stratif=="dn")]);Imm2<-Imm2[complete.cases(Imm2)]

        if((length(Imm1)>0)&(length(Imm2)>0)){
          t<-makeStatisticalTest(Imm1[!is.na(Imm1)],Imm2[!is.na(Imm2)],statistics)
          collect_stats[ct2[j],ct1[i]] <- ifelse(t$statistic>0,(-log10(t$p.value)),-(-log10(t$p.value)))
        }

      }
    }
  }


}

# Now for the grouped
strats<-c("Cell_Type_strat_1","Cell_Type_strat_2")
for(strat in strats){
  Proportions_grouped<-prop.table(table(meta_all[,which(colnames(meta_all) %in% c("Patient",strat))]),1)
  Proportions_grouped<-melt(Proportions_grouped)
  Proportions_grouped$Stratif<-Clin$Stratif[match(Proportions_grouped$Patient,Clin$Patient_ID)]
  colnames(Proportions_grouped)[2]<-"Cell_Type"

  runs<-list(list("Malignant",Malignant),
             list("Malignant",Stromal),
             list("Stromal",Malignant),
             list("Stromal",Stromal),
             list("CD8",Malignant),
             list("CD8",Stromal),
             list("CD4",Malignant),
             list("CD4",Stromal),
             list("DC",Malignant),
             list("DC",Stromal),
             list("Macrophages",Malignant),
             list("Macrophages",Stromal),
             list("T_cells",Malignant),
             list("T_cells",Stromal),
             list("Myeloid",Malignant),
             list("Myeloid",Stromal),
             list("B_cells",Malignant),
             list("B_cells",Stromal))

  keep<-c()
  for(k in 1:length(runs)){
    keep<-c(keep,ifelse(runs[[k]][[1]] %in% unique(Proportions_grouped$Cell_Type),k,NA))
  }
  keep<-keep[complete.cases(keep)]
  runs<-runs[keep]

  for(r in 1:length(runs)){
    ct1<-runs[[r]][[1]]
    ct2<-runs[[r]][[2]]
    for(j in 1:length(ct2)){
      tmp1<-Proportions_grouped[which(Proportions_grouped$Cell_Type==ct1),]
      tmp2<-Proportions[which(Proportions$Cell_Type==ct2[j]),]
      tmp2<-tmp2[match(tmp1$Patient,tmp2$Patient),]

      #tmp1$Ratio<-tmp1$value/tmp2$value
      tmp1$Ratio<-log2(tmp1$value/tmp2$value)

      if(length(which(tmp1$Ratio==Inf))>0){tmp1$Ratio[which(tmp1$Ratio==Inf)]<-NA}
      if(length(which(tmp1$Ratio==(-Inf)))>0){tmp1$Ratio[which(tmp1$Ratio==(-Inf))]<-NA}

      Imm1<-tmp1$Ratio[which(tmp1$Stratif=="up")];names(Imm1)<-as.character(tmp1$Patient[which(tmp1$Stratif=="up")])
      Imm2<-tmp1$Ratio[which(tmp1$Stratif=="dn")];names(Imm2)<-as.character(tmp1$Patient[which(tmp1$Stratif=="dn")])
      t<-makeStatisticalTest(Imm1[!is.na(Imm1)],Imm2[!is.na(Imm2)],statistics)

      collect_stats[ct2[j],ct1] <- ifelse(t$statistic>0,(-log10(t$p.value)),-(-log10(t$p.value)))
    }
    for(j in 1:length(ct2)){
      tmp1<-Proportions_grouped[which(Proportions_grouped$Cell_Type==ct1),]
      tmp2<-Proportions[which(Proportions$Cell_Type==ct2[j]),]
      tmp2<-tmp2[match(tmp1$Patient,tmp2$Patient),]

      #tmp1$Ratio<-tmp2$value/tmp1$value
      tmp1$Ratio<-log2(tmp2$value/tmp1$value)

      if(length(which(tmp1$Ratio==Inf))>0){tmp1$Ratio[which(tmp1$Ratio==Inf)]<-NA}
      if(length(which(tmp1$Ratio==(-Inf)))>0){tmp1$Ratio[which(tmp1$Ratio==(-Inf))]<-NA}

      Imm1<-tmp1$Ratio[which(tmp1$Stratif=="up")];names(Imm1)<-as.character(tmp1$Patient[which(tmp1$Stratif=="up")])
      Imm2<-tmp1$Ratio[which(tmp1$Stratif=="dn")];names(Imm2)<-as.character(tmp1$Patient[which(tmp1$Stratif=="dn")])
      t<-makeStatisticalTest(Imm1[!is.na(Imm1)],Imm2[!is.na(Imm2)],statistics)

      collect_stats[ct1,ct2[j]] <- ifelse(t$statistic>0,(-log10(t$p.value)),-(-log10(t$p.value)))
    }
  }



}


# Grouped between them
strat<-"Cell_Type_strat_1"
Proportions_gr1<-prop.table(table(meta_all[,which(colnames(meta_all) %in% c("Patient",strat))]),1)
Proportions_gr1<-melt(Proportions_gr1)
Proportions_gr1$Stratif<-Clin$Stratif[match(Proportions_gr1$Patient,Clin$Patient_ID)]
colnames(Proportions_gr1)[2]<-"Cell_Type"

strat<-"Cell_Type_strat_2"
Proportions_gr2<-prop.table(table(meta_all[,which(colnames(meta_all) %in% c("Patient",strat))]),1)
Proportions_gr2<-melt(Proportions_gr2)
Proportions_gr2$Stratif<-Clin$Stratif[match(Proportions_gr2$Patient,Clin$Patient_ID)]
colnames(Proportions_gr2)[2]<-"Cell_Type"

Proportions<-rbind(Proportions_gr1,Proportions_gr2)

for(i in 1:length(cell_types)){
  for(j in 1:length(cell_types)){
    if(cell_types[i] != cell_types[j]){

      tmp1<-Proportions[which(Proportions$Cell_Type==cell_types[i]),]
      tmp2<-Proportions[which(Proportions$Cell_Type==cell_types[j]),]
      tmp2<-tmp2[match(tmp1$Patient,tmp2$Patient),]

      #tmp1$Ratio<-tmp1$value/tmp2$value
      tmp1$Ratio<-log2(tmp1$value/tmp2$value)

      if(length(which(tmp1$Ratio==Inf))>0){tmp1$Ratio[which(tmp1$Ratio==Inf)]<-NA}
      if(length(which(tmp1$Ratio==(-Inf)))>0){tmp1$Ratio[which(tmp1$Ratio==(-Inf))]<-NA}

      Imm1<-tmp1$Ratio[which(tmp1$Stratif=="up")];names(Imm1)<-as.character(tmp1$Patient[which(tmp1$Stratif=="up")])
      Imm2<-tmp1$Ratio[which(tmp1$Stratif=="dn")];names(Imm2)<-as.character(tmp1$Patient[which(tmp1$Stratif=="dn")])
      t<-makeStatisticalTest(Imm1[!is.na(Imm1)],Imm2[!is.na(Imm2)],statistics)
      collect_stats[cell_types[j],cell_types[i]] <- ifelse(t$statistic>0,(-log10(t$p.value)),-(-log10(t$p.value)))

      #tmp1$Ratio<-tmp2$value/tmp1$value
      tmp1$Ratio<-log2(tmp2$value/tmp1$value)

      if(length(which(tmp1$Ratio==Inf))>0){tmp1$Ratio[which(tmp1$Ratio==Inf)]<-NA}
      if(length(which(tmp1$Ratio==(-Inf)))>0){tmp1$Ratio[which(tmp1$Ratio==(-Inf))]<-NA}

      Imm1<-tmp1$Ratio[which(tmp1$Stratif=="up")];names(Imm1)<-as.character(tmp1$Patient[which(tmp1$Stratif=="up")])
      Imm2<-tmp1$Ratio[which(tmp1$Stratif=="dn")];names(Imm2)<-as.character(tmp1$Patient[which(tmp1$Stratif=="dn")])
      t<-makeStatisticalTest(Imm1[!is.na(Imm1)],Imm2[!is.na(Imm2)],statistics)
      collect_stats[cell_types[i],cell_types[j]] <- ifelse(t$statistic>0,(-log10(t$p.value)),-(-log10(t$p.value)))

    }

  }
}

########################################################################

signif<-abs(collect_stats)
signif_annot<-signif
signif_annot[]<-lapply(signif_annot,as.character)

for(i in 1:ncol(signif_annot)){
  signif_annot[,i] <-" "
  signif_annot[which(signif[,i] > (-log10(0.05))),i]<-"*"
  signif_annot[which(signif[,i] > (-log10(0.01))),i]<-"**"
  signif_annot[which(signif[,i] > (-log10(0.001))),i]<-"***"
}

lim<-max(abs(collect_stats),na.rm=T)
if(!is.null(breaks)){
  lim<-breaks
}
col_fun = colorRamp2(c(-lim, 0, lim), c("#3A34DF","white","#DF3434"))
gaps_row<-c(rep("A",4),rep("B",length(Malignant)),rep("C",length(CD8_T_cells)),
            rep("D",length(NK_tgd)),rep("E",length(CD4_T_cells)),
            rep("F",length(DC_cells)),rep("G",length(Macro_cells)),rep("H",length(Myelo_non_Macro)),
            rep("I",length(B_cells)),
            rep("J",length(Stromal)),
            rep("K",3),rep("L",3),rep("M",2))
gaps_col<-gaps_row[-c(1:4)]

if(length(which(colnames(collect_stats)=="others"))>0){
  collect_stats<-collect_stats[,-which(colnames(collect_stats) == "others")]
}

rownames(collect_stats)[1:4]<-paste0("Out of ",rownames(collect_stats)[1:4])
cat<-c("Mal.","CD8 T cells","NK/Tgd","CD4 T cells","DC","Macrophages","Gran.","B cells","Strom.","T cells","Myelo.","Others")

heatmap<-Heatmap(as.matrix(collect_stats), name = "mat", na_col = "lightgray",
        col=col_fun,cluster_rows=F,
        #row_order = 1:nrow(collect_stats),
        column_order = 1:ncol(collect_stats),
        column_title = "Cell Type Proportion and Ratio",
        column_title_gp = gpar(fontsize = 15, fontface = "bold"),
        row_names_side = "left",column_names_side = "top",
        row_names_gp = gpar(fontsize = 8),column_names_gp = gpar(fontsize = 8),
        row_labels=gsub("_"," ",rownames(collect_stats)),
        column_labels=gsub("_"," ",colnames(collect_stats)),
        row_split=gaps_row,row_gap = unit(c(3,rep(1,11)), "mm"),
        column_split=gaps_col,column_gap = unit(c(rep(1,11)), "mm"),
        border = TRUE,
        cell_fun = function(j, i, x, y, w, h, col) {
          grid.text(signif_annot[i, j], x, y,gp = gpar(fontsize = 8))
        },
        row_title = c("Proportion",rep(" ",4),"Ratio between cell types",rep(" ",7)),
        row_title_gp = gpar(fontsize = 12),
        width = unit(14, "cm"), height = unit(16, "cm"),
        top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = "white"),
                                                            labels = cat,
                                                            labels_gp = gpar(col = "black", fontsize = 5))))

return(heatmap)
}


#' @export
CellType_Proportion_Boxplot <- function(Cohort_selection = c("ATATIL","NeoTIL","SOLTIL"),
                                        Tumor_Type_selection = NULL,
                                        Timepoint = c("T0"),
                                        Stratification,
                                        groups = list(),
                                        type = c("Proportion","Ratio"),
                                        cell_type = NULL,
                                        out_of = NULL,
                                        cell_type_ratio = NULL,
                                        col_var,
                                        excluded_patient = NULL,
                                        plot = T){

  require(ComplexHeatmap)
  require(limma)
  require(reshape2)
  require(crayon)
  require(circlize)
  require(usethis)
  require(scales)

  data("precomputed_celltype_frequencies", package = "CellProportionACT")

  # Define parameters to plot
  if(is.na(match(type,c("Proportion","Ratio")))){
    warning("The selected type should be either `Proportion` or `Ratio`")
    return()
  }
  if(length(which(Timepoint  %in% c("T0","T30")==F))>0){
    warning("The selected time-points should be `T0` (default; baseline), and/or `T30` (post-ACT)")
    return()
  }

  # Select Correct time-points
  if(!is.null(Timepoint)){
    Collect_Proportion<-Collect_Proportion[which(Collect_Proportion$Timepoint %in% Timepoint),]
  }

  if(type=="Proportion"){
    if(is.na(match(out_of,c("all_cells","CD45","T_Myelo_B","CD8_CD4_DC_Macro_B")))){
      warning("The selected value for the `out_of` parameter should be `all_cells`, `CD45`, `T_Myelo_B` or `CD8_CD4_DC_Macro_B`")
      return()
    }

    idx<-grep(paste0("^out_of_:",out_of,"$"),Collect_Proportion$Type)
    Collect_Proportion<-Collect_Proportion[idx,]

    if(is.na(match(cell_type,unique(Collect_Proportion$Cell_Type)))){
      warning(paste0("The selected value for the `cell_type` parameter should be one of :",paste(unique(Collect_Proportion$Cell_Type),collapse=", ")))
      return()
    }

    idx<-which(Collect_Proportion$Cell_Type %in% cell_type)
    Collect_Proportion<-Collect_Proportion[idx,]

    y_value<-paste0("Proportion out of ",out_of)
    title<-paste0(cell_type," out of ",out_of)
  }

  if(type=="Ratio"){

    idx<-grep("ratio",Collect_Proportion$Type)
    Collect_Proportion<-Collect_Proportion[idx,]
    Collect_Proportion$Type<-gsub("ratio:","",Collect_Proportion$Type)
    Collect_Proportion$Cell_Type1<-strsplit2(Collect_Proportion$Type,split="[:::]")[,1]
    Collect_Proportion$Cell_Type2<-strsplit2(Collect_Proportion$Type,split="[:::]")[,2]

    if(is.na(match(cell_type,unique(Collect_Proportion$Cell_Type1)))){
      warning(paste0("The selected value for the `cell_type` parameter should be one of :",paste(unique(Collect_Proportion$Cell_Type1),collapse=", ")))
      return()
    }

    if(is.na(match(cell_type_ratio,unique(Collect_Proportion$Cell_Type2)))){
      warning(paste0("The selected value for the `cell_type_ratio` parameter should be one of :",paste(unique(Collect_Proportion$Cell_Type2),collapse=", ")))
      return()
    }
    idx1<-which(Collect_Proportion$Cell_Type1 %in% cell_type)
    idx2<-which(Collect_Proportion$Cell_Type2 %in% cell_type_ratio)
    idx<-intersect(idx1,idx2)
    Collect_Proportion<-Collect_Proportion[idx,]

    y_value<-paste0("Ratio ",cell_type," / ",cell_type_ratio)
    title<-paste0("Ratio ",cell_type," / ",cell_type_ratio)
  }

  if(length(which(Cohort_selection  %in% c("ATATIL","NeoTIL","SOLTIL")==F))>0){
    warning("The selected cohorts should be `ATATIL`, `NeoTIL` and/or `SOLTIL`")
    return()
  }

  # Select only Patient in the selected cohorts
  Collect_Proportion$Cohort<-strsplit2(Collect_Proportion$Patient,split="_")[,2]
  Collect_Proportion<-Collect_Proportion[which(Collect_Proportion$Cohort %in% Cohort_selection),]

  idx<-match(Stratification,colnames(Clin))
  if(is.na(idx)){
    names<-colnames(Clin)[which(sapply(Clin, is.character)==T)]
    warning(paste0("Stratification factor not found in clinical data; available ones are: ",paste(names[-which(names %in% c("Patient_ID","ATATIL_ID","Use_all_cells","Use_CD45_cells"))],collapse=", ")
    ))
    return()
  }

  if(!is.na(idx)){
    Clin$Stratif<-Clin[,Stratification]

    if(length(which(unlist(groups)  %in% unique(Clin$Stratif)==F))>0){
      warning(paste0("The selected groups should be in : ",paste(unique(Clin$Stratif),collapse=", ")))
      return()
    }

  }

  # Select Cancer Type
  if(!is.null(Tumor_Type_selection)){

    if(length(which(Tumor_Type_selection  %in% Clin$Cancer_Type ==F))>0){
      warning("The selected cancer types should be `Melanoma`, `Lung` and/or `Other_Solid`")
      return()
    }

    patients <- Clin$Patient_ID[which(Clin$Cancer_Type %in% Tumor_Type_selection)]
    Clin<-Clin[which(Clin$Patient_ID %in% patients),]
  }

  # remove patient to exclude
  if(!is.null(excluded_patient)){
    idx<-as.vector(unlist(sapply(excluded_patient,function(x){grep(x,Clin$Patient_ID)})))
    if(length(idx)==0){
      warning(paste0("In the selected cohort, the patients that are available for exclusion are :",
                     paste(apply(strsplit2(Clin$Patient_ID,split="_")[,c(1,2)],1,function(x){paste(x,collapse="_")}),collapse=", ")))
      return()
    }

    if(length(idx)>0){
      Clin<-Clin[-idx,]
    }
  }

  # Remove potential duplicates
  Collect_Proportion$Patient_Cohort<-apply(strsplit2(Collect_Proportion$Patient,split="_")[,c(1,2)],1,function(x){paste(x,collapse="_")})
  if(length(which(duplicated(Collect_Proportion$Patient_Cohort)))>0){
    Collect_Proportion<-Collect_Proportion[-which(duplicated(Collect_Proportion$Patient_Cohort)),]
  }

  # Start collecting values for boxplots
  idx<-match(Collect_Proportion$Patient,Clin$Patient_ID)
  Collect_Proportion$Stratif<-Clin$Stratif[idx]
  Collect_Proportion$pch<-Clin$pch[idx]
  Collect_Proportion$col<-Clin$col[idx]
  Collect_Proportion<-Collect_Proportion[complete.cases(Collect_Proportion),]
  Collect_Proportion$x<-NA

  list_Imm<-list()
  max<-c()
  min<-c()
  for(j in 1:length(groups)){
    idx<-which(Collect_Proportion$Patient %in% Clin$Patient_ID[which(Clin$Stratif %in% groups[[j]])])
    tmp<-Collect_Proportion$value[idx]
    assign(paste0("Imm",j),tmp)
    list_Imm[[j]]<-tmp[!is.na(tmp)]
    max<-c(max,tmp[!is.na(tmp)])
    min<-c(min,tmp[!is.na(tmp)])
    Collect_Proportion$x[idx]<-j
  }
  max<-max(max)
  min<-min(min)
  range<-(max-min)

  Collect_Proportion<-Collect_Proportion[complete.cases(Collect_Proportion),]

  if(plot==T){
    b<-boxplot(list_Imm,at=1:length(groups),col=col_var,
               names=names(groups),main=title,ylab=y_value, xaxt = "n",whisklty=1,lwd=1.5,
               ylim=c(min,max+(0.3*range)), outline = FALSE,width=rep(0.3,length(groups)))
    axis(side = 1, at = seq_along(b$names), labels = b$names, tick = FALSE)

    for(i in 1:nrow(Collect_Proportion)){
      stripchart(Collect_Proportion$value[i],
                 vertical = TRUE, method = "jitter",jitter=0.2,at=Collect_Proportion$x[i],
                 pch = Collect_Proportion$pch[i],
                 bg=alpha(Collect_Proportion$col[i],1),col = "black",
                 add = TRUE,cex=1.4)
    }

    if(length(groups)==2){
      used_test<-"wilcox.test"
      res<-do.call(used_test,list(Imm1,Imm2))
      text(1.2,max+(0.05*range),labels = paste0("p-val: ",signif(res$p.value,digits = 2)),cex=0.8,pos=4)
    }


    if(length(groups)==3){
      adjusted<-F
      pval<-as.data.frame(pairwise.wilcox.test(Collect_Proportion$value, Collect_Proportion$Stratif,p.adjust.method = ifelse(adjusted,"fdr","none"))$p.value)
      pval<-signif(pval,digit=2)

      text(2,max+(0.3*range),colnames(pval)[1],cex=0.8)
      text(3,max+(0.3*range),colnames(pval)[2],cex=0.8)
      text(1,max+(0.23*range),rownames(pval)[1],cex=0.8)
      text(1,max+(0.16*range),rownames(pval)[2],cex=0.8)

      text(2,max+(0.23*range),pval[1,1],cex=0.7)
      text(2,max+(0.16*range),pval[2,1],cex=0.7)
      text(3,max+(0.16*range),pval[2,2],cex=0.7)

    }

    if(length(groups)==4){
      adjusted<-F
      pval<-as.data.frame(pairwise.wilcox.test(Collect_Proportion$value, Collect_Proportion$Stratif,p.adjust.method = ifelse(adjusted,"fdr","none"))$p.value)
      pval<-signif(pval,digit=2)

      text(2,max+(0.3*range),colnames(pval)[1],cex=0.8)
      text(3,max+(0.3*range),colnames(pval)[2],cex=0.8)
      text(4,max+(0.3*range),colnames(pval)[3],cex=0.8)
      text(1,max+(0.23*range),rownames(pval)[1],cex=0.8)
      text(1,max+(0.16*range),rownames(pval)[2],cex=0.8)
      text(1,max+(0.09*range),rownames(pval)[3],cex=0.8)

      text(2,max+(0.23*range),pval[1,1],cex=0.7)
      text(2,max+(0.16*range),pval[2,1],cex=0.7)
      text(2,max+(0.09*range),pval[3,1],cex=0.7)
      text(3,max+(0.16*range),pval[2,2],cex=0.7)
      text(3,max+(0.09*range),pval[3,2],cex=0.7)
      text(4,max+(0.09*range),pval[3,3],cex=0.7)

    }
  }

  return(Collect_Proportion)
}
