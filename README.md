# CellProportionR

## Overview

**CellProportionR** is an R package designed for analyzing and visualizing cell type proportions in single-cell datasets. The package provides functions for computing cell proportions, generating heatmaps and boxplots, and performing statistical comparisons between different groups.

## Installation

To install **CellProportionR**, use the following command in R:

```r
# Install from GitHub
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
devtools::install_github("dbarras2/CellProportionR")
```

## Functions

### 1. CellType\_Proportion\_Heatmap

**Description:** Generation of a heatmap displaying statistical comparisons of Cell Type proportions and ratios between cell type proportions between two groups in a single-cell annotated object. 

#### Example after editing
<img width="730" alt="Heatmap_example" src="https://github.com/user-attachments/assets/b4f3a359-10c9-4cf4-9ca3-351fe030547a" />

#### Usage

```r
data("example_data", package = "CellProportionR")
Proportion_Heatmap <- CellType_Proportion_Heatmap(single_cell_data = single_cell_data,
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
Proportion_Heatmap$heatmap
pvalues <- Proportion_Heatmap$statistics
```

#### Arguments

- `single_cell_data`: A data.frame (typically the meta.data from a Seurat object) where rows are individual cells and columns contain informations about cell types, their stratifications, the sample origin and the sample information needed to perform the comparison and for potential subsetting of the data. The cell type information has to be encoded in columns with the `Cell_Type_Strat` pattern followed by a numeric (e.g. Cell_Type_Strat1, Cell_Type_Strat2, Cell_Type_Strat3 and Cell_Type_Strat4) where the lower stratification (Cell_Type_Strat1) encodes the less granular cell types (e.g. `immune` and `non_immune`) and the highest stratification (Cell_Type_Strat4) encodes for the highest granularity (e.g. `exhausted_T_cell`, `CD4_TRegs` or `Pericytes`). The highest granularity stratification has to be a factor with the levels encoding the order of the cell types as they will appear in the heatmap.
- `sample_colname`: A character indicating which column name of the single_cell_data object corresponds to the column with the sample identifiers.
- `subset_data`: if provided, should be a list indicating which data to subset/use (such as list(`Technology`=c("5prime"),`Patient`=c(`Patient1`,`Patient2`))) where each element is named with the same name than the column in `single_cell_data` and each elements is a vector containing the identifiers contained in these columns and that allow to select the cells to be analyzed.
- `stratification`: A character mentioning which stratification factor to use for splitting the data into groups. This character has to be one of the column name of the single_cell_data object.
- `group_up`: A vector or character value indicating which group to use in the selected stratification. This group will be compared against the `group_dn` group. The cell types enriched in this group will appear in red (unless other color is specific in `group_up_col`) in the heatmap.
- `group_dn`: A vector or character value indicating which group to use in the selected stratification. This group will be compared against the `group_up` group. The cell types enriched in this group will appear in blue (unless other color is specific in `group_dn_col`) in the heatmap.
- `group_up_col`: If provided, a color that will be used to highlight cell types enriched in the `group_up`.
- `group_dn_col`: If provided, a color that will be used to highlight cell types enriched in the `group_dn`.
- `statistics`: Which statistical test to use. Should be either "wilcoxon" or "ttest". Default set to wilcoxon.
- `include_ratio`: A boolean indicating if ratio between cell type proportions should be added to the heatmap. If set to FALSE, only cell type proportions will be shown. Default set to TRUE.
- `stratification_names`: A vector that can used to supply the names of the stratification so that they are more meaningful that `Cell_Type_Strat1`, `Cell_Type_Strat2` etc... . It can be something like c(`Immune_NotImmune`,`Tcell_Bcell_Malignant_Stromal`,`CD4_CD8_Malignant_CAFs_Endo`).
- `pval_range`: A numerical setting the maxima of the -log10 p-values in the color code. That will force the range of -log10 p-values (encoded by the color code) to be contain within this range. If set to 3, the color code will range between -3 and 3.
- `stratification_for_ratio`: If provided, should be a list indicating which data to subset/use (such as list(`Cell_Type_Strat1`=c("Immune")) ) where each element is named with the same name than the column in `single_cell_data` and each elements is a vector containing the identifiers contained in these columns and that allow to select the cells to be analyzed for comouting the ratio of between cell proportions.

#### Returns

a list containing a heatmap (`heatmap`) and the statistics p-value table used to make the heatmap (`statistics`). The output can be plotted directly by calling the function or alternatively saved into a variable and plot using the draw() function. The values in the heatmap are showing directional log10 of the p-values. The red color-code indicates enrichment in the `group_up` group while blue color indicates enrichment in the `group_dn` group. Significance is defined, without adjustment for multiple comparisons, by `*` = p-value < 0.05, `**` = p-value < 0.01 and `***` = p-value < 0.001.

---

### 2. Compute\_Proportions\_Ratios

**Description:** Computes the proportions and ratios of different cell types within a single-cell dataset.

#### Usage

```r
data("example_data", package = "CellProportionR")
Compute_Proportions_Ratios(single_cell_data = single_cell_data,
                           sample_colname = "Patient",
                           stratification = "Category",
                           group_up = "group1",
                           group_dn = "group2")
```

#### Arguments

- `single_cell_data`: A data.frame (typically the meta.data from a Seurat object) where rows are individual cells and columns contain informations about cell types, their stratifications, the sample origin and the sample information needed to perform the comparison and for potential subsetting of the data. The cell type information has to be encoded in columns with the `Cell_Type_Strat` pattern followed by a numeric (e.g. Cell_Type_Strat1, Cell_Type_Strat2, Cell_Type_Strat3 and Cell_Type_Strat4) where the lower stratification (Cell_Type_Strat1) encodes the less granular cell types (e.g. `immune` and `non_immune`) and the highest stratification (Cell_Type_Strat4) encodes for the highest granularity (e.g. `exhausted_T_cell`, `CD4_TRegs` or `Pericytes`).
- `sample_colname`: A character indicating which column name of the single_cell_data object corresponds to the column with the sample identifiers.
- `stratification`: A character mentioning which stratification factor to use for splitting the data into groups. This character has to be one of the column name of the single_cell_data object.
- `group_up`: A vector or character value indicating which group to use in the selected stratification. This group will be compared against the `group_dn` group. The cell types enriched in this group will appear in red (unless other color is specific in `group_up_col`) in the heatmap.
- `group_dn`: A vector or character value indicating which group to use in the selected stratification. This group will be compared against the `group_up` group. The cell types enriched in this group will appear in blue (unless other color is specific in `group_dn_col`) in the heatmap.
- `stratification_for_ratio`: If provided, should be a list indicating which data to subset/use (such as list(`Cell_Type_Strat1`=c("Immune")) ) where each element is named with the same name than the column in `single_cell_data` and each elements is a vector containing the identifiers contained in these columns and that allow to select the cells to be analyzed for comouting the ratio of between cell proportions.
  
#### Returns

A data.frame containing the values for all cell proportion and ratios per sample.

---

### 3. CellType\_Proportion\_Boxplot

**Description:** Generates a boxplot comparing cell type proportions or ratios across groups.

#### Usage

```r
data("example_data", package = "CellProportionR")
Proportion_object <- Compute_Proportions_Ratios(single_cell_data = single_cell_data,
                                                sample_colname = "Patient",
                                                stratification = "Category",
                                                group_up = "group1",
                                                group_dn = "group2")

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
```

#### Arguments

- `cell_proportion_object`: A data.frame as output by the `Compute_Proportions_Ratios` function.
- `subset_data`: If provided, should be a list indicating which data to subset/use (such as list(`Technology`=c("5prime"),`Patient`=c(`Patient1`,`Patient2`))) where each element is named with the same name than the column in the `cell_proportion_object` object. Important, don't forget to annotate the `cell_proportion_object` object if the subset information is not already present.
- `stratification`: A character mentioning which stratification factor to use for splitting the data into groups. This character has to be one of the column name of `cell_proportion_object`.
- `groups`: A list indicating all the groups that will be used to make the boxplots. Each group will represented by a box. Each element of the list can be a simple character or a vector indicating groups to form (e.g "Responders"=c("CR","PR")). The name of each element will be used as the name of the category for plotting.
- `type`: A character indicating whether to plot `Proportion` or `Ratio`.
- `cell_type`: A character indicating which cell type to report. If type="Proportion", then the cell_type indicates the proportion of this cell type will be computed relative to the indicated cell population (see `out_of` below). If type="Ratio", then the cell_type will be the cell type that will be divided by the second cell type indicated by the `cell_type_ratio` parameter (see below). Of note, bulk cell types (parental cell types haveing more than one daughter cell type) were being added the "_bulk" extension to their name to avoid any misinterpretation with more granular/daughter cell type bearing the same name.
- `out_of`: If type="Proportion", `out_of` should indicated to which cell population the frequency of cell_type should be computed. If the proportion of CD8_Pex out of T cells has to be computed, then cell_type="CD8_Pex" and out_of="T_cells". Any values that corresponds to any levels of stratification is allowed as long as it was computable. Of note, bulk cell types (parental cell types haveing more than one daughter cell type) were being added the "_bulk" extension to their name to avoid any misinterpretation with more granular/daughter cell type bearing the same name.
- `cell_type_ratio`: If type="Ratio", `cell_type_ratio` indicates the cell type to normalize to. So the ratio would be the ratio of cell_type divided by cell_type_ratio. Importantly, the log2 ratio is reported, so positive values indicates that the proportion of cell_type is higher than the proportion of cell_type_ratio. If the ratio between CD8_Pex and CD4_CXCL13 has to be computed, then cell_type="CD8_Pex" and cell_type_ratio="CD4_CXCL13".
- `col_var`: A vector indicating which colors should be used for the boxes. Should be the same length as the length of the `groups` parameter.

#### Returns

A boxplot comparing cell proportions or ratios across different groups.

---

## Dependencies

- `ComplexHeatmap`
- `limma`
- `reshape2`
- `dplyr`
- `ggplot2`

## Author

Developed by **David Barras**. Contributions and feedback are welcome!

