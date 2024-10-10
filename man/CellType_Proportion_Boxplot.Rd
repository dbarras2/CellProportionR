\name{CellType_Proportion_Boxplot}
\alias{CellType_Proportion_Boxplot}
\title{Creation of a boxplot displaying Cell Type proportions or ratio in ACT cohort}
\description{The function will plot the cell type frequencies or ratio in the indicated cohort(s) by splitting for the chosen clionical factor. The cell type proportion can be shown according to several levels of resolution ranging from the overall stoechiometry to their proportion within their cell lineage. The ratio between different cell types can also be reported. The function supports cancer type selection within the cohort and exclusion of specific patients.}
\usage{CellType_Proportion_Boxplot(Cohort_selection = c("ATATIL","NeoTIL","SOLTIL"),
                                        Tumor_Type_selection = NULL,
                                        Timepoint = c("T0"),
                                        Stratification,
                                        groups = list(),
                                        type = c("Proportion","Ratio"),
                                        cell_type = NULL,
                                        out_of = NULL,
                                        cell_type_ratio = NULL,
                                        col_var,
                                        excluded_patient = NULL)}

\arguments{
  \item{Cohort_selection}{a vector indicating which ACT cohort to interrogate. Accepted values are `ATATIL`, `NeoTIL` or `SOLTIL`}
  \item{Tumor_Type_selection}{if provided, should be a vector indicating which tumor type(s) to focus the analysis on. Accepted values are `Melanoma`, `Lung` and/or `Other_Solid` but depends on the cohort used. Please run the function with a random values to see the possibilities offered by the selected cohort}
  \item{Timepoint}{a character or vector indicating which tumoral timepoint to use. Accepted values are `T0` for baseline tumors and `T30` for 30 days post-ACT (Adoptive cell therapy)}
  \item{Stratification}{a character mentioning which clinical category to use for splitting the data into groups. Please run the function with a random values to see the possibilities offered by the selected cohort}
  \item{groups}{a list indicating all the groups that will be used to make the boxplots. Each group will represented by a box. Each element of the list can be a simple character or a vector indicasting groups to make (e.g "Rs"=c("CR","PR")). The name of each element will be used as the name of the category for plotting. Please run the function with a random values to see the possibilities offered by the selected cohort}
  \item{type}{a character indicating whether to plot `Proportion` or `Ratio`}
  \item{cell_type}{a character indicating which cell type to report. If type="Proportion", then the cell_type indicates the proportion of this cell type will be computed relative to the indicated cell population (see `out_of` below).  If type="Ratio", then the cell_type will be the cell type that will be divided by the second cell type indicated by the `cell_type_ratio` parameter (see below)}
  \item{out_of}{If type="Proportion", `out_of` should indicated to which cell population the frequency of cell_type should be computed. Allowed values are `all_cells`,`CD45`, `T_Myelo_B` or `CD8_CD4_DC_Macro_B`. If the proportion of CD8_Pex out of T cells has to be computed, then cell_type="CD8_Pex" and out_of="T_Myelo_B". If the percentage of CD4_CXCL13 out of CD4 T cells has to be computed, then cell_type="CD4_CXCL13" and out_of="CD8_CD4_DC_Macro_B"}
  \item{cell_type_ratio}{If type="Ratio", `cell_type_ratio` indicates the cell type to normalize to. So the ratio would be the ratio of cell_type divided by cell_type_ratio. Importantly, the log2 ratio is reported, so positive values indicates that the proportion of cell_type is higher than the proportion of cell_type_ratio. If the ratio between CD8_Pex and CD4_CXCL13 has to be computed, then cell_type="CD8_Pex" and cell_type_ratio="CD4_CXCL13"}
  \item{col_var}{a vector indicating which colors should be used for the boxes. Should be the same length as the length of the `groups` parameter}
  \item{excluded_patient}{a character or vector names of patient(s) to exclude from the analysis. Please run the function with a random values to see the possibilities offered by the selected cohort}
  \item{plot}{a logical indicating whether to plot the boxplot or not. if set to FALSE, the return value will be only the data.frame containing the data}
}
\details{}
\value{Running the function will plot the boxplot in the graphic device and will return as value, a data.frame containg the values used to make the boxplot. p-values in the plot are computed using non-adjusted Wilcoxon rank-sum test}
\references{https://www.biorxiv.org/content/10.1101/2022.12.23.519261v1}
\author{David Barras}
\note{}
\seealso{}
\examples{CellType_Proportion_Boxplot(Cohort_selection = c("ATATIL","NeoTIL"),
                                      Tumor_Type_selection = "Melanoma",
                                      Stratification = "Best_Overall_Response",
                                      groups = list("CR"="CR","PR"="PR","SD"="SD","PD"="PD"),
                                      type = "Proportion",cell_type = "Tgd",out_of = "T_Myelo_B",
                                      col_var = c("darkred","red","lightblue","darkgray"))}

