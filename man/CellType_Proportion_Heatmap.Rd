\name{CellType_Proportion_Heatmap}
\alias{CellType_Proportion_Heatmap}
\title{Creation of a heatmap displaying Cell Type proportions and ratio in ACT cohort}
\description{The function will interrogate the indicated cohort(s) for any potential cell type enrichment in the two groups provided. The cell type proportion are shown according to several levels of resolution ranging from the overall stoechiometry to their proportion within their cell lineage. The reported heatmap also shows ratio between cell type proportions. The function supports cancer type selection within the cohort and exclusion of specific patients.}
\usage{CellType_Proportion_Heatmap(Cohort_selection = c("ATATIL","NeoTIL","SOLTIL"),
                                        Tumor_Type_selection = NULL,
                                        Timepoint = c("T0"),
                                        Stratification,
                                        group_up,
                                        group_dn,
                                        statistics = "wilcoxon",
                                        excluded_patient = NULL)}

\arguments{
  \item{Cohort_selection}{a vector indicating which ACT cohort to interrogate. Accepted values are `ATATIL`, `NeoTIL` or `SOLTIL`}
  \item{Tumor_Type_selection}{if provided, should be a vector indicating which tumor type(s) to focus the analysis on. Accepted values are `Melanoma`, `Lung` and/or `Other_Solid` but depends on the cohort used. Please run the function with a random values to see the possibilities offered by the selected cohort}
    \item{Timepoint}{a character or vector indicating which tumoral timepoint to use. Accepted values are `T0` for baseline tumors and `T30` for 30 days post-ACT (Adoptive cell therapy)}
  \item{Stratification}{a character mentioning which clinical category to use for splitting the data into groups. Please run the function with a random values to see the possibilities offered by the selected cohort}
  \item{group_up}{a vector or character value indicating which group to use in the selected stratification. This group will be compared against the `group_dn` group The cell types enriched in this group will appear in red in the heatmap}
  \item{group_dn}{a vector or character value indicating which group to use in the selected stratification. This group will be compared against the `group_up` group The cell types enriched in this group will appear in blue in the heatmap}
  \item{statistics}{which statistical test to use. Should be either "wilcoxon" or "ttest"}
  \item{excluded_patient}{a character or vector names of patient(s) to exclude from the analysis. Please run the function with a random values to see the possibilities offered by the selected cohort}
}
\details{}
\value{a heatmap from the complexHeatmap package. The output can be plotted directly by calling the function or alternatively saved into a variable and plot using the draw() function. The values in the heatmap are showing directional log10 of the p-values. The red color-code indicates enrichment in the `group_up` group while blue color indicates enrichment in the `group_dn` group. Significance is defined, without adjustment for multiple comparisons, by `*` = p-value < 0.05, `**` = p-value < 0.01 and `***` = p-value < 0.001.}
\references{https://www.biorxiv.org/content/10.1101/2022.12.23.519261v1}
\author{David Barras}
\note{}
\seealso{}
\examples{CellType_Proportion_Heatmap(Cohort_selection = c("ATATIL","NeoTIL"),
                                      Tumor_Type_selection = NULL,
                                      Stratification = "Best_Overall_Response",
                                      group_up = c("CR","PR"), group_dn = c("SD","PD"),
                                      statistics = "wilcoxon",excluded_patient = NULL)}

