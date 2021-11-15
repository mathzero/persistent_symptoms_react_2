
#' First clear the environment of variables
rm(list=ls(all=TRUE))
setwd("E:/Group/react2_study5/report_phases_combined/projects/long_covid_resubmission")
outpath <- "E:/Group/react2_study5/report_phases_combined/projects/long_covid_resubmission/output/"
figpath <- "E:/Group/react2_study5/report_phases_combined/projects/long_covid_resubmission/plots/"

#' Source any functions from the local file
source("E:/Group/functions/model_table.R")
source("E:/Group/functions/load_packages.R")
source("E:/Group/functions/model_maker.R")
source("E:/Group/functions/cats_and_covs.R")
source("E:/Group/functions/save_pheatmap.R")
source("E:/Group/functions/hamming_distance.R")
source("E:/Group/functions/longCovidPrevTables.R")
source("E:/Group/functions/LassoSatabilitySelection.R")
source("E:/Group/functions/relevel_all_categorical_vars.R")

#' Pull in packages needed
package.list <- c("prevalence","mgcv","knitr","MASS","kableExtra","table1","dplyr","factoextra",
                  "glmnet","fastDummies",
                  "tidyr","forcats", "cluster", "fpc", "mclust", "pheatmap","FactoMineR", 
                  "NbClust","parallel","ICoLour",
                  "ggplot2","gdata","ggsci", "RColorBrewer", "tidyverse", "lubridate", 
                  "egg", "circlize","ggnetwork",
                  "poLCA","snow","randomForest","ranger", "jtools", "moreparty", 
                  "ComplexHeatmap","blockcluster",
                  "readr","ggthemes", "questionr", "gridExtra", "future.apply", 
                  "foreach","doParallel","catboost",
                  "tidytext", "quanteda", "widyr", "igraph", "ggraph", "patchwork", "OverReact")
load_packages(package.list)



# Run data prep  ----------------------------------------------------------

source("E:/Group/react2_study5/report_phases_combined/projects/long_covid_resubmission/code/00_functions.R")
source("E:/Group/react2_study5/report_phases_combined/projects/long_covid_resubmission/code/00_bits_and_pieces.R")
source("E:/Group/react2_study5/report_phases_combined/projects/functions/univariate_models.R")


# Load data ---------------------------------------------------------------

dfRes <- readRDS("E:/Group/react2_study5/saved_objects/rep_all_r123456_validated_newwrangle.rds")
#load data key
data_key=read_csv("E:/Group/react2_study5/report_phases_combined/projects/long_covid_2021/long_covid_data_key.csv")


# Get Ids to exclude
ids_varselection <- readRDS("E:/Group/react2_study5/saved_objects/long_covid_clustering/ids_varselection.rds")

# Load cluster data
# save all clustering data
cluster_dat <- readRDS("E:/Group/react2_study5/saved_objects/long_covid_clustering/all_clustering_data.rds")
createMySubfolder("clustering")


# Main analysis
ids_5 <- cluster_dat$r3_5_main_analysis$cluster_obj$cluster_results$pamobject$clustering %>% names()
clusters_5 <- cluster_dat$r3_5_main_analysis$cluster_obj$cluster_results$pamobject$clustering
roworder=cluster_dat$r3_5_main_analysis$roworder

dfRes$lc_84 <- dfRes$lc_84_all_29

# dfRes <- makeVarsFactors(dfRes = dfRes) # convert all character variables to factors for modelling
moddat_5 <- prepMyData(dat = dfRes %>% filter(round%in%3:5, u_passcode %in% ids_5), t0 = F, min_followup = NULL, 
                     abresult_filter = NULL, 
                     covidafilter = c(1,2,3),
                     forclustering = T,
                     lc_names_12_weeks = data_key$symptom_duration_84days_binary[data_key$one_of_29==1])


moddat_5_L1 <- moddat_5[clusters_5==1,]
moddat_5_L2 <- moddat_5[clusters_5==2,]

# Cooccurrence mat
#### T0
cooccurmat_L1 <-crossprod(as.matrix(moddat_5_L1))
marginal_counts_L1=diag(cooccurmat_L1)
diag(cooccurmat_L1) <- 0
rownames(cooccurmat_L1) <- colnames(cooccurmat_L1) <- data_key$short_name[data_key$one_of_29==1]
cooccurmat_L2 <-crossprod(as.matrix(moddat_5_L2))
marginal_counts_L2=diag(cooccurmat_L2)
diag(cooccurmat_L2) <- 0
rownames(cooccurmat_L2) <- colnames(cooccurmat_L2) <- data_key$short_name[data_key$one_of_29==1]
maxval <- max(max(cooccurmat_L1,na.rm=T),max(cooccurmat_L2,na.rm=T))
col_fun=colorRamp2(breaks = c(0,maxval/2,maxval), colors = c("white",myCols[[2]],myCols[[1]]))


### Create heatmap
my.complex.heatmap.cooccur.5.L1 <- Heatmap(cooccurmat_L1,
                                         name = "Symptom \ncooccurrence",
                                         row_gap = unit(5,"mm"),
                                         column_order = roworder,
                                         # row_split = t0_cluster_cols$cluster_results$pamobject$clustering,
                                         # column_split = t0_cluster_cols$cluster_results$pamobject$clustering,
                                         column_gap = unit(5,"mm"),
                                         border = T,
                                         left_annotation = rowAnnotation(Sum =
                                                                            anno_barplot(marginal_counts_L1, 
                                                                                         gp = gpar(fill = myCols[[1]]),
                                                                          border=F,
                                                                          axis_param = list(direction = "reverse"))),
                                         # row_split = 3,
                                         # column_split = 3,
                                         # split = 3,
                                         row_order = roworder,
                                         column_title = "Cluster L1",
                                         cluster_rows = F,
                                         cluster_columns = F,
                                         rect_gp = gpar(col = "white", lwd = 1),
                                         cell_fun = function(j,i,x,y,width, height, fill){
                                           grid.text(sprintf("%.0f", cooccurmat_L1[i,j]), x, y, gp=gpar(fontsize=6,
                                                                                                        col = 
                                                                                                          ifelse(cooccurmat_L1[i,j]>maxval/3 |
                                                                                                                   cooccurmat_L1[i,j] == 0,
                                                                                                                 "white",
                                                                                                                 "black")))
                                         },
                                         col = col_fun)

my.complex.heatmap.cooccur.5.L1


### L2 cluster ### 



### Create heatmap
my.complex.heatmap.cooccur.5.L2 <- Heatmap(cooccurmat_L2,
                                           name = "Symptom \ncooccurrence",
                                           row_gap = unit(5,"mm"),
                                           column_order = roworder,
                                           # row_split = t0_cluster_cols$cluster_results$pamobject$clustering,
                                           # column_split = t0_cluster_cols$cluster_results$pamobject$clustering,
                                           column_gap = unit(5,"mm"),
                                           border = T,
                                           left_annotation = rowAnnotation(Sum =
                                                                             anno_barplot(marginal_counts_L2, 
                                                                                          gp = gpar(fill = myCols[[1]]),
                                                                                          border=F,
                                                                                          axis_param = list(direction = "reverse"))),
                                           # row_split = 3,
                                           # column_split = 3,
                                           # split = 3,
                                           row_order = roworder,
                                           column_title = "Cluster L2",
                                           cluster_rows = F,
                                           cluster_columns = F,
                                           rect_gp = gpar(col = "white", lwd = 1),
                                           cell_fun = function(j,i,x,y,width, height, fill){
                                             grid.text(sprintf("%.0f", cooccurmat_L2[i,j]), x, y, gp=gpar(fontsize=6,
                                                                                                          col = 
                                                                                                            ifelse(cooccurmat_L2[i,j]>maxval/3 |
                                                                                                                     cooccurmat_L2[i,j] == 0,
                                                                                                                   "white",
                                                                                                                   "black")))
                                           },
                                           col = col_fun)

my.complex.heatmap.cooccur.5.L2




# COmbine plots
combined_heatmap <- draw((my.complex.heatmap.cooccur.5.L1+my.complex.heatmap.cooccur.5.L2), ht_gap =unit(1,"cm"))

combined_heatmap




png(paste0(figpath,"cluster_symptom_cooccurrence_r3_5.png"), width = 20,
    height = 10, units = "in", res = 300)
combined_heatmap
dev.off()





# Repeat for replication analysis -----------------------------------------

# dfRes <- makeVarsFactors(dfRes = dfRes) # convert all character variables to factors for modelling
moddat_5 <- prepMyData(dat = dfRes %>% filter(round%in%3:5, u_passcode %in% ids_5), t0 = F, min_followup = NULL, 
                       abresult_filter = NULL, 
                       covidafilter = c(1,2,3),
                       forclustering = T,
                       lc_names_12_weeks = data_key$symptom_duration_84days_binary[data_key$one_of_29==1])


moddat_5_L1 <- moddat_5[clusters_5==1,]
moddat_5_L2 <- moddat_5[clusters_5==2,]

# Cooccurrence mat
#### T0
cooccurmat_L1 <-crossprod(as.matrix(moddat_5_L1))
marginal_counts_L1=diag(cooccurmat_L1)
diag(cooccurmat_L1) <- 0
rownames(cooccurmat_L1) <- colnames(cooccurmat_L1) <- data_key$short_name[data_key$one_of_29==1]
cooccurmat_L2 <-crossprod(as.matrix(moddat_5_L2))
marginal_counts_L2=diag(cooccurmat_L2)
diag(cooccurmat_L2) <- 0
rownames(cooccurmat_L2) <- colnames(cooccurmat_L2) <- data_key$short_name[data_key$one_of_29==1]
maxval <- max(max(cooccurmat_L1,na.rm=T),max(cooccurmat_L2,na.rm=T))
col_fun=colorRamp2(breaks = c(0,maxval/2,maxval), colors = c("white",myCols[[2]],myCols[[1]]))


### Create heatmap
my.complex.heatmap.cooccur.5.L1 <- Heatmap(cooccurmat_L1,
                                           name = "Symptom \ncooccurrence",
                                           row_gap = unit(5,"mm"),
                                           column_order = roworder,
                                           # row_split = t0_cluster_cols$cluster_results$pamobject$clustering,
                                           # column_split = t0_cluster_cols$cluster_results$pamobject$clustering,
                                           column_gap = unit(5,"mm"),
                                           border = T,
                                           left_annotation = rowAnnotation(Sum =
                                                                             anno_barplot(marginal_counts_L1, 
                                                                                          gp = gpar(fill = myCols[[1]]),
                                                                                          border=F,
                                                                                          axis_param = list(direction = "reverse"))),
                                           # row_split = 3,
                                           # column_split = 3,
                                           # split = 3,
                                           row_order = roworder,
                                           column_title = "Cluster L1",
                                           cluster_rows = F,
                                           cluster_columns = F,
                                           rect_gp = gpar(col = "white", lwd = 1),
                                           cell_fun = function(j,i,x,y,width, height, fill){
                                             grid.text(sprintf("%.0f", cooccurmat_L1[i,j]), x, y, gp=gpar(fontsize=6,
                                                                                                          col = 
                                                                                                            ifelse(cooccurmat_L1[i,j]>maxval/3 |
                                                                                                                     cooccurmat_L1[i,j] == 0,
                                                                                                                   "white",
                                                                                                                   "black")))
                                           },
                                           col = col_fun)

my.complex.heatmap.cooccur.5.L1


### L2 cluster ### 



### Create heatmap
my.complex.heatmap.cooccur.5.L2 <- Heatmap(cooccurmat_L2,
                                           name = "Symptom \ncooccurrence",
                                           row_gap = unit(5,"mm"),
                                           column_order = roworder,
                                           # row_split = t0_cluster_cols$cluster_results$pamobject$clustering,
                                           # column_split = t0_cluster_cols$cluster_results$pamobject$clustering,
                                           column_gap = unit(5,"mm"),
                                           border = T,
                                           left_annotation = rowAnnotation(Sum =
                                                                             anno_barplot(marginal_counts_L2, 
                                                                                          gp = gpar(fill = myCols[[1]]),
                                                                                          border=F,
                                                                                          axis_param = list(direction = "reverse"))),
                                           # row_split = 3,
                                           # column_split = 3,
                                           # split = 3,
                                           row_order = roworder,
                                           column_title = "Cluster L2",
                                           cluster_rows = F,
                                           cluster_columns = F,
                                           rect_gp = gpar(col = "white", lwd = 1),
                                           cell_fun = function(j,i,x,y,width, height, fill){
                                             grid.text(sprintf("%.0f", cooccurmat_L2[i,j]), x, y, gp=gpar(fontsize=6,
                                                                                                          col = 
                                                                                                            ifelse(cooccurmat_L2[i,j]>maxval/3 |
                                                                                                                     cooccurmat_L2[i,j] == 0,
                                                                                                                   "white",
                                                                                                                   "black")))
                                           },
                                           col = col_fun)

my.complex.heatmap.cooccur.5.L2




# COmbine plots
combined_heatmap <- draw((my.complex.heatmap.cooccur.5.L1+my.complex.heatmap.cooccur.5.L2), ht_gap =unit(1,"cm"))

combined_heatmap




png(paste0(figpath,"cluster_symptom_cooccurrence_r3_5.png"), width = 20,
    height = 10, units = "in", res = 300)
combined_heatmap
dev.off()
