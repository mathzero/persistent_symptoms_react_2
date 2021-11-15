##### ANalysing Long COVID Clusters ####

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
                  "NbClust","parallel",
                  "ggplot2","gdata","ggsci", "RColorBrewer", "tidyverse", "lubridate", 
                  "egg", "circlize",
                  "poLCA","snow","randomForest","ranger", "jtools", "moreparty", 
                  "ComplexHeatmap","blockcluster",
                  "readr","ggthemes", "questionr", "gridExtra", "future.apply", 
                  "foreach","doParallel",
                  "tidytext", "quanteda", "widyr", "igraph", "ggraph", "patchwork", "OverReact")
load_packages(package.list)



# Run data prep  ----------------------------------------------------------

source("E:/Group/react2_study5/report_phases_combined/projects/long_covid_resubmission/code/00_functions.R")
# source("E:/Group/react2_study5/report_phases_combined/projects/functions/univariate_models.R")

# Load data ---------------------------------------------------------------

dfRes <- readRDS("E:/Group/react2_study5/saved_objects/rep_all_r123456_validated_newwrangle.rds")
#load data key
data_key=read_csv("E:/Group/react2_study5/report_phases_combined/projects/long_covid_2021/long_covid_data_key.csv")
# # Uncomment this if not running from the cockpit
# myStudyPop="SuspectedCOVID"
# 
# # select study population
# studyPop=chooseMyStudyPop(subset = myStudyPop)

# Create folders
createMySubfolder("clustering")

# Uncomment to load all results without re-running
# load("E:/Group/react2_study5/saved_objects/long_covid_clustering/")

# Prep data
dfRes$lc_84 <- dfRes$lc_84_all_29
vars_29 <- data_key$symptom_duration_84days_binary[data_key$one_of_29==1]
vars_reduced <- setdiff(unique(data_key$symptom_duration_84days_grouped_binary[data_key$one_of_29==1]),
                        c("lc_grouped_84_chills","lc_grouped_84_heavy_arms_legs" ))

dat_t0 <- prepMyData(dat = dfRes %>% filter(round %in% 3:5), t0 = T, min_followup = 84, 
                          abresult_filter = NULL,forclustering = T, covidafilter = c(1,2,3),
                          lc_names_12_weeks = unique(data_key$symptom_duration_84days_binary[data_key$one_of_29==1]))
dat_t12_3_5 <- prepMyData(dat = dfRes %>% filter(round %in% 3:5), t0 = F, min_followup = 84, 
                          abresult_filter = NULL,forclustering = T, covidafilter = c(1,2,3),
                          lc_names_12_weeks = unique(data_key$symptom_duration_84days_binary[data_key$one_of_29==1]))

# For parallel clustering analysis
dat_t12_3_5_rep <- prepMyData(dat = dfRes %>% filter(round %in% 3:5), t0 = F, min_followup = 84, 
                          abresult_filter = NULL,forclustering = T, covidafilter = c(1,2,3),
                          lc_names_12_weeks = vars_reduced)
dat_t12_3_5_rep <- dat_t12_3_5_rep[(rowSums(dat_t12_3_5_rep) >0),]


dat_t12_6 <- prepMyData(dat = dfRes %>% filter(round %in% 6, sample_type == "MAIN"), t0 = F, min_followup = 84, 
                          abresult_filter = NULL,forclustering = T, covidafilter = c(1,2,3),
                          lc_names_12_weeks = vars_reduced)





# table(dat_t12_3_5$lc_84_nausea_vomiting,dat_t12_3_5$lc_84_diarrhoea)

### Get numbers for test-train split
nobs=nrow(dat_t0)
n_varsel=floor(0.3*nobs)
n_analysis=nobs-n_varsel

### DO a test-train split
if(nrow(dat_t0)> 63000){
  set.seed(123)
  dat_t0_samp <- dat_t0 %>% sample_n(size = n_analysis)
  dat_t12_samp <- dat_t12_3_5[(rownames(dat_t12_3_5) %in% rownames(dat_t0_samp)),]
  dat_t12_3_5_rep_samp <- dat_t12_3_5_rep[(rownames(dat_t12_3_5_rep) %in% rownames(dat_t0_samp)),]
  
}else{
  dat_t0_samp = dat_t0
  dat_t12_samp <- dat_t12_3_5
  dat_t12_3_5_rep_samp <- dat_t12_3_5_rep
  
}

### save obs for variable selection
ids_varselection=setdiff(rownames(dat_t0), rownames(dat_t0_samp))

# # Run cluster analysis t12
t12_cluster_rows_r3_5=runClusterAnalysis(dat_t12_samp, stability = T,sampsize = 0.2, num_bootstraps = 50)

# Replication
t12_cluster_rows_r3_5_rep=runClusterAnalysis(dat_t12_3_5_rep_samp, stability = T,sampsize = 0.2, num_bootstraps = 50,
                                             symptom_names = vars_reduced,force_k = 2)
t12_cluster_rows_r6_rep=runClusterAnalysis(dat_t12_6, stability = T,sampsize = 0.2, num_bootstraps = 50,
                                           symptom_names = vars_reduced,force_k = 2)

# what are the mediods?
t12_cluster_rows_r3_5_rep$cluster_results$pamobject$medoids
t12_cluster_rows_r6_rep$cluster_results$pamobject$medoids



# Save silhouette plots
OverReact::saveREACTplot(p = t12_cluster_rows_r3_5$silhouette_plot,figpath = figpath,
                         filename = "t12_3_5_silhouette_plot",width = 5,height = 4)
OverReact::saveREACTplot(p = t12_cluster_rows_r3_5_rep$silhouette_plot,figpath = figpath,
                         filename = "t12_cluster_rows_r3_5_rep_silhouette_plot",width = 5,height = 4)
OverReact::saveREACTplot(p = t12_cluster_rows_r6_rep$silhouette_plot,figpath = figpath,
                         filename = "t12_cluster_rows_r6_rep_silhouette_plot",width = 5,height = 4)



# No tiredness ------------------------------------------------------------

# get new dfs
dat_t12_3_5_rep_samp_notired <- dat_t12_3_5_rep_samp %>% select(-lc_grouped_84_tiredness)
dat_t12_3_5_rep_samp_notired <- dat_t12_3_5_rep_samp_notired[rowSums(dat_t12_3_5_rep_samp_notired)!=0,]
dat_t12_6_notired <- dat_t12_6 %>% select(-lc_grouped_84_tiredness)
dat_t12_6_notired <- dat_t12_6_notired[rowSums(dat_t12_6_notired)!=0,]

# Replication
t12_cluster_rows_r3_5_rep_notiredness=runClusterAnalysis(dat_t12_3_5_rep_samp_notired, 
                                                         stability = F,sampsize = 0.2, num_bootstraps = 20,
                                                         symptom_names = vars_reduced,force_k = 2)
t12_cluster_rows_r6_rep_notiredness=runClusterAnalysis(dat_t12_6_notired,
                                                       stability = F,sampsize = 0.2, num_bootstraps = 20,
                                                       symptom_names = vars_reduced,force_k = 2)

t12_cluster_rows_r3_5_rep_notiredness$cluster_results$pamobject$medoids
t12_cluster_rows_r6_rep_notiredness$cluster_results$pamobject$medoids





# Replication
t12_cluster_rows_r3_5_rep_notiredness_all_symptomatic=runClusterAnalysis(dat_t12_3_5_rep_samp_notired, 
                                                         stability = F,sampsize = 0.2, num_bootstraps = 20,
                                                         symptom_names = vars_reduced,force_k = 2)
t12_cluster_rows_r6_rep_notiredness_all_symptomatic=runClusterAnalysis(dat_t12_6_notired,
                                                       stability = F,sampsize = 0.2, num_bootstraps = 20,
                                                       symptom_names = vars_reduced,force_k = 2)

t12_cluster_rows_r3_5_rep_notiredness_all_symptomatic$cluster_results$pamobject$medoids
t12_cluster_rows_r6_rep_notiredness_all_symptomatic$cluster_results$pamobject$medoids






# Plot heatmaps -----------------------------------------------------------


figname="t12_cluster_rows_r3_5"
# Create heatmaps
t12_heatmap_r3_5_dummy=createMyHeatmap(rowclust=t12_cluster_rows_r3_5,
                                     cluster_prefix="_",
                                 cluster_rows = F,
                                     figname="Rounds 3-5",
                                     roworder = NULL,
                                     # rowsplit =  t0_cluster_cols$cluster_results$pamobject$clustering,
                                     save_heatmap=F,
                                     symptoms = unique(data_key$short_name[data_key$one_of_29==1])
                                     # symptom_indices=symptom_subset_indexes_datadriven
)
t12_heatmap_r3_5_dummy$heatmap
# Create heatmaps
t12_heatmap_r3_5=createMyHeatmap(rowclust=t12_cluster_rows_r3_5,
                                 cluster_prefix="_",
                                 cluster_rows = F,
                                 figname="Rounds 3-5",
                                 roworder = t12_heatmap_r3_5_dummy$result_df[order(-t12_heatmap_r3_5_dummy$result_df$`Cluster _2\nn=4441`),]$symptom,
                                 # rowsplit =  t0_cluster_cols$cluster_results$pamobject$clustering,
                                 save_heatmap=F,
                                 symptoms = unique(data_key$short_name[data_key$one_of_29==1])
                                 # symptom_indices=symptom_subset_indexes_datadriven
)

t12_heatmap_r3_5$heatmap


png(paste0(figpath,figname,".png"), 
    width = 5,
    height = 7, units = "in", res = 300)
t12_heatmap_r3_5$heatmap
dev.off()





# PLot replication heatmaps -----------------------------------------------



figname="t12_cluster_rows_r3_5_rep"

# Create heatmaps
t12_heatmap_r3_5_rep_dummy=createMyHeatmap(rowclust=t12_cluster_rows_r3_5_rep,
                            cluster_prefix="R3-5 L_",
                            figname="Rounds 3-5",
                            roworder = NULL,
                            # rowsplit =  t0_cluster_cols$cluster_results$pamobject$clustering,
                            save_heatmap=F,
                            symptoms = setdiff(unique(data_key$short_name_grouped[data_key$one_of_29==1]),
                                               c("Chills","Heavy arms/legs"  ))
                            # symptom_indices=symptom_subset_indexes_datadriven
)
t12_heatmap_r3_5_rep_dummy$heatmap

# Create heatmaps
t12_heatmap_r3_5_rep=createMyHeatmap(rowclust=t12_cluster_rows_r3_5_rep,
                                     cluster_prefix="R3-5 L_",
                                     figname="Rounds 3-5",
                                     # roworder = NULL,
                                     roworder = t12_heatmap_r3_5_rep_dummy$result_df[order(-t12_heatmap_r3_5_rep_dummy$result_df$`Cluster R3-5 L_2\nn=4677`),]$symptom,
                                     # rowsplit =  t0_cluster_cols$cluster_results$pamobject$clustering,
                                     save_heatmap=F,
                                     symptoms = setdiff(unique(data_key$short_name_grouped[data_key$one_of_29==1]),
                                                        c("Chills","Heavy arms/legs"  ))
                                     # symptom_indices=symptom_subset_indexes_datadriven
)
t12_heatmap_r3_5_rep$heatmap
png(paste0(figpath,figname,".png"), 
    width = 5,
    height = 7, units = "in", res = 300)
t12_heatmap_r3_5_rep$heatmap
dev.off()


# Create heatmap R6
figname="t12_cluster_rows_r6_rep"
t12_heatmap_r6_rep=createMyHeatmap(rowclust=t12_cluster_rows_r6_rep,
                                   # cluster_rows = T,show_row_dend = F,
                                     cluster_prefix="R6 L_",reverse_cluster_order = T,
                                     figname="Round 6",
                                   roworder = t12_heatmap_r3_5_rep_dummy$result_df[order(-t12_heatmap_r3_5_rep_dummy$result_df$`Cluster R3-5 L_2\nn=4677`),]$symptom,
                                   # rowsplit =  t0_cluster_cols$cluster_results$pamobject$clustering,
                                     save_heatmap=F,
                                     symptoms = setdiff(unique(data_key$short_name_grouped[data_key$one_of_29==1]),
                                                        c("Chills","Heavy arms/legs"  ))
                                     # symptom_indices=symptom_subset_indexes_datadriven
)
t12_heatmap_r6_rep$heatmap
png(paste0(figpath,figname,".png"), 
    width = 5,
    height = 7, units = "in", res = 300)
t12_heatmap_r6_rep$heatmap
dev.off()

# combined heatmap
combined_heatmap_rep <- draw((t12_heatmap_r3_5_rep$heatmap+t12_heatmap_r6_rep$heatmap), ht_gap =unit(1,"cm"))
png(paste0(figpath,"combined_rep_heatmap.png"), 
    width = 8,
    height = 7, units = "in", res = 300)
combined_heatmap_rep
dev.off()




# Not tiredness -----------------------------------------------------------

figname="t12_cluster_rows_r3_5_rep_notiredness"

# Create heatmaps
t12_heatmap_r3_5_rep_notiredness_dummy=createMyHeatmap(rowclust=t12_cluster_rows_r3_5_rep_notiredness,
                                     cluster_prefix="R3-5 L_",
                                     figname="Rounds 3-5",
                                     roworder = NULL,reverse_cluster_order = F,
                                     # rowsplit =  t0_cluster_cols$cluster_results$pamobject$clustering,
                                     save_heatmap=F,
                                     symptoms = setdiff(unique(data_key$short_name_grouped[data_key$one_of_29==1]),
                                                        c("Chills","Heavy arms/legs","Tiredness" ))
                                     # symptom_indices=symptom_subset_indexes_datadriven
)

# Create heatmaps
t12_heatmap_r3_5_rep_notiredness=createMyHeatmap(rowclust=t12_cluster_rows_r3_5_rep_notiredness,
                                                 cluster_prefix="R3-5 L_",
                                                 figname="Rounds 3-5",reverse_cluster_order = F,
                                                 roworder = t12_heatmap_r3_5_rep_notiredness_dummy$result_df[order(-t12_heatmap_r3_5_rep_notiredness_dummy$result_df[,3]),]$symptom,
                                                 # rowsplit =  t0_cluster_cols$cluster_results$pamobject$clustering,
                                                 save_heatmap=F,
                                                 symptoms = setdiff(unique(data_key$short_name_grouped[data_key$one_of_29==1]),
                                                                    c("Chills","Heavy arms/legs","Tiredness" ))
                                                 # symptom_indices=symptom_subset_indexes_datadriven
)

t12_heatmap_r3_5_rep_notiredness$heatmap

# Create heatmap R6
figname="t12_cluster_rows_r6_rep"
t12_heatmap_r6_rep_notiredness=createMyHeatmap(rowclust=t12_cluster_rows_r6_rep_notiredness,
                                   # cluster_rows = T,show_row_dend = F,
                                   cluster_prefix="R6 L_",
                                   figname="Round 6",reverse_cluster_order = F,
                                   # roworder = roworder,
                                   roworder = t12_heatmap_r3_5_rep_notiredness_dummy$result_df[order(-t12_heatmap_r3_5_rep_notiredness_dummy$result_df[,3]),]$symptom,
                                   save_heatmap=F,
                                   symptoms = setdiff(unique(data_key$short_name_grouped[data_key$one_of_29==1]),
                                                      c("Chills","Heavy arms/legs","Tiredness"))
                                   # symptom_indices=symptom_subset_indexes_datadriven
)
t12_heatmap_r6_rep_notiredness$heatmap

t12_heatmap_r3_5_rep_notiredness$cluster_results$pamobject$medoids
t12_cluster_rows_r6_rep_notiredness$cluster_results$pamobject$medoids

# combined heatmap
combined_heatmap_notiredness <- draw((t12_heatmap_r3_5_rep_notiredness$heatmap+t12_heatmap_r6_rep_notiredness$heatmap), ht_gap =unit(1,"cm"))
combined_heatmap_notiredness
png(paste0(figpath,"combined_rep_heatmap_notiredness",".png"), 
    width = 8,
    height = 7, units = "in", res = 300)
combined_heatmap_notiredness
dev.off()



t12_cluster_rows_r3_5_rep_notiredness$silhouette_plot
t12_cluster_rows_r6_rep_notiredness$silhouette_plot

t12_cluster_rows_r3_5_rep_notiredness_all_symptomatic$silhouette_plot
t12_cluster_rows_r6_rep_notiredness_all_symptomatic$silhouette_plot



as.data.frame(cbind(tiredness=dat_t12_3_5_rep_samp[names(t12_cluster_rows_r3_5_rep_notiredness_all_symptomatic$cluster_results$pamobject$clustering),]$lc_grouped_84_tiredness,
      cluster=t12_cluster_rows_r3_5_rep_notiredness_all_symptomatic$cluster_results$pamobject$clustering)) %>% 
  group_by(tiredness,cluster) %>% 
  summarise(n=n())





# Save output -------------------------------------------------------------

### Save clusters and ids
df_clusts_t12 <- data.frame(u_passcode=rownames(dat_t12_samp), 
                            cluster_t12=unlist(t12_cluster_rows_r3_5$cluster_results$pamobject$clustering))

OverReact::saveREACTtable(tab = df_clusts_t12,
                          outpath = "E:/Group/react2_study5/saved_objects/long_covid_clustering", 
                          filename = "t12_cluster_assignment_3_5")
OverReact::saveREACT(ids_varselection,
                     outpath = "E:/Group/react2_study5/saved_objects/long_covid_clustering/", 
                     filename = "ids_varselection")

# save all clustering data
OverReact::saveREACT(file = list(r3_5_main_analysis=list(cluster_obj=t12_cluster_rows_r3_5,
                                                         roworder=t12_heatmap_r3_5_dummy$result_df[order(-t12_heatmap_r3_5_dummy$result_df$`Cluster _2\nn=4441`),]$symptom),
                                 r3_5_replication=list(cluster_obj=t12_cluster_rows_r3_5_rep,
                                                       roworder=t12_heatmap_r3_5_rep_dummy$result_df[order(-t12_heatmap_r3_5_rep_dummy$result_df$`Cluster R3-5 L_2\nn=4677`),]$symptom),
                                 
                                 r3_6_replication=list(cluster_obj=t12_cluster_rows_r6_rep,
                                                       roworder=t12_heatmap_r3_5_rep_dummy$result_df[order(-t12_heatmap_r3_5_rep_dummy$result_df$`Cluster R3-5 L_2\nn=4677`),]$symptom)),
                     figpath ="E:/Group/react2_study5/saved_objects/long_covid_clustering/",
                     outpath = "E:/Group/react2_study5/saved_objects/long_covid_clustering/",
                     filename = "all_clustering_data")


save.image("E:/Group/react2_study5/saved_objects/long_covid_clustering/clustering_results.rdata")

