
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


# define outcome
outcome = "lc_84_all_29"

# Create folders
createMySubfolder("networks")
# createMySubfolder(outcome)
dfRes$u_passcode
# Prepare data
dfRes$lc_84 <- dfRes$lc_84_all_29
# dfRes <- makeVarsFactors(dfRes = dfRes) # convert all character variables to factors for modelling
moddat <- prepMyData(dat = dfRes %>% filter(round%in%3:5, !u_passcode %in% ids_varselection), t0 = F, min_followup = 84, 
                     abresult_filter = NULL, 
                     covidafilter = c(1,2,3),forclustering = T,
                     lc_names_12_weeks = data_key$symptom_duration_84days_binary[data_key$one_of_29==1])
vars_reduced <- setdiff(unique(data_key$symptom_duration_84days_grouped_binary[data_key$one_of_29==1]),
                        c("lc_grouped_84_chills","lc_grouped_84_heavy_arms_legs" ))
moddat_6 <- prepMyData(dat = dfRes %>% filter(round%in%6, sample_type=="MAIN"), t0 = F, min_followup = 84, 
                       abresult_filter = NULL, 
                       covidafilter = c(1,2,3),forclustering = T,
                       lc_names_12_weeks = vars_reduced)

# Networks ----------------------------------------------------------------
cooccurmat <-crossprod(as.matrix(moddat))
colnames(cooccurmat) <- rownames(cooccurmat) <- data_key$short_name[data_key$one_of_29==1]
n <- cooccurmat %>% graph_from_adjacency_matrix(mode="undirected", diag=F, weighted=T)
V(n)$group <- cluster_optimal(n) %>% membership() %>% as.character()
V(n)$degree <-as.numeric(diag(cooccurmat))

cl <-cluster_optimal(n)
weights <- ifelse(crossing(cl,n),1,15)

#plot
p_network <- ggplot(ggnetwork(n, layout=igraph::with_fr(weights = weights)),
       aes(x,y, xend=xend, yend=yend)) +
  geom_edges(aes(alpha=weight), size=2) +
  geom_nodelabel(aes(label = name, size=degree, color = group)) +
  scale_alpha_continuous(guide=F,range = c(0,1))+
  scale_color_imperial(palette = "extended", guide=F)+
  scale_size_continuous(range = c(1,2), guide=F) +
  theme_blank() +
  coord_cartesian(clip="off")
  # theme(plot.margin = unit(c(1,1,1,1),"in"))
p_network

OverReact::saveREACTplot(p = p_network,figpath = figpath,
                         filename = "network_symptoms_r3_5",width = 4,height=4)




# Round 6 -----------------------------------------------------------------

cooccurmat <-crossprod(as.matrix(moddat_6))
colnames(cooccurmat) <- rownames(cooccurmat) <- setdiff(unique(data_key$short_name_grouped[data_key$one_of_29==1]),
                                                        c("Chills","Heavy arms/legs"))
n <- cooccurmat %>% graph_from_adjacency_matrix(mode="undirected", diag=F, weighted=T)
V(n)$group <- cluster_optimal(n) %>% membership() %>% as.character()
V(n)$degree <-as.numeric(diag(cooccurmat))


cl <-cluster_optimal(n)
weights <- ifelse(crossing(cl,n),1,15)

#plot
p_network_6 <- ggplot(ggnetwork(n, layout=igraph::with_fr(weights = weights)),
                    aes(x,y, xend=xend, yend=yend)) +
  geom_edges(aes(alpha=weight), size=2) +
  geom_nodelabel(aes(label = name, size=degree, color = group)) +
  scale_alpha_continuous(guide=F,range = c(0,1))+
  scale_color_imperial(palette = "red", guide=F, reverse=T)+
  scale_size_continuous(range = c(1,2), guide=F) +
  theme_blank() +
  coord_cartesian(clip="off")
# theme(plot.margin = unit(c(1,1,1,1),"in"))
p_network_6

OverReact::saveREACTplot(p = p_network_6,figpath = figpath,
                         filename = "network_symptoms_r6",width = 4,height=4)






