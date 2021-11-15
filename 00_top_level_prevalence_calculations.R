
#' First clear the environment of variables
rm(list=ls(all=TRUE))
setwd("E:/Group/react2_study5/report_phases_combined/projects/long_covid_resubmission/")
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
                  "egg", "circlize","datapasta",
                  "poLCA","snow","randomForest","ranger", "jtools", "moreparty", 
                  "ComplexHeatmap","blockcluster",
                  "readr","ggthemes", "questionr", "gridExtra", "future.apply", "foreach",
                  "doParallel","SHAPforxgboost","catboost","ModelMetrics","shapper",
                  "tidytext", "quanteda", "widyr", "igraph", "ggraph", "patchwork", "OverReact")
load_packages(package.list)



# Run data prep  ----------------------------------------------------------

source("E:/Group/react2_study5/report_phases_combined/projects/long_covid_2021/code/00_functions.R")
source("E:/Group/react2_study5/report_phases_combined/projects/functions/univariate_models.R")
source("E:/Group/react2_study5/report_phases_combined/projects/long_covid_resubmission/code/00_bits_and_pieces.R")

# Create folders
createMySubfolder("top_level_stats")

# Load data ---------------------------------------------------------------

# df_clusters_t12 <- readRDS("E:/Group/react2_study5/report_phases_combined/projects/long_covid_2021/output/SuspectedCOVID//t12_cluster_assignment.rds")
# df_clusters_t0 <- readRDS("E:/Group/react2_study5/report_phases_combined/projects/long_covid_2021/output/SuspectedCOVID/t0_cluster_assignment.rds")
dfRes <- dfRes <- readRDS("E:/Group/react2_study5/saved_objects/rep_all_r123456_validated_newwrangle.rds")
data_key=read_csv("E:/Group/react2_study5/report_phases_combined/projects/long_covid_2021/long_covid_data_key.csv")

dfRes$dummy="All participants"
dfRes$symptomatic <- case_when(dfRes$symptomatic  ==1 ~ "Symptomatic",
                               dfRes$symptomatic ==0 ~ "Asymptomatic",
                               TRUE ~NA_character_)
dfRes$res <- case_when(dfRes$res  ==1 ~ "Positive",
                       dfRes$res ==0 ~ "Negative",
                       TRUE ~NA_character_)

dfRes$lc_84_all_15_symptom_count <- pmax(as.numeric(rowSums(dfRes[,data_key$symptom_duration_84days_binary[data_key$reduced_symptom_set == 1]], na.rm= T)),
                                         as.numeric(rowSums(dfRes[,unique(data_key$symptom_duration_84days_grouped_binary[data_key$reduced_symptom_set == 1])], 
                                                            na.rm= T)))
dfRes$lc_84_all_15 <- as.numeric(dfRes$lc_84_all_15_symptom_count >0)
dfRes$lc_84_all_15_2 <- as.numeric(dfRes$lc_84_all_15_symptom_count >1)
dfRes$lc_84_all_15_3 <- as.numeric(dfRes$lc_84_all_15_symptom_count >2)
table(dfRes$round,dfRes$lc_84_all_15)
# Prepare data
# dfRes <- makeVarsFactors(dfRes = dfRes) # convert all character variables to factors for modelling
moddat <- prepMyData(dat = dfRes %>% filter(round %in% 3:6), t0 = T, min_followup = 84, 
                     abresult_filter = NULL,forclustering = F, covidafilter = c(1,2,3))




# Functions ---------------------------------------------------------------



### Function to create table
createMyTable <- function(mydat=dfRes,result_var = "lc_84_all_29",mycovs){
  
  
  ### Calculate prevalences 
  lc_prevs_full_cohort_weighted <- OverReact::makeTablesNew(dat = mydat, result_var = result_var, 
                                                            covariates = mycovs,sens = 1,spec = 1,
                                                            weights = "wt_antibody", for_report = T,
                                                            cov_name_list = cov_name_list,separator = ",",
                                                            output_list = F,sf = 3)
  lc_prevs_full_cohort_weighted <- lc_prevs_full_cohort_weighted %>% 
    rename(Prevalence_weighted = Prevalence)
  
  
  ### Calculate prevalences 
  lc_prevs_full_cohort_unweighted <- OverReact::makeTablesNew(dat = mydat, result_var = result_var, 
                                                              covariates = mycovs,sens = 1,spec = 1,
                                                              weights = NULL, separator = ",",for_report = T,
                                                              cov_name_list = cov_name_list,output_list = F,
                                                              sf = 3)
  lc_prevs_full_cohort_unweighted <- lc_prevs_full_cohort_unweighted %>% 
    rename(Prevalence_unweighted = Prevalence)
  
  lc_prevs_full_cohort_weighted$Prevalence_adjusted
  
  ### Combine tables into one 
  lc_prevs_full_cohort <- lc_prevs_full_cohort_unweighted %>% dplyr::select(-Prevalence_adjusted) %>% 
    left_join(lc_prevs_full_cohort_weighted %>% 
                dplyr::select(-Prevalence_adjusted, -Positive, -Total), 
              by = c("Variable", "Category"))  
    lc_prevs_full_cohort <- lc_prevs_full_cohort[(lc_prevs_full_cohort$Category  != 0),]
    
    # Replace 1s with Variable name
    lc_prevs_full_cohort[lc_prevs_full_cohort$Category==1,"Category"] <- 
      lc_prevs_full_cohort[lc_prevs_full_cohort$Category==1,"Variable"]
    
  return(lc_prevs_full_cohort)
}

### Combined table
createTableOne <- function(mydat=dfRes,result_var = "lc_84_all_29"){
  
  result_var_2=paste0(result_var,"_2")
  result_var_3=paste0(result_var,"_3")
  
  prevs_lc_84 <- createMyTable(mydat=mydat,result_var = result_var,mycovs = mycovs)
  # OverReact::saveREACTtable(tab = prevs_lc_84, outpath = outpath, filename = paste0("prevs_",result_var))
  
  
  # Long COVID any time (2/3 symptoms) -----------------------------------------------------
  
  prevs_lc_84_2 <- createMyTable(mydat=mydat,result_var = result_var_2,mycovs = mycovs)
  # OverReact::saveREACTtable(tab = prevs_lc_84_2, outpath = outpath, filename = paste0("prevs_",result_var_2))
  prevs_lc_84_3 <- createMyTable(mydat=mydat,result_var = result_var_3,mycovs = mycovs)
  # OverReact::saveREACTtable(tab = prevs_lc_84_3, outpath = outpath, filename = paste0("prevs_",result_var_3))
  
  ############## Master symptom table with 1/2/3 symptoms (weighted)
  prevs_lc_84_123 <-  cbind(Variable=prevs_lc_84$Variable,
                            Category=prevs_lc_84$Category,
                            total=prevs_lc_84$Total,
                            num_1_symp_weighted=prevs_lc_84$Positive,
                            num_2_symp_weighted=prevs_lc_84_2$Positive,
                            num_3_symp_weighted=prevs_lc_84_3$Positive,
                            prev_1_symp_weighted=prevs_lc_84$Prevalence_weighted,
                            prev_2_symp_weighted=prevs_lc_84_2$Prevalence_weighted,
                            prev_3_symp_weighted=prevs_lc_84_3$Prevalence_weighted) %>% as.data.frame()
  colnames(prevs_lc_84_123)[4:9] <- c("Number with 1 or more symptoms (weighted)",
                                      "Number with 2 or more symptoms (weighted)",
                                      "Number with 3 or more symptoms (weighted)",
                                      "Prevalence 1 or more symptoms (weighted)",
                                      "Prevalence 2 or more symptoms (weighted)",
                                      "Prevalence 3 or more symptoms (weighted)")
  # OverReact::saveREACTtable(tab = prevs_lc_84_123, outpath = outpath, filename = paste0("prevs_",result_var,"_123_weighted"))
  
  
  
  ############## Master symptom table with 1/2/3 symptoms (unweighted)
  prevs_lc_84_123_uw <-  cbind(Variable=prevs_lc_84$Variable,
                               Category=prevs_lc_84$Category,
                               total=prevs_lc_84$Total,
                               num_1_symp_unweighted=prevs_lc_84$Positive,
                               num_2_symp_unweighted=prevs_lc_84_2$Positive,
                               num_3_symp_unweighted=prevs_lc_84_3$Positive,
                               prev_1_symp_unweighted=prevs_lc_84$Prevalence_unweighted,
                               prev_2_symp_unweighted=prevs_lc_84_2$Prevalence_unweighted,
                               prev_3_symp_unweighted=prevs_lc_84_3$Prevalence_unweighted) %>% as.data.frame()
  colnames(prevs_lc_84_123_uw)[4:9] <- c("Number with 1 or more symptoms",
                                         "Number with 2 or more symptoms",
                                         "Number with 3 or more symptoms",
                                         "Prevalence 1 or more symptoms (unweighted)",
                                         "Prevalence 2 or more symptoms (unweighted)",
                                         "Prevalence 3 or more symptoms (unweighted)")
  # OverReact::saveREACTtable(tab = prevs_lc_84_123_uw, outpath = outpath, filename = paste0("prevs_",result_var,"_123_weighted"))
  
  
  
  prevs_lc_84_123_uw$Category <- factor(prevs_lc_84_123_uw$Category,levels = tableorder)
  prevs_lc_84_123_uw <- prevs_lc_84_123_uw[order(prevs_lc_84_123_uw$Category),]
  prevs_lc_84_123$Category <- factor(prevs_lc_84_123$Category,levels = tableorder)
  prevs_lc_84_123 <- prevs_lc_84_123[order(prevs_lc_84_123$Category),]
  
  # remove care home worker
  prevs_lc_84_123_uw <- prevs_lc_84_123_uw[!(prevs_lc_84_123_uw$Category=="No" & prevs_lc_84_123_uw$Variable =="Healthcare or care home worker"),]
  prevs_lc_84_123 <- prevs_lc_84_123[!(prevs_lc_84_123$Category=="No" & prevs_lc_84_123$Variable =="Healthcare or care home worker"),]
  
  # Switch smokers back round
  df_order <- data.frame(Variable=unique(prevs_lc_84_123$Variable))
  df_order$num=1:nrow(df_order)
  prevs_lc_84_123 <- prevs_lc_84_123 %>% left_join(df_order) %>% arrange(num) %>% select(-num)
  prevs_lc_84_123_uw <- prevs_lc_84_123_uw %>% left_join(df_order) %>% arrange(num) %>% select(-num)
  
  
  
  return(list(weighted=prevs_lc_84_123,
              unweighted=prevs_lc_84_123_uw))
}


# Run table ones ----------------------------------------------------------
# mycovs=c("symptomatic",mycovs)


### Among full 29
tab1_lc84_3_5 <- createTableOne(mydat = dfRes %>% filter(round%in%3:5),result_var = "lc_84_all_29")
tab1_lc84_6 <- createTableOne(mydat = dfRes %>% filter(round==6,sample_type =="MAIN"),result_var = "lc_84_all_29")
tab1_lc84_6_extended <- createTableOne(mydat = dfRes %>% filter(round==6,sample_type =="MAIN"),result_var = "lc_84_all_37")

# save
OverReact::saveREACTtable(tab = tab1_lc84_3_5$weighted, outpath = outpath, filename = "tab1_lc84_weighted_r3_5",save_rds = F)
OverReact::saveREACTtable(tab = tab1_lc84_6$weighted, outpath = outpath, filename = "tab1_lc84_weighted_r6",save_rds = F)
OverReact::saveREACTtable(tab = tab1_lc84_6_extended$weighted, outpath = outpath, filename = "tab1_lc84_weighted_r6_extended_symtoms",save_rds = F)


# Same plot among symptomatics only ---------------------------------------

### Among full 29
tab1_lc84_3_5_symp <- createTableOne(mydat = moddat %>% filter(round%in%3:5),result_var = "lc_84_all_29")
tab1_lc84_6_symp <- createTableOne(mydat = moddat %>% filter(round==6,sample_type =="MAIN"),result_var = "lc_84_all_29")
tab1_lc84_6_extended_symp <- createTableOne(mydat = moddat %>% filter(round==6,sample_type =="MAIN"),result_var = "lc_84_all_37")

# save
OverReact::saveREACTtable(tab = tab1_lc84_3_5_symp$unweighted, outpath = outpath, filename = "tab1_symptomatic_lc84_unweighted_r3_5",save_rds = F)
OverReact::saveREACTtable(tab = tab1_lc84_6_symp$unweighted, outpath = outpath, filename = "tab1_symptomatic_lc84_unweighted_r6",save_rds = F)
OverReact::saveREACTtable(tab = tab1_lc84_6_extended_symp$unweighted, outpath = outpath, filename = "tab1_lc84_symptomatic_unweighted_r6_extended_symtoms",save_rds = F)



# With reduced symptom set ------------------------------------------------


### Among full 29
tab1_lc84_3_5_reduced <- createTableOne(mydat = dfRes %>% filter(round%in%3:5),result_var = "lc_84_all_15")
tab1_lc84_6_reduced <- createTableOne(mydat = dfRes %>% filter(round==6,sample_type =="MAIN"),result_var = "lc_84_all_15")
OverReact::saveREACTtable(tab = tab1_lc84_3_5_reduced$weighted, outpath = outpath, filename = "tab1_lc84_weighted_r3_5_reduced_symptom_set",save_rds = F)
OverReact::saveREACTtable(tab = tab1_lc84_6_reduced$weighted, outpath = outpath, filename = "tab1_lc84_weighted_r6_reduced_symptom_set",save_rds = F)


# Among symptomatics with reduced symptom set -----------------------------



### Among full 29
tab1_lc84_3_5_symp <- createTableOne(mydat = moddat %>% filter(round%in%3:5),result_var = "lc_84_all_15")
tab1_lc84_6_symp <- createTableOne(mydat = moddat %>% filter(round==6,sample_type =="MAIN"),result_var = "lc_84_all_15")

# save
OverReact::saveREACTtable(tab = tab1_lc84_3_5_symp$unweighted, outpath = outpath, filename = "tab1_symptomatic_lc84_unweighted_r3_5_reduced_symptom_set",save_rds = F)
OverReact::saveREACTtable(tab = tab1_lc84_6_symp$unweighted, outpath = outpath, filename = "tab1_symptomatic_lc84_unweighted_r6_reduced_symptom_set",save_rds = F)




# Top level prevalence of reported COVID ----------------------------------
dfRes <- dfRes %>% mutate(prev_covid=case_when(covida %in% 1:3 & covidc_cat != "No symptoms" ~ 1,
                                               TRUE ~ 0),
                          study_pop = case_when(round %in% 3:5 ~ "Rounds 3-5",
                                                sample_type=="MAIN" ~ "Round 6",
                                                T ~ NA_character_))


dfRes$covidc_cat %>% unique()

table(dfRes$study_pop,dfRes$round, exclude = "none")

prevs_prev_covid_weight <- OverReact::makeTablesNew(dat = dfRes,result_var = "prev_covid",
                         covariates = "study_pop",
                         sens = 1,cov_name_list = NULL,
                         spec = 1, weights = "wt_antibody",output_list = F)


dfRes %>% dplyr::filter(!is.na(all_of(mycovs[1:12])))

table(complete.cases(dfRes[dfRes$prev_covid == 1,mycovs[1:12]]),dfRes[dfRes$prev_covid == 1,]$study_pop)



