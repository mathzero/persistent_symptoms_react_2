
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
                  "readr","ggthemes", "questionr", "gridExtra", "future.apply", "foreach",
                  "doParallel","SHAPforxgboost","catboost","ModelMetrics","shapper",
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
createMySubfolder("EDA")

### Filter data to study cohort & Create new variable with asymptomatic 
mydat <- prepMyData(dat = dfRes %>% filter(round %in% 3:5), t0 = T, min_followup =150, abresult_filter = NULL,
                    forclustering = F,covidafilter = c(1,2,3))

symptom_duration_list <- data_key$symptom_duration_continuous[data_key$one_of_29==1]
mydat[,symptom_duration_list]



plotSymtomPersistenceOverTime <- function(mydat,myvar = "age_group_named", symptomatic = F,varname="Age", 
                                          symptom_count_var = "lc_84_all_29_symptom_count",
                                          symptom_durations=symptom_duration_list,
                                          ndays=150, use_subset=F){
  
  if(symptomatic){
    insert = "symptomatic "
  }else{
    insert = NULL
    
  }
  if(use_subset){
    symptom_durations_plot=symptom_durations[c(1,9, 10, 13, 14, 15, 17, 18, 21, 22, 23, 24, 25, 26, 27)]
  }else{
    symptom_durations_plot=symptom_durations
  }
  
  
  nobs=nrow(mydat)
  # uniques <- unique(unlist(mydat[myvar])) %>% as.character()
  if(class(unlist(mydat[,myvar]))!= "factor"){
    mydat[,myvar] <- as.factor(unlist(mydat[,myvar]))
  }
  uniques <- levels(unlist(mydat[,myvar]))
  uniques <- uniques[!is.na(uniques)]
  symptom_persistence_df_list <- list()

  for (x in 1:length(uniques)){
    mydat$myvar = mydat[,myvar]
    
    dfLongCovSubset <- mydat %>% filter(myvar==uniques[[x]])
    
    num_in_age <- nrow(dfLongCovSubset)
    
    symptom_persistence_df <- data.frame(day_count=1:ndays,
                                         nobs=NA,
                                         one_or_more=NA, 
                                         two_or_more=NA, 
                                         three_or_more=NA, 
                                         four_or_more=NA,
                                         five_or_more = NA)
    
    for (i in 1:ndays){
      count_censored <- sum(dfLongCovSubset$test_date - dfLongCovSubset$covidsta < i, na.rm = T)
      count_uncensored <- sum(dfLongCovSubset$test_date - dfLongCovSubset$covidsta >= i, na.rm = T) + sum(dfLongCovSubset$covidc_cat == "No symptoms", na.rm = T)
      
      sum_symps <- rowSums(dfLongCovSubset[,symptom_durations_plot] > i, na.rm=T)
      
      symptom_persistence_df[[i,2]] <- num_in_age
      symptom_persistence_df[[i,3]] <- as.numeric(sum(sum_symps>0, na.rm=T))
      symptom_persistence_df[[i,4]] <- as.numeric(sum(sum_symps>1, na.rm=T))
      symptom_persistence_df[[i,5]] <- as.numeric(sum(sum_symps>2, na.rm=T))
      symptom_persistence_df[[i,6]] <- as.numeric(sum(sum_symps>3, na.rm=T))
      symptom_persistence_df[[i,7]] <- as.numeric(sum(sum_symps>4, na.rm=T))
      
      
    }
    
    ### pivot longer
    symptom_persistence_df_long <- symptom_persistence_df %>% pivot_longer(cols=3:ncol(symptom_persistence_df))
    
    ### add CIs
    symptom_persistence_df_long$lower= propCI(x = pull(symptom_persistence_df_long,"value"),
                                              n = pull(symptom_persistence_df_long,"nobs"), 
                                              method="wilson")$lower * 100
    
    symptom_persistence_df_long$p= propCI(x = pull(symptom_persistence_df_long,"value"),
                                          n = pull(symptom_persistence_df_long,"nobs"), 
                                          method="wilson")$p * 100
    symptom_persistence_df_long$upper= propCI(x = pull(symptom_persistence_df_long,"value"),
                                              n = pull(symptom_persistence_df_long,"nobs"), 
                                              method="wilson")$upper * 100
    
    
    ### add age group
    symptom_persistence_df_long$group <- rep(uniques[[x]],nrow(symptom_persistence_df_long))
    
    symptom_persistence_df_list[[x]] <-symptom_persistence_df_long
  }
  
  ### Bind all together
  symptom_persistence_df_age <- do.call(rbind, symptom_persistence_df_list)

  highlight_df <- symptom_persistence_df_age %>% filter(day_count %in% c(28,84)) %>% 
    mutate(name=stringr::str_to_sentence(gsub("_", " ", name)),
           name = factor(name, levels = unique(name)), 
           value = 100*value/nobs,
           mygroup = factor(group, levels = uniques),
           four_weeks = ifelse(day_count == 28, 1,0),
           twelve_weeks = ifelse(day_count == 84, 1,0))
  
  ### plot as a proportion of symptomatic
  p <- symptom_persistence_df_age %>% 
    mutate(name=stringr::str_to_sentence(gsub("_", " ", name)),
           name = factor(name, levels = unique(name)), 
           mygroup = factor(group, levels = uniques), 
           value = 100*value/nobs,
           four_weeks = ifelse(day_count == 28, 1,0),
           twelve_weeks = ifelse(day_count == 84, 1,0)) %>%
    ggplot(aes(x=day_count, y=value, col = name, fill = name)) +
    geom_line(stat="identity") +
    geom_ribbon(aes(ymin = lower, ymax=upper, fill = name), alpha =0.3, colour = NA,show.legend = F) +
    geom_point(data=highlight_df, 
               aes(x=day_count, y = value,  col = name, 
                   fill = name), size=3,
               show.legend = F) +
    theme_clean() +
    ylim(0,100) +
    facet_wrap(.~mygroup, nrow=1) +
    scale_color_manual(values = myCols)+
    scale_fill_manual(values = myCols) + 
    geom_vline(xintercept = 28, linetype = "dashed", col = "grey60") +
    geom_vline(xintercept = 84, linetype = "dashed", col = "grey60") +
    labs(x="Days since symptom onset", y="Percentage of respondents", col = "Number of symptoms",
         title = varname,
         subtitle = paste0("n=",nobs)) + 
    theme_adjust +
    theme(plot.background  = element_blank(),
          panel.background  = element_rect(color=NA),
          legend.position = "bottom",
          legend.text.align = 0,
          legend.background = element_rect(color=NA))
    
  
  # ### Get nobs plot
  # p2=symptom_persistence_df_age %>%
  #   mutate(nobadjust=nobs/5) %>% 
  #   ggplot(aes(x=day_count, y = nobadjust)) +
  #   geom_col(fill=myCols[2], alpha=0.8) + 
  #   facet_wrap(.~group, nrow=1) +
  #   theme_bw() +
  #   labs(x="Days since symptom onset",y="Number of \nrespondents \nin denominator") + theme_adjust
  # p2
  p
  # p_comb=p/p2+plot_layout(guides="collect", heights = c(3,1))
  # p_comb
  return(p)
}




# Create plots ------------------------------------------------------------
mydat$sex_mw <- case_when(mydat$sex=="Male" ~ "Men",
                          mydat$sex=="Female" ~ "Women",
                          TRUE ~ NA_character_)

mydat$age_group_four_cat <- factor(mydat$age_group_four_cat, levels=c("Under 60","60-69", "70-79","80+"))

# create plots
age_plot <- plotSymtomPersistenceOverTime(mydat, myvar = "age_group_four_cat",varname = "Age")
sex_plot <- plotSymtomPersistenceOverTime(mydat, myvar = "sex_mw",varname = "Sex")
ethnic_plot <- plotSymtomPersistenceOverTime(mydat, myvar = "ethnic_new",varname = "Ethnicity")
smoke_plot <- plotSymtomPersistenceOverTime(mydat, myvar = "smokenow",varname = "Smoking")
bmi_plot <- plotSymtomPersistenceOverTime(mydat, myvar = "bmi_cat",varname = "Adiposity")

# save plots
OverReact::saveREACTplot(p = age_plot,figpath = figpath,filename = "symptom_persistence_age",width = 14,height = 6)
OverReact::saveREACTplot(p = sex_plot,figpath = figpath,filename = "symptom_persistence_sex",width = 10,height = 6)
OverReact::saveREACTplot(p = ethnic_plot,figpath = figpath,filename = "symptom_persistence_ethnicity",width = 16,height = 6)
OverReact::saveREACTplot(p = smoke_plot,figpath = figpath,filename = "symptom_persistence_smoking",width = 12,height = 6)
OverReact::saveREACTplot(p = bmi_plot,figpath = figpath,filename = "symptom_persistence_bmi",width = 14,height = 6)


# create panel plot
layout <- "
AA#####
BBBB###
CCC####
DDD####
EEEEE##
"
p_comb <- ((sex_plot+age_plot +smoke_plot+bmi_plot+ethnic_plot)) + 
  plot_layout(guides = "collect", design=layout) & 
  theme(legend.position = "bottom",
        plot.background  = element_blank()
        # plot.margin = unit(0.3,"cm")
        )
p_comb
OverReact::saveREACTplot(p = p_comb,figpath = figpath,filename = "symptom_persistence_panel",width = 12,height = 16)
