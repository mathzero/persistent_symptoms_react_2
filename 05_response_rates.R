
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
source("E:/Group/functions/OR_unconcatenator.R")
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
surveys_issued=read_csv("E:/Group/react2_study5/report_phases_combined/projects/long_covid_resubmission/surveys_issued.csv")
createMySubfolder("response_rates")


dfRes <- dfRes %>% mutate(age_group_newcats = case_when(age >= 85 ~ "85+",
                                                        age >= 75 ~ "75 to 84",
                                                        age >= 65 ~ "65 to 74",
                                                        age >= 55 ~ "55 to 64",
                                                        age >= 45 ~ "45 to 54",
                                                        age >= 35 ~ "35 to 44",
                                                        age >= 25 ~ "25 to 34",
                                                        age >= 18 ~ "18 to 24",
                                                        T~NA_character_))


generateResponseRateTab <- function(var, myround =myround,truevar, mydat){
  df <- mydat %>% 
    filter(round %in% myround, !is.na(.[[var]])) %>% 
    group_by(.[[var]]) %>% 
    summarise(n=n())
  tot=sum(df$n)
  names(df)[1] <- "Category"
  df$surveys_issued <- as.numeric(unlist(surveys_issued[surveys_issued$Variable == var,paste0("r",myround)]))
  df$survey_responses <- df$n
  df$Variable <- truevar
  df <- df %>% select(Variable,Category, surveys_issued, survey_responses)
  df$Category <- as.character(df$Category)
  df$survey_non_responders <- df$surveys_issued-df$survey_responses
  df$percent_responders <- round(100*df$survey_responses/sum(df$survey_responses),1)
  df$percent_non_responders <- round(100*df$survey_non_responders/sum(df$survey_non_responders),1)
  chiSq <- chisq.test(x=c(df$survey_responses,df$survey_non_responders))
  df$Chi_2_p_val <-   chiSq$p.value
  return(df)
  }



generateFullResponseRateTable <- function(surveys_issued,mydat, myround =3){
  vars <-  unique(surveys_issued$Variable)
  truevars <- c("Sex", "Age", "IMD decile")
  res_list <- list()
  for (i in 1:length(vars)){
    print(i)
    res_list[[i]] <- generateResponseRateTab(var = vars[[i]],myround=myround, truevar=truevars[[i]],mydat = mydat)
    
    }
  res <- bind_rows(res_list)
  colnames(res) <- c("Variable", "Category", "Surveys issued", "Surveys returned","Surveys not returned", 
                     "Percent of responses", "Percent of non-responses", "Chi2 p-value")
  return(res)
  }
dfRes_edit$sample_type %>% unique()
### Create variable for all 
surveys_issued$r35 <- surveys_issued$r3+surveys_issued$r4+surveys_issued$r5
dfRes_edit <- dfRes %>% filter(!sample_type %in% c("Boost_6574","Boost_5564"))
dfRes_edit$round <- case_when(dfRes_edit$round %in% c(3:5)~35,
                         TRUE ~ dfRes_edit$round)

# Generate tabs
tab_3_5 <- generateFullResponseRateTable(surveys_issued=surveys_issued,mydat = dfRes_edit,myround = 35)
tab_6 <- generateFullResponseRateTable(surveys_issued=surveys_issued,mydat = dfRes_edit,myround = 6)

# Save 
OverReact::saveREACTtable(tab = tab_3_5,outpath = outpath, filename = "response_rates_tab_r_3_5",save_rds = F)
OverReact::saveREACTtable(tab = tab_6,outpath = outpath, filename = "response_rates_tab_r_6",save_rds = F)


for(x in 1:6){
  tab <- generateFullResponseRateTable(surveys_issued=surveys_issued,
                                       mydat = dfRes %>% filter(!sample_type %in% c("Boost_6574","Boost_5564")),
                                       myround = x)
  OverReact::saveREACTtable(tab = tab,outpath = outpath, filename = paste0("response_rates_tab_r",x),save_rds = F)
  
}
