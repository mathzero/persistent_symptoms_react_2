
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

# Get Ids to exclude
ids_varselection <- readRDS("E:/Group/react2_study5/saved_objects/long_covid_clustering/ids_varselection.rds")


# define outcome
outcome = "lc_84_all_29"


# Prepare data
dfRes$lc_84 <- dfRes$lc_84_all_29
# dfRes <- makeVarsFactors(dfRes = dfRes) # convert all character variables to factors for modelling
moddat <- prepMyData(dat = dfRes %>% filter(round%in%3:5, !u_passcode %in% ids_varselection), t0 = F, 
                     min_followup = NULL, 
                     abresult_filter = NULL, 
                     covidafilter = NULL,
                     forclustering = F,
                     lc_names_12_weeks = data_key$symptom_duration_84days_binary[data_key$one_of_29==1])
vars_reduced <- setdiff(unique(data_key$symptom_duration_84days_grouped_binary[data_key$one_of_29==1]),
                        c("lc_grouped_84_chills","lc_grouped_84_heavy_arms_legs" ))
moddat_6 <- prepMyData(dat = dfRes %>% filter(round%in%6, sample_type=="MAIN"), t0 = F, min_followup = NULL, 
                       abresult_filter = NULL, 
                       covidafilter = NULL,forclustering = F,
                       lc_names_12_weeks = vars_reduced)



# Choose which data to work on --------------------------------------------------------------


setwd("E:/Group/react2_study5/report_phases_combined/projects/long_covid_resubmission")
outpath <- "E:/Group/react2_study5/report_phases_combined/projects/long_covid_resubmission/output/"
figpath <- "E:/Group/react2_study5/report_phases_combined/projects/long_covid_resubmission/plots/"

# Create folders
createMySubfolder("age_trajectory")

# Define which rounds we're working on
whichround <- "round_6"
# quick ifelse to switch data and assign subfolder
createMySubfolder(whichround)
if(whichround=="rounds_3_5"){
  mydat=moddat
}else{
  mydat=moddat_6
}

### change sex_mw to factor
mydat$sex_mw <- as.factor(case_when(mydat$sex == "Male" ~ "Men",
                                    mydat$sex == "Female" ~ "Women",
                                     TRUE ~ NA_character_))



# Run models --------------------------------------------------------------


### Define function to plot results of spline model
makeSplinePlot <- function(mod, dat, yaxisvar = "Symptomaticness"){
  
  ### Extract predictions to plot
  newdata=expand.grid(c("Men", "Women"), seq(20,80,1))
  preds <- predict(object = mod, type = "response", se.fit = T, newdata=
                     data.frame(sex_mw=newdata$Var1,
                                age=newdata$Var2))
  dat <- dat %>% filter(!is.na(sex_mw) ,!is.na(age_group_named))
  preddat <- cbind(newdata, preds=unlist(preds$fit), preds_se=unlist(preds$se.fit))
  preddat$lower = preddat$preds -1.96*preddat$preds_se
  preddat$upper = preddat$preds +1.96*preddat$preds_se
  names(preddat)[c(1:2)] <- c("sex_mw","age")
  preddat %>% 
    ggplot(aes(x=age, y=(preds), col = factor(sex_mw)) ) +
    geom_line(linetype = "dashed")+
    scale_fill_manual(values = myCols) +
    xlim(18,80) +
    # ylim(0,0.8) +
    scale_colour_manual(values = myCols) +
    geom_ribbon(aes(ymin = lower, 
                    ymax = upper,
                    colour =factor(sex_mw),
                    fill =factor(sex_mw)), 
                alpha = 0.1) +
    # facet_wrap(.~age_group_four_cat, ncol =4) +
    theme_bw() +
    labs(fill="", col="", x="Age",y= paste0("Modelled probability of", yaxisvar))
  
}


# A. Generate plots from GAM models ------------------------------------------


# 1. Probability of symptoms at any time ----------------------------------


### adjusted for age and gender, clinical vulnerability and history of COVID
f <- as.formula("symptomatic ~s(as.numeric(age), by = sex_mw) + sex_mw")
mod_gam_symp <- gam(f, data = mydat,family = binomial(link = "logit"),
                    method = "REML"
                    # ,weights=mydat$wt_antibody
)
summary(mod_gam_symp)
### Create OR table
mod_tab_symp <- makeORTable(mod_gam_symp, ref_level=NULL)

p1 <- makeSplinePlot(mod = mod_gam_symp,dat = mydat,yaxisvar = "symptoms at any time")



# 2. Probability of Long COVID  -------------------------------------------




### adjusted for age and gender, clinical vulnerability and history of COVID
f <- as.formula("lc_84 ~s(as.numeric(age), by = sex_mw) + sex_mw")
mod_gam <- gam(f, data = mydat,family = binomial(link = "logit"), 
               method = "REML"
               # ,weights=mydat$wt_antibody
)
summary(mod_gam)
### Create OR table
mod_tab <- makeORTable(mod_gam, ref_level=NULL)

### plot 
p2 <- makeSplinePlot(mod = mod_gam, 
                     dat = mydat,
                     yaxisvar = "persistent symptoms for 12 weeks")
p2



# 3. Probability of long COVID, conditional on being symptomatic -------------

### adjusted for age and gender, clinical vulnerability and history of COVID
f <- as.formula("lc_84 ~s(as.numeric(age), by = sex_mw) + sex_mw")
mod_gam_symp_only <- gam(f, data = mydat %>%  filter(symptomatic == 1),
                         family = binomial(link = "logit"), method = "REML"
                         # ,
                         # weights=mydat[mydat$symptomatic==1,]$wt_antibody
)

### plot 
p3 <- makeSplinePlot(mod = mod_gam_symp_only,dat = mydat%>%  filter(symptomatic == 1),yaxisvar = "persistent symptoms for 12 weeks")




lc_84_comb <-   p1+ggtitle("Probability of experiencing symptoms at any time") + 
  p2+ggtitle("Probability of persistent symptoms (12 weeks)") +
  p3+ggtitle("Probability of persistent symptoms (12 weeks), \nconditional on being symptomatic") +
  plot_layout(guides="collect")

lc_84_comb

### Save plot
OverReact::saveREACTplot(p = lc_84_comb,figpath = figpath,filename = "lc_84_age_sex_mw_spline_plot",
                         width = 15, height=6)

# ### Save plot as file for joining with real data
# OverReact::saveREACT(file = lc_84_comb,outpath = outpath, filename = "lc_84_age_sex_mw_spline_plot")
# 



# B. Generate equivalent plots with real data --------------------------------


vars_of_interest <- c( "age_group_named")


# 1. Probability of symptoms at any time ----------------------------------

# Get table of symptomatics
prevs_symptoms <- OverReact::stratifiedTables(dat = mydat,result_var = "symptomatic", 
                                              strat_var="sex_mw",
                                              covariates = vars_of_interest,
                                              cov_name_list = cov_name_list,
                                              sens = 1,spec = 1, output_list = F,
                                              weights = "wt_antibody", 
                                              percent = T)  %>% 
  pivot_longer(cols=c("Men","Women"))%>% rename(sex_mw=name)
prevs_symptoms <- cbind(prevs_symptoms,OR_unconcat(prevs_symptoms$value,
                                                   lbracket = "(", rbracket = ")", 
                                                   separator = "-"))



### Create plot
dodger=position_dodge(width=0.4)

p_symptoms <- prevs_symptoms %>% 
  ggplot(aes(x=Category, y=OR_sep, col=factor(sex_mw)))+
  geom_errorbar(aes(ymin=CI_lower, ymax=CI_upper), width=0.4, position=dodger) + 
  geom_point(position=dodger) + 
  scale_colour_manual(values = myCols[c(1,2)],guide=guide_legend(reverse=F)) +
  theme_bw() +
  labs(fill="", col="", x="Age",
       y="Weighted percentage of respondents"
       # subtitle = "Percentage of respondents who experienced symptoms"
  )

p_symptoms


# 2. Probability of Long COVID  -------------------------------------------

# Get table of lcs
prevs_lc_all <- OverReact::stratifiedTables(dat = mydat,result_var = "lc_84", 
                                            covariates = vars_of_interest,
                                            strat_var="sex_mw",
                                            
                                            cov_name_list = cov_name_list,
                                            sens = 1,spec = 1, output_list = F,
                                            weights = "wt_antibody", 
                                            percent = T) %>% 
  pivot_longer(cols=c("Men","Women"))%>% rename(sex_mw=name)
prevs_lc_all <- cbind(prevs_lc_all,OR_unconcat(prevs_lc_all$value,
                                               lbracket = "(", rbracket = ")", 
                                               separator = "-"))



### Plot
p_lc_all <- prevs_lc_all %>% 
  ggplot(aes(x=Category, y=OR_sep, col=factor(sex_mw)))+
  geom_errorbar(aes(ymin=CI_lower, ymax=CI_upper), width=0.4, position=dodger) + 
  geom_point(position=dodger) + 
  scale_colour_manual(values = myCols[c(1,2)],guide=guide_legend(reverse=F)) +
  theme_bw() +
  labs(fill="", col="", x="Age",
       y="Weighted percentage of respondents"
       # subtitle = "Percentage of respondents who experienced symptoms"
  )

p_lc_all



# 3. Probability of long COVID, conditional on being symptomatic -------------


# Get table of lcs conditional on being symptomatic
prevs_lc_symptomatic <- OverReact::stratifiedTables(dat = mydat %>% filter(symptomatic==1),
                                                    result_var = "lc_84", 
                                                    covariates = vars_of_interest,
                                                    strat_var="sex_mw",
                                                    
                                                    cov_name_list = cov_name_list,
                                                    sens = 1,spec = 1, output_list = F,
                                                    weights = "wt_antibody", 
                                                    percent = T) %>% 
  pivot_longer(cols=c("Men","Women")) %>% rename(sex_mw=name)
prevs_lc_symptomatic <- cbind(prevs_lc_symptomatic,OR_unconcat(prevs_lc_symptomatic$value,
                                                               lbracket = "(", rbracket = ")", 
                                                               separator = "-"))


### Plot
p_lc_symptomatic <- prevs_lc_symptomatic %>% 
  ggplot(aes(x=Category, y=OR_sep, col=factor(sex_mw)))+
  geom_errorbar(aes(ymin=CI_lower, ymax=CI_upper), width=0.4, position=dodger) + 
  geom_point(position=dodger) + 
  scale_colour_manual(values = myCols[c(1,2)],guide=guide_legend(reverse=F)) +
  theme_bw() +
  labs(fill="", col="", x="Age",
       y="Weighted percentage of respondents"
       # subtitle = "Percentage of respondents who experienced symptoms"
  )



p_lc_symptomatic




### Combine all plots
p_comb_realdata=p_symptoms+p_lc_all+ p_lc_symptomatic + plot_layout(guides="collect")
p_comb_realdata


# Combine all plots into one giant plot -----------------------------------

p_six_panel <- lc_84_comb/p_comb_realdata + plot_layout(guides="collect")
p_six_panel

OverReact::saveREACTplot(p = p_six_panel,figpath = figpath,filename = "lc_six_panel",
                         width = 19,height = 12
)

