
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

# define outcome
outcome = "lc_84_all_29"

# Create folders
createMySubfolder("prediction")
createMySubfolder(outcome)

# Prepare data
dfRes <- makeVarsFactors(dfRes = dfRes) # convert all character variables to factors for modelling
moddat <- prepMyData(dat = dfRes %>% filter(round%in%3:5), t0 = T, min_followup = 84, 
                     abresult_filter = NULL, 
                     covidafilter = c(1,2,3),forclustering = F)
moddat_6 <- prepMyData(dat = dfRes %>% filter(round%in%6), t0 = T, min_followup = 84, 
                     abresult_filter = NULL, 
                     covidafilter = c(1,2,3),forclustering = F)

# Add hospitalised to my covs
mycovs <- c(mycovs,"hospitalised_covid")
# remove health conditions
mycovs <- setdiff(mycovs,healthnames)


# Run univariate analysis
p_univ <- runUnivariate(mydat = moddat,
                        outcome_var =  outcome,
                        negatives ="vs_all",
                        outpath_new = outpath,figpath_new = figpath,
                        myheight = 10,mywidth = 8)
p_univ
p_univ_r6 <- runUnivariate(mydat = moddat_6,outcome_var =  outcome,
                        negatives ="vs_all_r6",
                        outpath_new = outpath,
                        figpath_new = figpath,
                        myheight = 10,mywidth = 8)
p_univ_r6



# CATboost modelling -------------------------------------------------------

my_predictors=c("age","bmi","imd_score","sex","ethnic_new",
                "smokenow","vapenow","hospitalised_covid",
                "gross_household_cat",
                "work1_healthcare_or_carehome_worker"
)
df_mod=moddat %>% select(outcome,my_predictors)
# i <- 6
# replace colnames with informative oness
for(i in 1:length(colnames(df_mod))){
  print(colnames(df_mod)[[i]])
  if(colnames(df_mod)[[i]] %in% predictors_df$variable_name){
    featnm=colnames(df_mod)[[i]]
    featdesc=predictors_df[predictors_df$variable_name==featnm,]$variable_desc
    colnames(df_mod)[[i]]=featdesc
  }
}

### add dummies for ethnicity
df_mod <- df_mod %>% dummy_columns("Ethnicity",remove_selected_columns = T, ignore_na = T)

### Create test train split

# create holdout data set
set.seed(123)
split <- caret::createDataPartition(pull(df_mod,outcome),times = 1,p = 0.7, list =F)
splitindex <- 1:nrow(df_mod) %in% split
table(splitindex)
train=df_mod[splitindex,]
test=df_mod[!splitindex,]

# split training again for anti-overfitting pool
split2 <- caret::createDataPartition(pull(train,outcome),times = 1,p = 0.8, list =F)
splitindex2 <- 1:nrow(train) %in% split2
train_learn <- train[splitindex2,]
train_valid=train[!splitindex2,]



# Set parameters
ncores=6
my_params=list(thread_count=ncores,
               loss_function="Logloss",
               eval_metric = "Logloss",
               iterations =10^5,
               early_stopping_rounds=50,
               border_count=254,
               depth =7,
               learning_rate = 0.1)


X_learn=train_learn %>% select(-outcome)
X_valid=train_valid %>% select(-outcome)
X_test=test %>% select(-outcome)

# create learn and test pools
learn_pool = catboost.load_pool(data = X_learn,label = pull(train_learn,outcome), cat_features = 3:(ncol(X_learn)-1))
test_pool = catboost.load_pool(data = X_valid,label = pull(train_valid,outcome), cat_features = 3:(ncol(X_valid)-1))
validation_pool = catboost.load_pool(data =X_test,label = pull(test,outcome), cat_features = 3:(ncol(X_test)-1))


# train model
cb_model = catboost.train(learn_pool = learn_pool,test_pool = test_pool, params = my_params)



####################  ####################  ####################  ####################
####################  ####################  ####################  ####################

basemodel <- getFeatureImpAndSHAP(cb_model = cb_model,validation_pool = validation_pool,learn_pool=learn_pool,modelvars = "covariates",
                                  X = test %>% select(-outcome),
                                  y = pull(test,outcome),
                                  nperm = 1000,
                                  myheight = 6.5,bandwidth_scaler = 15,
                                  train_learn = train_learn,
                                  outcome = outcome,
                                  parallelise = T,
                                  nclust = 40)
basemodel$permutation_varimps_boxplot

# ROC analysis ------------------------------------------------------------


### Prediction evaluation ###

### FULL MODEL ###
preds <- catboost.predict(model = cb_model, pool = validation_pool,verbose = T,prediction_type = "Probability")
auc_ci_mod <- pROC::ci.auc(pull(test,outcome),preds)
auc_mod <- pROC::auc(pull(test,outcome),preds)
roc_mod <- pROC::roc(pull(test,outcome),preds)
ci.sp.obj <- pROC::ci.sp(roc_mod,sensitivities=seq(0,1,0.01), boot.n=500)

### Evaluate and plot
png(filename = paste0(figpath,"roc_",outcome,".png"), width = 9,height = 9,units = "in",res = 300)
plot(roc_mod, print.auc=F, grid = T, legacy.axes=T)
plot(ci.sp.obj, type="shape", col =  rgb(0.1,0.1,0.1, alpha = 0.2))
# plot(ci.sp.obj_symps, type = "shape", col = rgb(1,0.1,0.1, alpha = 0.2),
#      lwd = 0.01,
#      lty = 1)
# 
# plot(ci.sp.obj_clust, type = "shape", col = rgb(0.1,0.1,1, alpha = 0.2),
#      lwd = 0.01,
#      lty = 1)

lines(pROC::roc(pull(test,outcome),preds),lty=1, col = rgb(0.1,0.1,0.1, alpha = 1))
# lines(pROC::roc(pull(test,outcome),preds_symps),lty=1, col = rgb(1,0.1,0.1, alpha = 1))
# lines(pROC::roc(pull(test,outcome),preds_clust),lty=1, col = rgb(0.1,0.1,1, alpha = 1))

legend("bottomright", inset = c(0,0.05),lty=c(1,1,1), bty="n",
       col = c(col = rgb(0.1,0.1,0.1, alpha = 1), 
               rgb(0.1,0.1,1, alpha = 1),rgb(1,0.1,0.1, alpha = 1)),
       legend=c(paste0("AUC = ",
                       round(auc_ci_mod[2],2)," [",
                       round(auc_ci_mod[1],2),"-",
                       round(auc_ci_mod[3],2),"]")
                
                # paste0("Covariates + acute symptom cluster: AUC = ",
                #        round(auc_ci_mod_clust[2],2)," [",
                #        round(auc_ci_mod_clust[1],2),"-",
                #        round(auc_ci_mod_clust[3],2),"]"),
                # paste0("Covariates + acute symptoms: AUC = ",
                #        round(auc_ci_mod_symps[2],2)," [",
                #        round(auc_ci_mod_symps[1],2),"-",
                #        round(auc_ci_mod_symps[3],2),"]"))
       ))

dev.off()




# Combine into panel ------------------------------------------------------

p_or_varimp_shapmean <- p_univ + (basemodel$permutation_varimps_boxplot/ basemodel$shap_mean ) + plot_layout(widths=c(3,2))
OverReact::saveREACTplot(p = p_or_varimp_shapmean, figpath = figpath,
                         filename = paste0("OR_permute_varimps_shapmean_compare_",outcome),
                         width =14, height=9)

p_or_varimp_spacer <- p_univ + (basemodel$permutation_varimps_boxplot/ plot_spacer()) + plot_layout(widths=c(3,2))
OverReact::saveREACTplot(p = p_or_varimp_spacer, figpath = figpath,
                         filename = paste0("OR_permute_varimps_spacer_compare_",outcome),
                         width =14, height=9)


# ORs vs varimps
p_panel <- plot_spacer()/ 
  (p_univ ) / 
  plot_spacer() /
  ((basemodel$permutation_varimps_boxplot )) + plot_layout(heights=c(0.1,3.5,0.1,1))


p_panel

OverReact::saveREACTplot(p = p_panel, figpath = figpath,
                         filename = paste0("OR_permute_varimps_compare_",outcome),
                         width =10, height=13)


# ORs vs SHAP local
p_panel_shaplocal <- p_univ + basemodel$shaps + plot_layout(widths=c(2,1))
p_panel_shaplocal

OverReact::saveREACTplot(p = p_panel_shaplocal, figpath = figpath,
                         filename = paste0("panel_OR_SHAP_local_compare_",outcome),
                         width =14, height=9)

# ORs vs SHAP mean



## First general global mean shap plot
p_shap_mean <- basemodel$shaps$data %>% distinct(variable, .keep_all = T) %>% 
  ggplot(aes(x=mean_value, y=reorder(variable,mean_value))) +
  geom_col(fill = myCols[[2]])+
  labs(x="Mean SHAP value \n(variable importance in multivariable model)", y= "Predictor") +
  theme_bw()


p_panel_mean <- p_univ + p_shap_mean + plot_layout(widths=c(2,1))
p_panel_mean

OverReact::saveREACTplot(p = p_panel_mean, figpath = figpath,
                         filename = paste0("panel_OR_SHAP_global_compare_",outcome),
                         width =14, height=9)




# CATboost replication -------------------------------------------------------
createMySubfolder("replication")

my_predictors=c("age","bmi","imd_score","sex","ethnic_new",
                "smokenow","vapenow","hospitalised_covid",
                "gross_household_cat",
                "work1_healthcare_or_carehome_worker"
)
df_mod=moddat_6 %>% select(outcome,my_predictors)
# i <- 6
# replace colnames with informative oness
for(i in 1:length(colnames(df_mod))){
  print(colnames(df_mod)[[i]])
  if(colnames(df_mod)[[i]] %in% predictors_df$variable_name){
    featnm=colnames(df_mod)[[i]]
    featdesc=predictors_df[predictors_df$variable_name==featnm,]$variable_desc
    colnames(df_mod)[[i]]=featdesc
  }
}

### add dummies for ethnicity
df_mod <- df_mod %>% dummy_columns("Ethnicity",remove_selected_columns = T, ignore_na = T)

### Create test train split

# create holdout data set
set.seed(123)
split <- caret::createDataPartition(pull(df_mod,outcome),times = 1,p = 0.7, list =F)
splitindex <- 1:nrow(df_mod) %in% split
table(splitindex)
train=df_mod[splitindex,]
test=df_mod[!splitindex,]

# split training again for anti-overfitting pool
split2 <- caret::createDataPartition(pull(train,outcome),times = 1,p = 0.8, list =F)
splitindex2 <- 1:nrow(train) %in% split2
train_learn <- train[splitindex2,]
train_valid=train[!splitindex2,]



# Set parameters
ncores=6
my_params=list(thread_count=ncores,
               loss_function="Logloss",
               eval_metric = "Logloss",
               iterations =10^5,
               early_stopping_rounds=50,
               border_count=254,
               depth =7,
               learning_rate = 0.1)


X_learn=train_learn %>% select(-outcome)
X_valid=train_valid %>% select(-outcome)
X_test=test %>% select(-outcome)

# create learn and test pools
learn_pool = catboost.load_pool(data = X_learn,label = pull(train_learn,outcome), cat_features = 3:(ncol(X_learn)-1))
test_pool = catboost.load_pool(data = X_valid,label = pull(train_valid,outcome), cat_features = 3:(ncol(X_valid)-1))
validation_pool = catboost.load_pool(data =X_test,label = pull(test,outcome), cat_features = 3:(ncol(X_test)-1))


# train model
cb_model = catboost.train(learn_pool = learn_pool,test_pool = test_pool, params = my_params)



####################  ####################  ####################  ####################
####################  ####################  ####################  ####################

basemodel <- getFeatureImpAndSHAP(cb_model = cb_model,validation_pool = validation_pool,learn_pool=learn_pool,modelvars = "covariates",
                                  X = test %>% select(-outcome),
                                  y = pull(test,outcome),
                                  nperm = 1000,
                                  myheight = 6.5,bandwidth_scaler = 15,
                                  train_learn = train_learn,
                                  outcome = outcome,
                                  parallelise = T,
                                  nclust = 40)
basemodel$permutation_varimps_boxplot

# ROC analysis ------------------------------------------------------------


### Prediction evaluation ###

### FULL MODEL ###
preds <- catboost.predict(model = cb_model, pool = validation_pool,verbose = T,prediction_type = "Probability")
auc_ci_mod <- pROC::ci.auc(pull(test,outcome),preds)
auc_mod <- pROC::auc(pull(test,outcome),preds)
roc_mod <- pROC::roc(pull(test,outcome),preds)
ci.sp.obj <- pROC::ci.sp(roc_mod,sensitivities=seq(0,1,0.01), boot.n=500)

### Evaluate and plot
png(filename = paste0(figpath,"roc_",outcome,".png"), width = 9,height = 9,units = "in",res = 300)
plot(roc_mod, print.auc=F, grid = T, legacy.axes=T)
plot(ci.sp.obj, type="shape", col =  rgb(0.1,0.1,0.1, alpha = 0.2))
# plot(ci.sp.obj_symps, type = "shape", col = rgb(1,0.1,0.1, alpha = 0.2),
#      lwd = 0.01,
#      lty = 1)
# 
# plot(ci.sp.obj_clust, type = "shape", col = rgb(0.1,0.1,1, alpha = 0.2),
#      lwd = 0.01,
#      lty = 1)

lines(pROC::roc(pull(test,outcome),preds),lty=1, col = rgb(0.1,0.1,0.1, alpha = 1))
# lines(pROC::roc(pull(test,outcome),preds_symps),lty=1, col = rgb(1,0.1,0.1, alpha = 1))
# lines(pROC::roc(pull(test,outcome),preds_clust),lty=1, col = rgb(0.1,0.1,1, alpha = 1))

legend("bottomright", inset = c(0,0.05),lty=c(1,1,1), bty="n",
       col = c(col = rgb(0.1,0.1,0.1, alpha = 1), 
               rgb(0.1,0.1,1, alpha = 1),rgb(1,0.1,0.1, alpha = 1)),
       legend=c(paste0("AUC = ",
                       round(auc_ci_mod[2],2)," [",
                       round(auc_ci_mod[1],2),"-",
                       round(auc_ci_mod[3],2),"]")
                
                # paste0("Covariates + acute symptom cluster: AUC = ",
                #        round(auc_ci_mod_clust[2],2)," [",
                #        round(auc_ci_mod_clust[1],2),"-",
                #        round(auc_ci_mod_clust[3],2),"]"),
                # paste0("Covariates + acute symptoms: AUC = ",
                #        round(auc_ci_mod_symps[2],2)," [",
                #        round(auc_ci_mod_symps[1],2),"-",
                #        round(auc_ci_mod_symps[3],2),"]"))
       ))

dev.off()




# Combine into panel ------------------------------------------------------

p_or_varimp_shapmean <- p_univ_r6 + (basemodel$permutation_varimps_boxplot/ basemodel$shap_mean ) + plot_layout(widths=c(3,2))
OverReact::saveREACTplot(p = p_or_varimp_shapmean, figpath = figpath,
                         filename = paste0("OR_permute_varimps_shapmean_compare_",outcome),
                         width =14, height=9)

p_or_varimp_spacer <- p_univ_r6 + (basemodel$permutation_varimps_boxplot/ plot_spacer()) + plot_layout(widths=c(3,2))
OverReact::saveREACTplot(p = p_or_varimp_spacer, figpath = figpath,
                         filename = paste0("OR_permute_varimps_spacer_compare_",outcome),
                         width =14, height=9)


# ORs vs varimps
p_panel <- plot_spacer()/ 
  (p_univ_r6 ) / 
  plot_spacer() /
  ((basemodel$permutation_varimps_boxplot )) + plot_layout(heights=c(0.1,3.5,0.1,1))


p_panel

OverReact::saveREACTplot(p = p_panel, figpath = figpath,
                         filename = paste0("OR_permute_varimps_compare_",outcome),
                         width =10, height=13)


# ORs vs SHAP local
p_panel_shaplocal <- p_univ_r6 + basemodel$shaps + plot_layout(widths=c(2,1))
p_panel_shaplocal

OverReact::saveREACTplot(p = p_panel_shaplocal, figpath = figpath,
                         filename = paste0("panel_OR_SHAP_local_compare_",outcome),
                         width =14, height=9)

# ORs vs SHAP mean



## First general global mean shap plot
p_shap_mean <- basemodel$shaps$data %>% distinct(variable, .keep_all = T) %>% 
  ggplot(aes(x=mean_value, y=reorder(variable,mean_value))) +
  geom_col(fill = myCols[[2]])+
  labs(x="Mean SHAP value \n(variable importance in multivariable model)", y= "Predictor") +
  theme_bw()


p_panel_mean <- p_univ_r6 + p_shap_mean + plot_layout(widths=c(2,1))
p_panel_mean

OverReact::saveREACTplot(p = p_panel_mean, figpath = figpath,
                         filename = paste0("panel_OR_SHAP_global_compare_",outcome),
                         width =14, height=9)



# Reduced symptom set -----------------------------------------------------
setwd("E:/Group/react2_study5/report_phases_combined/projects/long_covid_resubmission")
outpath <- "E:/Group/react2_study5/report_phases_combined/projects/long_covid_resubmission/output/"
figpath <- "E:/Group/react2_study5/report_phases_combined/projects/long_covid_resubmission/plots/"

# define outcome
outcome = "lc_84_all_15"

# Create folders
createMySubfolder("prediction")
createMySubfolder(outcome)


dfRes$lc_84_all_15_symptom_count <- pmax(as.numeric(rowSums(dfRes[,data_key$symptom_duration_84days_binary[data_key$reduced_symptom_set == 1]], na.rm= T)),
                                         as.numeric(rowSums(dfRes[,unique(data_key$symptom_duration_84days_grouped_binary[data_key$reduced_symptom_set == 1])], 
                                                            na.rm= T)))
dfRes$lc_84_all_15 <- as.numeric(dfRes$lc_84_all_15_symptom_count >0)
dfRes$lc_84_all_15_2 <- as.numeric(dfRes$lc_84_all_15_symptom_count >1)
dfRes$lc_84_all_15_3 <- as.numeric(dfRes$lc_84_all_15_symptom_count >2)
table(dfRes$round,dfRes$lc_84_all_15)
moddat <- prepMyData(dat = dfRes %>% filter(round%in%3:5), t0 = T, min_followup = 84, 
                     abresult_filter = NULL, 
                     covidafilter = c(1,2,3),forclustering = F)
# Run univariate

# Run univariate analysis
p_univ_15 <- runUnivariate(mydat = moddat,
                        outcome_var =  outcome,
                        negatives ="vs_all_reduced",
                        outpath_new = outpath,figpath_new = figpath,
                        myheight = 10,mywidth = 7)
p_univ_15



# CATboost modelling -------------------------------------------------------

my_predictors=c("age","bmi","imd_score","sex","ethnic_new",
                "smokenow","vapenow","hospitalised_covid",
                "gross_household_cat",
                "work1_healthcare_or_carehome_worker"
)
df_mod=moddat %>% select(outcome,my_predictors)
# i <- 6
# replace colnames with informative oness
for(i in 1:length(colnames(df_mod))){
  print(colnames(df_mod)[[i]])
  if(colnames(df_mod)[[i]] %in% predictors_df$variable_name){
    featnm=colnames(df_mod)[[i]]
    featdesc=predictors_df[predictors_df$variable_name==featnm,]$variable_desc
    colnames(df_mod)[[i]]=featdesc
  }
}

### add dummies for ethnicity
df_mod <- df_mod %>% dummy_columns("Ethnicity",remove_selected_columns = T, ignore_na = T)

### Create test train split

# create holdout data set
set.seed(123)
split <- caret::createDataPartition(pull(df_mod,outcome),times = 1,p = 0.7, list =F)
splitindex <- 1:nrow(df_mod) %in% split
table(splitindex)
train=df_mod[splitindex,]
test=df_mod[!splitindex,]

# split training again for anti-overfitting pool
split2 <- caret::createDataPartition(pull(train,outcome),times = 1,p = 0.8, list =F)
splitindex2 <- 1:nrow(train) %in% split2
train_learn <- train[splitindex2,]
train_valid=train[!splitindex2,]



# Set parameters
ncores=6
my_params=list(thread_count=ncores,
               loss_function="Logloss",
               eval_metric = "Logloss",
               iterations =10^5,
               early_stopping_rounds=50,
               border_count=254,
               depth =7,
               learning_rate = 0.1)


X_learn=train_learn %>% select(-outcome)
X_valid=train_valid %>% select(-outcome)
X_test=test %>% select(-outcome)

# create learn and test pools
learn_pool = catboost.load_pool(data = X_learn,label = pull(train_learn,outcome), cat_features = 3:(ncol(X_learn)-1))
test_pool = catboost.load_pool(data = X_valid,label = pull(train_valid,outcome), cat_features = 3:(ncol(X_valid)-1))
validation_pool = catboost.load_pool(data =X_test,label = pull(test,outcome), cat_features = 3:(ncol(X_test)-1))


# train model
cb_model = catboost.train(learn_pool = learn_pool,test_pool = test_pool, params = my_params)



####################  ####################  ####################  ####################
####################  ####################  ####################  ####################

basemodel <- getFeatureImpAndSHAP(cb_model = cb_model,validation_pool = validation_pool,learn_pool=learn_pool,modelvars = "covariates",
                                  X = test %>% select(-outcome),
                                  y = pull(test,outcome),
                                  nperm = 1000,
                                  myheight = 6.5,bandwidth_scaler = 15,
                                  train_learn = train_learn,
                                  outcome = outcome,
                                  parallelise = T,
                                  nclust = 40)
basemodel$permutation_varimps_boxplot

# ROC analysis ------------------------------------------------------------


### Prediction evaluation ###

### FULL MODEL ###
preds <- catboost.predict(model = cb_model, pool = validation_pool,verbose = T,prediction_type = "Probability")
auc_ci_mod <- pROC::ci.auc(pull(test,outcome),preds)
auc_mod <- pROC::auc(pull(test,outcome),preds)
roc_mod <- pROC::roc(pull(test,outcome),preds)
ci.sp.obj <- pROC::ci.sp(roc_mod,sensitivities=seq(0,1,0.01), boot.n=500)

### Evaluate and plot
png(filename = paste0(figpath,"roc_",outcome,".png"), width = 9,height = 9,units = "in",res = 300)
plot(roc_mod, print.auc=F, grid = T, legacy.axes=T)
plot(ci.sp.obj, type="shape", col =  rgb(0.1,0.1,0.1, alpha = 0.2))
# plot(ci.sp.obj_symps, type = "shape", col = rgb(1,0.1,0.1, alpha = 0.2),
#      lwd = 0.01,
#      lty = 1)
# 
# plot(ci.sp.obj_clust, type = "shape", col = rgb(0.1,0.1,1, alpha = 0.2),
#      lwd = 0.01,
#      lty = 1)

lines(pROC::roc(pull(test,outcome),preds),lty=1, col = rgb(0.1,0.1,0.1, alpha = 1))
# lines(pROC::roc(pull(test,outcome),preds_symps),lty=1, col = rgb(1,0.1,0.1, alpha = 1))
# lines(pROC::roc(pull(test,outcome),preds_clust),lty=1, col = rgb(0.1,0.1,1, alpha = 1))

legend("bottomright", inset = c(0,0.05),lty=c(1,1,1), bty="n",
       col = c(col = rgb(0.1,0.1,0.1, alpha = 1), 
               rgb(0.1,0.1,1, alpha = 1),rgb(1,0.1,0.1, alpha = 1)),
       legend=c(paste0("AUC = ",
                       round(auc_ci_mod[2],2)," [",
                       round(auc_ci_mod[1],2),"-",
                       round(auc_ci_mod[3],2),"]")
                
                # paste0("Covariates + acute symptom cluster: AUC = ",
                #        round(auc_ci_mod_clust[2],2)," [",
                #        round(auc_ci_mod_clust[1],2),"-",
                #        round(auc_ci_mod_clust[3],2),"]"),
                # paste0("Covariates + acute symptoms: AUC = ",
                #        round(auc_ci_mod_symps[2],2)," [",
                #        round(auc_ci_mod_symps[1],2),"-",
                #        round(auc_ci_mod_symps[3],2),"]"))
       ))

dev.off()




# Combine into panel ------------------------------------------------------

p_or_varimp_shapmean <- p_univ_15 + (basemodel$permutation_varimps_boxplot/ basemodel$shap_mean ) + plot_layout(widths=c(3,2))
OverReact::saveREACTplot(p = p_or_varimp_shapmean, figpath = figpath,
                         filename = paste0("OR_permute_varimps_shapmean_compare_",outcome),
                         width =14, height=9)

p_or_varimp_spacer <- p_univ_15 + (basemodel$permutation_varimps_boxplot/ plot_spacer()) + plot_layout(widths=c(3,2))
OverReact::saveREACTplot(p = p_or_varimp_spacer, figpath = figpath,
                         filename = paste0("OR_permute_varimps_spacer_compare_",outcome),
                         width =14, height=9)


# ORs vs varimps
p_panel <- plot_spacer()/ 
  (p_univ_15 ) / 
  plot_spacer() /
  ((basemodel$permutation_varimps_boxplot )) + plot_layout(heights=c(0.1,3.5,0.1,1))


p_panel

OverReact::saveREACTplot(p = p_panel, figpath = figpath,
                         filename = paste0("OR_permute_varimps_compare_",outcome),
                         width =10, height=13)


# ORs vs SHAP local
p_panel_shaplocal <- p_univ_15 + basemodel$shaps + plot_layout(widths=c(2,1))
p_panel_shaplocal

OverReact::saveREACTplot(p = p_panel_shaplocal, figpath = figpath,
                         filename = paste0("panel_OR_SHAP_local_compare_",outcome),
                         width =14, height=9)

# ORs vs SHAP mean



## First general global mean shap plot
p_shap_mean <- basemodel$shaps$data %>% distinct(variable, .keep_all = T) %>% 
  ggplot(aes(x=mean_value, y=reorder(variable,mean_value))) +
  geom_col(fill = myCols[[2]])+
  labs(x="Mean SHAP value \n(variable importance in multivariable model)", y= "Predictor") +
  theme_bw()


p_panel_mean <- p_univ_15 + p_shap_mean + plot_layout(widths=c(2,1))
p_panel_mean

OverReact::saveREACTplot(p = p_panel_mean, figpath = figpath,
                         filename = paste0("panel_OR_SHAP_global_compare_",outcome),
                         width =14, height=9)



# Checking vaccination breakthroughs --------------------------------------


dfRes %>% 
  mutate(breakthrough_oneVax=as.numeric(as.numeric(as.Date(vaccinefirst, format = "%m/%d/%Y")))<
           as.numeric(as.Date(covidsta, format = "%m/%d/%Y")),
         breakthrough_twoVax=as.numeric(as.numeric(as.Date(vaccinesecond, format = "%m/%d/%Y")))<as.numeric(as.Date(covidsta, format = "%m/%d/%Y"))) %>%
  # mutate(second_vax_infection_days=(as.Date(vaccinesecond)-as.Date(covidsta))) %>% 
  # filter(!is.na(vaccinefirst), symptomatic==1) %>%
  group_by(round) %>% 
  summarise(n=n(),
            # n_12_weeks = sum(follow_up_time >= 84, na.rm=T),
            n_breakthrough_1vax = sum(breakthrough_oneVax, na.rm=T),
            n_breakthrough_2vax = sum(breakthrough_twoVax, na.rm=T),
            n_breakthrough_plus12weeks = sum(breakthrough_twoVax == 1 & follow_up_time >=84, na.rm=T))



