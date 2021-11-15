
#' First clear the environment of variables
rm(list=ls(all=TRUE))
setwd("E:/Group/react2_study5/report_phases_combined/projects/long_covid_resubmission//")
outpath <- "E:/Group/react2_study5/report_phases_combined/projects/long_covid_resubmission/output/"
figpath <- "E:/Group/react2_study5/report_phases_combined/projects/long_covid_resubmission/plots/"

source("E:/Group/functions/load_packages.R", local = T)
source("E:/Group/functions/full_join_multiple_datasets.R", local = T)
source("E:/Group/functions/wrangle_cleaner_functions.R", local = T)
source("E:/Group/functions/cats_and_covs.R", local = T)
source("E:/Group/react2_study5/report_phases_combined/projects/delta_symptom_prediction/code/00_bits_and_pieces.R", local = T)
source("E:/Group/react2_study5/report_phases_combined/projects/delta_symptom_prediction/code/00_functions.R", local = T)
source("E:/Group/react2_study5/report_phases_combined/projects/symptom_prediction_children/code/00_bits_and_pieces.R", local = T)
source("E:/Group/react2_study5/report_phases_combined/projects/long_covid_resubmission/code/00_functions.R")


#' Pull in packages needed
package.list <- c("prevalence","mgcv","knitr","MASS","kableExtra","table1","dplyr","factoextra","tableone","networkD3",
                  "tidyr","forcats", "cluster", "fpc", "mclust", "pheatmap","FactoMineR", "NbClust","clValid","plotly",
                  "ggplot2","ggsci", "RColorBrewer", "tidyverse", "lubridate", "egg", "poLCA", "Rcpp","xml2","splitstackshape",
                  "fs", "later", "promises","proxy","dendextend", "ComplexHeatmap","circlize","doSNOW","ClusterR","htmlwidgets",
                  "readr","ggthemes", "questionr", "gridExtra", "foreach", "doParallel","withr","patchwork",
                  "purrr", "httr", "htmltools","ggalluvial","datapasta","xgboost","SHAPforxgboost", "focus","ICoLour"
)

load_packages(package.list)


# Import REACT-1 data -----------------------------------------------------

dfRes <- readRDS("E:/Group/react2_study5/saved_objects/react_1_saved_objects/react1_r1_to_r14.rds")
dfRes <- dfRes %>% filter(round!=1)
as.Date(min(dfRes$d_col, na.rm=T), origin="1970-01-01")
as.Date(max(dfRes$d_col, na.rm=T), origin="1970-01-01")
createMySubfolder("react1_background_symptom_prevalence")


# Symptomaticness by age group --------------------------------------------
xtab_symptomatic_by_age <- OverReact::crossTabMulti(dat = dfRes %>% filter(estbinres==1, age >=18),
                                                    rowvar_list = "age_group_char",
                                                    colvar = "symptomatic")
# reorder rows
xtab_symptomatic_by_age <- xtab_symptomatic_by_age[c(4,5,1,6,7,8),c(1,2,5,3,4)]

OverReact::saveREACTtable(tab = xtab_symptomatic_by_age,outpath = outpath,
                          filename = "xtab_symptomatics_by_age_react1")



# 11 days or more ---------------------------------------------------------
dfRes$symptfirst_
# Clean up all sympfirst variables
dfRes <- dfRes %>%
  mutate_at(vars(matches("symptfirst_")), binaryCleaner_1_0 ) 

# Loop adding new variables for new 11+ day vars
for (i in 1:26){
  print(i)
  newvar <- case_when(dfRes$symptst ==11 & dfRes[,sympnames_type_df$first_symptom_code[[i]]] == 1 ~ 1,
                      TRUE ~ 0)
  newvarname = paste0(sympnames_type_df$first_symptom_code[[i]], "_11plus_days")
  dfRes[,newvarname] <- newvar
}

## add 11plus days var to symptom df
sympnames_type_df$first_symptom_code_11plus <- paste0(sympnames_type_df$first_symptom_code, "_11plus_days")

dfRes$symptom_count_11plus_days <- as.numeric(rowSums(dfRes[,sympnames_type_df$first_symptom_code_11plus]))
dfRes$any_symptom_11plus_days <- as.numeric(dfRes$symptom_count_11plus_days>=1)
cov_name_list[sympnames_type_df$first_symptom_code_11plus] <- sympnames_type_df$symptom
cov_name_list$any_symptom_11plus_days <- "Any symptom"
# dir.create(path = paste0(figpath,"11_day_plus_symptoms/"))

# Function to create plot


# Analysis of prevalence of each symptom in test negatives by round -------
createMyPlot <- function(dat, filename, 
                         rowvar_list = c("any_symptom",
                                         sympnames_type_df$symptom_code),
                         rowvar_cats = c("Any symptom",sympnames_type_df$symptom_type)){
  
  xtab_symps_negs_all <- OverReact::crossTabMulti(dat =dat,
                                                  rowvar_list = rowvar_list,
                                                  colvar = "round",
                                                  cov_names = cov_name_list,confint = F, include_percentages = T,
                                                  rowwise_precentages = F)
  
  xtab_symps_negs_all <- xtab_symps_negs_all %>% filter(Category ==1)
  
  # replace with percentages only
  xtab_symps_negs_all[,3:15] <- apply(xtab_symps_negs_all[,3:15],MARGIN = 2,FUN = xtabPercentageExtractor)
  
  # Drop unneccesary columns
  xtab_symps_negs_all <- xtab_symps_negs_all %>% select(-Category, -`Sum / mean(SD)`)
  
  # Add symptom type
  xtab_symps_negs_all$Symptom_cat=rowvar_cats
  
  # rename cols by round
  # names(xtab_symps_negs_all)[2:13] <- round_date_midpoints
  
  # Get max prevalence
  maxprev=ceiling(max(xtab_symps_negs_all[,2:14]))
  
  # Plot
  myplots <- list()
  # palette_list <- c("pink", "green", "red", "blue","orange_green_imperial_blue","brick_navy_teal")
  types_list <- unique(rowvar_cats)
  for(i in 1:length(types_list)){
    type <- types_list[[i]]
    p <- xtab_symps_negs_all %>% pivot_longer(cols = 2:14) %>% 
      filter(Symptom_cat==type) %>% 
      rename(round = name) %>% 
      mutate(round=as.numeric(round),
             Variable=factor(Variable, levels =unique(Variable))) %>% 
      left_join(round_date_midpoints_df) %>% 
      ggplot(aes(x=round_date, y=value, col = Variable)) +
      geom_line() +
      geom_point() +
      scale_color_imperial(palette = "orange_green_imperial_blue") +
      theme_bw() +
      theme_adjust +
      ylim(c(0,maxprev)) +
      theme(axis.text.x = element_text(angle=45, vjust =0.5, hjust =0.5))+
      facet_wrap(.~Symptom_cat, ncol = 1)+
      labs(x="Date", y= "Percentage with symptom", col = "")
    p
    if(i<length(types_list)){
      p <- p+theme(axis.text.x = element_blank(),
                   axis.title.x = element_blank()) 
      
    }
    myplots[[i]] <- p
  }
  
  
  
  p_comb <- wrap_plots(myplots,ncol = 1)& theme(legend.justification = "left")
  OverReact::saveREACTplot(p = p_comb,figpath = figpath,filename = filename,
                           width = 11, height =15
  )
  return(p_comb)
}





# Run in adults
p_neg_11day_adults <- createMyPlot(dat = dfRes %>% filter(estbinres==0, age >=18),
                                   rowvar_list = c("any_symptom_11plus_days",
                                                   sympnames_type_df$first_symptom_code_11plus),
                                   filename = "11_day_plus_symptoms/test_negatives_11day_symptoms_age18_plus")
p_pos_11day_adults <- createMyPlot(dat = dfRes %>% filter(estbinres==1, age >=18),
                                   rowvar_list = c("any_symptom_11plus_days",
                                                   sympnames_type_df$first_symptom_code_11plus),
                                   filename = "11_day_plus_symptoms/test_positives_11day_symptoms_age18_plus")



# Overall plot with CIs ---------------------------------------------------

dfRes$dummy="All"
# PLot of % symptomatic by round among negatives ------------------------------------------
dfRes$round <- as.character(dfRes$round)
xtab_symptomatic_adults_negs=OverReact::makeTablesNew(dat = dfRes %>% filter(age>=18, estbinres==0),
                                                      result_var = "any_symptom_11plus_days",covariates = "round",
                                                      cov_name_list =cov_name_list,weights = "wt_antigen",
                                                      output_list = F,sens = 1,spec = 1)
# Add ages
xtab_symptomatic_adults_negs$date <- round_date_midpoints

# combine and rename
xtab_symptomatic_comb_negs=xtab_symptomatic_adults_negs %>% 
  rename(Positive = Total,
         Symptomatic = Positive,
         Prop_symptomatic = Prevalence,
         Round = Category)


# Plot
pd=position_dodge(width=20)
maxdate="08-10-2021"

p_prop_negatives=xtab_symptomatic_comb_negs %>% 
  mutate(Round=factor(Round, levels=2:13),
         date=as.Date(date,format = "%d-%m-%Y")) %>% 
  ggplot(aes(x=date, y= Prop_symptomatic)) + 
  geom_point(position = pd) +
  geom_errorbar(aes(x=date, y= Prop_symptomatic,ymin =Lower, ymax = Upper), 
                width=10,
                position = pd) +
  scale_x_date(date_breaks = "2 month",date_labels = "%B \n%Y") +
  theme_bw() + 
  theme_adjust +
  scale_color_manual(values = myCols[c(5,1)]) +
  # annotate("rect",xmin=0.5, xmax=6.5, ymin=0, ymax=Inf, alpha=0.1) +
  annotate("rect",xmin=as.Date("23-12-2020",format = "%d-%m-%Y"), 
           xmax=as.Date("9-5-2021",format = "%d-%m-%Y"),
           ymin=0, ymax=Inf, alpha=0.1) +
  annotate("rect",xmin=as.Date("9-5-2021",format = "%d-%m-%Y"), 
           xmax=as.Date(maxdate,format = "%d-%m-%Y"), 
           ymin=0, ymax=Inf, alpha=0.17)+
  annotate("text",x=as.Date("23-08-2020",format = "%d-%m-%Y"), y=5, label = "Wild type") +
  annotate("text",x=as.Date("13-02-2021",format = "%d-%m-%Y"), y=5, label = "Alpha")  +
  annotate("text",x=as.Date("03-07-2021",format = "%d-%m-%Y"), y=5, label = "Delta")  +
  labs(x="Date", y="Weighted prevalence of persistent symptoms \nfor 11+ days among PCR negatives in REACT-1", col = "")
p_prop_negatives

OverReact::saveREACTplot(p = p_prop_negatives,figpath = figpath,filename = "background_symptom_prevalence_11plus_days_weighted",
                         width = 7, height =5
)

