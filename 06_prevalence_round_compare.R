
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
source("E:/Group/functions/prev_difference_test.R")

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
                  "foreach","doParallel","catboost","textTools",
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
data_key$short_name <- paste0(data_key$short_name,c("","*","*","*","**","","**",
                                                    rep("", 3),
                                                    rep("***", 2),
                                                    rep("", 4),
                                                    rep("****", 2),
                                                    rep("", 21)
                                                    ))

# Get Ids to exclude
ids_varselection <- readRDS("E:/Group/react2_study5/saved_objects/long_covid_clustering/ids_varselection.rds")


# define outcome
outcome = "lc_84_all_29"

createMySubfolder("prevalence_compare_rounds")


# Main analysis -----------------------------------------------------------


# Add variable that shows whether data is in main or valudation data
dfRes <- dfRes %>% mutate(data_set = case_when(sample_type == "MAIN" ~ "Replication (May 2021)",
                                               round %in% 3:5 ~ "Main analysis (Sep 2020-Feb 2021)",
                                               TRUE ~ NA_character_
))

# Get symptom codes and names
symptom_codes <- c(data_key$symptom_duration_84days_binary[which(data_key$one_of_37==1)])
symptom_names <- c(data_key$short_name[which(data_key$one_of_37==1)])



# Build prevalence table
prevs_symptoms <- OverReact::crossTabMulti(dat = dfRes %>% filter(follow_up_time >= 84, 
                                                                  covida>=2, 
                                                                  round >= 3,
                                                                  !sample_type%in% c("Boost_6574","Boost_5564")),
                                           rowvar_list =  symptom_codes,
                                           colvar = "round", 
                                           cov_names = NULL,
                                           include_percentages = T,rowwise_precentages = F)

# Wrangle
prevs_symptoms <- prevs_symptoms %>% filter(Category==1)
prevs_symptoms$Variable <- symptom_names

# save
OverReact::saveREACTtable(tab = prevs_symptoms,outpath = outpath,filename = "lc_symptom_prevalences_by_round")


# Prevs for analysis / validation compare ---------------------------------



# Build prevalence table
prevs_symptoms <- OverReact::crossTabMulti(dat = dfRes %>% filter(follow_up_time >= 84, 
                                                                  covida>=2, 
                                                                  round >= 3,
                                                                  !sample_type%in% c("Boost_6574","Boost_5564")),
                                           rowvar_list = symptom_codes,
                                           colvar = "data_set", 
                                           cov_names = NULL,
                                           # weights = "wt_antibody",
                                           include_percentages = T,rowwise_precentages = F)

# Wrangle
prevs_symptoms <- prevs_symptoms %>% filter(Category==1)
prevs_symptoms$Variable <- symptom_names

# save
OverReact::saveREACTtable(tab = prevs_symptoms,outpath = outpath,filename = "lc_symptom_prevalences_by_main_validation")






# Plot --------------------------------------------------------------------


#Mutate data
dfRes_long <-dfRes %>% filter(follow_up_time >= 84, 
                                   covida>=2, 
                                   !is.na(data_set)) %>% 
  select(all_of(symptom_codes),data_set) %>% 
  pivot_longer(cols = symptom_codes) %>% 
  mutate(value = case_when(data_set =="Main analysis (Sep 2020-Feb 2021)" & 
                             name=="lc_84_thrombosis" ~ NA_real_,
                           TRUE ~ value))



# run stratified table
prev_tabs <- OverReact::stratifiedTables(dat = dfRes_long,
                                         result_var = "value",
                                         covariates = "name",
                                         strat_var = "data_set",
                                         cov_name_list = cov_name_list,
                                         sens = 1,spec = 1)

# # Replace names 
# # prev_tabs$Category <- str_sub(prev_tabs$Category,start = 7,end = 100)
# prev_tabs <-prev_tabs %>%  left_join(data_key, by =c("Category" = "short_name_combined"))
# 
# prev_tabs$Category <- prev_tabs$R6_Symptom_description
# prev_tabs <- prev_tabs[,1:4]

# Convert to plottable table
prev_tabs_plot <- makeStratifiedPrevalenceTablePlottable(tab = prev_tabs)




# Add proportion change ---------------------------------------------------

# Build prevalence table
prevs_symptoms_raw <- OverReact::crossTabMulti(dat = dfRes %>% filter(follow_up_time >= 84, 
                                                                           covida>=2, 
                                                                           round >= 3,
                                                                           !sample_type%in% c("Boost_6574","Boost_5564")) %>% 
                                                 mutate(lc_84_thrombosis = case_when(round != 6 ~ NA_real_,
                                                                                     TRUE ~ lc_84_thrombosis)),
                                               rowvar_list = symptom_codes,
                                               colvar = "data_set", 
                                               cov_names = NULL,
                                               include_percentages = F,rowwise_precentages = F)

prevs_symptoms_raw_wide <- prevs_symptoms_raw %>% select(-`Sum / mean(SD)`) %>% 
  pivot_wider(names_from = Category,
              values_from = c(`Main analysis (Sep 2020-Feb 2021)`, `Replication (May 2021)`))

diff_cis <- getPvalforPrevDiff(n1 = prevs_symptoms_raw_wide$`Main analysis (Sep 2020-Feb 2021)_0`+prevs_symptoms_raw_wide$`Main analysis (Sep 2020-Feb 2021)_1`,
                               npos1 = prevs_symptoms_raw_wide$`Main analysis (Sep 2020-Feb 2021)_1`,
                               n2 = prevs_symptoms_raw_wide$`Replication (May 2021)_0`+ prevs_symptoms_raw_wide$`Replication (May 2021)_1`,
                               npos2 = prevs_symptoms_raw_wide$`Replication (May 2021)_1`,
                               sens = 1,spec = 1
) %>% as.data.frame()

diff_cis$symptom <- symptom_codes
diff_cis$symptom_name <- symptom_names


# Bind with original df to make plotting easier
prev_tabs_plot <- prev_tabs_plot %>% left_join(diff_cis, by = c("Category"="symptom"))



# Now main plots ----------------------------------------------------------
prev_tabs_plot$name <- factor(prev_tabs_plot$name,levels = c("Main analysis (Sep 2020-Feb 2021)",
                                                             "Replication (May 2021)"))

prev_tabs_plot$significant_delta = case_when(prev_tabs_plot$P_val <= 0.05 ~ "Significant",
                                             TRUE ~ "Not significant")

prev_tabs_plot$delta = factor(case_when(prev_tabs_plot$Difference_in_prevalence>=0 & prev_tabs_plot$P_val <= 0.05~"Increase",
                                 prev_tabs_plot$Difference_in_prevalence<0 & prev_tabs_plot$P_val <= 0.05~"Decrease",
                                 TRUE ~ "No significant change"),
                              levels = c("Decrease","No significant change","Increase"))

# figure out a variable for which rounds certain symptoms are included in
prev_tabs_plot$symptom_present <- factor(case_when(!is.nan(prev_tabs_plot$P_val) ~ "All rounds",
                                            prev_tabs_plot$Category %in% 
                                              c("lc_84_chills","lc_84_heavy_arms_legs") ~ "Rounds 3-5 only",
                                            TRUE ~ "Round 6 only"),
                                         levels=c("All rounds","Rounds 3-5 only","Round 6 only"))

prev_tabs_plot$name <- factor(prev_tabs_plot$name, levels=c("Main analysis (Sep 2020-Feb 2021)","Replication (May 2021)"))

# Plot
dodger = position_dodge2(width=0.9,reverse = T)

p_prev <- prev_tabs_plot %>% ggplot(aes(x=reorder(symptom_name,-Difference_in_prevalence),
                                        y=prev, col = name, 
                                        fill = name)) +
  geom_col(position = dodger, width=0.7, col = "black", size=0.1) +
  geom_errorbar(aes(ymin=lower, ymax = upper), size=0.3, width=0.7,position = dodger,color = "grey20") +
  theme_bw() +
  coord_flip() +
  scale_fill_manual(values = c(myCols[c(3)], "grey90"))+
  scale_color_manual(values = myCols[c(2,5)])+
  labs(x="", y="Prevalence of each symptom \nfor 12 weeks or more", fill ="", col = "")+
  ggforce::facet_col(.~symptom_present, scales = "free_y",space = "free") +
  theme(legend.position = "top",
        legend.text.align = 0,
        legend.justification  = "left",
        axis.title.x = element_text(size=9))+ 
  theme_adjust +
  guides(fill=guide_legend(nrow=2, byrow = T))
p_prev

p_diff <- prev_tabs_plot %>% filter(name=="Replication (May 2021)") %>% 
  ggplot(aes(x=(reorder(symptom_name,-Difference_in_prevalence)),y=Difference_in_prevalence, 
             fill = delta, alpha=significant_delta)) +
  geom_col(col = "black", size=0.1) +
  geom_errorbar(aes(ymin=Lower_CI, ymax = Upper_CI), size=0.4, width=0.5,color = "grey20") +
  theme_bw() +
  coord_flip() +
  scale_alpha_manual(values =c(0.5,1), guide= F)+
  scale_fill_manual(values = c(myCols[[4]],myCols[[5]],myCols[[7]])) +
  labs(x="", y="Change in symptom prevalence \nfrom rounds 3-5 to round 6", 
       fill ="", 
       col = "") +
  ggforce::facet_col(.~symptom_present, scales = "free_y",space = "free") +
  theme(axis.text.y = element_blank(),
        legend.position = "top",
        legend.text.align = 0,
        legend.justification = "left",
        axis.title.x = element_text(size=9)) + 
  theme_adjust +
  guides(fill=guide_legend(nrow=3, byrow = T))



p_comb <- p_prev+ p_diff +plot_layout(widths=c(4,3))
p_comb


# Now you can save the plot as both png and pdf simultaneously
OverReact::saveREACTplot(p = p_comb,figpath = figpath,filename = "symptom_prevalence_and_delta",
                         width = 8, height = 10)


