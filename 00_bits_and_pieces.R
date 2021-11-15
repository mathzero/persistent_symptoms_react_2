# Miscellaneous bits and pieces

cov_name_list$imd_quintile_cat
symptom_durations <- paste0("symptom_duration_",1:29)
acute_symptoms <- names(sympnames_list)[c(1:26,28:30)]

### Create ordered sympnames list for loop
sympnames_list_reordered <- sympnames_list[c(1:26,28:30)]

healthnames <- c(paste0("healtha_0",c(1:4,6:9)),
                 paste0("healtha_1",0:7))
predictors_df <- data.frame(variable_name=c("age","height_cm","bmi", "bmi_cat","hh_size",
                                            names(cov_name_list[c(1:6,8:10,13:18,20)]),
                                            healthnames,acute_symptoms,"hospitalised_covid",
                                            "cluster","work1_healthcare_or_carehome_worker","imd_quintile_cat",
                                            "imd_score","smokenow","vapenow"),
                            variable_desc=c("Age","Height (cm)","BMI","BMI (categorical)","Household size",
                                            unlist(cov_name_list[c(1:6,8:10,13:18,20)]),
                                            health_conditions[1:16],
                                            as.character(unlist(sympnames_list[c(1:26,28:30)])),
                                            "Hospitalised with COVID",
                                            "Acute symptom cluster",
                                            "Healthcare or care home worker",
                                            "Index of multiple deprivation (IMD) quintile",
                                            "Index of multiple deprivation (IMD)",
                                            "Smoking", "Vaping"

                            ),
                            variable_type = c(rep("Biological",4),
                                              rep("Sociodeographic",1),
                                              rep("Biological",4),
                                              rep("Sociodeographic",5),
                                              rep("COVID history",3),
                                              rep("Sociodeographic",4),
                                              rep("Health",16),
                                              rep("Acute symptoms",29),
                                              "COVID history","COVID history",
                                              rep("Sociodeographic",5)
                                              
                            )
)

### Save as a list
my_new_list <- as.list(predictors_df$variable_desc)
names(my_new_list) <- predictors_df$variable_name


### Subsets of symptoms, per Helen's research
symptom_subset_indexes <- c(1,2,3,4,9,10,
                            13,14,15,16,
                            17,18,19,20,
                            21,22,23,24,26)

symptom_subset<- sympnames_list[symptom_subset_indexes]


### get indexes of selected symptoms (this is from univariate analysis)
symptom_subset_indexes_datadriven <- c(1, 9, 10, 13, 14, 15, 17, 18, 21, 22, 23, 24, 25, 26, 28)
symptom_subset_datadriven<- sympnames_list[symptom_subset_indexes_datadriven]



## Long covid binary symptom names
lc_names_12_weeks <- c(paste0("covidsym_12_week_0",1:9),paste0("covidsym_12_week_",10:29))


### Covs for tables and figs
mycovs <- c("dummy", "sex", "age_group_named","ethnic_new","bmi_cat", 
            "work1_healthcare_or_carehome_worker", "imd_quintile_cat",
            "smokenow","vapenow","covidc_cat","covid_severity",
            # "hospitalised_covid", 
            "gross_household_cat", 
            "res",
             healthnames)
cov_name_list[healthnames] <-health_conditions[c(1:16)]
cov_name_list$work1_healthcare_or_carehome_worker <- "Healthcare or care home worker"
# cov_name_list$res = "Lateral flow test result (antibody)"
cov_name_list$hospitalised_covid <- "Hospitalised with COVID"



tableorder <-c("Symptomatic", "All participants",  "Male", "Female", "18-24", "25-34", "35-44", "45-54", "55-64", "65-74", "74+",
               "Asian", "Black", "Mixed", "Other", "White", 
               "Underweight", "Normal weight", "Overweight","Obese", 
               "Yes", "No", "1 - most deprived", "2", "3", "4", "5 - least deprived", 
               "Current cigarette smoker", "Not current cigarette smoker", 
               "Current vaper", "Not current vaper", "Prefer not to say", 
               "Lateral flow test result (antibody)",
               "No symptoms", "Mild symptoms", "Moderate symptoms", "Severe symptoms", 
               "No medical attention sought", "Sought medical attention from pharmacist / by phone (NHS 111/GP)", 
               "Visited GP/walk-in centre/A&E", "Hospital admission",
               "0-14,999", "15,000-49,999", "50,000-149,999", ">150,000",
               "Organ transplant recipient", "Diabetes (type I or II)", "Heart disease or heart problems", 
               "Hypertension", "Stroke", "Kidney disease", "Liver disease", "Anemia", "Asthma", 
               "Other lung condition", "Cancer", "Neurological condition", "Immunocompromised*", 
               "Depression", "Anxiety", "Psychiatric disorder (other)")
