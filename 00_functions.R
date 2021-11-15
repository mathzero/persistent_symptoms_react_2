
# Create data frame of symptoms -------------------------------------------

symptom_df <- data.frame(symp_name =names(sympnames_list[c(1:26, 28:30)]), 
                         symp_duration_name =paste0("symptom_duration_",1:29),
                         symp_desc =unlist(sympnames_list[c(1:26, 28:30)]),
                         col_cluster=NA,
                         symp_type=c(rep("Gastrointestinal",4),
                                     rep("Flu / fever",4),
                                     rep("Anosmia / ageusia",2),
                                     rep("Upper respiratory",2),
                                     rep("Flu / fever",2),
                                     rep("Upper respiratory",4),
                                     rep("Flu / fever",2),
                                     rep("Fatigue",3),
                                     rep("Flu / fever",3),
                                     rep("Dermatological",3)))



# Create new folder and change outpaths -----------------------------------

createMySubfolder <- function(subfolderName="example"){
  
  outpath<<-paste0(outpath,subfolderName,"/")
  figpath<<-paste0(figpath,subfolderName,"/")
  
  dir.create(outpath,showWarnings = F)
  dir.create(figpath,showWarnings = F)
  
}


# Prepare data for clustering ---------------------------------------------

prepMyData <- function(dat, t0 = T, min_followup = 84, abresult_filter = 1, covidafilter=NULL,
                       forclustering=T, reduced_symptom_set=NULL,
                       lc_names_0_weeks = c(paste0("covidsym_0",1:9),paste0("covidsym_",10:26),paste0("covidsym_",28:30)),
                       lc_names_12_weeks = c(paste0("covidsym_12_week_0",1:9),paste0("covidsym_12_week_",10:29))){
  
  
  
  
  subset_index_15 <- reduced_symptom_set
  
  if(!is.null(reduced_symptom_set)){
    lc_names_0_weeks <- lc_names_0_weeks[subset_index_15]
    lc_names_12_weeks <- lc_names_12_weeks[subset_index_15]
  }
  

  if(!is.null(min_followup)){
    dat <- dat %>% filter(test_date - covidsta >=min_followup) # minimum 12 weeks' follow-up
  }
  
  
  if(!is.null(abresult_filter)){
    dat <- dat %>% filter(res %in% abresult_filter) 
  }
  
  if(!is.null(covidafilter)){
    dat <- dat %>% filter(covida %in% covidafilter) 
  }
  
  if(forclustering){
    if(t0){
      dat <- dat %>% 
        dplyr::select(all_of(lc_names_0_weeks), u_passcode) %>% 
        column_to_rownames("u_passcode")
    }else{
      dat <- dat %>% 
        filter(lc_84==1) %>% 
        dplyr::select(all_of(lc_names_12_weeks), u_passcode) %>% 
        column_to_rownames("u_passcode")
    }
  }
  
  return(dat)
  
}


### Choose study population
chooseMyStudyPop <- function(subset="SuspectedCOVID"){
  out=list()
  if(subset=="SuspectedCOVID"){
    out$outpath="SuspectedCOVID"
    out$abresult_filter = NULL
    out$covidafilter = c(1,2,3)
    out$sampsize_t0 = 0.05
    out$sampsize_t12 = 0.1
    out$reduced_symptom_set = NULL
    
  }else if(subset=="AbPos"){
    out$outpath="AbPos"
    out$abresult_filter = 1
    out$covidafilter = NULL
    out$sampsize_t0 = 0.1
    out$sampsize_t12 = 0.2
    out$reduced_symptom_set = NULL
    
  }else if(subset=="subset15"){
    out$outpath="subset15"
    out$abresult_filter = NULL
    out$covidafilter = c(1,2,3)
    out$sampsize_t0 = 0.05
    out$sampsize_t12 = 0.1
    out$reduced_symptom_set = c(1, 9, 10, 13, 14, 15, 17, 18, 21, 22, 23, 24, 25, 26, 27)
  }
  return(out)
}




# Create Heatmap plot -----------------------------------------------------

createMyHeatmap <- function(rowclust, cluster_prefix="A",figname="example",cluster_rows=F,show_row_dend=F,
                            roworder,rowsplit, save_heatmap=T, symptoms=sympnames[c(1:26, 28:30)],
                            reverse_cluster_order=F){
  
  
  clus_counts <- table(unlist(rowclust$cluster_results$pamobject$clustering))
  result_df <- data.frame(symptom=symptoms)
  sums_df <- data.frame(symptom=symptoms)
  
  for(i in 1:length(clus_counts)){
    if(reverse_cluster_order){
      x=3-i
    }else{
      x=i
    }
    result_df <- cbind(result_df,round(100*colMeans(rowclust$data[rowclust$cluster_results$pamobject$clustering == x,]),1))
    sums_df <- cbind(sums_df,round(colSums(rowclust$data[rowclust$cluster_results$pamobject$clustering == x,]),1))
    colnames(result_df)[1+i]=paste0("Cluster ",cluster_prefix,i,"\nn=",clus_counts[[x]])
  }
  sums_df$total <- rowSums(sums_df[,2:ncol(sums_df)])
  
  
  
  
  ### Add rownames
  rownames(result_df) <- symptoms
  
  col_fun=colorRamp2(breaks = c(0,50,100), colors = c("white",myCols[[2]],myCols[[1]]))
  ### Create heatmap
  my.complex.heatmap.prevs <- Heatmap(result_df[,2:ncol(result_df)],
                                      name = paste0("Symptom \nprevalence \n",figname),
                                      row_gap = unit(5,"mm"),
                                      left_annotation = rowAnnotation(
                                        Total =anno_barplot(sums_df$total, gp=gpar(fill=myCols[[1]]),border=F,
                                                            axis_param = list(direction = "reverse"))),
                                      column_order = colnames(result_df)[c(2:ncol(result_df))],
                                      row_order = roworder,
                                      # row_split = rowsplit,
                                      cluster_row_slices = F,
                                      cluster_rows=cluster_rows,
                                      show_row_dend = show_row_dend,
                                      heatmap_legend_param = list(at=c(0,50,100), col_fun = col_fun),
                                      
                                      # column_split = t0_cluster_cols$cluster_results$pamobject$clustering,
                                      column_gap = unit(5,"mm"),
                                      rect_gp = gpar(col = "white", lwd = 1),
                                      cell_fun = function(j,i,x,y,width, height, fill){
                                        grid.text(sprintf("%.1f", result_df[,2:ncol(result_df)][i,j]), 
                                                  x, y, gp=gpar(fontsize=10, col = 
                                                                  ifelse(result_df[,2:ncol(result_df)][i,j]>50,"white","black")))
                                      },
                                      col = col_fun
                                      )
  
  if(save_heatmap){
    my.complex.heatmap.prevs
    png(paste0(figpath,figname,".png"), 
        width = 5,
        height = 9, units = "in", res = 300)
    my.complex.heatmap.prevs
    dev.off()
  }
 return(list(heatmap=my.complex.heatmap.prevs,
             result_df=result_df))
}



# Long COVID wrangling ---------------------------------------------------------------


### Get symptom count at X days
xDaySymptomCount <- function(mydat,ndays){
  lc_symp_durations <- paste0("symptom_duration_",1:29)
  symp_count_12_weeks <- rowSums(mydat[,lc_symp_durations] >= ndays, na.rm = T)
  return(symp_count_12_weeks)
}



  
  



#' This function takes either a list of long COVID symptom names or a count of symtoms and
#' returns a vector of the number of days that each person in the data has experienced 
#' at least one of the named symptoms, or at least num_of_symps symptoms


engineerMultipleSymptomDaycounts <- function(dat,names_of_symps=NULL, num_of_symps=NULL,
                                             symptom_df){
  if(!is.null(names_of_symps)){
    selected_symptom_durations <- symptom_df[names_of_symps,]$symp_duration_name
    symptomatic_at_infection=(rowSums(dat[,names_of_symps],na.rm = T)>0)
    sympdurs <- matrixStats::rowMaxs(as.matrix(dat[,selected_symptom_durations]), na.rm = T)
    sympdurs[symptomatic_at_infection] <- pmax(sympdurs[symptomatic_at_infection],1)
    sympdurs <- continuousCleaner(sympdurs,floornumber = 0)
  }else if(!is.null(num_of_symps)){
    sympdurs<- apply(dat[,symptom_df$symp_duration_name],1,Rfast::nth, k=num_of_symps, descending = T)
  }else{
    print("Please enter something for names_of_symps or num_of_symps")
  }
  return(sympdurs)
}



#' This function takes a day count and returns a vector of binaries indicating whether each 
#' person still had numsymps (or more) symptoms after numdays days.

engineerLongCOVIDBinary <- function(dat,
                                    numdays = 28, 
                                    numsymps=1,
                                    symptom_durations=paste0("symptom_duration_",1:29)){
  res <- rowSums(dat[,symptom_durations]>=numdays, na.rm = T)
  # allnas <- rowSums(is.na(dat[,symptom_durations])) == length(symptom_durations)
  # res[allnas] <- NA_real_
  out=as.numeric(res>=numsymps)
  # out[allnas] <- NA_real_
  return(out)
}

### Function to calculate rates of 4 and 12 week covid by a certain variable ###

longCOVIDtable <- function(dat, var = "age_group_broad", levs = c("25-44","18-24","45-64", "65+"), 
                           long_cov_names = long_sympnames[1:29]){
  
  if(is.null(var)){
    num_pos = nrow(dat)
    num_pos_symp = sum(!is.na(dat$covidc))
    
    ### get varnames
    covid_yesnos <- dat %>% select(covidsym_01:covidsym_26,covidsym_28:covidsym_30) %>% colnames()
    longcovid_yesnos <- dat %>% select(longcovida_1:longcovida_26,covidsym_28:covidsym_30) %>% colnames()
    
    long_covid_tab <- data.frame(symptom=long_sympnames[1:29], symptomatic = NA, symptomatic_4_weeks = NA, 
                                 symptomatic_12_weeks = NA, mean_symptom_duration = NA)
    
    ### Loop to generate results
    for (i in 1:length(long_cov_names)){
      long_covid_tab[i,1] <- long_cov_names[[i]]
      long_covid_tab[i,2] <- sum(dat[,covid_yesnos[[i]]] == 1 & !is.na(dat[,covid_yesnos[[i]]]) & 
                                   dat$covidc %in% c(2:4))
      long_covid_tab[i,3] <- sum(dat[,symptom_durations[[i]]] >= 28 & !is.na(dat[,symptom_durations[[i]]]))
      long_covid_tab[i,4] <- sum(dat[,symptom_durations[[i]]] >= 84 & !is.na(dat[,symptom_durations[[i]]]))
      long_covid_tab[i,5] <- round(colMeans(dat[!is.na(dat[,symptom_durations[[i]]]),symptom_durations[[i]]]),1)
    }
    
    ### Loop to add %s
    for (i in 1:length(long_cov_names)){
      long_covid_tab[i,2] <- paste0(as.numeric(long_covid_tab[i,2]), " (", round(100* as.numeric(long_covid_tab[i,2]) / num_pos,1), "%)")
      long_covid_tab[i,3] <- paste0(as.numeric(long_covid_tab[i,3]), " (", round(100*as.numeric(long_covid_tab[i,3]) / num_pos,1), "%)")
      long_covid_tab[i,4] <- paste0(as.numeric(long_covid_tab[i,4]), " (", round(100*as.numeric(long_covid_tab[i,4]) / num_pos,1), "%)")
    }
  }
  
  
  else{
    
    
    
    numlevels <- length(as.character(sort(unique(unlist(pull(dat, var))))))
    alllevels <- as.character(sort(unique(unlist(pull(dat, var)))))
    
    ### get numbers of positive
    numposlist <- list()
    for(i in 1:numlevels){
      numposlist[[i]] <- sum(dat[,var] == alllevels[[i]], na.rm=T)
      names(numposlist)[[i]] <-alllevels[[i]]
    }
    
    
    ### Get numbers of symptomatic
    numposlistSymp <- list()
    for(i in 1:numlevels){
      numposlistSymp[[i]] <- sum(dat[,var] == alllevels[[i]] & !is.na(dat$covidc) & 
                                   dat$covidc %in% c(2:4), na.rm=T)
      names(numposlistSymp)[[i]] <-alllevels[[i]]
      
    }
    
    ### create data frame for results
    long_covid_tab <- data.frame(matrix(ncol = length(alllevels)*3 + 1,nrow = 0,dimnames = list(NULL,c("Symptom", paste0("Symptomatic_", alllevels), 
                                                                                                       paste0("Four_week_symptomatic_", alllevels), 
                                                                                                       paste0("Twelve_week_symptomatic_", alllevels)))))
    
    
    ### get varnames
    covid_yesnos <- dat %>% select(covidsym_01:covidsym_26,covidsym_28:covidsym_30) %>% colnames()
    longcovid_yesnos <- dat %>% select(longcovida_1:longcovida_26,covidsym_28:covidsym_30) %>% colnames()
    
    ### Loop to generate results
    for (i in 1:length(long_cov_names)){
      for(j in 1:length(levs)){
        k=j+1
        long_covid_tab[i,k] <- sum(dat[,covid_yesnos[[i]]] == 1 & !is.na(dat[,covid_yesnos[[i]]]) &pull(dat, var) == alllevels[[j]], na.rm=T)
        long_covid_tab[i,k] <- paste0(as.numeric(long_covid_tab[i,k]), " (", round(100* as.numeric(long_covid_tab[i,k]) / numposlist[[j]],1), "%)")
        
      }
      
      for(j in 1:length(levs)){
        k <- length(levs) + j + 1
        long_covid_tab[i,k] <- sum(dat[,symptom_durations[[i]]] >= 28 & !is.na(dat[,symptom_durations[[i]]]) & pull(dat, var) == alllevels[[j]], na.rm=T)
        long_covid_tab[i,k] <- paste0(as.numeric(long_covid_tab[i,k]), " (", round(100* as.numeric(long_covid_tab[i,k]) / numposlist[[j]],1), "%)")
        
      }
      
      for(j in 1:length(levs)){
        k <- (2*length(levs)) + j + 1
        long_covid_tab[i,k] <- sum(dat[,symptom_durations[[i]]] >= 84 & !is.na(dat[,symptom_durations[[i]]]) & pull(dat, var) == alllevels[[j]], na.rm=T)
        long_covid_tab[i,k] <- paste0(as.numeric(long_covid_tab[i,k]), " (", round(100* as.numeric(long_covid_tab[i,k]) / numposlist[[j]],1), "%)")
      }
    }
    
    ### add symtom name
    long_covid_tab$Symptom <- long_cov_names
  }
  
  return(long_covid_tab)
  
}


### Function to get table of prevalence and  symtomaticsness

longCOVIDtableSimple <- function(dat, var = "age_group_broad", levs = c("25-44","18-24","45-64", "65+"), 
                                 long_cov_names = long_sympnames,
                                 show_props_of_symptomatics=F){
  
  if(is.null(levs)){
    levs <- unique(pull(dat,var))[!is.na(unique(pull(dat,var)))]
  }
  
  # Create new 4 and 12 week symptom count variables
  dat$symp_count_4_weeks <- rowSums(dat[,symptom_durations[1:29]] >= 28, na.rm = TRUE)
  dat$symp_count_12_weeks <- rowSums(dat[,symptom_durations[1:29]] >= 84, na.rm = TRUE)
  
  if(is.null(var)){
    df_results <- data.frame(category = "Full_cohort", 
                             Antibody_positive=nrow(dat), 
                             Symptomatic_at_infection = sum(!is.na(dat$covidc) & 
                                                              dat$covidc %in% c(2:4)), 
                             Symptomatic_4_weeks = sum(pull(dat,symp_count_4_weeks) > 0), 
                             Symptomatic_12_weeks = sum(pull(dat,symp_count_12_weeks) > 0))
    if(show_props_of_symptomatics){
      # df_results$Symptomatic_at_infection <- paste0(df_results$Symptomatic_at_infection, " (",round(100*df_results$Symptomatic_at_infection / nrow(dat), 0),"%)")
      numasymp = sum(dat$covidc_cat == "No symptoms", na.rm=T)
      denom_4_weeks <-  nrow(dat[dat$test_date - dat$covidsta >= 28,]) - numasymp
      denom_12_weeks <-  nrow(dat[dat$test_date - dat$covidsta >= 84,]) - numasymp
      df_results$Symptomatic_4_weeks <- paste0(df_results$Symptomatic_4_weeks, " (",round(100*df_results$Symptomatic_4_weeks / denom_4_weeks, 1),"%)")
      df_results$Symptomatic_12_weeks <- paste0(df_results$Symptomatic_12_weeks, " (",round(100*df_results$Symptomatic_12_weeks /denom_12_weeks, 1),"%)")
      
    }else{
      numasymp = sum(dat$covidc_cat == "No symptoms", na.rm=T)
      denom_4_weeks <-  nrow(dat[dat$test_date - dat$covidsta >= 28,])
      denom_12_weeks <-  nrow(dat[dat$test_date - dat$covidsta >= 84,])
      df_results$Symptomatic_4_weeks <- paste0(df_results$Symptomatic_4_weeks, " (",round(100*df_results$Symptomatic_4_weeks / denom_4_weeks, 1),"%)")
      df_results$Symptomatic_12_weeks <- paste0(df_results$Symptomatic_12_weeks, " (",round(100*df_results$Symptomatic_12_weeks /denom_12_weeks, 1),"%)")
    }
    
  }
  
  else{
    
    numlevels <- length(as.character(sort(unique(unlist(pull(dat, var))))))
    alllevels <- as.character(sort(unique(unlist(pull(dat, var)))))
    
    df_results <- data.frame(category = alllevels, Antibody_positive=NA, 
                             Symptomatic_at_infection = NA, 
                             Symptomatic_4_weeks = NA,  
                             Symptomatic_12_weeks = NA )
    longcovid_durations <- dat %>% select(symptom_durations[1:29]) %>% colnames()
    
    
    ### get numbers of positive
    numposlist <- list()
    for(i in 1:numlevels){
      numposlist[[i]] <- sum(dat[,var] == alllevels[[i]],na.rm = T)
    }
    names(numposlist) <-alllevels 
    
    ### Get numbers of sympotomatic
    numposlistSymp <- list()
    for(i in 1:numlevels){
      numposlistSymp[[i]] <- sum(dat[,var] == alllevels[[i]] & !is.na(dat$covidc) & 
                                   dat$covidc %in% c(2:4), na.rm=T)
    }
    names(numposlistSymp) <-alllevels 
    
    ### Get numbers of sympotomatic at 4 weeks
    numposlistSymp_4week <- list()
    for(i in 1:numlevels){
      numposlistSymp_4week[[i]] <- sum(dat[,var] == alllevels[[i]] & pull(dat,symp_count_4_weeks) > 0, na.rm=T)
    }
    names(numposlistSymp_4week) <-alllevels 
    
    ### Get numbers of followed up at 4 weeks
    denom_4_weeks <- list()
    for(i in 1:numlevels){
      if(show_props_of_symptomatics){
        asymp=0
      }else{
        asymp <- sum(dat[,var] == alllevels[[i]] & dat$covidc_cat  == "No symptoms", na.rm=T)
      }
      denom_4_weeks[[i]] <- sum(dat[,var] == alllevels[[i]] & dat$test_date - dat$covidsta >= 28, na.rm=T)+ asymp
    }
    names(denom_4_weeks) <-alllevels 
    
    ### Get numbers of sympotomatic at 12 weeks
    numposlistSymp_12week <- list()
    for(i in 1:numlevels){
      numposlistSymp_12week[[i]] <- sum(dat[,var] == alllevels[[i]] & 
                                          pull(dat,symp_count_12_weeks) > 0, na.rm=T)
    }
    names(numposlistSymp_12week) <-alllevels 
    
    ### Get numbers of followed up at 412 weeks
    denom_12_weeks <- list()
    for(i in 1:numlevels){
      if(show_props_of_symptomatics){
        asymp=0
      }else{
        asymp <- sum(dat[,var] == alllevels[[i]] & dat$covidc_cat  == "No symptoms", na.rm=T)
      }      
      denom_12_weeks[[i]] <- sum(dat[,var] == alllevels[[i]] & dat$test_date - dat$covidsta >= 84, na.rm=T) + asymp
    }
    names(denom_12_weeks) <-alllevels 
    
    ### Enter results into data frame
    for(i in 1:numlevels){
      df_results[i,2] <- numposlist[[i]]
      if(!show_props_of_symptomatics){
        df_results[i,3] <- paste0(numposlistSymp[[i]]," (", round(100*numposlistSymp[[i]]/numposlist[[i]],1),"%)")
      }else{
        df_results[i,3] <- numposlistSymp[[i]]
      }
      df_results[i,4] <- paste0(numposlistSymp_4week[[i]]," (", round(100*numposlistSymp_4week[[i]]/denom_4_weeks[[i]],1),"%)")
      df_results[i,5] <- paste0(numposlistSymp_12week[[i]]," (", round(100*numposlistSymp_12week[[i]]/denom_12_weeks[[i]],1),"%)")
      
    }
  }
  
  return(df_results)
}




# Get symptom persistence over time for all symptoms ----------------------

### Define function to get symptom persistence
getSymptomPersistence <- function(dat){
  symptom_persistence_df <- data.frame(matrix(data = NA, nrow = 300, ncol = 31,
                                              dimnames = list(NULL, c("day_count", "nobs",
                                                                      long_sympnames[c(1:26,28:30)]))))
  symptom_persistence_df$day_count <- 1:300
  num_in_total <- nrow(dat)
  for (i in 1:300){
    count_censored <- sum(dat$test_date - dat$covidsta <= i, na.rm = T)
    symptom_persistence_df[[i,2]] <- num_in_total- count_censored
    
    for (y in 1:29){
      
      symptom_persistence_df[[i,y+2]] <- sum(pull(dat, symptom_durations[[y]]) >= i, na.rm = T)
      
    }
    
  }
  return(symptom_persistence_df)
}








# Function to run full cluster analysis in one go -------------------------



runClusterAnalysis <- function(datmat,stability=T,sampsize=0.1, 
                               num_bootstraps=20,
                               k.max=10, 
                               force_k=NULL,
                               symptom_names=sympnames[c(1:26,28:30)]){
  
  #### FPC package
  dist.mat.t0 <- hammingDist(X = datmat %>% t() %>% as.matrix(), Y=NULL)
  
  ### Silhouette method for optimal num clusters
  nclust_plot_sil <- fviz_nbclust(x = datmat,
                                  FUNcluster = cluster::clara, 
                                  method = "silhouette",
                                  diss = dist.mat.t0,
                                  k.max = k.max,
                                  verbose = T,
                                  metric = "manhattan",
                                  samples = 50,
                                  sampsize = ceiling(sampsize*nrow(datmat)),
                                  pamLike = TRUE)
  
  # plot
  nclust_plot_sil + theme_bw()
  
  
  # Get  optimal k
  if(is.null(force_k)){
    k_opt <- which.max(nclust_plot_sil$data$y)
  }else{
    k_opt <- force_k
  }
  
  # run final clustering
  fpc.clara <- fpc::pamk(data = datmat, 
                         krange = k_opt, 
                         criterion = "asw",
                         usepam = F,
                         critout = T, 
                         seed = 123, 
                         diss = dist.mat.t0,
                         metric = "manhattan",
                         samples = 100,
                         sampsize = ceiling(sampsize*nrow(datmat)),
                         pamLike = TRUE
  )
  
  
  ### Now clusterboot
  boot.obj <- fpc::clusterboot(data = datmat,
                               B = num_bootstraps, 
                               distances = F,
                               # datatomatrix = F,
                               bootmethod = "boot",
                               bscompare = T,
                               # subtuning =ceiling(0.8*nrow(datmat)),
                               clustermethod = claraCBI,
                               seed = 123,
                               usepam=F,
                               count = T,
                               metric = "manhattan",
                               k= k_opt,
                               samples = 50,
                               sampsize = ceiling(sampsize*nrow(datmat)),
                               pamLike = TRUE
  )
  
  ### Get by-cluster stability
  avgJaccard=boot.obj$bootmean
  instability=boot.obj$bootbrd/num_bootstraps
  clusters=c(1:k_opt)
  bootstrap_eval=cbind(clusters,avgJaccard,instability)
  
  
  if(stability){
    ### APN stability
    stab <- matrix(0,nrow=ncol(datmat), ncol=4, 
                   dimnames = list(NULL,c("APN","AD","ADM", "FOM")))
    
    for(del in 1:ncol(datmat)){
      print(del)
      dat <- datmat[,-del]
      fpc.tmp <- fpc::pamk(data = dat, 
                           krange =  k_opt, 
                           criterion = "asw",
                           usepam = F,critout = T, seed = 123, diss = F,
                           metric = "manhattan",
                           samples = 50,
                           sampsize = ceiling(sampsize*nrow(dat)),
                           pamLike = TRUE
      )
      stab[del,] <-   clValid::stability(mat = datmat,del = del,
                                         cluster = fpc.clara$pamobject$clustering,
                                         method = "manhattan", 
                                         clusterDel = fpc.tmp$pamobject$clustering
      )
    }
    
    stab <- as.data.frame(stab)
    stab$symptom <- symptom_names
    #' APN = average proportion of non-overlap - % of obs that change cluster
    #' AD - average distance ave dist between obs in same cluster in full vs del
    #' ADM - average distance between means (doesn't apply for non-Euclidean)
    #' FOM - figure of merit
    
    p1_APN <- stab %>% ggplot(aes(x=reorder(symptom,APN),y=APN)) +
      geom_col() +
      theme_bw() +
      coord_flip() +
      labs(x= "Symptom", y="Average proportion of non-overlap")
    
    p2_FOM <- stab %>% ggplot(aes(x=reorder(symptom,FOM),y=FOM)) +
      geom_col() +
      theme_bw() +
      coord_flip() +
      labs(x= "Symptom", y="Figure of Merit")
    
    p3_APN_FOM <- stab %>% ggplot(aes(x=APN,y=FOM)) +
      geom_point() +
      theme_bw() +
      coord_flip() +
      labs(x= "APN", y="Figure of Merit")
    
    return(list(data=datmat,
                silhouette_plot=nclust_plot_sil + theme_bw(),
                # distance_mat=dist.mat.t0,
                cluster_results=fpc.clara,
                stability_results=stab,
                stability_plots=list(apn=p1_APN,
                                     fom=p2_FOM,
                                     apn_fom=p3_APN_FOM),
                bootstrap_stability=list(boot.obj,
                                         bootstrap_eval)))
  }else{
    return(list(data=datmat,
                silhouette_plot=nclust_plot_sil + theme_bw(),
                # distance_mat=dist.mat.t0,
                cluster_results=fpc.clara,
                bootstrap_stability=list(boot.obj,
                                         bootstrap_eval)))
  }
  
  
}





# Make all vars categorical -----------------------------------------------


makeVarsFactors <- function(dfRes){
  dfRes$region_named <- as.factor(dfRes$region_named)
  dfRes <- within(dfRes, region_named <- relevel(region_named,ref="South East"))
  dfRes$hh_size_cat <- as.factor(dfRes$hh_size_cat)
  dfRes <- within(dfRes, hh_size_cat <- relevel(hh_size_cat,ref="2"))
  dfRes$carehome <- as.factor(dfRes$carehome)
  dfRes <- within(dfRes, carehome <- relevel(carehome,ref="No"))
  dfRes$edu_cat <- as.factor(dfRes$edu_cat)
  dfRes <- within(dfRes, edu_cat <- relevel(edu_cat,ref="GCSE"))
  dfRes$shielding <- as.factor(dfRes$shielding)
  dfRes <- within(dfRes, shielding <- relevel(shielding,ref="No")) 
  dfRes$smokenow <- as.factor(dfRes$smokenow)
  dfRes <- within(dfRes, smokenow <- relevel(smokenow,ref="Not current cigarette smoker"))
  dfRes$vapenow <- as.factor(dfRes$vapenow)
  dfRes <- within(dfRes, vapenow <- relevel(vapenow,ref="Not current vaper"))
  dfRes$sex <- as.factor(dfRes$sex)
  dfRes <- within(dfRes, sex <- relevel(sex,ref="Male"))
  dfRes$ethnic_new <- as.factor(dfRes$ethnic_new)
  dfRes <- within(dfRes, ethnic_new <- relevel(ethnic_new,ref="White"))
  dfRes$work1_healthcare_or_carehome_worker <- as.factor(dfRes$work1_healthcare_or_carehome_worker)
  dfRes <- within(dfRes, work1_healthcare_or_carehome_worker <- relevel(work1_healthcare_or_carehome_worker,ref="No"))
  dfRes$gross_household_cat <- as.factor(dfRes$gross_household_cat)
  dfRes <- within(dfRes, gross_household_cat <- relevel(gross_household_cat,ref="15,000-49,999"))
  dfRes$imd_quintile_cat <- as.factor(dfRes$imd_quintile_cat)
  dfRes <- within(dfRes, imd_quintile_cat <- relevel(imd_quintile_cat,ref="3"))
  dfRes$clin_vulnerable <- as.factor(dfRes$clin_vulnerable)
  dfRes <- within(dfRes, clin_vulnerable <- relevel(clin_vulnerable,ref="No"))
  dfRes$covid_severity <- as.factor(dfRes$covid_severity)
  dfRes <- within(dfRes, covid_severity <- relevel(covid_severity,ref="No medical attention sought"))
  dfRes$covida_cat <- as.factor(dfRes$covida_cat)
  dfRes <- within(dfRes, covida_cat <- relevel(covida_cat,ref="No COVID"))
  dfRes$covidc_cat <- as.factor(dfRes$covidc_cat)
  dfRes <- within(dfRes, covidc_cat <- relevel(covidc_cat,ref="No symptoms"))
  dfRes$age_group_named <- as.factor(dfRes$age_group_named)
  dfRes <- within(dfRes, age_group_named <- relevel(age_group_named,ref="55-64"))
  dfRes <- within(dfRes, age_group_four_cat <- relevel(age_group_four_cat,ref="Under 60"))
  dfRes$age_group_paul <- as.factor(dfRes$age_group_paul)
  dfRes <- within(dfRes, age_group_paul <- relevel(age_group_paul,ref="18-29"))
  dfRes$bmi_cat <- as.factor(dfRes$bmi_cat)
  dfRes <- within(dfRes, bmi_cat <- relevel(bmi_cat,ref="Normal weight"))
  dfRes$hospitalised_covid <- as.factor(dfRes$hospitalised_covid)
  dfRes <- within(dfRes, hospitalised_covid <- relevel(hospitalised_covid,ref="No"))
  dfRes <- dfRes %>% 
    mutate(across(c(paste0("healtha_0",c(1:9)),
                    paste0("healtha_1",0:7)),as_factor))
  return(dfRes)
}







# Prediction functions ----------------------------------------------------



## Univariate --------------------------------------------------------------


runUnivariate <- function(mydat,outcome_var = "lc_84",negatives ="vs_all", outpath_new,
                          figpath_new,mywidth=10,myheight=6){
  
  ### Relevel caovidc_cat as there are no people with no symptoms
  mydat$covidc_cat <- as.factor(mydat$covidc_cat)
  mydat <- within(mydat, covidc_cat <- relevel(covidc_cat,ref="Mild symptoms"))
  mydat$covidc_cat <- drop.levels(mydat$covidc_cat)
  
  adjustments=c("sex", "age_group_named","ethnic_new","bmi_cat", 
            "work1_healthcare_or_carehome_worker", "imd_quintile_cat",
            "smokenow","vapenow","covid_severity",
            "gross_household_cat")
  newcovs=mycovs[mycovs!="dummy"]
  
  
  
  if(negatives=="vs_all"){
    mydat[,outcome_var] <- as.numeric(pull(mydat,outcome_var) ==1 & !is.na(pull(mydat,outcome_var)))
  }else{
    mydat <- mydat[!is.na(pull(mydat,outcome_var)),]
  }
  
  lc_84_cluster_1_univ <- OverReact::ModelMakerMulti(dat = mydat,
                                                     list_of_variables_of_interest = newcovs,
                                                     outcome = outcome_var,
                                                     joint_adjustment_vars = adjustments,
                                                     cov_name_list = (cov_name_list)
  )
  
  
  
  lc_84_cluster_1_ORs <- 	lc_84_cluster_1_univ$df_output[lc_84_cluster_1_univ$df_output$Level != "0 [reference]",]
  write_csv(x = lc_84_cluster_1_ORs, file = paste0(outpath_new, "univ_ORs_",outcome_var,"_",negatives,".csv"))
  
  
  
  # Forest plot -------------------------------------------------------------
  dodger = position_dodge2(width=0.8,reverse = T)
  
  # dodger = position_dodge(width=0.8)
  order_levels <- c("Underweight" ,
                    "Normal weight [reference]", 
                    "Overweight",
                    "Obese",                                                       
                    
                    "Female" , 
                    "Male [reference]",
                    "18-24",
                    "25-34",
                    "35-44",
                    "45-54",
                    "55-64 [reference]" ,
                    "65-74",
                    "74+",
                    "0-14,999",
                    "15,000-49,999 [reference]" ,
                    "50,000-149,999",
                    ">150,000",
                    "1 - most deprived" , "2","3 [reference]","4","5 - least deprived",
                    "No qualification","Other","Post-GCSE qualification","GCSE [reference]",
                    "Degree or above",
                    "Mild symptoms [reference]","Moderate symptoms", "Severe symptoms",
                    "Not current cigarette smoker [reference]",
                    "Current cigarette smoker",
                    "Not current vaper [reference]",
                    "Current vaper",
                    "White [reference]", 
                    "Asian"  , "Black", "Mixed",
                    "No [reference]", "Yes"
  )
  
  
  ###Switch names of comorbidities
  
  cov_name_list$imd_quintile_cat
  lc_84_cluster_1_univ$plot_output[lc_84_cluster_1_univ$plot_output$Level==1,"Level"] <- lc_84_cluster_1_univ$plot_output[lc_84_cluster_1_univ$plot_output$Level==1,"predictor"]

  
  lc_84_cluster_1_univ$plot_output[lc_84_cluster_1_univ$plot_output$Level==lc_84_cluster_1_univ$plot_output$predictor &
                                     lc_84_cluster_1_univ$plot_output$Level !="Current vaper","predictor"] <- "Comorbidities"
  
  
  ### Add to level order
  order_levels_2 <- c(order_levels,
                      predictors_df[predictors_df$variable_type=="Health",]$variable_desc)

  
  myplot_84_all <- lc_84_cluster_1_univ$plot_output %>% 
    filter(!Level %in% c("Prefer not to say","Household size", "Anxiety", "Depression", 
                         "Psychiatric disorder"),
           adjustment %in% c(3,11),
           predictor %in% c("Ethnicity",
                            # "Severity of COVID symptoms",
                            "Adiposity",
                            "Sex",
                            "Age group",
                            "Gross household income",
                            "Index of multiple deprivation (IMD) quintile",                            "Current smoker",
                            "Current vaper",
                            "Healthcare or care home worker",
                            "Hospitalised with COVID"
           )) %>% 
    mutate(Level = factor(Level, levels = order_levels_2),
           model_type = ifelse(adjustment==3, "Adjusted on \nage and sex", "Mutually \nadjusted"),
           predictor = factor(predictor, levels = c("Sex",
                                                    "Age group",
                                                    "Ethnicity",
                                                    "Adiposity",
                                                    "Current smoker",
                                                    "Current vaper",
                                                    "Gross household income",
                                                    "Index of multiple deprivation (IMD) quintile",
                                                    "Healthcare or care home worker",
                                                    "Hospitalised with COVID"))) %>% 
    ggplot(aes(x=(Level), y = OR, col = model_type)) +
    geom_point(position = dodger) +
    geom_errorbar(aes(ymin = Lower, ymax = Upper), position = dodger, width = 0.7) +
    coord_flip() +
    # scale_x_discrete(labels=function(x) str_wrap(x, width=33)) +
    # scale_x_reordered() +
    geom_hline(yintercept = 1, linetype = "dashed", col = "grey60" )+
    scale_color_manual(values = myCols) +
    scale_y_log10()+
    theme_bw() +
    # ylim(c(0.3,1.8))+
    
    theme_adjust +
    ggforce::facet_col(.~predictor,  scales = "free_y", space = "free") +
    labs(col = "Model", x="", y="Odds ratio")
  
  
  myplot_84_all
  OverReact::saveREACTplot(p = myplot_84_all,figpath = figpath_new,
                           filename = paste0("ORs_all_", outcome_var,"_",negatives),
                           width = mywidth,height = myheight)
  
  return(myplot_84_all)
  
}



## Multivariate ------------------------------------------------------------



PermuteCatboost <- function(mymodel,X, y, nperm = 100, lossmeasure="auc"){
  predictors <- colnames(X)
  perm_df= data.frame(matrix(data=NA, nrow=nperm, 
                             ncol = ncol(X), 
                             dimnames = list(1:nperm,predictors)))
  pool_unperm = catboost.load_pool(data =X,label = y, cat_features = 3:ncol(X))
  preds_full=catboost.predict(model = mymodel, pool = pool_unperm,verbose = T,prediction_type = "Probability")
  
  logloss_full= ModelMetrics::logLoss(as.vector(y),as.vector(preds_full))
  AUC_full=ModelMetrics::auc(as.vector(y),as.vector(preds_full))
  if(lossmeasure=="auc"){
    loss_full=AUC_full
  }else if(lossmeasure=="logloss"){
    loss_full=logloss_full
  }else{
    print("Please set lossmeasure to 'auc' or 'logloss'" )
  }
  
  for(i in 1:length(predictors)){
    print(paste("Permuting", predictors[[i]]))
    for(j in 1:nperm){
      print(paste0("Permutation ", j, " complete"))
      newX=X
      newX[,predictors[[i]]]=gtools::permute((unlist(X[,predictors[[i]]])))
      pool_perm = catboost.load_pool(data =newX,label = y, cat_features = 3:ncol(X))
      preds_new=catboost.predict(model = mymodel, pool = pool_perm,verbose = T,prediction_type = "Probability")
      logloss_new= ModelMetrics::logLoss(as.vector(y),as.vector(preds_new))
      auc_new= ModelMetrics::auc(as.vector(y),as.vector(preds_new))
      if(lossmeasure=="auc"){
        loss_new=auc_new
      }else if(lossmeasure=="logloss"){
        loss_new=logloss_new
      }
      perm_df[j,i] = loss_new
    }
  }
  
  colnames(perm_df) <- predictors
  perm_df_long <- perm_df %>% tidyr::pivot_longer(cols = predictors)
  if(lossmeasure=="auc"){
    perm_df_long$value=loss_full-perm_df_long$value
  }else if(lossmeasure=="logloss"){
    perm_df_long$value=perm_df_long$value-loss_full
  }
  lossname=paste0("mean_",lossmeasure,"_permute")
  lossname_sd=paste0("SD_",lossmeasure,"_permute")
  
  results.df = data.frame(predictor=predictors, 
                          mean_loss_permute = NA_real_,
                          SD_loss_permute = NA_real_,
                          prop_permutes_greater_than_unpermute=NA_real_,
                          mean_loss_permute_delta = NA_real_)
  
  results.df$mean_loss_permute <- colMeans(perm_df, na.rm = T)
  results.df$SD_loss_permute <- sapply(perm_df, sd)
  
  ### Calculate proportion of permuted runs that had a superior log loss to the unpermuted
  results.df$prop_permutes_greater_than_unpermute <- colSums(perm_df>=loss_full) /nperm
  
  ### Get deltas
  results.df$mean_loss_permute_delta =abs(loss_full-results.df$mean_loss_permute)
  results.df$mean_loss_permute_delta_lower <-abs(results.df$mean_loss_permute_delta- 1.96*(
    results.df$SD_loss_permute) / sqrt(nperm))
  results.df$mean_loss_permute_delta_upper <- abs(results.df$mean_loss_permute_delta+ 1.96*(
    results.df$SD_loss_permute) / sqrt(nperm))
  
  return(list(summary_stats=results.df,
              full_results=perm_df_long))
}



# Permute catboost parallel -----------------------------------------------


PermuteCatboostParallel <- function(mymodel,X, y, nperm = 100, lossmeasure="auc", 
                                    parallelise=T,nclust=20){
  predictors <- colnames(X)
  # perm_df= data.frame(matrix(data=NA, nrow=nperm, 
  #                            ncol = ncol(X), 
  #                            dimnames = list(1:nperm,predictors)))
  pool_unperm = catboost.load_pool(data =X,label = y, cat_features = 3:ncol(X))
  preds_full=catboost.predict(model = mymodel, pool = pool_unperm,verbose = T,prediction_type = "Probability")
  
  logloss_full= ModelMetrics::logLoss(as.vector(y),as.vector(preds_full))
  AUC_full=ModelMetrics::auc(as.vector(y),as.vector(preds_full))
  if(lossmeasure=="auc"){
    loss_full=AUC_full
  }else if(lossmeasure=="logloss"){
    loss_full=logloss_full
  }else{
    print("Please set lossmeasure to 'auc' or 'logloss'" )
  }
  if(parallelise){
    ### parallelised loop for getting stability coefficients
    cl <- parallel::makeCluster(nclust)
    doParallel::registerDoParallel(cl)
    foreach::getDoParWorkers()
    perm_df <- foreach(i = 1:length(predictors), .combine = cbind, 
                       .packages = c("caret", "pROC", "ModelMetrics", 
                                     "catboost","gtools")) %dopar% {
                                       perm_vec=(1:nperm)
                                       print(paste("Permuting", predictors[[i]]))
                                       for(j in 1:nperm){
                                         # print(paste0("Permutation ", j, " complete"))
                                         newX=X
                                         newX[,predictors[[i]]]=gtools::permute((unlist(X[,predictors[[i]]])))
                                         pool_perm = catboost.load_pool(data =newX,label = y, cat_features = 3:ncol(X))
                                         preds_new=catboost.predict(model = mymodel, pool = pool_perm,verbose = T,prediction_type = "Probability")
                                         logloss_new= ModelMetrics::logLoss(as.vector(y),as.vector(preds_new))
                                         auc_new= ModelMetrics::auc(as.vector(y),as.vector(preds_new))
                                         if(lossmeasure=="auc"){
                                           loss_new=auc_new
                                         }else if(lossmeasure=="logloss"){
                                           loss_new=logloss_new
                                         }
                                         perm_vec[j] = loss_new
                                       }                   
                                       
                                       return(perm_vec)
                                       
                                     }
    stopCluster(cl)
  }else{
    # Unparallelised
    for(i in 1:length(predictors)){
      print(paste("Permuting", predictors[[i]]))
      for(j in 1:nperm){
        print(paste0("Permutation ", j, " complete"))
        newX=X
        newX[,predictors[[i]]]=gtools::permute((unlist(X[,predictors[[i]]])))
        pool_perm = catboost.load_pool(data =newX,label = y, cat_features = 3:ncol(X))
        preds_new=catboost.predict(model = mymodel, pool = pool_perm,verbose = T,prediction_type = "Probability")
        logloss_new= ModelMetrics::logLoss(as.vector(y),as.vector(preds_new))
        auc_new= ModelMetrics::auc(as.vector(y),as.vector(preds_new))
        if(lossmeasure=="auc"){
          loss_new=auc_new
        }else if(lossmeasure=="logloss"){
          loss_new=logloss_new
        }
        perm_df[j,i] = loss_new
      }
    }
  }
  
  # wrangle output data frame
  perm_df <- perm_df %>% as.data.frame()
  colnames(perm_df) <- predictors
  perm_df_long <- perm_df %>% tidyr::pivot_longer(cols = predictors)
  if(lossmeasure=="auc"){
    perm_df_long$value=loss_full-perm_df_long$value
  }else if(lossmeasure=="logloss"){
    perm_df_long$value=perm_df_long$value-loss_full
  }
  lossname=paste0("mean_",lossmeasure,"_permute")
  lossname_sd=paste0("SD_",lossmeasure,"_permute")
  
  results.df = data.frame(predictor=predictors, 
                          mean_loss_permute = NA_real_,
                          SD_loss_permute = NA_real_,
                          prop_permutes_greater_than_unpermute=NA_real_,
                          mean_loss_permute_delta = NA_real_)
  
  results.df$mean_loss_permute <- colMeans(perm_df, na.rm = T)
  results.df$SD_loss_permute <- sapply(perm_df, sd)
  
  ### Calculate proportion of permuted runs that had a superior log loss to the unpermuted
  results.df$prop_permutes_greater_than_unpermute <- colSums(perm_df>=loss_full) /nperm
  
  ### Get deltas
  results.df$mean_loss_permute_delta =abs(loss_full-results.df$mean_loss_permute)
  results.df$mean_loss_permute_delta_lower <-abs(results.df$mean_loss_permute_delta- 1.96*(
    results.df$SD_loss_permute) / sqrt(nperm))
  results.df$mean_loss_permute_delta_upper <- abs(results.df$mean_loss_permute_delta+ 1.96*(
    results.df$SD_loss_permute) / sqrt(nperm))
  
  return(list(summary_stats=results.df,
              full_results=perm_df_long))
}






# Get variable importance measures ----------------------------------------


### Define function to get variable importances and SHAPS
getFeatureImpAndSHAP <- function(cb_model, validation_pool, learn_pool,
                                 modelvars="covariates",myheight=5,
                                 train_learn,figpath_new=figpath,
                                 X,y, bandwidth_scaler=5,nperm=100,outcome,bonferroni=FALSE,
                                 parallelise=T,nclust=20){
  
  
  
  # get feature importances
  cb_varimps=catboost.get_feature_importance(model = cb_model, pool = validation_pool, type = "LossFunctionChange") %>% as.data.frame()
  cb_varimps$feature =rownames(cb_varimps)
  colnames(cb_varimps)[[1]] <- "Mean_delta_logloss"
  
  plot_varimps <- cb_varimps %>% 
    left_join(predictors_df, by = c("feature" = "variable_desc")) %>% 
    ggplot(aes(x=reorder(feature,Mean_delta_logloss), y=Mean_delta_logloss, 
               col=variable_type,fill=variable_type)) + 
    geom_col(width=0.03) +
    geom_point() +
    coord_flip() +
    scale_fill_manual(values = myCols) +
    scale_colour_manual(values = myCols) +
    theme_bw() +
    labs(x="Variable",y="Mean change in logloss \nafter variable removal", col = "Variable type",  fill = "Variable type")
  
  
  OverReact::saveREACTplot(p = plot_varimps, figpath = figpath_new,
                           filename = paste0("varimps_plot_",outcome,"_",modelvars),width =7, height=myheight)
  
  
  

  ### Get permutation varimps
  set.seed(12345)
  cb_varimps_permute=PermuteCatboostParallel(mymodel = cb_model, X = X,y = y,nperm = nperm,
                                             parallelise=parallelise,nclust=nclust) 
  
  # Get bonferroni adjusted pvalue
  bonfpval=0.05/ncol(X)
  if(bonferroni){
    mypval=bonfpval
  }else{
    mypval=0.05
  }
  
  cb_varimps_permute$summary_stats <- cb_varimps_permute$summary_stats %>% 
    mutate(selected=case_when(prop_permutes_greater_than_unpermute < mypval ~ "Selected",
                              TRUE ~ "Not selected"))
  # 
  # ### First plot, with simple columns and points
  # plot_varimps_permute <- cb_varimps_permute$summary_stats %>% 
  #   left_join(predictors_df, by = c("predictor" = "variable_desc")) %>% 
  #   ggplot(aes(x=reorder(predictor,mean_loss_permute_delta ), y=mean_loss_permute_delta, 
  #              col=selected,fill=selected)) + 
  #   geom_col(width=0.03) +
  #   geom_point() +
  #   coord_flip() +
  #   scale_fill_manual(values = myCols[c(5,3)]) +
  #   scale_colour_manual(values = myCols[c(5,3)]) +
  #   theme_bw() +
  #   labs(x="Variable",y="Mean change in AUC \nafter variable permutation", col = "Variable type",  fill = "Variable type")
  # 
  # plot_varimps_permute
  # 
  # # bandwidth_scaler=30
  # # Get bandwidth
  # bandwidth=mean(cb_varimps_permute$full_results$value)/bandwidth_scaler
  # 
  # ### seconds plot, with distributions
  # plot_varimps_permute_distribution <- cb_varimps_permute$full_results %>% 
  #   rename(predictor=name,
  #          loss_permute_delta=value) %>% 
  #   left_join(predictors_df, by = c("predictor" = "variable_desc")) %>% 
  #   left_join(cb_varimps_permute$summary_stats) %>% 
  #   arrange(prop_permutes_greater_than_unpermute,mean_loss_permute_delta) %>% 
  #   mutate(indx=row_number()) %>% 
  #   ggplot(aes(y=reorder(predictor,indx ), x=loss_permute_delta , 
  #              # height=loss_permute_delta,
  #              col=selected,fill=selected)) + 
  #   ggridges::geom_density_ridges2(alpha=0.7, rel_min_height=0.001,scale=1.1, bandwidth=bandwidth) +
  #   geom_vline(xintercept = 0, linetype="dashed", col="grey50") +
  #   scale_fill_manual(values = myCols[c(5,3)]) +
  #   scale_colour_manual(values = myCols[c(5,3)]) +
  #   theme_bw() +
  #   scale_y_discrete(labels=function(x) str_wrap(x, width=33)) +
  #   labs(y="Variable",x="Change in AUC \nafter variable permutation", 
  #        col = "Variable selected",  fill = "Variable selected")
  # 
  # plot_varimps_permute_distribution
  # 
  # OverReact::saveREACTplot(p = plot_varimps_permute_distribution, figpath = figpath_new,
  #                          filename = paste0("permutation_varimps_distributions_plot_",outcome,"_",modelvars),
  #                          width =6, height=myheight)
  # 
  # 
  # 
  
  
  ### seconds plot, with distributions
  plot_varimps_permute_boxplot <- cb_varimps_permute$full_results %>% 
    rename(predictor=name,
           loss_permute_delta=value) %>% 
    left_join(predictors_df, by = c("predictor" = "variable_desc")) %>% 
    left_join(cb_varimps_permute$summary_stats) %>% 
    arrange(selected,mean_loss_permute_delta) %>% 
    mutate(indx=row_number()) %>% 
    ggplot(aes(y=reorder(predictor,indx), x=loss_permute_delta , 
               # height=loss_permute_delta,
               col=selected)) + 
    geom_boxplot() +
    geom_vline(xintercept = 0, linetype="dashed", col="grey50") +
    scale_fill_manual(values = myCols[c(5,3)]) +
    scale_colour_manual(values = myCols[c(5,3)]) +
    theme_bw() +
    scale_y_discrete(labels=function(x) str_wrap(x, width=33)) +
    labs(y="Variable",x="Change in AUC \nafter variable permutation", 
         col = "",  fill = "")
  
  plot_varimps_permute_boxplot
  
  OverReact::saveREACTplot(p = plot_varimps_permute_boxplot, figpath = figpath_new,
                           filename = paste0("permutation_varimps_boxplot_",outcome,"_",modelvars),
                           width =8, height=myheight)
  
  
  
  ### Get feature interactions
  interacts_basemod <- catboost.get_feature_importance(model = cb_model,
                                                       pool = learn_pool,
                                                       type = "Interaction") %>% 
    as.data.frame()
  
  ### get varnames to replace indexes
  varnames_basemod <- X %>% colnames()
  varnames_basemod_df=data.frame(varnames=varnames_basemod,
                                 varindex=0:(length(varnames_basemod)-1))
  interacts_basemod$feature1_name <- NA
  interacts_basemod$feature2_name <- NA
  
  i <- 1
  for (i in 1:nrow(interacts_basemod)){
    indx <- interacts_basemod[i,]$feature1_index
    interacts_basemod[i,]$feature1_name <- varnames_basemod_df[indx+1,]$varnames
  }
  for (i in 1:nrow(interacts_basemod)){
    indx <- interacts_basemod[i,]$feature2_index
    interacts_basemod[i,]$feature2_name <- varnames_basemod_df[indx+1,]$varnames
  }
  
  ### Concatenate names
  interacts_basemod$features_comb <- paste0(interacts_basemod$feature1_name, " x ",
                                            interacts_basemod$feature2_name)
  
  ###plot
  ### seconds plot, with distributions
  plot_interacts <- interacts_basemod%>% 
    arrange(score) %>% 
    top_n(20) %>% 
    ggplot(aes(y=reorder(features_comb,score ), x=score )) + 
    geom_col(fill=myCols[2]) +
    theme_bw() +
    scale_y_discrete(labels=function(x) str_wrap(x, width=33)) +
    labs(y="Variable",x="Interaction strength", 
         col = "",  fill = "")
  
  plot_interacts
  OverReact::saveREACTplot(p = plot_interacts, figpath = figpath_new,
                           filename = paste0("interactions_plot",outcome,"_",modelvars),
                           width =6, height=6)
  
  
  
  # get SHAP values
  cb_shaps=catboost.get_feature_importance(model = cb_model, pool = learn_pool, type = "ShapValues") %>% as.data.frame()
  
  colnames(cb_shaps) <- c(cb_varimps$feature,"Bias")
  cb_shaps <- cb_shaps %>% select(-Bias)
  cb_shaps_long=cb_shaps %>% mutate(ID=row_number()) %>% 
    pivot_longer(-ID) %>% 
    rename(variable=name)
  cb_shaps_long_means=cb_shaps %>% mutate(ID=row_number()) %>% 
    pivot_longer(-ID) %>% 
    rename(variable=name) %>% 
    group_by(variable) %>% 
    summarise(mean_value=mean(abs(value)))
  cb_shaps_long <- left_join(cb_shaps_long,cb_shaps_long_means)
  
  # add values from original data
  ## quick function for 0/1 scaling
  scaleZeroOne <- function(x){
    (x-min(x,na.rm=T))/(max(x,na.rm=T)-min(x,na.rm=T))
  }
  
  cb_shaps_stdfvalues=train_learn %>% 
    select(-outcome) %>%
    # mutate(across(where(is.numeric),scale)) %>% 
    mutate_all(as.numeric) %>% 
    mutate_all(scaleZeroOne) %>% 
    mutate(ID=row_number()) %>% 
    pivot_longer(-ID) %>% 
    rename(variable=name,
           stdfvalue=value)
  
  cb_shaps_long=left_join(cb_shaps_long,cb_shaps_stdfvalues)
  
  
  # create plot
  x_bound <- max(abs(cb_shaps_long$value)) * 1.1
  label_format = "%.3f"
  plot1 <- ggplot(data = cb_shaps_long) + coord_flip(ylim = c(-x_bound, 
                                                              x_bound)) + 
    geom_hline(yintercept = 0) + ggforce::geom_sina(aes(x = reorder(variable,mean_value), 
                                                        y = value, color = stdfvalue), method = "counts", 
                                                    maxwidth = 0.7, alpha = 0.7) + 
    geom_text(data = unique(cb_shaps_long[, c("variable", "mean_value")]), aes(x =  reorder(variable,mean_value), y = -Inf, 
                                                                               label = sprintf(label_format, mean_value)),
              size = 3, alpha = 0.7, hjust = -0.2, fontface = "bold") + 
    scale_color_gradient(low = myCols[[1]], high = myCols[[3]], 
                         breaks = c(0, 1), labels = c(" Low", "High "), 
                         guide = guide_colorbar(barwidth = 12, barheight = 0.3)) + 
    theme_bw() + theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(), 
                       legend.position = "bottom", legend.title = element_text(size = 10), 
                       legend.text = element_text(size = 8), axis.title.x = element_text(size = 10)) + 
    # scale_x_discrete(limits = rev(levels(cb_shaps_long$variable)), 
    #                  labels = label.feature(rev(levels(cb_shaps_long$variable)))) + 
    labs(y = "SHAP value (impact on model output)", 
         x = "", color = "Feature value  ")
  
  plot1
 
  OverReact::saveREACTplot(p = plot1, figpath = figpath_new,
                           filename = paste0("shap_plot_",outcome,"_",modelvars),width =8, height=myheight)
  
  
  
  ### Final absolute shap value plot
  p_shap_mean <- cb_shaps_long %>% distinct(variable, .keep_all = T) %>% 
    ggplot(aes(x=mean_value, y=reorder(variable,mean_value))) +
    geom_col(fill = myCols[[2]])+
    labs(x="Mean SHAP value \n(variable importance in multivariable model)", y= "Predictor") +
    theme_bw()
  
  
  
  return(list(varimps=cb_varimps,
              varimps_plot=plot_varimps,
              plot_interacts=plot_interacts,
              permutation_varimps=cb_varimps_permute,
              # permutation_varimps_plot=plot_varimps_permute,
              # permutation_varimps_dist_plot=plot_varimps_permute_distribution,
              permutation_varimps_boxplot=plot_varimps_permute_boxplot,
              shaps=plot1,
              shap_mean=p_shap_mean))
}
