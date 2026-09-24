.libPaths(c("~/R_libs/4.5.0", .libPaths()))
library(TraMineR)
library(weights)
library(seqHMM)
library(R.matlab)
library(mice)

rm(list = ls())
gender<-1
split<-1
if (gender == 1) {
  setwd(paste0("~/male/split",split))
  print(paste('Process male, split',split))
} else {
  setwd(paste0("~/female/split",split))
  print(paste('Process male, split',split))
}

source("~/compress_paths_multi.R")

#Covariates guide cluster discovery 
covariates <- read.csv("covariates.csv",
                       na.strings = c("", "NA", "NaN", ".")) # subjects x covariates
str(covariates)

#Impute missing values in covariates
imp <- mice(
  covariates,
  m = 5,    
  maxit = 20,
  seed = 123
)
covariates_imp <- complete(imp)

#####
ns<-c(3,4,5) # number of states
nc<-c(2,3,4,5,6) # number of clusters

task_id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
print(paste("Number of Trial:",task_id))

set.seed(NULL)

# csv files showing cumulative LTC counts per category over each consecutive year from year 0 to year 90
# see disease.csv for an example
files <- c(
  "Cardiovascular.csv",
  "Eye.csv",
  "Gastrointestinal.csv",
  "Kidney.csv",
  "Liver.csv",
  "Mental.csv",
  "Metabolic.csv",
  "Musculoskeletal.csv",
  "Neurological.csv",
  "Reproductive.csv",
  "Respiratory.csv",
  "Death.csv"
)
#Death.csv MUST ALWAYS BE AT THE END OF THE LIST


data_list <- lapply(
  files,
  function(x) {
    message("Processing: ", x)
    dat <- read.csv(
      x,
      header = TRUE,
      stringsAsFactors = FALSE,
      na.strings = "NA" )
    
    # convert non-ID columns to character
    dat[-1] <- lapply(dat[-1], as.character)
    
    return(dat)
    
  }
)

names(data_list) <- tools::file_path_sans_ext(files)


seq_list <- list()

for (name in names(data_list)) {
  
  cat("Creating sequence list:", name, "\n")
  
  seq_list[[name]] <- seqdef(
    data_list[[name]][, -1],
    missing = NA
  )
  
}
#########################################

for (kk in 1:length(ns)) {
  
  for (jj in 1:length(nc)) {
    
    #Number of states
    Number_of_states <- ns[kk]  
    
    #Number of clusters 
    Number_of_clusters <- nc[jj]  
    
    # Construct output filename
    fname <- paste0(
      "hmm_results_", Number_of_states, "_states_",Number_of_clusters,"_clusters_cov_run",task_id,".mat")
    
    if (!file.exists(fname)) {
      
      print(paste0("Process Number of states ",Number_of_states,", Number of clusters ",Number_of_clusters))
      
      #Number of hidden states per class (cluster)
      states_per_class <- rep(Number_of_states, Number_of_clusters)
      
      message(sprintf("Building HMM with %s states", Number_of_states))
      message(sprintf("And %s clusters", length(states_per_class)))
      
      # Fit one HMM per intended cluster
      hmm_models <- vector("list", length(states_per_class))
      hmm_fits   <- vector("list", length(states_per_class))
      
      for (k in seq_along(states_per_class)) {
        
        hmm_models[[k]] <- build_hmm(
          observations = seq_list,
          n_states = states_per_class[k]
        )
        
        hmm_fits[[k]] <- fit_model(
          hmm_models[[k]],
          em_step = TRUE,
          local_step = TRUE,
          threads = 4
        )
      }
      
      initial_probs <- lapply(hmm_fits, function(x) x$model$initial_probs)
      
      transition_probs <- lapply(hmm_fits, function(x) x$model$transition_probs)
      
      emission_probs <- lapply(hmm_fits, function(x) x$model$emission_probs)
      
      model <- build_mhmm(
        observations      = seq_list,
        initial_probs     = initial_probs,
        transition_probs  = transition_probs,
        emission_probs    = emission_probs,
        formula = ~ Cov1 + Cov2 + Cov3 + Cov4 + Cov5 + Cov6 + Cov7 + Cov8 + Cov9 + Cov10 + 
          Cov11 + Cov12 + Cov13 + Cov14 + Cov15 + Cov16 + Cov17 + Cov18 + Cov19 + Cov20 + 
          Cov21 + Cov22 + Cov23 + Cov24 + Cov25 + Cov26 + Cov27 + Cov28 + Cov29 + Cov30 + 
          Cov31 + Cov32 + Cov33 + Cov34 + Cov35 + Cov36 + Cov37 + Cov38 + Cov39,  #covariates to guide clusters
        data = covariates_imp
      )
      
      fit <- fit_model(
        model,
        em_step = TRUE,
        local_step = TRUE,
        threads = 4
      )
      
      #Most likely path for each subject 
      message("Estimating hidden paths")
      paths <- hidden_paths(fit$model)
      
      # Compress trajectories
      message('Compressing paths')
      paths_compact <- compress_paths_multi(paths)
      paths_compact_mat <- as.matrix(paths_compact)
      storage.mode(paths_compact_mat) <- "double"
      
      #Cluster membership 
      message('Inferring cluster membership')
      pp <- posterior_cluster_probabilities(fit$model)
      pp <- unname(simplify2array(pp))
      
      message("Writing outputs to MATLAB mat files")
      # Extract quantities from fitted model
      transition_probs <- fit$model$transition_probs
      initial_probs <- fit$model$initial_probs
      emission_probs <- fit$model$emission_probs
      
      #### Save the estimated probs #####
      fname_probs<- paste0("hmm_results_", Number_of_states, "_states_",
                           Number_of_clusters,"_clusters_cov_run",task_id,"_probs.RData")
      save(emission_probs,transition_probs,initial_probs, file=fname_probs)
      ###
      
      transition_probs <- as.matrix(transition_probs)
      transition_probs <- unname(simplify2array(transition_probs))
      
      initial_probs <- as.matrix(initial_probs)
      initial_probs <- unname(simplify2array(transition_probs))
      
      for(i in seq_along(emission_probs)) {
        assign(
          paste0("Cluster_", i),
          emission_probs[[i]],
          envir = .GlobalEnv
        )
      }
      cluster_vars <- mget(ls(pattern = "^Cluster_[0-9]+$"))
      
      # Number of hidden states
      number_of_states <- fit$model$n_states
      
      logLik_value <- logLik(fit$model)
      AIC_value <- AIC(fit$model)
      BIC_value <- BIC(fit$model)
      
      # Write everything to one MAT file
      cluster_args <- setNames(
        mget(paste0("Cluster_", 1:Number_of_clusters)),
        paste0("Cluster_", 1:Number_of_clusters)
      )
      
      mat_args <- c(
        list(
          con = fname,
          transition_probs = transition_probs,
          initial_probs = initial_probs,
          pp = pp,
          paths_compact_mat = paths_compact_mat,
          logLik_value = as.numeric(logLik_value),
          AIC_value = AIC_value,
          BIC_value = BIC_value
        ),
        cluster_args
      )
      
      do.call(R.matlab::writeMat, mat_args)
      
      cat("Saved:", fname, "\n")
    
      rm(fit,model,paths,paths_compact,paths_compat_mat,emission_probs)
      
    }
  }
}

