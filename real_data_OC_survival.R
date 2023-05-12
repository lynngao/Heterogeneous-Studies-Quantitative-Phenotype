source("helper.R")
source(system.file("extdata", "patientselection.config",package="curatedOvarianData"))
sapply(ls(), function(x) if(!x %in% c("remove.samples", "duplicates")) print(get(x)))
source(system.file("extdata", "createEsetList.R", package = "curatedOvarianData"))


data(E.MTAB.386_eset) #Bentink
data(GSE30161_eset) #Ferriss
data(GSE17260_eset) #Yoshihara

d1 = GSE30161_eset
d3 = E.MTAB.386_eset
d2 = GSE17260_eset

ds1 = as.data.frame(t(exprs(GSE30161_eset))) #Ferriss
ds3 = as.data.frame(t(exprs(E.MTAB.386_eset))) #Bentink
ds2 = as.data.frame(t(exprs(GSE17260_eset))) #Yoshihara

trim_count_data <- function(features, data){
  for (i in 1:length(features)){
    if (features[i] %!in% colnames(data)){
      data = cbind(data, 0)
      colnames(data)[ncol(data)] = features[i]
    }
  }
  data = data[,colnames(data)%in%features]
  data = data[,sort(colnames(data))]
  return(as.data.frame(data))
}

features = intersect(intersect(colnames(ds1), colnames(ds2)), colnames(ds3))
ds1 = trim_count_data(features, ds1)
ds2 = trim_count_data(features, ds2)
ds3 = trim_count_data(features, ds3)

curated_ds <- function(ds, original_ds) {
  ds$tumorstage = original_ds$tumorstage
  ds$debulking = original_ds$debulking
  ds$time = original_ds$days_to_death
  ds$status = original_ds$vital_status
  ds$status[ds$status == "living"] = 0
  ds$status[ds$status == "deceased"] = 1
  #ds$status = factor(ds$status)
  ds$status = as.integer(ds$status)
  return(ds)
}

ds1 = na.omit(curated_ds(ds1, d1))
ds2 = na.omit(curated_ds(ds2, d2))
ds3 = na.omit(curated_ds(ds3, d3))

mod1 <- model.matrix(~as.factor(tumorstage) + as.factor(debulking), data=rbind(ds3, ds1))
batch1 <- as.factor(c(rep(0, nrow(ds3)), rep(1, nrow(ds1))))
combat_edata1 <- ComBat(dat=t(rbind(ds3, ds1)[,1:length(features)]), batch=batch1, mod=mod1)
#combat_edata1 <- ComBat(dat=t(rbind(ds3, ds1)[,1:length(features)]), batch=batch1)
mod2 <- model.matrix(~as.factor(tumorstage) + as.factor(debulking), data=rbind(ds3, ds2))
batch2 <- as.factor(c(rep(0, nrow(ds3)), rep(1, nrow(ds2))))
combat_edata2 <- ComBat(dat=t(rbind(ds3, ds2)[,1:length(features)]), batch=batch2, mod=mod2)
#combat_edata2 <- ComBat(dat=t(rbind(ds3, ds2)[,1:length(features)]), batch=batch2)

ds1_norm = as.data.frame(t(combat_edata1))[(nrow(ds3) + 1): nrow(t(combat_edata1)),]
ds2_norm = as.data.frame(t(combat_edata2))[(nrow(ds3) + 1): nrow(t(combat_edata2)),]
ds1_norm$tumorstage = ds1$tumorstage
ds1_norm$debulking = ds1$debulking
ds1_norm$time = ds1$time
ds1_norm$status = ds1$status
ds2_norm$tumorstage = ds2$tumorstage
ds2_norm$debulking = ds2$debulking
ds2_norm$time = ds2$time
ds2_norm$status = ds2$status

## Merge two training sets into one
merged_training = rbind(ds1, ds2)
merged_training_norm = rbind(ds1_norm, ds2_norm)


# Set a dataframe to save results
perf_df <- as.data.frame(matrix(NA, nrow=14, ncol=20))
colnames(perf_df) = c("Training1","Training2", "Training1_norm","Training2_norm","Merged", "Merged_norm",
                      "Avg", "n_Avg", "CS_Avg", "Reg_a", "Reg_s", "val_rmse", "LOSO_rmse", 
                      "Avg_norm", "n_Avg_norm", "CS_Avg_norm", "Reg_a_norm", "Reg_s_norm", "val_rmse_norm", "LOSO_rmse_norm")
rownames(perf_df) = c(1:10, "mean", "sd", "median", "IQR")

method = "coxboost"

for (i in 1:10) {
  ## Split testing set into validation and testing
  testing = ds3
  sample = sample.int(n=nrow(testing),size=floor(0.5*nrow(testing)),replace=F)
  val = testing[sample,]
  testing_val = testing[rownames(testing)%!in%rownames(val),]
  
  ## Obtain predictions from learner trained within each training set
  ml1 = CoxBoost(time=ds1$time,status=ds1$status,x=as.matrix(ds1[,1:length(features)]))
  ml2 = CoxBoost(time=ds2$time,status=ds2$status,x=as.matrix(ds2[,1:length(features)]))
  ml1_norm = CoxBoost(time=ds1_norm$time,status=ds1_norm$status,x=as.matrix(ds1_norm[,1:length(features)]))
  ml2_norm = CoxBoost(time=ds2_norm$time,status=ds2_norm$status,x=as.matrix(ds2_norm[,1:length(features)]))
  pred_prob1 = as.vector(predict(ml1, as.matrix(testing[,1:length(features)]), testing$time,testing$status))
  pred_prob2 = as.vector(predict(ml2, as.matrix(testing[,1:length(features)]), testing$time,testing$status))
  pred_sgbatch_res <- list(pred_prob1, pred_prob2)
  pred_prob1_norm = as.vector(predict(ml1_norm, as.matrix(testing[,1:length(features)]), testing$time,testing$status))
  pred_prob2_norm = as.vector(predict(ml2_norm, as.matrix(testing[,1:length(features)]), testing$time,testing$status))
  pred_sgbatch_res_norm  <- list(pred_prob1_norm, pred_prob2_norm)
  
  ## Prediction from combine training together (Merged)
  ml_merged = CoxBoost(time=merged_training$time,status=merged_training$status,x=as.matrix(merged_training[,1:length(features)]))
  pred_merged_res =  as.vector(predict(ml_merged, as.matrix(testing[,1:length(features)]), testing$time,testing$status))
  
  ## Prediction from training after batch adjustment (Merged combat)
  ml_merged_norm = CoxBoost(time=merged_training_norm$time,status=merged_training_norm$status,x=as.matrix(merged_training_norm[,1:length(features)]))
  pred_norm_res = as.vector(predict(ml_merged_norm, as.matrix(testing[,1:length(features)]), testing$time,testing$status))
  
  ## Prediction from training on testing_val
  pred_sgbatch_res_testing_val <- list(as.vector(predict(ml1, as.matrix(testing_val[,1:length(features)]), testing_val$time,testing_val$status)), 
                                       as.vector(predict(ml2, as.matrix(testing_val[,1:length(features)]), testing_val$time,testing_val$status)))
  names(pred_sgbatch_res_testing_val) <- paste0("ds", 1:2)
  pred_sgbatch_res_norm_testing_val  <- list(as.vector(predict(ml1_norm, as.matrix(testing_val[,1:length(features)]), testing_val$time,testing_val$status)), 
                                             as.vector(predict(ml2_norm, as.matrix(testing_val[,1:length(features)]), testing_val$time,testing_val$status)))
  names(pred_sgbatch_res_norm_testing_val) <- paste0("ds", 1:2)
  
  ## Prediction on validation dataset
  pred_sgbatch_res_val <- list(as.vector(predict(ml1, as.matrix(val[,1:length(features)]), val$time,val$status)), 
                               as.vector(predict(ml2, as.matrix(val[,1:length(features)]), val$time,val$status)))
  names(pred_sgbatch_res_val) <- paste0("ds", 1:2)
  pred_sgbatch_res_val_norm <- list(as.vector(predict(ml1_norm, as.matrix(val[,1:length(features)]), val$time,testing_val$status)), 
                                    as.vector(predict(ml2_norm, as.matrix(val[,1:length(features)]), val$time,testing_val$status)))
  names(pred_sgbatch_res_val_norm) <- paste0("ds", 1:2)
  
  ##  Aggregate with different weights
  pred_test_lst <- lapply(pred_sgbatch_res, function(tmp){return(tmp)})
  pred_mat <- do.call(cbind, pred_test_lst)
  pred_test_lst_norm <- lapply(pred_sgbatch_res_norm, function(tmp){return(tmp)})
  pred_mat_norm <- do.call(cbind, pred_test_lst_norm)
  
  pred_test_lst_testing_val <- lapply(pred_sgbatch_res_testing_val, function(tmp){return(tmp)})
  pred_mat_testing_val <- do.call(cbind, pred_test_lst_testing_val)
  pred_test_lst_norm_testing_val <- lapply(pred_sgbatch_res_norm_testing_val, function(tmp){return(tmp)})
  pred_mat_norm_testing_val <- do.call(cbind, pred_test_lst_norm_testing_val)
  
  # Avg: simple average
  pred_avg <- rowMeans(pred_mat)
  pred_avg_norm <- rowMeans(pred_mat_norm)
  # n-Avg: sample-size-weighted average
  pred_N_avg <- pred_mat %*% (as.matrix(rbind(nrow(ds1), nrow(ds2))) / (nrow(ds1) + nrow(ds2)))
  pred_N_avg_norm <- pred_mat_norm %*% (as.matrix(rbind(nrow(ds1_norm), nrow(ds2_norm))) / (nrow(ds1_norm) + nrow(ds2_norm)))
  
  # CS-Avg: replicability weights
  cs_zmat <- CS_zmatrix(n_batch=2, training = list(ds1, ds2), perf_name=NULL, method)
  cs_weights_seq <- CS_weight(cs_zmat)
  pred_cs_avg <- pred_mat %*% cs_weights_seq
  
  cs_zmat_norm <- CS_zmatrix(n_batch=2, training = list(ds1_norm, ds2_norm), perf_name=NULL, method)
  cs_weights_seq_norm <- CS_weight(cs_zmat_norm)
  pred_cs_avg_norm <- pred_mat_norm %*% cs_weights_seq_norm
  
  
  # Reg-a: use each function to predict on one study, bind predictions and do regression
  reg_ssl_res <- Reg_SSL_pred(n_batch=2, training = list(ds1, ds2), method)
  reg_a_beta <- Reg_a_weight(coef_mat=do.call(rbind, reg_ssl_res$coef), 
                             n_seq=sapply(list(ds1, ds2), nrow))
  pred_reg_a <- pred_mat %*% reg_a_beta
  
  reg_ssl_res_norm <- Reg_SSL_pred(n_batch=2, training = list(ds1_norm, ds2_norm), method)
  reg_a_beta_norm <- Reg_a_weight(coef_mat=do.call(rbind, reg_ssl_res_norm$coef), 
                                  n_seq=sapply(list(ds1_norm, ds2_norm), nrow))
  pred_reg_a_norm <- pred_mat_norm %*% reg_a_beta_norm
  
  
  # Reg-s:
  stacked_pred <- do.call(rbind, reg_ssl_res$pred)
  stacked_label = c(ds1$time, ds2$time)
  reg_s_beta <- pnnls(a=stacked_pred, b=as.numeric(stacked_label), sum=1)$x
  pred_reg_s <- pred_mat %*% reg_s_beta
  
  stacked_pred_norm <- do.call(rbind, reg_ssl_res_norm$pred)
  reg_s_beta_norm <- pnnls(a=stacked_pred_norm, b=as.numeric(stacked_label), sum=1)$x
  pred_reg_s_norm <- pred_mat_norm %*% reg_s_beta_norm
  
  
  # Weighting by validation set performance
  coef1 = get.cindex(time = val$time, censoring = val$status, predicted = pred_sgbatch_res_val$ds1)
  coef2 = get.cindex(time = val$time, censoring = val$status, predicted = pred_sgbatch_res_val$ds2)
  pred_val_auc = as.matrix(1/coef1*pred_mat_testing_val[,1] + 1/coef2*pred_mat_testing_val[,2])/(1/coef1+1/coef2)
  
  coef1_norm = get.cindex(time = val$time, censoring = val$status, predicted = pred_sgbatch_res_val_norm$ds1)
  coef2_norm = get.cindex(time = val$time, censoring = val$status, predicted = pred_sgbatch_res_val_norm$ds2)
  pred_val_auc_norm = as.matrix(1/coef1_norm*pred_mat_norm_testing_val[,1] + 1/coef2_norm*pred_mat_norm_testing_val[,2])/(1/coef1_norm+1/coef2_norm)
  
  
  # LOSO
  LOSO = as.data.frame(matrix(NA, ncol = 4, nrow = nrow(testing)))
  colnames(LOSO) = c("auc_ds1","auc_ds2","auc_ds1_norm","auc_ds2_norm")
  for (k in 1:nrow(LOSO)){
    LOSO[k,1] = 1/get.cindex(time = testing$time[-k], censoring = testing$status[-k], predicted = pred_prob1[-k])
    LOSO[k,2] = 1/get.cindex(time = testing$time[-k], censoring = testing$status[-k], predicted = pred_prob2[-k])
    LOSO[k,3] = 1/get.cindex(time = testing$time[-k], censoring = testing$status[-k], predicted = pred_prob1_norm[-k])
    LOSO[k,4] = 1/get.cindex(time = testing$time[-k], censoring = testing$status[-k], predicted = pred_prob2_norm[-k])
  }
  pred_LOSO = (LOSO[,1]*pred_prob1 + LOSO[,2]*pred_prob2)/ (LOSO[,1]+LOSO[,2])
  pred_LOSO_norm = (LOSO[,3]*pred_prob1_norm + LOSO[,4]*pred_prob2_norm)/ (LOSO[,3]+LOSO[,4])
  
  
  ####  Evaluate performance
  tst_scores <- c(list(Training1 = pred_prob1, Training2 = pred_prob2, Training1_norm = pred_prob1_norm, Training2_norm = pred_prob2_norm, Merged=pred_merged_res, Merged_norm=pred_norm_res, 
                       Avg=pred_avg, n_Avg=pred_N_avg, CS_Avg=pred_cs_avg, Reg_a=pred_reg_a, Reg_s=pred_reg_s, val_auc = pred_val_auc, LOSO_auc = pred_LOSO, 
                       Avg_norm=pred_avg_norm, n_Avg_norm=pred_N_avg_norm, CS_Avg_norm=pred_cs_avg_norm, Reg_a_norm=pred_reg_a_norm, Reg_s_norm=pred_reg_s_norm, val_auc_norm = pred_val_auc_norm, LOSO_auc_norm = pred_LOSO_norm)) 
  
  for (col in 1:20){
    if (col == 12 | col == 19){
      perf_df[i,col] <- get.cindex(time = testing_val$time, censoring = testing_val$status, predicted = tst_scores[[col]])
    }else{
      perf_df[i,col] <- get.cindex(time = testing$time, censoring = testing$status, predicted = tst_scores[[col]])
    } 
  }
  
}


# Get mean and sd
perf_df[11,] = colMeans(perf_df[1:10,])
perf_df[12,] = apply(perf_df[1:10,], 2, sd)
perf_df[13,] = apply(perf_df[1:10,], 2, median)
perf_df[14,] = apply(perf_df[1:10,], 2, IQR)

####  Output results
file_name <- sprintf('quantitative_real_data_survival_rmse_%s_combat_logRelAbun_covar.csv', 
                     gsub('.', '', method, fixed=T))
write.csv(perf_df, paste0(command_args[1],file_name))
