source("helper.R")

command_args <- commandArgs(trailingOnly=TRUE)
#command_args <- c(100,10,2,"rfr","log","conqur",
#"/Users/lynngao/Desktop/abundance_metadata/count_data/centrifuge_count_hannigan_ctr.rds",
#"/Users/lynngao/Desktop/abundance_metadata/count_data/centrifuge_count_yu_ctr.rds", 
#"/Users/lynngao/Desktop/abundance_metadata/count_data/centrifuge_count_feng_ctr.rds")
if(length(command_args)!=10){stop("Not enough input parameters!")} 
#population difference alpha(0,0.5,1), library size factor lambda(0.5,1,2), sample size, 
#genes, disease effect factor(1,2,5), dataset1, dataset2, saved table location
sample_size = as.numeric(command_args[1])
num_gene = as.numeric(command_args[2])
overlap_gene = as.numeric(command_args[3])
method = as.character(command_args[4]) #ml method
phenotype = as.character(command_args[5]) #type of phenotype
norm_method = as.character(command_args[6]) #combat or conqur

count_data1 <- readRDS(command_args[7])
count_data1 <- count_data1[,which(colSums(count_data1) != 0)]
count_data2 <- readRDS(command_args[8])
count_data2 <- count_data2[,which(colSums(count_data2) != 0)]
count_data3 <- readRDS(command_args[9])
count_data3 <- count_data3[,which(colSums(count_data3) != 0)]

#get fixed library size as median of all samples
library_size = 1000000

top_otu1 = colnames(count_data1)[order(apply(count_data1, 2, var), decreasing = TRUE)[1:1000]]
top_otu2 = colnames(count_data2)[order(apply(count_data2, 2, var), decreasing = TRUE)[1:1000]]
count_data1 = count_data1[,colnames(count_data1)%in%top_otu1]
count_data2 = count_data2[,colnames(count_data2)%in%top_otu2]

#get union of species
union_species = unique(c(top_otu1, top_otu2))
#union_species = unique(c(colnames(count_data1), colnames(count_data2)))
for (i in 1:length(union_species)){
  if (union_species[i] %!in% colnames(count_data1)){
    count_data1 = cbind(count_data1, 0)
    colnames(count_data1)[ncol(count_data1)] = union_species[i]  
  }
  if (union_species[i] %!in% colnames(count_data2)){
    count_data2 = cbind(count_data2, 0)
    colnames(count_data2)[ncol(count_data2)] = union_species[i]  
  }
  if (union_species[i] %!in% colnames(count_data3)){
    count_data3 = cbind(count_data3, 0)
    colnames(count_data3)[ncol(count_data3)] = union_species[i]  
  }
}
count_data1 = count_data1[,colnames(count_data1)%in%union_species]
count_data1 = count_data1[,sort(colnames(count_data1))]
count_data2 = count_data2[,colnames(count_data2)%in%union_species]
count_data2 = count_data2[,sort(colnames(count_data2))]
count_data3 = count_data3[,colnames(count_data3)%in%union_species]
count_data3 = count_data3[,sort(colnames(count_data3))]

#real probabilties
prob_data1 = as.vector(colMeans(count_data1)/sum(colMeans(count_data1)))
prob_data2 = as.vector(colMeans(count_data2)/sum(colMeans(count_data2)))
prob_data3 = as.vector(colMeans(count_data3)/sum(colMeans(count_data3)))


otu_list = intersect(which(prob_data1!=0), which(prob_data2!=0))
set.seed(1)
otu_selected = sample(otu_list, num_gene)
otu_list = otu_list[!otu_list %in% otu_selected]
set.seed(1)
otu_extra = sample(otu_list, 8)
if (overlap_gene == num_gene){
  otu_selected1 = otu_selected2 = otu_selected
}else{
  otu_selected1 = otu_selected
  otu_selected2 = c(otu_selected1[1:round(overlap_gene / 2)], otu_selected1[(round(num_gene / 2) + 1):(round(num_gene / 2) + round(overlap_gene / 2))])
  otu_selected2 = c(otu_selected2, otu_extra[1:(num_gene - overlap_gene)])
}

#change to relative abundance
log_relative_abundance <- function(ds){
  ds = ds/rowSums(ds)
  min = min(apply(ds[,1:ncol(ds)], 1, function(x) min(x[x>0])))
  ds[ds == 0] = min*0.65
  ds = log(ds)
  return(ds)
}

generate_samples <- function(prob, count_data, sample_size, otu_selected, coefs, phenotype, norm_method){
  simulation = rmultinom(sample_size, size = library_size, prob = prob)
  simulation_adjust = as.data.frame(t(simulation))
  colnames(simulation_adjust) = colnames(count_data)
  if (norm_method == "conqur"){
    simulation_adjust_count = simulation_adjust
  }
  simulation_adjust = log_relative_abundance(simulation_adjust)
  noise = rnorm(nrow(simulation_adjust))
  if (phenotype == "linear"){
    for (i in 1:sample_size){
      simulation_adjust$status[i] =  as.matrix(simulation_adjust[i,otu_selected]) %*% coefs + noise[i]
    }
  }
  if (phenotype == "quadratic"){
    for (i in 1:sample_size){
      simulation_adjust$status[i] =  1/10 * as.matrix(apply(simulation_adjust[i,otu_selected], c(1,2), function(x) x^2)) %*% coefs + noise[i]
    }
  }
  if (phenotype == "inverse"){
    for (i in 1:sample_size){
      simulation_adjust$status[i] = 1e+3 / as.matrix(simulation_adjust[i,otu_selected]) %*% coefs + noise[i]
    }
  }
  if (phenotype == "log"){
    for (i in 1:sample_size){
      simulation_adjust$status[i] = 1e+15 / (1 + exp(as.matrix(simulation_adjust[i,otu_selected]) %*% coefs))  + noise[i]
    }
  }
  if (phenotype == "sin"){
    for (i in 1:sample_size){
      simulation_adjust$status[i] = 1e-15 *sin(as.matrix(simulation_adjust[i,otu_selected]) %*% coefs)* exp(as.matrix(simulation_adjust[i,otu_selected]) %*% coefs) + noise[i]
    }
  }
  
  if (norm_method == "conqur"){
    return (list(simulation_adjust, simulation_adjust_count))
  }else{
    return(simulation_adjust)
  }
}

prob_v1 = prob_data1
prob_v2 = prob_data2
prob_v3 = prob_data3

# Set a dataframe to save results
perf_df <- as.data.frame(matrix(NA, nrow=104, ncol=20))
colnames(perf_df) = c("Training1","Training2", "Training1_norm","Training2_norm","Merged", "Merged_norm",
                      "Avg", "n_Avg", "CS_Avg", "Reg_a", "Reg_s", "val_rmse", "LOSO_rmse", 
                      "Avg_norm", "n_Avg_norm", "CS_Avg_norm", "Reg_a_norm", "Reg_s_norm", "val_rmse_norm", "LOSO_rmse_norm")
rownames(perf_df) = c(1:100, "mean", "sd","median", "IQR")

set.seed(101)
coefs1 <- c(runif(num_gene - round(num_gene/2), 3, 5), runif(round(num_gene/2), -5, -3))
if (overlap_gene == num_gene){
  coefs2 = coefs1
}else{
  coefs2 = c(coefs1[1:round(overlap_gene / 2)], coefs1[(round(num_gene / 2) + 1):(round(num_gene / 2) + round(overlap_gene / 2))])
  set.seed(101)
  coefs_extra = c(runif(round((num_gene - overlap_gene)/2), 3, 5), runif(round((num_gene - overlap_gene)/2), -5, -3))
  coefs2 = c(coefs2, coefs_extra)
}

# response1 = NA
# response2 = NA
# response3 = NA
# for (i in 1:100) {
#   set.seed(i)
#   simulation1 = as.data.frame(generate_samples(prob_v1, count_data1, sample_size, otu_selected1, coefs1, phenotype)) #training set1response1 = c(response1, simulation1$status)
#   response1 = c(response1, simulation1$status)
#   set.seed(i+100)
#   simulation2 = as.data.frame(generate_samples(prob_v2, count_data2, sample_size, otu_selected1, coefs1, phenotype)) #training set2
#   response2 = c(response2, simulation2$status)
#   set.seed(i)
#   simulation3 = as.data.frame(generate_samples(prob_v3, count_data3, sample_size, otu_selected2, coefs2, phenotype)) #testing set
#   response3 = c(response3, simulation3$status)
# }
# response1 = response1[-1]
# response2 = response2[-1]
# response3 = response3[-1]
# summary(c(response1,response2,response3))

for (i in 1:100) {
  if (norm_method == "conqur"){
    set.seed(i)
    temp = generate_samples(prob_v1, count_data1, sample_size, otu_selected1, coefs1, phenotype, norm_method)
    simulation1 = as.data.frame(temp[1]) #training set1
    simulation1_count = as.data.frame(temp[2]) #training set1 count data
    set.seed(i+100)
    temp = generate_samples(prob_v2, count_data2, sample_size, otu_selected1, coefs1, phenotype, norm_method)
    simulation2 = as.data.frame(temp[1]) #training set2
    simulation2_count = as.data.frame(temp[2]) #training set2 count data
    set.seed(i)
    temp = generate_samples(prob_v3, count_data3, sample_size, otu_selected2, coefs2, phenotype, norm_method)
    simulation3 = as.data.frame(temp[1]) #testing set
    simulation3_count = as.data.frame(temp[2]) #testing set count data
  }else{
    set.seed(i)
    simulation1 = as.data.frame(generate_samples(prob_v1, count_data1, sample_size, otu_selected1, coefs1, phenotype, norm_method)) #training set1
    set.seed(i+100)
    simulation2 = as.data.frame(generate_samples(prob_v2, count_data2, sample_size, otu_selected1, coefs1, phenotype, norm_method)) #training set2
    set.seed(i)
    simulation3 = as.data.frame(generate_samples(prob_v3, count_data3, sample_size, otu_selected2, coefs2, phenotype, norm_method)) #testing set
  }
  
  ## normalize two training sets
  if(norm_method == "combat"){
    simulation1_norm = combat_norm(simulation3, simulation1)
    simulation2_norm = combat_norm(simulation3, simulation2)
  }
  if(norm_method == "conqur"){
    simulation1_norm = conqur_norm(simulation3_count, simulation1_count)
    simulation2_norm = conqur_norm(simulation3_count, simulation2_count)
    #convert to log relative abundance
    simulation1_norm = log_relative_abundance(simulation1_norm)
    simulation2_norm = log_relative_abundance(simulation2_norm)
    #add y
    simulation1_norm$status = simulation1$status
    simulation2_norm$status = simulation2$status
  }
  ## Merge two training sets into one
  merged_training = rbind(simulation1, simulation2)
  merged_training_norm = rbind(simulation1_norm, simulation2_norm)
  
  ## Split testing set into validation and testing
  testing = simulation3
  sample = sample.int(n=nrow(testing),size=floor(0.5*nrow(testing)),replace=F)
  val = testing[sample,]
  testing_val = testing[rownames(testing)%!in%rownames(val),]
  
  ## Obtain predictions from learner trained within each training set
  ml1 = ml_model(simulation1, method)
  ml2 = ml_model(simulation2, method)
  ml1_norm = ml_model(simulation1_norm, method)
  ml2_norm = ml_model(simulation2_norm, method)
  pred_prob1 = predict(ml1, testing)
  pred_prob2 = predict(ml2, testing)
  pred_sgbatch_res <- list(pred_prob1, pred_prob2)
  pred_prob1_norm = predict(ml1_norm, testing)
  pred_prob2_norm = predict(ml2_norm, testing)
  pred_sgbatch_res_norm  <- list(pred_prob1_norm, pred_prob2_norm)
  
  ## Prediction from combine training together (Merged)
  ml_merged = ml_model(merged_training, method)
  pred_merged_res = predict(ml_merged, testing)
  
  ## Prediction from training after batch adjustment (Merged combat)
  ml_merged_norm = ml_model(merged_training_norm, method)
  pred_norm_res = predict(ml_merged_norm, testing)
  
  ## Prediction from training on testing_val
  pred_sgbatch_res_testing_val <- list(predict(ml1, testing_val), predict(ml2, testing_val))
  names(pred_sgbatch_res_testing_val) <- paste0("Simulation", 1:2)
  pred_sgbatch_res_norm_testing_val  <- list(predict(ml1_norm, testing_val), predict(ml2_norm, testing_val))
  names(pred_sgbatch_res_norm_testing_val) <- paste0("Simulation", 1:2)
  
  ## Prediction on validation dataset
  pred_sgbatch_res_val <- list(predict(ml1, val), predict(ml2, val))
  names(pred_sgbatch_res_val) <- paste0("Simulation", 1:2)
  pred_sgbatch_res_val_norm <- list(predict(ml1_norm, val), predict(ml2_norm, val))
  names(pred_sgbatch_res_val_norm) <- paste0("Simulation", 1:2)
  
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
  pred_N_avg <- pred_mat %*% (as.matrix(rbind(nrow(simulation1), nrow(simulation2))) / (nrow(simulation1) + nrow(simulation2)))
  pred_N_avg_norm <- pred_mat_norm %*% (as.matrix(rbind(nrow(simulation1_norm), nrow(simulation2_norm))) / (nrow(simulation1_norm) + nrow(simulation2_norm)))
  
  # CS-Avg: replicability weights
  cs_zmat <- CS_zmatrix(n_batch=2, training = list(simulation1, simulation2), perf_name="rmse", method)
  cs_weights_seq <- CS_weight(cs_zmat)
  pred_cs_avg <- pred_mat %*% cs_weights_seq
  
  cs_zmat_norm <- CS_zmatrix(n_batch=2, training = list(simulation1_norm, simulation2_norm), perf_name="rmse", method)
  cs_weights_seq_norm <- CS_weight(cs_zmat_norm)
  pred_cs_avg_norm <- pred_mat_norm %*% cs_weights_seq_norm
  
  
  # Reg-a: use each function to predict on one study, bind predictions and do regression
  reg_ssl_res <- Reg_SSL_pred(n_batch=2, training = list(simulation1, simulation2), method)
  reg_a_beta <- Reg_a_weight(coef_mat=do.call(rbind, reg_ssl_res$coef), 
                             n_seq=sapply(list(simulation1, simulation2), nrow))
  pred_reg_a <- pred_mat %*% reg_a_beta
  
  reg_ssl_res_norm <- Reg_SSL_pred(n_batch=2, training = list(simulation1_norm, simulation2_norm), method)
  reg_a_beta_norm <- Reg_a_weight(coef_mat=do.call(rbind, reg_ssl_res_norm$coef), 
                                    n_seq=sapply(list(simulation1_norm, simulation2_norm), nrow))
  pred_reg_a_norm <- pred_mat_norm %*% reg_a_beta_norm
  
  
  # Reg-s:
  stacked_pred <- do.call(rbind, reg_ssl_res$pred)
  stacked_label = c(simulation1$status, simulation2$status)
  reg_s_beta <- pnnls(a=stacked_pred, b=as.numeric(stacked_label), sum=1)$x
  #reg_s_beta <- reg_s_beta / sum(reg_s_beta)
  pred_reg_s <- pred_mat %*% reg_s_beta
  
  stacked_pred_norm <- do.call(rbind, reg_ssl_res_norm$pred)
  reg_s_beta_norm <- pnnls(a=stacked_pred, b=as.numeric(stacked_label), sum=1)$x
  #reg_s_beta_norm <- reg_s_beta_norm / sum(reg_s_beta_norm)
  pred_reg_s_norm <- pred_mat_norm %*% reg_s_beta_norm
  
  
  # Weighting by validation set performance
  coef1 = RMSE(val$status, pred_sgbatch_res_val$Simulation1)
  coef2 = RMSE(val$status, pred_sgbatch_res_val$Simulation2)
  pred_val_auc = as.matrix(1/coef1*pred_mat_testing_val[,1] + 1/coef2*pred_mat_testing_val[,2])/(1/coef1+1/coef2)
  
  coef1_norm = RMSE(val$status, pred_sgbatch_res_val_norm$Simulation1)
  coef2_norm = RMSE(val$status, pred_sgbatch_res_val_norm$Simulation2)
  pred_val_auc_norm = as.matrix(1/coef1_norm*pred_mat_norm_testing_val[,1] + 1/coef2_norm*pred_mat_norm_testing_val[,2])/(1/coef1_norm+1/coef2_norm)
  
  # LOSO
  LOSO = as.data.frame(matrix(NA, ncol = 4, nrow = nrow(testing)))
  colnames(LOSO) = c("auc_ds1","auc_ds2","auc_ds1_norm","auc_ds2_norm")
  for (k in 1:nrow(LOSO)){
    LOSO[k,1] = 1/RMSE(testing$status[-k], pred_prob1[-k])
    LOSO[k,2] = 1/RMSE(testing$status[-k], pred_prob2[-k])
    LOSO[k,3] = 1/RMSE(testing$status[-k], pred_prob1_norm[-k])
    LOSO[k,4] = 1/RMSE(testing$status[-k], pred_prob2_norm[-k])
  }
  pred_LOSO = (LOSO[,1]*pred_prob1 + LOSO[,2]*pred_prob2)/ (LOSO[,1]+LOSO[,2])
  pred_LOSO_norm = (LOSO[,3]*pred_prob1_norm + LOSO[,4]*pred_prob2_norm)/ (LOSO[,3]+LOSO[,4])
  
  ####  Evaluate performance
  tst_scores <- c(list(Training1 = pred_prob1, Training2 = pred_prob2, Training1_norm = pred_prob1_norm, Training2_norm = pred_prob2_norm, 
                       Merged=pred_merged_res, Merged_norm=pred_norm_res, 
                       Avg=pred_avg, n_Avg=pred_N_avg, CS_Avg=pred_cs_avg, Reg_a=pred_reg_a, Reg_s=pred_reg_s, val_auc = pred_val_auc, LOSO_auc = pred_LOSO, 
                       #rank_mean = rank_mat$mean, rank_geo_mean = rank_mat$geo_mean, rank_Stuart=rank_mat$Q, rank_RRA=rank_mat$rho, rank_Bayesian=rank_mat$Bayesian,
                       Avg_norm=pred_avg_norm, n_Avg_norm=pred_N_avg_norm, CS_Avg_norm=pred_cs_avg_norm, Reg_a_norm=pred_reg_a_norm, Reg_s_norm=pred_reg_s_norm, val_auc_norm = pred_val_auc_norm, LOSO_auc_norm = pred_LOSO_norm)) 
                       #rank_mean_norm = rank_mat_norm$mean, rank_geo_mean_norm = rank_mat_norm$geo_mean, rank_Stuart_norm=rank_mat_norm$Q, rank_RRA_norm=rank_mat_norm$rho, rank_Bayesian_norm=rank_mat_norm$Bayesian))
  
  perf_df[i,] <- sapply(tst_scores, function(preds){
    testing_auc <- RMSE(testing$status, as.vector(preds))
    return(testing_auc[1])})
  perf_df[i,12] <- RMSE(testing_val$status, as.vector(tst_scores$val_auc))
  perf_df[i,19] <- RMSE(testing_val$status, as.vector(tst_scores$val_auc_norm))
}

# Get mean and sd
perf_df[101,] = colMeans(perf_df[1:100,])
perf_df[102,] = apply(perf_df[1:100,], 2, sd)
perf_df[103,] = apply(perf_df[1:100,], 2, median)
perf_df[104,] = apply(perf_df[1:100,], 2, IQR)

####  Output results
file_name <- sprintf('scenario3_%s_rmse_%s_%s_s%s_g%s_overlap%s_test%s_logRelAbun.csv', 
                     gsub('.', '', method, fixed=T),
                     gsub('.', '', phenotype, fixed=T),
                     gsub('.', '', norm_method, fixed=T),
                     gsub('.', '', sample_size, fixed=T), gsub('.', '', num_gene, fixed=T), gsub('.', '', overlap_gene, fixed=T),
                     gsub('.', '', sub(".*_count_*(.*?)_ctr.*", "\\1", command_args[9]), fixed=T))
write.csv(perf_df, paste0(command_args[10],file_name))
