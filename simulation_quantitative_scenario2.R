source("helper.R")

####  Simulate data from same background population
command_args <- commandArgs(trailingOnly=TRUE)
#command_args <- c(100,10,"/Users/lynngao/Desktop/abundance_metadata/count_data/centrifuge_count_yu_ctr.rds",
                  #2,100,100,3,1,"rfr","linear","conqur")
if(length(command_args)!=12){stop("Not enough input parameters!")}
test_sample_size = as.numeric(command_args[1])
num_gene = as.numeric(command_args[2])
count_data1 <- readRDS(command_args[3]) #Asian population
phenotype = as.character(command_args[10]) #type of phenotype
norm_method = as.character(command_args[11]) #combat or conqur

top_otu1 = colnames(count_data1)[order(apply(count_data1, 2, var), decreasing = TRUE)[1:1000]]
count_data1 = count_data1[,colnames(count_data1)%in%top_otu1]
count_data1 = count_data1[,sort(colnames(count_data1))]

#get fixed library size as median of all samples
library_size = 1000000

#real probabilties
prob_data1 = as.vector(colMeans(count_data1)/sum(colMeans(count_data1)))

otu_list = which(prob_data1!=0)
set.seed(1)
otu_selected = sample(otu_list, num_gene)

#change to relative abundance
log_relative_abundance <- function(ds){
  ds = ds/rowSums(ds)
  min = min(apply(ds[,1:ncol(ds)], 1, function(x) min(x[x>0])))
  ds[ds == 0] = min*0.65
  ds = log(ds)
  return(ds)
}

generate_samples <- function(prob, count_data, sample_size, num_gene, otu_selected, coefs, phenotype, norm_method){
  prob = rdirichlet(sample_size, 100000 * prob)
  
  simulation = matrix(nrow=length(top_otu1))
  for (i in 1:nrow(prob)){
    set.seed(i)
    sim = rmultinom(1, size = library_size, prob = as.vector(prob[i,]))
    simulation = cbind(simulation, c(sim[,1]))
  }
  simulation = simulation[,-1]
  simulation_adjust = as.data.frame(t(simulation))
  colnames(simulation_adjust) = colnames(count_data)
  simulation_adjust = log_relative_abundance(simulation_adjust)
  noise = rnorm(nrow(simulation_adjust))
  if (phenotype == "linear"){
    for (i in 1:sample_size){
      simulation_adjust$status[i] =  10 * as.matrix(simulation_adjust[i,otu_selected]) %*% coefs + noise[i]
    }
  }
  if (phenotype == "quadratic"){
    for (i in 1:sample_size){
      simulation_adjust$status[i] =  1/10 * as.matrix(apply(simulation_adjust[i,otu_selected], c(1,2), function(x) x^2)) %*% coefs + noise[i]
    }
  }
  if (phenotype == "inverse"){
    for (i in 1:sample_size){
      simulation_adjust$status[i] = 10 / as.matrix(simulation_adjust[i,otu_selected]) %*% coefs + noise[i]
    }
  }
  if (phenotype == "log"){
    for (i in 1:sample_size){
      simulation_adjust$status[i] = 100 / (1 + exp(as.matrix(simulation_adjust[i,otu_selected]) %*% coefs))  + noise[i]
    }
  }
  if (phenotype == "sin"){
    for (i in 1:sample_size){
      simulation_adjust$status[i] = 100 * sin(as.matrix(simulation_adjust[i,otu_selected]) %*% coefs)* exp(as.matrix(simulation_adjust[i,otu_selected]) %*% coefs) + noise[i]
    }
  }
  
  return(simulation_adjust)
}

prob_v1 = prob_data1

####  Parameters 
#command_args <- commandArgs(trailingOnly=TRUE)  
#command_args <- c("20", "5", "1")
#if(length(command_args)!=3){stop("Not enough input parameters!")}

## Degree of batch effect (strength of signal)
N_batch <- as.numeric(command_args[4])
N_sample_size <- c(as.numeric(command_args[5]),as.numeric(command_args[6]))
max_batch_mean <- as.numeric(command_args[7])
max_batch_var <- as.numeric(command_args[8])
#N_sample_size <- as.numeric(command_args[1])   # number of samples per batch (case + control)
#max_batch_mean <- as.numeric(command_args[2]) 
# median 0.12, range(0, 12), recommended values 0 0.3 0.5 (additive)
#max_batch_var <- as.numeric(command_args[3]) 
# median 0.02, range(0.002, 0.33), recommended values 1 2 4 (multiplicative)
hyper_pars <- list(hyper_mu=seq(from=-max_batch_mean, to=max_batch_mean, length.out=N_batch),  
                   hyper_sd=sqrt(rep(0.01, N_batch)),
                   hyper_alpha=mv2ab(m=seq(from=1/max_batch_var, to=max_batch_var, length.out=N_batch), 
                                     v=rep(0.01, N_batch))$alpha, 
                   hyper_beta=mv2ab(m=seq(from=1/max_batch_var, to=max_batch_var, length.out=N_batch), 
                                    v=rep(0.01, N_batch))$beta)  
cat("\nBatch changes\n");
print(hyper_pars$hyper_mu);
print(hyper_pars$hyper_beta/(hyper_pars$hyper_alpha-1))
#print((hyper_pars$hyper_beta^2)/((hyper_pars$hyper_alpha-1)^2*(hyper_pars$hyper_alpha-2)))

## Pipeline
iterations <- 100
#exp_name <- sprintf('batchN%s_m%s_v%s', N_sample_size, 
#gsub('.', '', max_batch_mean, fixed=T), gsub('.', '', max_batch_var, fixed=T)) 
exp_name <- sprintf('batch_m%s_v%s', 
                    gsub('.', '', max_batch_mean, fixed=T), gsub('.', '', max_batch_var, fixed=T))
perf_measures <- "rmse"    
method = as.character(command_args[9]) #ml method
l_type <- method

####  Run pipeline
#c; #l_type = learner_types[2]
#start_time <- Sys.time()
perf_df <- as.data.frame(matrix(NA, nrow=104, ncol=21))
colnames(perf_df) = c("NoBatch","Batch","Batch1","Batch2", "Batch1_norm","Batch2_norm","Merged_norm",
                      "Avg", "n_Avg", "CS_Avg", "Reg_a", "Reg_s", "val_rmse", "LOSO_rmse", 
                      "Avg_norm", "n_Avg_norm", "CS_Avg_norm", "Reg_a_norm", "Reg_s_norm", "val_rmse_norm", "LOSO_rmse_norm")
rownames(perf_df) = c(1:100, "mean", "sd", "median", "IQR")

set.seed(101)
coefs <- c(runif(num_gene - round(num_gene/2), 3, 5), runif(round(num_gene/2), -5, -3))


# response1 = NA
# response2 = NA
# response3 = NA
# for (i in 1:100) {
#   set.seed(i)
#   simulation1 = as.data.frame(generate_samples(prob_v1, count_data1, N_sample_size[1], num_gene, otu_selected, coefs, phenotype)) #training set1
#   response1 = c(response1, simulation1$status)
#   set.seed(i+100)
#   simulation2 = as.data.frame(generate_samples(prob_v1, count_data1, N_sample_size[2], num_gene, otu_selected, coefs, phenotype)) #training set2
#   response2 = c(response2, simulation2$status)
#   set.seed(i+200)
#   simulation3 = as.data.frame(generate_samples(prob_v1, count_data1, test_sample_size, num_gene, otu_selected, coefs, phenotype))#testing set
#   response3 = c(response3, simulation3$status)
# }
# response1 = response1[-1]
# response2 = response2[-1]
# response3 = response3[-1]
# summary(c(response1,response2,response3))


for(i in 1:iterations){
  set.seed(i)
  simulation1 = as.data.frame(generate_samples(prob_v1, count_data1, N_sample_size[1], num_gene, otu_selected, coefs, phenotype, norm_method)) #training set1
  set.seed(i+100)
  simulation2 = as.data.frame(generate_samples(prob_v1, count_data1, N_sample_size[2], num_gene, otu_selected, coefs, phenotype, norm_method)) #training set2
  set.seed(i+200)
  simulation3 = as.data.frame(generate_samples(prob_v1, count_data1, test_sample_size, num_gene, otu_selected, coefs, phenotype, norm_method)) #testing set
  set.seed(i+300)
  validation = as.data.frame(generate_samples(prob_v1, count_data1, test_sample_size, num_gene, otu_selected, coefs, phenotype, norm_method)) #validation set
  
  # Define testing set and validation set
  testing = simulation3
  sample = sample.int(n=nrow(testing),size=floor(0.5*nrow(testing)),replace=F)
  val = validation
  testing_val = testing

  tst_scores_modlst <- cs_zmat_lst <- list()  # list of results for each model
  
  ## Subset training set in batches
  batch1 = simulation1
  batch2 = simulation2
  rownames(batch2) = (N_sample_size[1]+1):(N_sample_size[1]+N_sample_size[2])
  batches_ind = list(rownames(batch1), rownames(batch2))
  batch = c(rep(1,N_sample_size[1]), rep(2,N_sample_size[2]))
  
  ## Remove genes with only 0 values in any batch in current training set
  g1 = colnames(batch1[colSums(batch1[,-ncol(batch1)]) != 0])
  g2 = colnames(batch2[colSums(batch2[,-ncol(batch2)]) != 0])
  g_keep <- intersect(g1, g2)
  curr_test_expr <- testing[,g_keep]
  curr_val_expr <- val[,g_keep]
  batch1 = batch1[,g_keep]
  batch2 = batch2[,g_keep]
  curr_train_expr = rbind(batch1, batch2)
  
  ## Simulate batch effect 
  set.seed(i)
  sim_batch_res <- simBatch(curr_train_expr, N_sample_size, batches_ind, batch, hyper_pars)
  train_expr_batch <- sim_batch_res$new_dat
  batch1_simulated = train_expr_batch[1:N_sample_size[1],]
  batch1_sim = batch1_simulated
  batch2_simulated = train_expr_batch[(N_sample_size[1]+1):nrow(train_expr_batch),]
  batch2_sim = batch2_simulated
  
  ##normalization on batch and testing dataset (simulation3 as reference)
  if (norm_method == "combat"){
    batch1_norm <- combat_norm(simulation3, batch1_simulated)
    rownames(batch1_norm) = rownames(batch1_simulated)
    batch2_norm <- combat_norm(simulation3, batch2_simulated)
    rownames(batch2_norm) = rownames(batch2_simulated)
  }
  if (norm_method == "conqur"){
    #convert to count data by exp(data) * library size
    simulation3[,-ncol(simulation3)] = round(exp(simulation3[,-ncol(simulation3)]) * library_size)
    batch1_simulated[,-ncol(batch1_simulated)] = round(exp(batch1_simulated[,-ncol(batch1_simulated)]) * library_size)
    batch2_simulated[,-ncol(batch2_simulated)] = round(exp(batch2_simulated[,-ncol(batch2_simulated)]) * library_size)
    
    #normalize data by using testing as reference
    batch1_norm <- conqur_norm(simulation3, batch1_simulated)
    rownames(batch1_norm) = rownames(batch1_simulated)
    temp_status = batch1_norm$status
    batch1_norm = log_relative_abundance(batch1_norm[,-ncol(batch1_norm)])
    batch1_norm$status = temp_status
    
    #convert back to log relative abundance data
    batch2_norm <- conqur_norm(simulation3, batch2_simulated)
    rownames(batch2_norm) = rownames(batch2_simulated)
    temp_status = batch2_norm$status
    batch2_norm = log_relative_abundance(batch2_norm[,-ncol(batch2_norm)])
    batch2_norm$status = temp_status
  }
  
  train_expr_norm = rbind(batch1_norm, batch2_norm)
  train_expr_norm_norm <- train_expr_norm
  
  train_expr_norm <- as.data.frame(curr_train_expr)
  test_expr_norm <- curr_test_expr
  val_expr_norm <- curr_val_expr
  train_expr_batch_whole_norm <- train_expr_batch_norm <- as.data.frame(train_expr_batch)
  
  train_lst <- lapply(batches_ind, function(ind){train_expr_batch_norm[as.character(ind),]})
  train_lst_norm <- lapply(batches_ind, function(ind){train_expr_norm_norm[as.character(ind),]})
  
  ####  Training with single model
  print(sprintf("Simulation: %s, Model: %s", i, l_type))
  
  
  ## Prediction from original training to test, without batch effect
  ml_merged <- ml_model(train_expr_norm, method)
  pred_base_res <- predict(ml_merged, test_expr_norm)
  
  
  ## Prediction from training WITH batch effect to test
  ml_merged_batch <- ml_model(train_expr_batch_whole_norm, method)
  pred_batch_res <- predict(ml_merged_batch, test_expr_norm)
  
  
  ##  Prediction from training after batch adjustment (Merged)
  ml_merged_batch_norm <- ml_model(train_expr_norm_norm, method)
  pred_norm_res <- predict(ml_merged_batch_norm, test_expr_norm)
  
  
  ## Obtain predictions from learner trained within each batch
  ml1 <- ml_model(batch1_sim, method)
  ml2 <- ml_model(batch2_sim, method)
  ml1_norm <- ml_model(batch1_norm, method)
  ml2_norm <- ml_model(batch2_norm, method)
  pred_prob1 = predict(ml1, test_expr_norm)
  pred_prob2 = predict(ml2, test_expr_norm)
  pred_prob1_norm = predict(ml1_norm, test_expr_norm)
  pred_prob2_norm = predict(ml2_norm, test_expr_norm)
  pred_sgbatch_res <- list(pred_prob1, pred_prob2)
  pred_sgbatch_res_norm <- list(pred_prob1_norm, pred_prob2_norm)
  names(pred_sgbatch_res) <- paste0("Batch", 1:N_batch)
  names(pred_sgbatch_res_norm) <- paste0("Batch", 1:N_batch, "_norm")
  
  ## Prediction from training on testing_val
  pred_sgbatch_res_testing_val <- list(predict(ml1, testing_val), predict(ml2, testing_val))
  names(pred_sgbatch_res_testing_val) <- paste0("Batch", 1:N_batch)
  pred_sgbatch_res_norm_testing_val  <- list(predict(ml1_norm, testing_val), predict(ml2_norm, testing_val))
  names(pred_sgbatch_res_norm_testing_val) <- paste0("Batch", 1:N_batch, "_norm")
  
  ##prediction on validation dataset
  pred_prob1_val = predict(ml1, val_expr_norm)
  pred_prob2_val = predict(ml2, val_expr_norm)
  pred_prob1_val_norm = predict(ml1_norm, val_expr_norm)
  pred_prob2_val_norm = predict(ml2_norm, val_expr_norm)
  pred_sgbatch_res_val <- list(pred_prob1_val, pred_prob2_val)
  names(pred_sgbatch_res_val) <- paste0("Batch", 1:N_batch)
  pred_sgbatch_res_val_norm <- list(pred_prob1_val_norm, pred_prob2_val_norm)
  names(pred_sgbatch_res_val_norm) <- paste0("Batch", 1:N_batch, "_norm")
  
  
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
  pred_N_avg <- pred_mat %*% (as.matrix(sapply(batches_ind, length)) / nrow(curr_train_expr))
  pred_N_avg_norm <- pred_mat_norm %*% (as.matrix(sapply(batches_ind, length)) / nrow(curr_train_expr))
  
  
  # CS-Avg: replicability weights
  cs_zmat <- CS_zmatrix(n_batch=N_batch, training = train_lst, perf_name="rmse", method)
  cs_weights_seq <- CS_weight(cs_zmat)
  pred_cs_avg <- pred_mat %*% cs_weights_seq
  
  cs_zmat_norm <- CS_zmatrix(n_batch=N_batch, training = train_lst_norm, perf_name="rmse", method)
  cs_weights_seq_norm <- CS_weight(cs_zmat_norm)
  pred_cs_avg_norm <- pred_mat_norm %*% cs_weights_seq_norm
  
  # Reg-a: use each function to predict on one study, bind predictions and do regression
  reg_ssl_res <- Reg_SSL_pred(n_batch=N_batch, training = train_lst, method)
  reg_a_beta <- Reg_a_weight(coef_mat=do.call(rbind, reg_ssl_res$coef), 
                             n_seq=sapply(batches_ind, length))
  pred_reg_a <- pred_mat %*% reg_a_beta
  
  reg_ssl_res_norm <- Reg_SSL_pred(n_batch=N_batch, training = train_lst_norm, method)
  reg_a_beta_norm <- Reg_a_weight(coef_mat=do.call(rbind, reg_ssl_res_norm$coef), 
                                    n_seq=sapply(batches_ind, length))
  pred_reg_a_norm <- pred_mat_norm %*% reg_a_beta_norm
  
  
  # Reg-s:
  stacked_pred <- do.call(rbind, reg_ssl_res$pred)
  stacked_label = curr_train_expr$status
  reg_s_beta <- pnnls(a=stacked_pred, b=as.numeric(stacked_label), sum=1)$x
  #reg_s_beta <- reg_s_beta / sum(reg_s_beta)
  pred_reg_s <- pred_mat %*% reg_s_beta
  
  stacked_pred_norm <- do.call(rbind, reg_ssl_res_norm$pred)
  reg_s_beta_norm <- pnnls(a=stacked_pred_norm, b=as.numeric(stacked_label), sum=1)$x
  #reg_s_beta_norm <- reg_s_beta_norm / sum(reg_s_beta_norm)
  pred_reg_s_norm <- pred_mat_norm %*% reg_s_beta_norm
  
  
  # Weighting by validation set performance
  coef1 = RMSE(val_expr_norm$status, pred_sgbatch_res_val$Batch1)
  coef2 = RMSE(val_expr_norm$status, pred_sgbatch_res_val$Batch2)
  pred_val_auc = as.matrix(1/coef1*pred_mat_testing_val[,1] + 1/coef2*pred_mat_testing_val[,2])/(1/coef1+1/coef2)
  
  coef1_norm = RMSE(val_expr_norm$status, pred_sgbatch_res_val_norm$Batch1)
  coef2_norm = RMSE(val_expr_norm$status, pred_sgbatch_res_val_norm$Batch2)
  pred_val_auc_norm = as.matrix(1/coef1_norm*pred_mat_norm_testing_val[,1] + 1/coef2_norm*pred_mat_norm_testing_val[,2])/(1/coef1_norm+1/coef2_norm)
  
  
  # LOSO
  LOSO = as.data.frame(matrix(NA, ncol = N_batch*2, nrow = nrow(test_expr_norm)))
  colnames(LOSO) = c("auc_ds1","auc_ds2", "auc_ds1_norm","auc_ds2_norm")
  for (k in 1:nrow(LOSO)){
    LOSO[k,1] = 1/RMSE(testing$status[-k], pred_sgbatch_res$Batch1[-k])
    LOSO[k,2] = 1/RMSE(testing$status[-k], pred_sgbatch_res$Batch2[-k])
    LOSO[k,3] = 1/RMSE(testing$status[-k], pred_sgbatch_res_norm$Batch1[-k])
    LOSO[k,4] = 1/RMSE(testing$status[-k], pred_sgbatch_res_norm$Batch2[-k])
  }
  pred_LOSO = (LOSO[,1]*pred_prob1 + LOSO[,2]*pred_prob2)/ (LOSO[,1]+LOSO[,2])
  pred_LOSO_norm = (LOSO[,3]*pred_prob1_norm + LOSO[,4]*pred_prob2_norm)/ (LOSO[,3]+LOSO[,4])
  
  
  ####  Evaluate performance
  tst_scores <- c(list(NoBatch=pred_base_res, Batch=pred_batch_res), 
                  pred_test_lst, pred_test_lst_norm,
                  list(Merged_norm=pred_norm_res, Avg=pred_avg, n_Avg=pred_N_avg, CS_Avg=pred_cs_avg, Reg_a=pred_reg_a, Reg_s=pred_reg_s, val_auc = pred_val_auc, LOSO_auc = pred_LOSO,
                       Avg_norm=pred_avg_norm, n_Avg_norm=pred_N_avg_norm, CS_Avg_norm=pred_cs_avg_norm, Reg_a_norm=pred_reg_a_norm, Reg_s_norm=pred_reg_s_norm, val_auc_norm = pred_val_auc_norm, LOSO_auc_norm = pred_LOSO_norm))
  
  perf_df[i,] <- sapply(tst_scores, function(preds){
    testing_auc <- RMSE(testing$status, as.vector(preds))
    return(testing_auc[1])})
  perf_df[i,13] <- RMSE(testing_val$status, as.vector(tst_scores$val_auc))
  perf_df[i,20] <- RMSE(testing_val$status, as.vector(tst_scores$val_auc_norm))
}

# Get mean and sd
perf_df[101,] = colMeans(perf_df[1:100,])
perf_df[102,] = apply(perf_df[1:100,], 2, sd)
perf_df[103,] = apply(perf_df[1:100,], 2, median)
perf_df[104,] = apply(perf_df[1:100,], 2, IQR)
  
  ####  Output results
file_name <- sprintf('scenario2_%s_rmse_%s_%s_%s_logRelAbun_m%s_v%s.csv', 
                     gsub('.', '', method, fixed=T),
                     gsub('.', '', phenotype, fixed=T),
                     gsub('.', '', norm_method, fixed=T),
                     gsub('.', '', sub(".*_count_*(.*?)_ctr.*", "\\1", command_args[3]), fixed=T),
                     gsub('.', '', max_batch_mean, fixed=T),
                     gsub('.', '', max_batch_var, fixed=T))
write.csv(perf_df, paste0(command_args[12],file_name))

#end_time <- Sys.time()
#print(end_time - start_time)