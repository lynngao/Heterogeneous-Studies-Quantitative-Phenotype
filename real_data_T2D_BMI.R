library(ggplot2)
library(ggfortify)
library(ggpubr)
library(vegan)

source("helper.R")

# remember to change the directories
metadata_qin = read_csv("metadata_qin.csv")
metadata_karlsson = read_csv("metadata_karlsson.csv")
metadata_li = read_csv("metadata_li.csv")
metadata_san = read_csv("metadata_san.csv")
qin = read_csv("count_qin.csv")
karlsson = read_csv("count_karlsson.csv")
li = read_csv("count_li.csv")
san = read_csv("count_san.csv")
rownames(qin) = qin$...1
qin$...1 = NULL
rownames(karlsson) = karlsson$...1
karlsson$...1= NULL
rownames(li) = li$...1
li$...1 = NULL
rownames(san) = san$...1
san$...1 = NULL

#features = unique(c(colnames(qin), colnames(karlsson), colnames(li), colnames(san)))
features = intersect(intersect(intersect(colnames(qin), colnames(karlsson)), colnames(li)), colnames(san))

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

qin_filtered = trim_count_data(features, qin)
karlsson_filtered = trim_count_data(features, karlsson)
li_filtered = trim_count_data(features, li)
san_filtered = trim_count_data(features, san)

qin_filtered$status = metadata_qin$BMI
qin_filtered$condition = metadata_qin$study_condition
qin_filtered = qin_filtered[!is.na(qin_filtered$status),]
qin_filtered = qin_filtered[!is.na(qin_filtered$condition),]
covar_qin = qin_filtered$condition

karlsson_filtered$status = metadata_karlsson$BMI
karlsson_filtered$condition = metadata_karlsson$study_condition
karlsson_filtered = karlsson_filtered[!is.na(karlsson_filtered$status),]
karlsson_filtered = karlsson_filtered[!is.na(karlsson_filtered$condition),]
covar_karlsson = karlsson_filtered$condition

li_filtered$status = metadata_li$BMI
li_filtered$condition = metadata_li$study_condition
li_filtered = li_filtered[!is.na(li_filtered$status),]
li_filtered = li_filtered[!is.na(li_filtered$condition),]
li_filtered =  subset(li_filtered, metadata_li$study_condition == "T2D" | metadata_li$study_condition == "control")
covar_li = li_filtered$condition

san_filtered$status = metadata_san$BMI
san_filtered$condition = metadata_san$study_condition
san_filtered = san_filtered[!is.na(san_filtered$status),]
san_filtered = san_filtered[!is.na(san_filtered$condition),]
covar_san = san_filtered$condition
###############################################################################

#command_args <- commandArgs(trailingOnly=TRUE)
command_args <- c("rfr","combat","san")
method = as.character(command_args[1])
norm_method = as.character(command_args[2])
test_ds = as.character(command_args[3])


#change to relative abundance
log_relative_abundance <- function(ds){
  ds = ds/rowSums(ds)
  min = min(apply(ds[,1:ncol(ds)], 1, function(x) min(x[x>0])))
  ds[ds == 0] = min*0.65
  ds = log(ds)
  #ds = ds / rowSums(ds)
  return(ds)
}

assign_ds <- function(ds_orig, covar_orig){
  ds_count = ds_orig[,1:length(features)]
  ds = log_relative_abundance(ds_count)
  ds$status = ds_orig$status
  covar = covar_orig
  return(list(ds_count, ds, covar))
}

if (test_ds == "qin"){
  temp1 = assign_ds(karlsson_filtered, covar_karlsson)
  temp2 = assign_ds(li_filtered, covar_li)
  temp3 = assign_ds(san_filtered, covar_san)
  temp4 = assign_ds(qin_filtered, covar_qin)
}
if (test_ds == "karlsson"){
  temp1 = assign_ds(qin_filtered, covar_qin)
  temp2 = assign_ds(li_filtered, covar_li)
  temp3 = assign_ds(san_filtered, covar_san)
  temp4 = assign_ds(karlsson_filtered, covar_karlsson)
}
if (test_ds == "li"){
  temp1 = assign_ds(qin_filtered, covar_qin)
  temp2 = assign_ds(karlsson_filtered, covar_karlsson)
  temp3 = assign_ds(san_filtered, covar_san)
  temp4 = assign_ds(li_filtered, covar_li)
}

if (test_ds == "san"){
  temp1 = assign_ds(qin_filtered, covar_qin)
  temp2 = assign_ds(karlsson_filtered, covar_karlsson)
  temp3 = assign_ds(li_filtered, covar_li)
  temp4 = assign_ds(san_filtered, covar_san)
}


for (i in 1:4){
  nam_count <- paste("ds_count", i, sep = "")
  nam <- paste("ds", i, sep = "")
  nam_cov <- paste("covar", i, sep = "")
  assign(nam_count, get(paste("temp", i, sep = ""))[[1]])
  assign(nam, get(paste("temp", i, sep = ""))[[2]])
  assign(nam_cov, get(paste("temp", i, sep = ""))[[3]])
}


##ComBat normalization
combat_covar <- function(ref_ds, ds, ref_covar, covar){
  mod = as.data.frame(factor(c(ref_covar, covar)))
  mod = model.matrix(~.,data = mod)
  ds_norm = combat_norm(ref_ds, ds, mod=mod)
  return(ds_norm)
}

##normalization process
if (norm_method == "combat"){
    ds1_norm = combat_covar(ds4, ds1, covar4, covar1)
    ds2_norm = combat_covar(ds4, ds2, covar4, covar2)
    ds3_norm = combat_covar(ds4, ds3, covar4, covar3)
}
  
if (norm_method == "conqur"){
  batchid = factor(c(rep(1, nrow(ds1)), rep(2, nrow(ds2)), rep(3, nrow(ds3)), rep(0, nrow(ds4))))
  covar = factor(c(covar1, covar2, covar3, covar4))
  taxa = rbind(ds_count1, ds_count2, ds_count3, ds_count4)
  #taxa = taxa[,1:length(features)]
  taxa_corrected = ConQuR(tax_tab=taxa, batchid=batchid, covariates=covar, batch_ref="0")
  taxa_corrected = as.data.frame(taxa_corrected)
  taxa_corrected$batchid = batchid
  count_data1_norm <- taxa_corrected[taxa_corrected$batchid == 1,!colnames(taxa_corrected)%in%c("batchid")] #training 1
  count_data2_norm <- taxa_corrected[taxa_corrected$batchid == 2,!colnames(taxa_corrected)%in%c("batchid")] #training 2
  count_data3_norm <- taxa_corrected[taxa_corrected$batchid == 3,!colnames(taxa_corrected)%in%c("batchid")] #training 3
  # count_data1_norm = conqur_norm(karlsson_filtered, qin_filtered)
  # count_data2_norm = conqur_norm(karlsson_filtered, li_filtered)
  ds1_norm = log_relative_abundance(count_data1_norm)
  ds1_norm$status = ds1$status
  ds2_norm = log_relative_abundance(count_data2_norm)
  ds2_norm$status = ds2$status
  ds3_norm = log_relative_abundance(count_data3_norm)
  ds3_norm$status = ds3$status
} 

## Merge two training sets into one
merged_training = rbind(ds1, ds2, ds3)
merged_training_norm = rbind(ds1_norm, ds2_norm, ds3_norm)

#############PCA plots######################
##unnormalized
data_unnormalized = rbind(ds1[,-ncol(ds1)], ds2[,-ncol(ds2)], ds3[,-ncol(ds3)], ds4[,-ncol(ds4)])
data_normalized = rbind(ds1_norm[,-ncol(ds1_norm)], ds2_norm[,-ncol(ds2_norm)], ds3_norm[,-ncol(ds3_norm)], ds4[,-ncol(ds4)])
#data_normalized = rbind(count_data1_norm, count_data2_norm, count_data3_norm, ds_count4)
#data_normalized = data_normalized[ , which(apply(data_normalized, 2, var) != 0)]
# data_unnormalized = data_unnormalized / rowSums(data_unnormalized)
# data_normalized = data_normalized / rowSums(data_normalized)
pca_res <- prcomp(data_unnormalized, center=TRUE, scale.=TRUE)
pca_res1 <- prcomp(data_normalized, center=TRUE, scale.=TRUE)

if (test_ds == "qin"){
  dataset = c(rep("Karlsson",nrow(ds1)),rep("Li",nrow(ds2)), rep("Sankaranarayanan",nrow(ds3)), rep("Qin",nrow(ds4)))
}
if (test_ds == "karlsson"){
  dataset = c(rep("Qin",nrow(ds1)),rep("Li",nrow(ds2)), rep("Sankaranarayanan",nrow(ds3)), rep("Karlsson",nrow(ds4)))
}
if (test_ds == "li"){
  dataset = c(rep("Qin",nrow(ds1)),rep("Karlsson",nrow(ds2)), rep("Sankaranarayanan",nrow(ds3)), rep("Li",nrow(ds4)))
}
if (test_ds == "san"){
  dataset = c(rep("Qin",nrow(ds1)), rep("Karlsson",nrow(ds2)), rep("Li",nrow(ds3)), rep("Sankaranarayanan",nrow(ds4)))
}

data_unnormalized$dataset = dataset
data_normalized$dataset = dataset

a  = autoplot(pca_res, data = data_unnormalized , colour = "dataset")  + 
  stat_ellipse(aes(colour = dataset)) +
  ggtitle("Unnormalized") +
  theme_classic() +
  theme(legend.title = element_blank())
#ggsave('plots_quantitative/PCA_T2D.jpg', width = 15, height = 15, units = "cm", dpi=600, bg="white")

#combat
b  = autoplot(pca_res1, data = data_normalized , colour = "dataset")  + 
  stat_ellipse(aes(colour = dataset)) +
  ggtitle("ComBat normalized") +
  theme_classic() +
  theme(legend.title = element_blank())

#conqur
c  = autoplot(pca_res1, data = data_normalized , colour = "dataset")  + 
  stat_ellipse(aes(colour = dataset)) +
  ggtitle("ConQuR normalized") + xlim(c(-0.05, 0.1)) + ylim(c(-0.15, 0.15)) +
  theme_classic() +
  theme(legend.title = element_blank())

plot1 = ggarrange(a, b, c, labels = c("A","B", "C"),
                  nrow = 1, ncol = 3, common.legend = TRUE, hjust = -0.35)
#annotate_figure(plot1, top = text_grob("PCA: T2D datasets, No covariate in normalization", face = "bold", size = 17, color = "red"))
#ggsave('/Users/lynngao/Desktop/final_results/PCA_T2D_no_covar.jpg', width = 50, height = 25, units = "cm", dpi=600, bg="white")
annotate_figure(plot1, top = text_grob("Test dataset: Sankaranarayanan", face = "bold", size = 17, color = "red"))
ggsave('/Users/lynngao/Desktop/final_results/PCA_T2D_covar_test_Sankaranarayanan.jpg', width = 50, height = 25, units = "cm", dpi=600, bg="white")
###########################################
# Set a dataframe to save results
perf_df <- as.data.frame(matrix(NA, nrow=34, ncol=22))
colnames(perf_df) = c("Training1","Training2", "Training3", "Training1_norm","Training2_norm","Training3_norm", "Merged", "Merged_norm",
                      "Avg", "n_Avg", "CS_Avg", "Reg_a", "Reg_s", "val_rmse", "LOSO_rmse", 
                      "Avg_norm", "n_Avg_norm", "CS_Avg_norm", "Reg_a_norm", "Reg_s_norm", "val_rmse_norm", "LOSO_rmse_norm")
rownames(perf_df) = c(1:30, "mean", "sd", "median", "IQR")

for (i in 1:30) {
  ## Split testing set into validation and testing
  testing = ds4
  sample = sample.int(n=nrow(testing),size=floor(0.5*nrow(testing)),replace=F)
  val = testing[sample,]
  testing_val = testing[rownames(testing)%!in%rownames(val),]
  
  ## Obtain predictions from learner trained within each training set
  ml1 = ml_model(ds1, method)
  ml2 = ml_model(ds2, method)
  ml3 = ml_model(ds3, method)
  ml1_norm = ml_model(ds1_norm, method)
  ml2_norm = ml_model(ds2_norm, method)
  ml3_norm = ml_model(ds3_norm, method)
  pred_prob1 = predict(ml1, testing)
  pred_prob2 = predict(ml2, testing)
  pred_prob3 = predict(ml3, testing)
  pred_sgbatch_res <- list(pred_prob1, pred_prob2, pred_prob3)
  pred_prob1_norm = predict(ml1_norm, testing)
  pred_prob2_norm = predict(ml2_norm, testing)
  pred_prob3_norm = predict(ml3_norm, testing)
  pred_sgbatch_res_norm  <- list(pred_prob1_norm, pred_prob2_norm, pred_prob3_norm)
  
  ## Prediction from combine training together (Merged)
  ml_merged = ml_model(merged_training, method)
  pred_merged_res = predict(ml_merged, testing)
  
  ## Prediction from training after batch adjustment (Merged combat)
  ml_merged_norm = ml_model(merged_training_norm, method)
  pred_norm_res = predict(ml_merged_norm, testing)
  
  ## Prediction from training on testing_val
  pred_sgbatch_res_testing_val <- list(predict(ml1, testing_val), predict(ml2, testing_val), predict(ml3, testing_val))
  names(pred_sgbatch_res_testing_val) <- paste0("ds", 1:3)
  pred_sgbatch_res_norm_testing_val  <- list(predict(ml1_norm, testing_val), predict(ml2_norm, testing_val), predict(ml3_norm, testing_val))
  names(pred_sgbatch_res_norm_testing_val) <- paste0("ds", 1:3)
  
  ## Prediction on validation dataset
  pred_sgbatch_res_val <- list(predict(ml1, val), predict(ml2, val), predict(ml3, val))
  names(pred_sgbatch_res_val) <- paste0("ds", 1:3)
  pred_sgbatch_res_val_norm <- list(predict(ml1_norm, val), predict(ml2_norm, val), predict(ml3_norm, val))
  names(pred_sgbatch_res_val_norm) <- paste0("ds", 1:3)
  
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
  pred_N_avg <- pred_mat %*% (as.matrix(rbind(nrow(ds1), nrow(ds2), nrow(ds3))) / (nrow(ds1) + nrow(ds2) + nrow(ds3)))
  pred_N_avg_norm <- pred_mat_norm %*% (as.matrix(rbind(nrow(ds1_norm), nrow(ds2_norm), nrow(ds3_norm))) / (nrow(ds1_norm) + nrow(ds2_norm) + nrow(ds3_norm)))
  
  # CS-Avg: replicability weights
  cs_zmat <- CS_zmatrix(n_batch=3, training = list(ds1, ds2, ds3), perf_name="rmse", method)
  cs_weights_seq <- CS_weight(cs_zmat)
  pred_cs_avg <- pred_mat %*% cs_weights_seq
  
  cs_zmat_norm <- CS_zmatrix(n_batch=3, training = list(ds1_norm, ds2_norm, ds3_norm), perf_name="rmse", method)
  cs_weights_seq_norm <- CS_weight(cs_zmat_norm)
  pred_cs_avg_norm <- pred_mat_norm %*% cs_weights_seq_norm
  
  
  # Reg-a: use each function to predict on one study, bind predictions and do regression
  reg_ssl_res <- Reg_SSL_pred(n_batch=3, training = list(ds1, ds2, ds3), method)
  reg_a_beta <- Reg_a_weight(coef_mat=do.call(rbind, reg_ssl_res$coef), 
                             n_seq=sapply(list(ds1, ds2, ds3), nrow))
  pred_reg_a <- pred_mat %*% reg_a_beta
  
  reg_ssl_res_norm <- Reg_SSL_pred(n_batch=3, training = list(ds1_norm, ds2_norm, ds3_norm), method)
  reg_a_beta_norm <- Reg_a_weight(coef_mat=do.call(rbind, reg_ssl_res_norm$coef), 
                                  n_seq=sapply(list(ds1_norm, ds2_norm, ds3_norm), nrow))
  pred_reg_a_norm <- pred_mat_norm %*% reg_a_beta_norm
  
  
  # Reg-s:
  stacked_pred <- do.call(rbind, reg_ssl_res$pred)
  stacked_label = c(ds1$status, ds2$status, ds3$status)
  reg_s_beta <- pnnls(a=stacked_pred, b=as.numeric(stacked_label), sum=1)$x
  pred_reg_s <- pred_mat %*% reg_s_beta
  
  stacked_pred_norm <- do.call(rbind, reg_ssl_res_norm$pred)
  reg_s_beta_norm <- pnnls(a=stacked_pred_norm, b=as.numeric(stacked_label), sum=1)$x
  pred_reg_s_norm <- pred_mat_norm %*% reg_s_beta_norm
  
  
  # Weighting by validation set performance
  coef1 = RMSE(val$status, pred_sgbatch_res_val$ds1)
  coef2 = RMSE(val$status, pred_sgbatch_res_val$ds2)
  coef3 = RMSE(val$status, pred_sgbatch_res_val$ds3)
  pred_val_auc = as.matrix(1/coef1*pred_mat_testing_val[,1] + 1/coef2*pred_mat_testing_val[,2] + 1/coef3*pred_mat_testing_val[,3])/(1/coef1+1/coef2+1/coef3)
  
  coef1_norm = RMSE(val$status, pred_sgbatch_res_val_norm$ds1)
  coef2_norm = RMSE(val$status, pred_sgbatch_res_val_norm$ds2)
  coef3_norm = RMSE(val$status, pred_sgbatch_res_val_norm$ds3)
  pred_val_auc_norm = as.matrix(1/coef1_norm*pred_mat_norm_testing_val[,1] + 1/coef2_norm*pred_mat_norm_testing_val[,2] + 1/coef3_norm*pred_mat_norm_testing_val[,3])/(1/coef1_norm+1/coef2_norm+1/coef3_norm)
  
  
  # LOSO
  LOSO = as.data.frame(matrix(NA, ncol = 6, nrow = nrow(testing)))
  for (k in 1:nrow(LOSO)){
    LOSO[k,1] = 1/RMSE(testing$status[-k], pred_prob1[-k])
    LOSO[k,2] = 1/RMSE(testing$status[-k], pred_prob2[-k])
    LOSO[k,3] = 1/RMSE(testing$status[-k], pred_prob3[-k])
    LOSO[k,4] = 1/RMSE(testing$status[-k], pred_prob1_norm[-k])
    LOSO[k,5] = 1/RMSE(testing$status[-k], pred_prob2_norm[-k])
    LOSO[k,6] = 1/RMSE(testing$status[-k], pred_prob3_norm[-k])
  }
  pred_LOSO = (LOSO[,1]*pred_prob1 + LOSO[,2]*pred_prob2 + LOSO[,3]*pred_prob3)/ (LOSO[,1]+LOSO[,2]+LOSO[,3])
  pred_LOSO_norm = (LOSO[,4]*pred_prob1_norm + LOSO[,5]*pred_prob2_norm + + LOSO[,6]*pred_prob3_norm)/ (LOSO[,4]+LOSO[,5]+LOSO[,6])
  
  
  ####  Evaluate performance
  tst_scores <- c(list(Training1 = pred_prob1, Training2 = pred_prob2, Training3 = pred_prob3, 
                       Training1_norm = pred_prob1_norm, Training2_norm = pred_prob2_norm, Training3_norm = pred_prob3_norm,
                       Merged=pred_merged_res, Merged_norm=pred_norm_res, 
                       Avg=pred_avg, n_Avg=pred_N_avg, CS_Avg=pred_cs_avg, Reg_a=pred_reg_a, Reg_s=pred_reg_s, val_auc = pred_val_auc, LOSO_auc = pred_LOSO, 
                       Avg_norm=pred_avg_norm, n_Avg_norm=pred_N_avg_norm, CS_Avg_norm=pred_cs_avg_norm, Reg_a_norm=pred_reg_a_norm, Reg_s_norm=pred_reg_s_norm, val_auc_norm = pred_val_auc_norm, LOSO_auc_norm = pred_LOSO_norm)) 
  
  perf_df[i,] <- sapply(tst_scores, function(preds){
    testing_auc <- RMSE(testing$status, as.vector(preds))
    return(testing_auc[1])})
  perf_df[i,14] <- RMSE(testing_val$status, as.vector(tst_scores$val_auc))
  perf_df[i,21] <- RMSE(testing_val$status, as.vector(tst_scores$val_auc_norm))
}

# Get mean and sd
perf_df[31,] = colMeans(perf_df[1:30,])
perf_df[32,] = apply(perf_df[1:30,], 2, sd)
perf_df[33,] = apply(perf_df[1:30,], 2, median)
perf_df[34,] = apply(perf_df[1:30,], 2, IQR)

####  Output results
file_name <- sprintf('quantitative_T2D_data_bmi_rmse_%s_%s_test_%s_logRelAbun_covar.csv', 
                     gsub('.', '', method, fixed=T),
                     gsub('.', '', norm_method, fixed=T),
                     gsub('.', '', test_ds, fixed=T))
write.csv(perf_df, paste0(command_args[3],file_name))
