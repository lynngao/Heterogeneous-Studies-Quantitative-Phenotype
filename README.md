# Mitigating Heterogeneity Effects in Microbiome-based Quantitative Phenotype Prediction: A Comprehensive Workflow for Integrating Multiple Studies with Batch Normalization

This repository contains R scripts for running the analysis described in our manuscript. 
The R scripts can be divided into 2 basic components:

#### 1. Simulation studies

- [simulation_quantitative_scenario1.R](https://github.com/lynngao/Heterogeneous-Studies-Quantitative-Phenotype/blob/main/simulation_quantitative_scenario1.R): Simulations for Scenario 1 (Different background distributions of genomic features in populations).<br/>

Use the following command to run the script:<br/>

***Rscript simulation_quantitative_scenario1.R alpha sample_size num_related_gene ml_model phenotype normalization_method training_background_data1 training_background_data2 test_background_data result_directory***<br/>

**alpha:** Population difference factor between training and test data in the range of 0 and 1. 0 means no population difference while 1 means largest difference.<br/>
**sample_size:** Number of samples you want to simulate. We fixed the sample size of training and test datasets to be the same.<br/>
**num_related_gene:** Number of phenotype related genes. We chose 10 genes and the script will automatically use the first 10 genes in the data as phenotype related.<br/>
**ml_model:** Choose between 'rfr' or 'lasso' for random forests regression or logistic regression with L1 regularization. In our analysis we used 'rfr'. You can also implement your own ml models in helper.R file.<br/>
**phenotype:** Choose between 'linear', 'inverse', quadratic' and 'log'.<br/>
**normalization_method:** Choose between 'combat' or 'conqur'.<br/>
**training_background_data1:** The directory of first training background data (count).<br/>
**training_background_data2:** The directory of second training background data (count).<br/>
**test_background_data:** The directory of test background data (count).<br/>
**result_directory:** The directory of the RMSE table to be saved.<br/>


- [simulation_quantitative_scenario2.R](https://github.com/lynngao/Heterogeneous-Studies-Quantitative-Phenotype/blob/main/simulation_quantitative_scenario2.R): Simulations for Scenario 2 (Different batch effects in studies with the same background distribution of genomic features in a population).<br/>

Use the following command to run the script:<br/>

***Rscript simulation_quantitative_scenario2.R test_sample_size num_related_gene background_data num_batches batch1_size batch2_size mean_change variance_change ml_model phenotype normalization_method result_directory***<br/>

**test_sample_size:** Number of test samples you want to simulate.<br/>
**num_related_gene:** Number of phenotype related genes. We chose 10 genes and the script will automatically use the first 10 genes in the data as phenotype related.<br/>
**background_data:** The directory of background data (count).<br/>
**num_batches:** 2 batches.<br/>
**batch1_size:** Number of samples in batch 1.<br/>
**batch2_size:** Number of samples in batch 2.<br/>
**mean_change:** Severity levels for the effect on the mean. In our study we chose among 0, 3 and 5.<br/>
**variance_change:** Severity levels for the effect on the variance. In our study we chose among 0, 1 and 2.<br/>
**ml_model:** Choose between 'rfr' or 'lasso' for random forests regression or logistic regression with L1 regularization. In our analysis we used 'rfr'. You can also implement your own ml models in helper.R file.<br/>
**phenotype:** Choose between 'linear', 'inverse', quadratic' and 'log'.<br/>
**normalization_method:** Choose between 'combat' or 'conqur'.<br/>
**result_directory:** The directory of the RMSE table to be saved.<br/>


- [simulation_quantitative_scenario3.R](https://github.com/lynngao/Heterogeneous-Studies-Quantitative-Phenotype/blob/main/simulation_quantitative_scenario3.R): Simulations for Scenario 3 (Different phenotype models in different studies).<br/>

Use the following command to run the script:<br/>

***Rscript simulation_quantitative_scenario3.R sample_size num_related_gene number_overlap_gene ml_model phenotype normalization_method training_background_data1 training_background_data2 test_background_data result_directory***<br/>

**sample_size:** Number of samples you want to simulate. We fixed the sample size of training and test datasets to be the same.<br/>
**num_related_gene:** Number of phenotype related genes. We chose 10 genes and the script will automatically use the first 10 genes in the data as phenotype related.<br/>
**num_overlap_gene:** Number of overlapping phenotype related genes in training and test model. We chose between 2 to 10 genes. 10 genes means the training and test have same phenotype models.<br/>
**ml_model:** Choose between 'rfr' or 'lasso' for random forests regression or logistic regression with L1 regularization. In our analysis we used 'rfr'. You can also implement your own ml models in helper.R file.<br/>
**phenotype:** Choose between 'linear', 'inverse', quadratic' and 'log'.<br/>
**normalization_method:** Choose between 'combat' or 'conqur'.<br/>
**training_background_data1:** The directory of first training background data (count).<br/>
**training_background_data2:** The directory of second training background data (count).<br/>
**test_background_data:** The directory of test background data (count).<br/>
**result_directory:** The directory of the AUC table to be saved.<br/>

- [helper.R](https://github.com/lynngao/Heterogeneous-Studies-Quantitative-Phenotype/blob/main/helper.R): helper file contains functions for both simulations and real data applications.

#### 2. Real data applications
- [real_data_T2D_BMI.R](https://github.com/lynngao/Heterogeneous-Studies-Quantitative-Phenotype/blob/main/real_data_T2D_BMI.R): Applications on body mass index (BMI) prediction from Type-2 diabetes (T2D) metagenomic datasets.<br/>
- [real_data_OC_survival.R](https://github.com/lynngao/Heterogeneous-Studies-Quantitative-Phenotype/blob/main/real_data_OC_survival.R): Applications on survival prediction from ovarian cancer gene expression datasets.<br/>
