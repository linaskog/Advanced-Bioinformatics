# Modelling

## Data Preparation

### Response Variable

The Log2 Fold Change between the Day 0 and Day 14 sgRNA abundance counts for the BC-3 cell line was computed using the DESeq2 library. The gene-level average of Log2 Fold Change was subtracted from this measure to account for differences in sgRNA abundance between genes. The raw Log2 Fold Change was stored for downstream analysis, and the centered Log2 Fold Change as the measure of sgRNA efficacy.

### Predictor Variables

Using the guide sequence composition, position-specific nucleotide information and GC-content variables were extracted. From the annotation experiment, the strand, chromosome, region type and gene overlap counts were selected as predictor variables. The ChIP-Seq peaks from the genome mapping was selected as a final variable.

All variable columns were converted to numerical representations for modelling. The script for response and predictor variable preparation is `R/09_data_preparation.R` and the computed variables can be found in `data/processed/variables.rds` .

## Linear Model Training

A base R linear model was trained on a training set of the variables. The first model was trained using one-hot encoded nucleotides of the guides, GC-content, strand, chromosome, region type, gene overlap count and H3K27ac peak information. This model is presented in the summary below:

``` text
Residual standard error: 0.3255 on 41459 degrees of freedom Multiple R-squared:  0.04302,   Adjusted R-squared:  0.04152  F-statistic: 28.67 on 65 and 41459 DF,  p-value: < 2.2e-16
```

In lieu of computing new model variables, cross-terms between the one-hot encoded guide sequence were considered. Two additional models were trained including cross-terms between pairs and triplets of adjacent position variables. The summary of the model containing the original model variables and the cross-terms between adjacent pairs of position variables is shown below.

``` text
Residual standard error: 0.3236 on 41288 degrees of freedom Multiple R-squared:  0.05838,   Adjusted R-squared:  0.053  F-statistic: 10.85 on 236 and 41288 DF,  p-value: < 2.2e-16
```

This model showed a decrease in residual standard error and an increase in R-squared compared to the previous model. A third model was also trained, including the variables of the second model along with cross-terms of triplets of adjacent position variables. The summary of this model is presented below:

``` text
Residual standard error: 0.3222 on 40640 degrees of freedom Multiple R-squared:  0.08071,   Adjusted R-squared:  0.06071  F-statistic: 4.036 on 884 and 40640 DF,  p-value: < 2.2e-16
```

This model shows a decrease in residual squared error and an increase in R-squared compared to the other models.

### Model Evaluation

All three models were then evaluated on the held out test data set, computing the R-Squared, Root Mean Squared Error, Mean Absolute Error and Spearman Correlation between the model predicted and observed centered Log2 Fold Change. The result of this evaluation is shown in the table below.

|  |  |  |  |  |  |  |  |  |
|--------|--------|--------|--------|--------|--------|--------|--------|--------|
| **Model** | **R2Train** | **R2Test** | **RMSETrain** | **RMSETest** | **MAETrain** | **MAETest** | **SpearmanTrain** | **SpearmanTest** |
| Linear Model (No cross terms | 0.043022517 | 0.041210384 | 0.3252493 | 0.3278084 | 0.2375401 | 0.2383944 | 0.20165180 | 0.20123585 |
| LM Dinucleotide Cross-Terms | 0.058379492 | 0.046440428 | 0.3226290 | 0.3269131 | 0.2360659 | 0.2380683 | 0.23250990 | 0.21631314 |
| LM Trinucleotide Cross-Terms | 0.080706956 | 0.036782904 | 0.3187810 | 0.3285644 | 0.2342924 | 0.2402606 | 0.26666786 | 0.20433756 |

Both models with cross-terms showed a decrease in training error compared to the first model, with the third model achieving the lowest training error/highest explained variance. However, when evaluated using the held-out test set, the third model was shown to have less predictive accuracy than the other two methods. These results indicate that the increase in training R-Squared was caused due to over-fitting. The second model performed the best overall with the highest test R-Squared and correlation, and lowest test RMSE and MAE values. The script used to train and evaluate the linear models can be found in `R/10_linear_model.R` and the second model was saved `models/xgb_model.ubj`.

The fit of the models was visualized plotting the observations against the predictions made by each model for both training and testing data. The script for generating these plots can be found in `R/10.1_linear_model_fit_plots.R` plots can be found under "results/models/". These plots were shown to be quite similar for the three models.

### Model Variables

The regression coefficients of the models were analyzed in an attempt to simplify the models. Filtering out variables with Pr(\>\|t\|) \> 0.10 for each of the three linear models. The models were retrained and evaluated and it was shown that removing the less-significant variables increased both training and generalization errors of the models. The script for these analyses can be found in `R/10.2_lm_sign_coeff_train.R` and `R/10.3_lm_evaluate_retrained_models.R.`

The second model, including cross-terms between pairs of adjacent position-variables was retained modelling.

## Extreme Gradient Boosting Model

XGBoost (Extreme Gradient Boosting) was used to model gene-centered Log2 Fold Change in sgRNA abundance. The training data was used to fit the model and for cross-validation during hyperparameter optimizations, whereas the testing data was only used for the evaluation of the final model.

### Hyperparameter Selection

Values of hyperparameters eta, gamma, max_depth, subsample, colsample_bytree, and min_child_weight were selected using 5-fold cross-validation on the training data. A small set of parameter values was investigated for each parameter and value producing the lowest cross-validation RMSE (CV-RMSE) was selected. The number of boosting rounds was also determined during cross-validation using early stopping. For each parameter configuration up to 1,000 boosting rounds were evaluated, stopping training when RMSE failed to improve for 30 consecutive rounds. The parameters were selected sequentially in the previously stated order, retaining the selected parameter value from the previous optimization. The optimization procedure can be found in `R/12_xgb_model_optimization.R`

### Final Model Evaluation 

The final XGBoost model was trained on the complete training data set using the selected hyperparameters. The held-out test data set was then used to obtain an independent evaluation of model performance and the results are shown in the table below.

| rounds | trainRMSE | testRMSE | trainMAE | testMAE | trainR2 | testR2 | trainSpearman | testSpearman |
|--------|--------|--------|--------|--------|--------|--------|---------|--------|
|    251 | 0.3216164 | 0.3263604 | 0.234871 | 0.2370772 | 0.06428116 | 0.04966186 | 0.2423636 | 0.2198846 |

The script for the final model evaluation can be found in `R/13_xgb_model_evaluation.R`.

### SHAP 

SHaplet Additive exPlanation was used to investigate how sequence features contributed to the predictions of the final model. SHAP analysis was performed on the final model after hyperparameter optimization. A SHAP-plot was generated and saved in `results/models/shap_plot_xgb_model.pdf` and the script for producing the SHAP-plot can be found in `R/13_xgb_model_evaluation.R`

## Final Model Exploration

The final models were evaluated in script `R/14_final_model_exploration.R`. The XGBoost Model predictions were plotted against the Linear Model predictions in scatter plots, stored in `results/models`. The Spearman correlation between the two model predictions was also computed, along with the correlation between absolute model prediction and observed difference. The resulting correlations are shown in the table below.

|  |  |  |  |
|-------------|-------------|------------------|-------------------|
| **Correlation Test** | **Correlation Train** | **Correlation Residual Test** | **Correlation Residual Train** |
| 0.894 | 0.891 | 0.966 | 0.99 |

### Variables

The five most important variables for either model were extracted and are presented below. The script for extracting these variables can be found in `R/14.1_summarize_results.R`

|                 |             |
|-----------------|-------------|
| **Linear**      | **XGBoost** |
| pos19_G         | pos20_G     |
| pos20_G         | pos20_C     |
| pos19_G:pos20_C | pos17_T     |
| pos19_C         | pos18_C     |
| pos18_C         | pos19_T     |

Guanine at position 20 is an important feature in both models, as are several other position-specific features near the PAM sequence.

### Essential Gene Investigation

The Log2 Fold Change was averaged per gene, and bottom 25% LFC scores was used to select "essential" genes. The essential genes were stored for re-use in `data/processed/essential_genes.rds`. Additionally, gene level averages of centered LFC, the linear model predicted centered LFC and the residuals between these LFCs were computed. The distribution of scores was plotted in histograms, coloring columns by gene groups. It was shown that essential and non-essential genes had highly similar distributions of these scores. The distribution plots were saved to `results/models` and the script can be found in `R/10.4_linear_model_genes.R`. This procedure was repeated for the XGBoost model and the script can be found in `R/13.1_xgb_model_genes.R`, here too the distributions of scores were highly similar for the two groups of genes.

### Computing Classification Metrics

The predicted scores for the training and testing data set were compared to the actual gene centered Log2 Fold Change score. Guides were marked as having "high" activity ( a negative score) or "low" activity (a positive score). Confusion matrices were computed independently for the training and testing data set to be able to compare training and generalization error. Classification metrics were computed for the linear model predictions in the scripts `R/10.5_lm_classification.R` and `R/10.6_lm_class_cont.R`. This procedure was repeated for the XGBoost model and the scripts can be found in `R/13.2_lm_classification.R` and `R/13.3_lm_class_cont.R`. A summarisation of these analyses was computed in `R/14.1_summarize_results.R` and the following table generated.

|  |  |  |  |  |  |  |
|----|----|----|----|----|----|----|
| **Data** | **Accuracy** | **Precision** | **Sensitivity** | **Specificity** | **NPV** | **F1 Score** |
| LMTrain | 0.579 | 0.576 | 0.569 | 0.589 | 0.582 | 0.572 |
| LMTest | 0.575 | 0.573 | 0.563 | 0.586 | 0.577 | 0.568 |
| XGBoostTrain | 0.582 | 0.579 | 0.573 | 0.591 | 0.586 | 0.576 |
| XGBoostTest | 0.578 | 0.576 | 0.566 | 0.589 | 0.580 | 0.571 |

The classification metrics computed for the either model are quite similar, however quite modest as either model predictions have ~58% accuracy. The metrics computed for training and testing data are similar for both models and all metrics. 
