# Advanced-Bioinformatics

Improving CRISPR KO Screening

## Data Preprocessing

Following the quality control, the data points from GeCKO library B were extracted resulting in 58,028 sgRNAs. The data points corresponding to non-targeting control sgRNAs were removed along with the three sgRNAs with mismatched gene identifiers between the data sheets. The resulting data set contains 57,025 sgRNAs.

Initially only the sgRNA counts from the BC3 cell line were used. The log2 Fold Change between day 0 and day 14 was computed using the DESeq2 library. These values were then used as the dependent variable to estimate sgRNA efficacy.

To model the sgRNA efficacy, the guide sequence composition was used. The overall GC-content was computed as the ratio of Guanine (G) and Cytosine (C) in the guide. Position-wise nucleotide composition was also used as input for the model.

This data was subseqeuntly divided into training data (80%) and testing data (20%).

The data preprocessing procedure can be found in: "R/03_data_preprocessing.R"

## Model

XGBoost (Extreme Gradient Boosting) was used to model log2FoldChange of sgRNA counts by sequence composition. XGBoost is a scalable gradient boosting framework. The model builds sequential decision trees, with each tree attempting to improve the predictions of the previous trees.

The response variable was the log2 fold change in sgRNA abundance and the predictor variables the sequence-composition features generated during pre-processing of data.

The training data was used to fit the model and for cross-validation during hyperparameter optimizations, whereas the testing data was only used for the evaluation of the final model.

## Hyperparameter Optimization

Hyperparameters eta, gamma, max_depth, subsample, colsample_bytree, and min_child_weight were optimized. Hyperparameter optimization was performed using 5-fold cross-validation on the training data. The number of boosting rounds was also determined during cross-validation using early stopping.

For each parameter configuration up to 1,000 boosting rounds were evaluated, stopping training when RMSE failed to improve for 30 consecutive rounds. The parameters were optimized sequentially in the previously stated order, retaining the selected parameter value from the previous optimization.

The cross-validation Root Mean Squared Error (CV RMSE) was used to select parameter values, where a lower value indicates a lower prediction error.

### eta

eta or learning rate sets the step size shrinkage used in update to prevent overfitting. The default value for eta for a Tree Booster model is 0.3 and the parameter has a range of 0 to 1. For the optimization three values larger than 0.3 (0.4, 0.5 and 0.6) and four values smaller than 0.3 (0.15, 0.1, 0.05 and 0.01) were investigated.

Reported for each tested value is the cross-validation RMSE and Mean Absolute Error (MAE).

|         |            |             |           |
|---------|------------|-------------|-----------|
| **eta** | **rounds** | **cvRMSE**  | **cvMAE** |
| 0.01    | 403        | 0.4782638   | 0.3282582 |
| 0.05    | 102        | 0.4781006   | 0.3284748 |
| 0.10    | 59         | 0.4782608   | 0.3289320 |
| 0.15    | 36         | 0.4786662   | 0.3288769 |
| 0.30    | 12         | 0.4798945   | 0.3301756 |
| 0.40    | 7          | 0.4809617   | 0.3306126 |
| 0.50    | 5          | 0.4814286   | 0.3308764 |
| 0.60    | 4          | 0.4815781   | 0.3310022 |

Out of the tested values, eta = 0.05 produced the lowest CV RMSE and was used for subsequent optimizations.

### gamma

gamma or min_split_loss refers to the minimum loss reduction required to make a split. A larger value makes the model more conservative. The default value for a Tree Booster is 0 and the range of the parameter is 0 to infinity. Gamma values 0, 0.05, 0.1, 0.15, 0.2 and 0.3 were investigated.

|           |            |            |           |
|-----------|------------|------------|-----------|
| **gamma** | **rounds** | **cvRMSE** | **cvMAE** |
| 0.00      | 102        | 0.4781006  | 0.3284748 |
| 0.05      | 103        | 0.4780736  | 0.3284346 |
| 0.10      | 102        | 0.4781087  | 0.3283958 |
| 0.15      | 102        | 0.4780777  | 0.3284907 |
| 0.20      | 102        | 0.4780952  | 0.3284844 |
| 0.30      | 98         | 0.4781303  | 0.3284101 |

Out of the tested values, gamma = 0.05, generated the lowest CV RMSE and was selected for subsequent optimizations.

### max_depth

The max_depth hyperparameter controls the maximum depth of a tree. Increasing this value causes the model to become more complex and more likely to overfit the data. The default value for max_depth for a Tree Booster is 6, with the range 0 to infinity. Investigated values of max_depth were 3, 4, 5, 6, 7, 8 and 9.

|               |            |            |           |
|---------------|------------|------------|-----------|
| **max_depth** | **rounds** | **cvRMSE** | **cvMAE** |
| 3             | 332        | 0.4782452  | 0.3282694 |
| 4             | 202        | 0.4781503  | 0.3282752 |
| 5             | 144        | 0.4780084  | 0.3282320 |
| 6             | 103        | 0.4780736  | 0.3284346 |
| 7             | 93         | 0.4780518  | 0.3290080 |
| 8             | 74         | 0.4783910  | 0.3292971 |
| 9             | 60         | 0.4797385  | 0.3306223 |

From the tested values, max_depth = 5 generated the lowest CV RMSE and was selected for subsequent optimizations.

### subsample

This parameter controls the subsample ratio of the training data. A value of 0.5 means that the model would randomly sample half the data prior to growing trees. The default value for a Tree Booster is subsample = 1 with a range of 0 to 1. Selected for investigation were values 1, 0.9, 0.8, 0.7, 0.6 and 0.5.

|               |            |            |           |
|---------------|------------|------------|-----------|
| **subsample** | **rounds** | **cvRMSE** | **cvMAE** |
| 1.0           | 144        | 0.4780084  | 0.3282320 |
| 0.9           | 152        | 0.4775343  | 0.3279694 |
| 0.8           | 139        | 0.4775003  | 0.3279005 |
| 0.7           | 160        | 0.4774903  | 0.3279962 |
| 0.6           | 170        | 0.4771984  | 0.3279112 |
| 0.5           | 114        | 0.4776754  | 0.3279129 |

From the tested values, subsample = 0.6 generated the lowest CV RMSE and was selected for subsequent optimizations.

### colsample_bytree

The colsample_bytree parameter is the subsample ratio of columns when constructing each tree. The default value for a Tree Booster is 1 with a range of 0 to 1. Selected for investigation were values 1, 0.9, 0.8, 0.7, 0.6 and 0.5.

|                      |            |            |           |
|----------------------|------------|------------|-----------|
| **colsample_bytree** | **rounds** | **cvRMSE** | **cvMAE** |
| 1.0                  | 170        | 0.4771984  | 0.3279112 |
| 0.9                  | 192        | 0.4773508  | 0.3279299 |
| 0.8                  | 132        | 0.4776795  | 0.3279114 |
| 0.7                  | 191        | 0.4773535  | 0.3278682 |
| 0.6                  | 192        | 0.4772528  | 0.3278777 |
| 0.5                  | 192        | 0.4770440  | 0.3276524 |

From the tested values, colsample_bytree = 0.5 generated the lowest CV RMSE and was selected for subsequent optimizations.

### min_child_weight

This parameter determines the minimum sum of instance weight needed in a child. For a Tree Booster the default value i 1 with a range of 0 to infinity, where 0 indicates no limit. Selected for investigation were values 0, 0.5, 1, 2, 5 and 10.

|                      |            |            |           |
|----------------------|------------|------------|-----------|
| **min_child_weight** | **rounds** | **cvRMSE** | **cvMAE** |
| 0.0                  | 192        | 0.4770440  | 0.3276524 |
| 0.5                  | 192        | 0.4770440  | 0.3276524 |
| 1.0                  | 192        | 0.4770440  | 0.3276524 |
| 2.0                  | 188        | 0.4771280  | 0.3276454 |
| 5.0                  | 192        | 0.4770498  | 0.3276026 |
| 10.0                 | 162        | 0.4773224  | 0.3276340 |

From the tested values, min_child_weight 0, 0.5 and 1 generated identical CV RMSE values.

## Selected Hyperparameters

The sequential optimization resulted in the following selected parameters:

| Parameter        | Selected Value |
|------------------|----------------|
| eta              | 0.05           |
| gamma            | 0.05           |
| max_depth        | 5              |
| subsample        | 0.6            |
| colsample_bytree | 0.5            |
| min_child_weight | 0              |

After selecting these parameter values, a final cross-validation identified 192 boosting rounds.

The optimization procedure can be found in: "R/04_model_optimization.R"

## Final Model Evaluation

The final XGBoost model was trained on the complete training dataset using the selected hyperparameters and 192 boosting rounds. The held-out test dataset was then used to obtain an independent evaluation of model performance.

|  |  |  |  |  |  |  |
|--------|----------|----------|----------|---------|--------------|--------------|
| **rounds** | **trainRMSE** | **testRMSE** | **trainMAE** | **testMAE** | **trainSpearman** | **testSpearman** |
| 192 | 0.4633858 | 0.483437 | 0.318533 | 0.330785 | 0.3373936 | 0.2719628 |

The test-set RMSE was 0.483 compared with the training-set RMSE of 0.463, and the corresponding MAE values were 0.331 and 0.319 for the test- and training-sets. Furthermore a Spearman correlation between observed and predicted log2 fold changes was 0.272 for the test-set and 0.337 for the training-set.

The model showed positive association between predicted and observed sgRNA efficacy, while prediction errors for the test data were higher than those for the training data.

The final model evaluation can be found in: "R/05_final_model.R"
