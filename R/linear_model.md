# Linear Model

## Model Training

A base R linear model was trained on a training set of the variables. Initially the centered Log2 Fold Change was modeled using only the GC-content of the guide sequence. This initial model did not perform well with an R-Squared of only 0.001837, but was kept as a baseline for evaluating the other models.

``` text
Residual standard error: 0.3322 on 41523 degrees of freedom
Multiple R-squared:  0.001837,  Adjusted R-squared:  0.001813 
F-statistic: 76.43 on 1 and 41523 DF,  p-value: < 2.2e-16
```

A new model was trained using one-hot encoded nucleotides of the guides, GC-content, strand, chromosome, region type, gene overlap count and H3K27ac peak information. This model is presented in the summary below:

``` text
Residual standard error: 0.3255 on 41459 degrees of freedom
Multiple R-squared:  0.04302,   Adjusted R-squared:  0.04152 
F-statistic: 28.67 on 65 and 41459 DF,  p-value: < 2.2e-16
```

This model shows a large increase in predictive ability compared to the baseline (x23 increase in R-Squared). In lieu of computing new model variables, cross-terms between the one-hot encoded guide sequence were considered. Two additional models were trained including cross-terms between pairs and triplets of adjacent position variables.

The summary of the model containing the original model variables and the cross-terms between adjacent pairs of position variables.

``` text
Residual standard error: 0.3236 on 41288 degrees of freedom
Multiple R-squared:  0.05838,   Adjusted R-squared:  0.053 
F-statistic: 10.85 on 236 and 41288 DF,  p-value: < 2.2e-16
```

This model showed a decrease in residual standard error and a 32x increase in R-squared compared to the baseline. A third model was also trained, including the variables of the second model along with cross-terms of triplets of adjacent position variables. The summary of this model is presented below:

``` text
Residual standard error: 0.3222 on 40640 degrees of freedom
Multiple R-squared:  0.08071,   Adjusted R-squared:  0.06071 
F-statistic: 4.036 on 884 and 40640 DF,  p-value: < 2.2e-16
```

This model shows a decrease in residual squared error and a 44x increase in R-squared compared to the baseline.

## Model Evaluation

All three models were then evaluated on the held out test data set, computing the R-Squared, Root Mean Squared Error, Mean Absolute Error and Spearman Correlation between the model predicted and observed centered Log2 Fold Change. The result of this evaluation is shown in the table below.

|  |  |  |  |  |  |  |  |  |
|----|----|----|----|----|----|----|----|----|
| **Model** | **R2Train** | **R2Test** | **RMSETrain** | **RMSETest** | **MAETrain** | **MAETest** | **SpearmanTrain** | **SpearmanTest** |
| Baseline (only GC) | 0.001837295 | 0.001798662 | 0.3321744 | 0.3344780 | 0.2411301 | 0.2417564 | 0.04075234 | 0.04586712 |
| Linear Model (No cross terms | 0.043022517 | 0.041210384 | 0.3252493 | 0.3278084 | 0.2375401 | 0.2383944 | 0.20165180 | 0.20123585 |
| LM Dinucleotide Cross-Terms | 0.058379492 | 0.046440428 | 0.3226290 | 0.3269131 | 0.2360659 | 0.2380683 | 0.23250990 | 0.21631314 |
| LM Trinucleotide Cross-Terms | 0.080706956 | 0.036782904 | 0.3187810 | 0.3285644 | 0.2342924 | 0.2402606 | 0.26666786 | 0.20433756 |

All three models showed increase in predictive accuracy compared to the baseline in all metrics. However, the third model was shown to have less predictive accuracy for the held-out test set than the other two methods. The R-Squared and correlation computed on the training data set for this model were the highest overall, and the RMSE and MAE were the lowest. Inveresly, the R-Squared was the lower, and the RMSE and MAE higher than the other two models. These results indicate that the increase in trainining R-Squared was caused due to over-fitting.

The second model performed the best overall with the highest test R-Squared and Spearman, and lowest test RMSE and MAE.

The fit of the models was visualized plotting the observations against the predictions made by each model for both training and testing data. These plots can be found under "results/models/". These plots were shown to be quite similar for the three models.

## Model Variables

Next the model variables for each of the three linear models were investigated. The coefficients of the models were extracted and ranked by significance level Pr(\>\|t\|). The top variables for the three models are shown in the table below.

|                  |                         |                                   |
|------------------|-------------------------|-----------------------------------|
| **Linear Model** | **Linear Model + Pair** | **Linear Model + Pair + Triplet** |
| *pos20_G*        | *pos19_G*               | *pos20_G*                         |
| pos18_C          | *pos20_G*               | *pos19_G*                         |
| *pos19_G*        | pos19_G:pos20_C         | pos19_G:pos20_C                   |
| pos17_G          | pos19_C                 | (Intercept)                       |
| pos19_C          | pos18_C                 | H3K27ac_peak_overlap              |
| pos17_C          | (Intercept)             | pos18_C:pos20_G                   |
| pos17_A          | pos17_G                 | pos18_C                           |
| pos19_A          | pos17_C                 | pos19_C:pos20_G                   |
| pos20_C          | pos18_G                 | pos1_C:pos2_C                     |
| pos10_A          | pos18_A                 | pos19_C                           |

The three models share a lot of significant variables. A Guanine in positions 20 and 19 is shown as highly significant for all three models. Indeed majority of the significant variables are position variables in position 17-20, and Guanines or Cytosines. The second and third model also share a significant cross-term, pos19_G to pos20_C. The H3K27ac peak overlap variable is also shown to be highly significant.

The least 10 significant variables were also extracted for each of the three models and are shown in the table below.

|  |  |  |
|----|----|----|
| **Linear Model** | **Linear Model + Pair** | **Linear Model + Pair + Triplet** |
| pos1_A | pos6_C:pos7_A | pos9_A:pos10_A |
| pos7_C | pos1_C:pos2_C | pos7_A:pos8_G |
| gene_overlap_count | pos5_G:pos6_C | pos11_G:pos13_C |
| pos8_C | gene_overlap_count | pos18_C:pos19_G:pos20_A |
| pos16_G | pos10_G:pos11_A | pos4_G:pos6_C |
| pos3_A | pos5_G:pos6_G | pos5_C:pos7_C |
| pos2_C | pos3_A:pos4_C | pos1_G:pos2_G:pos3_A |
| pos9_C | pos16_C:pos17_G | pos3_A:pos4_C |
| pos6_A | pos16_A:pos17_A | pos11_G:pos12_C:pos13_A |
| region_type | pos4_C:pos5_A | pos8_C:pos9_G:pos10_G |

The variables of least significance are shown to be position variables in the first half of the guide and various cross-terms between position varialbes. Gene overlap count is also shown as an insignificant variable.

### Retraining Models with Significant Variables

In order to justify the increasing complexity of adding cross-terms to the linear models, a retraining was performed only using significant variables in the intial models. Arbitrarily a significance level of 0.1 was used to determine which variables would be used to train the new models.

The new models used only 35, 94 and 135 of each respective models original set of variables. The models were trained resulting in the following summaries.

**Linear Model:**

``` text
Residual standard error: 0.3255 on 41489 degrees of freedom
Multiple R-squared:  0.04243,   Adjusted R-squared:  0.04162 
F-statistic: 52.53 on 35 and 41489 DF,  p-value: < 2.2e-16
```

**Linear Model + Pair Cross-Terms:**

``` text
Residual standard error: 0.3239 on 41430 degrees of freedom
Multiple R-squared:  0.05326,   Adjusted R-squared:  0.05111 
F-statistic:  24.8 on 94 and 41430 DF,  p-value: < 2.2e-16
```

**Linear Model + Pair + Triplet Cross-Terms**

``` text
Residual standard error: 0.3247 on 41389 degrees of freedom
Multiple R-squared:  0.04953,   Adjusted R-squared:  0.04643 
F-statistic: 15.98 on 135 and 41389 DF,  p-value: < 2.2e-16
```

Overall these models showed a slight decrease in predictive accuracy compared to the original models.

|  |  |  |  |  |  |  |  |  |
|----|----|----|----|----|----|----|----|----|
| **Model** | **R2Train** | **R2Test** | **RMSETrain** | **RMSETest** | **MAETrain** | **MAETest** | **SpearmanTrain** | **SpearmanTest** |
| Linear Model | 0.04243236 | 0.04081382 | 0.3253495 | 0.3278762 | 0.2375673 | 0.2384014 | 0.2005482 | 0.2002227 |
| LM + Dinucleotide Cross-Terms | 0.05326290 | 0.04509475 | 0.3235044 | 0.3271437 | 0.2365310 | 0.2382005 | 0.2220417 | 0.2100304 |
| LM + Di+ Trinucleotide Cross-Terms | 0.04952743 | 0.03827260 | 0.3241420 | 0.3283102 | 0.2369924 | 0.2384865 | 0.2104659 | 0.1959135 |

The evaluation fo the retrained models further confirms that the predictive accuracy of the models had decreased.

The coefficients of these models were investigated to detect possible discrepancies in variable significance. The top 10 most significant variables for each of the retrained models is shown below.

|                  |                         |                                   |
|------------------|-------------------------|-----------------------------------|
| **Linear Model** | **Linear Model + Pair** | **Linear Model + Pair + Triplet** |
| (Intercept)      | (Intercept)             | (Intercept)                       |
| pos20_G          | pos20_G                 | pos19_G                           |
| pos19_G          | pos19_G                 | pos18_C                           |
| pos18_C          | pos18_C                 | pos20_G                           |
| pos17_G          | pos17_C                 | pos17_G                           |
| pos19_C          | pos19_C                 | pos19_C                           |
| pos17_C          | pos17_G                 | pos19_G:pos20_C                   |
| pos17_A          | pos17_A                 | pos17_C                           |
| pos19_A          | pos19_A                 | pos18_G                           |
| pos20_C          | pos19_G:pos20_C         | pos17_G:pos18_C:pos19_C           |

The significant variables are highly similar for the original and retrained models, with Guanines in positions 20 and 19 showing the highest significance still.
