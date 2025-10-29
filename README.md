# SKRtools: Distributions analysis

This package offers simple tools to study complex distributions (non Gaussian) and their dynamic: 
- The computation of the individual moments (mean, variance, skewness and kurtosis),
- The Skewness-Kurtosis Relationship (SKR) analysis and related regression parameters (Slope, Y-intercept, R², distance from predicted SKR, distance from reference SKR)

The SKR framework consider the complexity found in distribution shapes (non Gaussian) and analyses similarities in these shapes. 

Installing the package
---------------------------------
The development version on Github can be installed using the devtools package:

devtools::install_github("RondNath/SKRtools")


Step 1 - Function: Compute the individual moments of distributions - "Moments"
---
Compute the mean, the variance, the skewness and the kurtosis for groups of values (related by explanatory variables) with the possibility to add weighting factors


Step 2 - Function: SKR analysis of the distributions - "SKRparam"
---
-   Analyse the SKR for a group of distributions (related by explanatory variables)
-   Extract the parameters:
    - Slope
    - Y-intercept  
    - R²
    - distance from predited distribution family (Slope and Y-intercept predicted, index of stability in the shapes of the distributions)
    - distance from a reference distribution family (default skew-uniform family: slope = 1; intercept = 1.86, the skew-uniform is a potential maximum evenness so interpreted as an index of evenness in the shapes of the distributions)
    - CV of distance from a predited distribution family (Slope and Y-intercept predicted)
    - CV of distance from a reference distribution family
