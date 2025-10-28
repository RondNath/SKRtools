# SKRtools: Distributions analysis

This package offers simple tools to study complex distributions and their dynamic: 
- the computation of the individual moments (mean, variance, skewness and kurtosis),
- a Skewness-Kurtosis Relationship (SKR) analysis and related regression parameters (Slope, Y-intercept, R², distance from predicted SKR, distance from reference SKR)

The SKR framework consider the complexity found in distribution shapes (non Gaussian) and analyse similarities in these shapes. 

Installing the package
---------------------------------
The development version on Github can be installed using the devtools package:

devtools::install_github("RondNath/SKRtools")


Step 1 - Function: Compute the individual moments of distributions - "Moments"
---
Compute the mean, the variance, the skewness and the kurtosis


Step 2 - Function: SKR analysis of the distributions - "SKRparam"
---
-   Calculate the SKR for the distributions group per statistic factors 
-   Extract the parameters:
    - Slope
    - Y-intercept  
    - R²
    - distance from predited distribution family (Slope and Y-intercept predicted)
    - distance from a reference distribution family (default skew-uniform family: slope = 1; intercept = 1.86)
    - CV of distance from a predited distribution family (Slope and Y-intercept predicted)
    - CV of distance from a reference distribution family
