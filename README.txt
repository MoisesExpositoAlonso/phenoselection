################################################################################
##### R script implementing selection gradients from Lande & Arnold 1983  ######
################################################################################
###################      by Moises Exposito-Alonso           ###################
################################################################################

This script is free to use, modify and share, but it comes with no warranty.

In total, there are seven equations implemented which deal with linear, disruptive, blancing and correlative selection, and, if heritabilities and genetic correlations are provided, response to selection is calculated. 

Specifically from the paper of Lande & Arnold 1983, linear gradients and coefficients are implementations of the equation 4 and 6c. Response to selection is the implementation of equation 6b. Quadratic coefficients and gradients are implementations of the equations 13b and 14a, respectively. To infer changes in phenotypic variance and additive variances due to disruptive or balancing selection, equations 15a and 12. 

This scripts handles missing values by removing them, then only complete information is used. If specific treatments of missing values are required, for instance inputing from mean or including and indicator variable to treat it as "invisible fraction" (Moorad and Wade 2013). Whatever the transformation or inputation technique is, it must be done prior to the analyses by the user. 

Additionally, there is a step previos to calculation of gradients that is normalization of phenotypes, substracting mean and deviding by standard deviation. This allows to compute variance standarized selection gradients wich express the gain in fitnes by standard deviation of phenotypes, something that allows comparisons across phenotypes and with other studies (Kingsolver et al. 2001). Regarding fitness, equations also transforme it to relative fitness.




See cited papers:

Lande, R., and S. J. Arnold. 1983. The measurement of selection on correlated characters. Evolution (N. Y). 37:1210–1226.

Kingsolver, J. G., H. E. Hoekstra, J. M. Hoekstra, D. Berrigan, S. N. Vignieri, C. E. Hill, A. Hoang, P. Gibert, and P. Beerli. 2001. The strength of phenotypic selection in natural populations. Am. Nat. 157:245–61.

Moorad, J. a, and M. J. Wade. 2013. Selection gradients, the opportunity for selection, and the coefficient of determination. Am. Nat. 181:291–300.

