##########################################################################################
##### R script implementing bivariate selection gradients from Lande & Arnold 1983  ######
##########################################################################################
###############         by Moises Exposito-Alonso            #############################
###############      moisesexpositoalonso@gmail.com          #############################
##########################################################################################

These scripts are free to use, modify and share, as long as this repo is acknowledged. It comes with no warranty, however, emails with bugs or questions are wellcomed.

This code was used for the analyses in Exposito-Alonso et al. (201X). It is meant for 2 phenotypes and a measure of fitness; although can be extended to a different number of phenotypes with little tuning.

In total, there are seven equations implemented dealing with linear, disruptive, blancing and correlative selection, and, if heritabilities and genetic correlations are provided, response to selection is calculated. Significance values are calculated using bootstrapping.

Specifically, the equations from the paper of Lande & Arnold 1983 are: Equation 4 and 6c for linear gradients and coefficients. Equation 6b for response to selection. Equations 13b and 14a for quadratic coefficients and gradients, respectively. Equations 15a and 12 are used to infer changes in phenotypic variance and additive variances due to disruptive or balancing selection. 

This script handles missing values by removing them. Only complete information is used. If specific treatments of missing values are required, e.g. inputing from mean or dummy variable to infer the "invisible fraction" (Moorad and Wade 2013), those must be done prior to the analyses using my scripts. 

Additionally, I recommend to do a prior step before the calculation of gradients, normalization of phenotypes (substracting mean and deviding by standard deviation). This allows to compute variance standarized selection gradients wich express the gain in fitnes by standard deviation units of phenotypes, something that allows comparisons across phenotypes and with other studies (Kingsolver et al. 2001). NO need of fitness transformation to relative fitness, it is already implemented in the scripts.

Bibliography

Exposito-Alonso, M., Brennan, A., Alonso-Blanco, C., Picó, F.X. 201X. Local and temporal adaptation of life history in Arabidopsis thalianan. (to be submitted) 

Lande, R., and S. J. Arnold. 1983. The measurement of selection on correlated characters. Evolution (N. Y). 37:1210–1226.

Kingsolver, J. G., H. E. Hoekstra, J. M. Hoekstra, D. Berrigan, S. N. Vignieri, C. E. Hill, A. Hoang, P. Gibert, and P. Beerli. 2001. The strength of phenotypic selection in natural populations. Am. Nat. 157:245–61.

Moorad, J. a, and M. J. Wade. 2013. Selection gradients, the opportunity for selection, and the coefficient of determination. Am. Nat. 181:291–300.

