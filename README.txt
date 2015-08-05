######################################################################################################
### phenoselection, R script implementing phenotypic selection gradients from Lande & Arnold 1983  ###
######################################################################################################
#######################           by Moises Exposito-Alonso        ###################################

This script is free to use and share, but come with no warrante.

In total, there are seven equations that deal with linear, disruptive,blancing and correlative selection,
and if heritabilities and genetic correlations are provided, response to selection is calculated. Specifically,
linear gradients and coefficients are implementations of the equation 4 and 6c. Response to selection is the
implementation of equation 6b. Quadratic coefficients and gradients are implementations of the equations 13b and 14a,
respectively. To understand changes in phenotypic variance and additive variances due to disruptive or balancing
selection, equations 15a and 12. It can handle missing values, but they are just removed, then only complete
information is used and if specific treatments of missing values are required, should be transformed or imputed a
priori by the user. Although it can be modified, our equations also have a previous step of normalization of
phenotypes, substracting mean and deviding by standard deviation, to allow comparisons across phenotypes and with
other studies. Regarding fitness, equations also transforme it to relative fitness.





