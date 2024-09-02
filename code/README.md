# Code

This folder contains all the necessary code create the dyadic data and run Bayesian logistic regressions and Mantel tests.


## Data

The datasets can be found in the `data` folder.

- The data to extract nearly-private rare allele similarity was obtained from the GitHub folder referred to in Fontsere et al. 2022 (https://github.com/kuhlwilm/rareCAGA).


## Code
- `Get_nepra_proportions.R`: code to generate dyadic nearly-private rare allele similarity proportions (data see above).
- `Create_dyadic_data.R`: code to generate the m.dyads data file including all dyadic variables in the data folder.
- `BRMS_models.R`: code to run Bayesian logistic regressions including models testing for robustness.
- `Mantel_tests.R`: code to run Mantel tests.


