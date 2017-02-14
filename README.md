# Calibration and Sharpness Examples from Gneiting+07

The script examples.py and associated plotting outputs implement the probabilistic model calibration and sharpness examples introduced in Section 2.2 of Gneiting et al. 2007 (DOI: 10.1111/j.1467-9868.2007.00587.x)

**Note:** The results from the script for most of the forecasters appear to match the published results from Gneiting+07, except for Hammil's forecaster.  The results here for Hammil have significantly higher PIT coverage, less sharpness, and a lower CRPS than Gneiting+07 report.  The reason for the discrepancy is unclear, but presumably due to some coding error for Hammil's forecaster in this script.

# Discussion of Results

## Basic diagnostics

As a first diagnostic, we can simply compare the observation in each simulation with a random draw from the forecaster's predictive distribution.  Presumably, a good forecaster will generate predictions strongly correlated with the observations.  For example, the ideal forecaster achieves correlation $\rho=0.5$.  Interestingly, the climatological forecaster's predictions are uncorrelated with the observations, despite being probabilistically and marginally calibrated.  The sign-based forecaster's predictions are, by design, negatively correlated with the observations, yet still achieve exceedence and marginal calibration.

![][diag_scatter]

[diag_scatter]: http://172.31.51.217/nsanders/Gneiting07-probability-examples/raw/dbd69d57f2130a2a5c08bdf6e858d4bc5a089cd8/diagnostic_scatter.png

This diagnostic plots the marginal distribution of the observations and each example distribution of predictions.  Agreement between the marginal distributions tests for marginal calibration.  Note that Gneiting et al. propose using the CDF rather than a probability distribution histogram chart for this purpose (Section 3.2).  This diagram should indicates marginal calibration for all examples except the mean-biased, unfocused, and (more subtly) Hammil's forecaster.

![][diag_marginals]

[diag_marginals]: http://172.31.51.217/nsanders/Gneiting07-probability-examples/raw/dbd69d57f2130a2a5c08bdf6e858d4bc5a089cd8/diagnostic_marginals.png

## Marginal Calibration Test

Gneiting recommend a stronger diagnostic plot based on the difference between the CDFs of the predictive distributions and observations.  This should indicate the same results.

![][diag_marginal_calibration]

[diag_marginal_calibration]: http://172.31.51.217/nsanders/Gneiting07-probability-examples/raw/f3ebd2aacd2e522746cddb5dd88da1738e03e59d/diagnostic_marginal_calibration.png

## Probabilistic Calibration Test

The Probability Integral Transform (PIT) is a 'traditional' diagnostic for assessing calibration.  Uniformity, or flatness, of the PIT indicates probabilistic calibration

![][diag_PIT]

[diag_PIT]: http://172.31.51.217/nsanders/Gneiting07-probability-examples/raw/dbd69d57f2130a2a5c08bdf6e858d4bc5a089cd8/diagnostic_PIT.png

There is no emperical test for exceedence calibration.

## Sharpness Test

In addition to assessing these forms of calibration, Gneiting et al. suggest evaluating forecasters on the basis of the sharpness of their predictions.  This violin plot compares the interval width distributions (distributed across the simulations) for each of the example forecasters, with the ideal forecaster of course representing optimum sharpness.

![][diag_sharp]

[diag_sharp]: http://172.31.51.217/nsanders/Gneiting07-probability-examples/raw/f3ebd2aacd2e522746cddb5dd88da1738e03e59d/diagnostic_sharpness_violin.png

## Brier Scores and CRPS

The Brier score plot shows the continuous ranked probability score (CRPS) as a function of the threshold value.  The integral of these curves is the aggregate CRPS score.  Curves that lie lower at a given threshold value have a superior mix of sharpness and calibration in that regime.

![][diag_brier]

[diag_brier]: http://172.31.51.217/nsanders/Gneiting07-probability-examples/raw/f3ebd2aacd2e522746cddb5dd88da1738e03e59d/diagnostic_Brier.png

Finally, the aggregate CRPS itself can be plotted as Gneiting's ultimate meta-score for predictive performance of a forecaster.  Note that the sign based forecaster, which produced predictions that were actually negatively correlated with the observations, is penalized the most by CRPS.

![][diag_crps]

[diag_crps]: http://172.31.51.217/nsanders/Gneiting07-probability-examples/raw/f3ebd2aacd2e522746cddb5dd88da1738e03e59d/diagnostic_crps.png

Note that the CRPS score picks winners among the forecasters fairly consistent with the correlations we looked at in the very first scatterplot diagnostic.  Since the predictive correlation is a point estimate and does not take sharpness into account, it might be expected that the sharpness score primarily explains the differences between the two metrics.  Among these examples, that seems not to be the case.  In particular, the mean-based forecaster has a relatively high (poor) CRPS score despite quite a high (good) correlation, and yet it is among the sharper distributions.

![][diag_crps_cor]

[diag_crps_cor]: http://172.31.51.217/nsanders/Gneiting07-probability-examples/raw/b3e1ef387780061a899e1267f2d8d38688ebf7c0/diagnostic_crps_vs_cor.png



