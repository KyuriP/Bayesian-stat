
## BAYESIAN ASSIGNMENT README FILE
__________________________________________________________

## AUTHOR: KYURI PARK (student number: 5439043)
## This file provides the specification of "BayesianAssignment_KyuriPark.zip" file.
__________________________________________________________

## CONTENTS
- DATA file: Birthweight.csv
- REPORT file: Report_KP.pdf
- CODE file: Code_KP.R
- README file: Readme_KP.txt
__________________________________________________________

## DATA DESCRIPTION
Birthweight.csv contains the information on 47 new-born babies and their parents.
It is a data frame with 47 cases and 16 variables.
In the main analysis, 3 variables are used (i.e., "Birthweight": birth weight of babies, "mage": mother's age, "mnocig": number of cigarettes smoked by mothers per day).
In line 9-19 of the code file, you can find how the data is loaded and processed.
__________________________________________________________

## REPORT DESCRIPTION
Report_KP.pdf consists of 5 sections.

1) INTRODUCTION: 
It contains the detailed description of the data, and specification of research questions.(page 1)

2) METHODS: 
It contains the explanation on the model specification, and the following analysis procedure.(page 2)

3) RESULTS: 
It consists of 6 sub-sections.
   3-1) Convergence assessment: includes trace plot, density plot, 
        autocorrelation plot, Gelman-Rubin statistics, and MC error.(page 2-3)
   3-2) Model assumption check: tests homoscedasticity and 
        normality of residuals using posterior predictive check.(page 3-4)
   3-3) Parameter estimates: contains the result of analysis including mean, median, 
        standard deviation, credible intervals, and highest density intervals.(page 4)
   3-4) Model comparison: compares WAIC values for 4 different models.(page 4-5)
   3-5) Hypothesis evaluation: contains Bayes factors and 
        PMP of two informative hypotheses + Hu (unconstrained hypothesis).(page 5)
   3-6) Bayesian causal inference: presents the estimate 
        of ACE (average causal effect) of mother's smoking.(page 5-6)

4) CONCLUSION: 
It contains the summary of the results and comparison of frequentist and Bayesian approach.(page 6)

5) REFERENCES (page 6)
__________________________________________________________

## CODE DESCRIPTION
Code_KP.R consists of 9 sections. 
All code is necessary in order to obtain the results that are presented in the report. 
In total, 531 lines are used: 334 lines are the code and the rest 197 lines are comment + empty space.

1) Preparation (line 1-19; corresponds to INTRODUCTION section of the report):
It contains the code for loading the required libraries, read the data, and process the data.

2) Gibbs Sampler with MH step (line 21-104; corresponds to METHODS section of the report): 
It contains mini-functions that sample each of the parameters, 
a big sampling function that loops over the mini-functions, 
and a function that implements multiple chains. 
Basically, calling the last function ("multiple_chains") would 
run the Gibbs sampler with MH-step incorporated for b1.
   2-1) Mini function of Gibbs sampling for tau (line 27-32)
   2-2) Mini function of Gibbs sampling for b0 (line 34-39)
   2-3) Mini function of MH sampling for b1 (line 41-52)
   2-4) Mini function of Gibbs sampling for b1 (line 54-59)
   2-5) Mini function of Gibbs sampling for b2 (line 61-66)
   2-7) Big sampling function looping through mini functions (line 68-82)
   2-8) Function that implement multiple chains (line 84-104)

3) Perform MCMC Sampling: runs the sampler function (line 106-124; corresponds to RESULTS section of the report):
   3-1) Specify initial values (line 109-110)
   3-2) Specify parameters for prior distribution (line 112-114)
   3-3) Obtain the samples from the posteriors (line 116-124)

4) Assess Convergence (line 126-167; corresponds to RESULTS 3-1 section of the report):
   4-1) Plotting trace plot, density plot, and autocorrelation plot (line 129-145)
   4-2) Compute Gelman-Rubin statistic (line 146-164)
   4-3) Compute MC error (line 165-167)

5) Posterior Predictive Check (line 169-244; corresponds to RESULTS 3-2 section of the report):
   5-1) Generate the simulated data (line 175-188)
   5-2) Compute the simulated and observed residuals (line 190-200)
   5-3) Test homoscedasticity (line 202-219)
   5-4) Test normality of residuals (line 221-244)

6) Parameter Estimates (line 246-275; corresponds to RESULTS 3-3 section of the report):
   6-1) Compute the mean (line 250)
   6-2) Compute the median (line 252)
   6-3) Compute the standard deviation (line 254)
   6-4) Compute the 95% credible intervals (line 256)
   6-5) Compute the 95% highest density intervals (line 258-275)

7) Model Comparison (WAIC) (line 277-408; corresponds to RESULTS 3-4 section of the report):
   7-1) Obtain the posterior using the sampler for Model 1 (line 288-316)
   7-2) Obtain the posterior using the sampler for Model 2 (line 319-348)
   7-3) Obtain the posterior using the sampler for Model 3 (line 351-380)
   7-4) Define the function that computes WAIC ("getWAIC") (line 383-401)
   7-5) Compute WAIC for all models using the "getWAIC" function (line 404-408)

8) Hypothesis Evaluation (Bayes Factor) (line 410-487; corresponds to RESULTS 3-5 section of the report):
   8-1) Obtain the posterior for the standardized data (line 414-424)
   8-2) Specify the means, and covariance matrix for the parameters of interest (line 432-436)
   8-3) Compute the Bayes factors for H1 (line 440-460)
   8-4) Compute the Bayes factors for H2 (line 462-482)
   8-5) Compute posterior model probabilities (PMP) for H1, H2, and Hu (line 484-487)

9) Bayesian Causal Inference (IPW) (line 489-531; corresponds to RESULTS 3-6 section of the report):
   9-1) Run Bayesian logistic regression (line 498-501)
   9-2) Compute the posterior of propensity scores (line 506-508)
   9-3) Estimate the ACE using inverse propensity weighting (line 512-523)
   9-4) Compute the mean and 95% credible interval of ACE (line 528-531)




