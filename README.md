# Malaria-Dry-Season-Carriage

Data and code for the article “Evidence for exposure dependent carriage of malaria parasites across the dry season: modelling analysis of longitudinal data”. 

There are two folders, one for simulated data and one for previously described data from a longitudinal study in Mali [1,2]. 

 

“Mali data” folder: 

Mali-2012-data.xlsx: data from a longitudinal study in Mali previously described by Tran et al. (2013) [1] and Portugal et al. (2017) [2]. 

Mali-data-analysis.R: analysis of the Mali data including  

age distribution with status of infection (Fig. S12), 

Kaplan-Meier curves for the time from enrollment to infection (Fig. 7a, Fig. S13 and Fig. S15) and clinical malaria (Fig. 7b) and for time from first follow-up to infection (Fig. S16a, Fig. S17 and Fig. S18) and clinical malaria (Fig. S16b) and 

Cox Proportional Hazards models (Table S3, Table S4 and Table S5). 

Note that the Kaplan-Meier curves of the time from enrollment to clinical malaria for carriers and non-carriers (Fig. 7b) is a reproduction of data published by Portugal et al. (2017) [2] who have shown that carriers are more immune to clinical malaria than non-carriers. 

 

“Simulated data” folder: 

“Data” folder: data simulated using “Model_simulations.m”: 

Data-1.mat: simulation of 1,000 individuals with heterogeneous and random FOI and age 

Data-1Ind-low-FOI.mat, Data-1Ind-med-FOI.mat and Data-1Ind-high-FOI.mat: simulated data of a single individuals with low, intermediate or high FOI (Fig. 2, Fig. S1 and Fig. S3) 

Data-2-car-vs-non-car.mat: subset of the simulation of 100,000 individuals for the comparison of the number of infectious bites for carriers and non-carriers in the previous transmission season for individuals with random ages and different FOIs (Fig. 4a and Fig. S5) and the time since the last infection for carriers and non-carriers (Fig. 4b and Fig. S6) 

Simulation of 1,000 ind (homogeneous population):  

These files were too large to be uploaded here. The stochastic simulations can be repeated with the provided code. 

Simulation of 1,000 individuals with the same Force Of Infection (FOI, bites per day during the wet season) for 10 different FOIs. The different FOIs are noted in the README file. 

Simulation of 100,000 ind:  

These files were too large to be uploaded here. The stochastic simulations can be repeated with the provided code. 

100 files each with data of 1,000 simulated individuals (simulated for 20 years and with 10 different FOIs) 

Survival curve data: survival curve data for Data-1.mat and Data-3-x.mat used to visualize survival curves in R with the script Simulated_data_analysis.R. 

Matlab functions: 

intrahost.m: deterministic intrahost dynamics in the model simulations 

KM.m: visualization of data as Kaplan-Meier curves 

person.m: simulates an individual from birth to a certain age and returns parasite concentration and immunity at the specified age 

ttnb.m: computes the time to the next random infectious mosquito bite 

Model_simulations.m: script for all model simulations and plots including: 

Simulation of 1 individual with seasonal biting rate from birth to age 20: visualization of parasite concentration and immunity over time and duration of infections and peak parasite concentration by age (Fig. 2, Fig. S1 and Fig. S2) 

Simulation of 1,000 individuals with random mosquito biting rate and age for 1 year (from the end of the dry season to the end of the next dry season): time to first infection, biting rate distribution, parasite distribution by age, age distribution by carrier status, Parasite Multiplication Rate (PMR) distribution and immunity (Fig. S8), mean immunity by FOI and age heatmap and fraction of carriers by age and exposure (Fig. S9) and characterization of carriers after 1 year of model simulations (Fig. S11). 

Simulation of 100,000 individuals with 10 different biting rates from birth to age 20: age of first parasite carriage over the dry season heatmap (Fig. S3) and boxplot (Fig. S4). 

Simulated_data_analysis.R: script to analyze Data-2-car-vs-non-car.mat and comparing carriers and non-carriers previous exposure, i.e. number of bites in the previous transmission season (Fig. 4a and Fig. S5) and the time since the last infectious bite (Fig. 4b and Fig. S6), change of the status of infection at the end of the dry season after simulation of 1,000 individuals with a heterogeneous FOI for one year (Fig. S10), survival curves for the comparison of homogeneous and heterogeneous infection risk (Fig. 5g and h, Fig. S8a) and time to first infection by FOI for 1,000 individuals with the same FOI (Fig. S7 and p-values in Table S2). 

 

References: 

[1] Tran, T.M., et al., An Intensive Longitudinal Cohort Study of Malian Children and Adults Reveals No Evidence of Acquired Immunity to Plasmodium falciparum Infection. Clinical Infectious Diseases, 2013. 57(1): p. 40-7. 

[2] Portugal, S., et al., Treatment of Chronic Asymptomatic Plasmodium falciparum Infection Does Not Increase the Risk of Clinical Malaria Upon Reinfection. Clinical Infectious Diseases, 2017. 64(5): p. 645-53. 
