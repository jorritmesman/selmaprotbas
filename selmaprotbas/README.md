# SELMAPROTBAS 
Starts a branch of SELMA-model at GIT version:

	SHA1 ID: 59ff624f3c9a224917ceb62c290bd0d34257dbc2
	Author: Jorn Bruggeman <jorn@bolding-bruggeman.com>  2019-09-05 09:57:15
	Committer: Jorn Bruggeman <jorn@bolding-bruggeman.com>  2019-09-05 09:57:15
	Parent: 30c190621cd07cea8ab55fca921904afe13e662c (re-activated friction by flat bottom)
	Child:  0000000000000000000000000000000000000000 (Local uncommitted changes, not checked in to index)
	Branches: lake, remotes/old-origin/lake, remotes/origin/lake
	Follows: v5.3
The model has changed from being based on nitrogen to now base all calculations on carbon, has fixed issues
with calculations of chlorophyll, maintains the nutrient balance in case of different ratios between groups
in the model, and other changes as described in the supplementary material of Mesman et al. (2022), 
"Future  effects of storms on nutrient distribution and phytoplankton dynamics".

* SELMA (see below)
* PROTBAS - PROTECH based phytplankton model (Markensten & Pierson)
* PROTECH - Phytoplankton RespOnses To Environmental CHange  (Reynolds et al. 2001)
* DOMCAST - *See Clayer et al., 2021, doi:10.1029/2021JG006359* (PROGNOS project, Raoul-Marie Couture, José-Luis Guerrero, François Clayer)

Adding PROTECH phytoplankton behaviour to SELMA lake branch within FABM-GOTM framework

## SELMA - Simple EcoLogical Model for the Aquatic 
## Selmaprotbas = SELMA + PROTECH

See description of development history and contributors in file selmaprotbas.F90. 
