# Selmaprotbas
Adaptation to the SELMA model, coupled to the Framework for Aquatic Biogeochemical Models (FABM)

This model describes temporal dynamics of inorganic nutrients, oxygen, phytoplankton, and zooplankton. Coupled to a physical model, it
does this for multiple layers or cells in the water column. The SELMA model which it was based on
("github.com/fabm-model/fabm/tree/master/src/models/selma") was itself derived from the ERGOM model (Ecological ReGional Ocean Model,
"https://ergom.net/"). See the comments at the start of ./selmaprotbas/selmaprotbas.F90 for contributions and references. 

## Model functionality
The Selmaprotbas model describes dynamics of:
- Phosphorus, nitrogen, and silicon (sediment, detritus, and dissolved)
- Oxygen concentration
- An unrestricted amount of phytoplankton and zooplankton groups

Simulated processes include:
- Sedimentation and resuspension of sediment
- Mineralisation, bio-resuspension, and burial
- (De-)Nitrification
- Phosphorus adsorption to iron (influenced by oxygen)
- Temperature-dependent phytoplankton growth rates
- Nutrient and light limitation in phytoplankton
- Nitrogen-fixation and buoyancy regulation
- Zooplankton grazing (potentially varying preference for prey)

Nutrient ratios in phyto- and zooplankton are kept constant, but may vary in detritus or the sediment. 

All parameters can be accessed in text (YAML) format. 

### What was changed from the SELMA model?
#### v1.0
- The model is now carbon-based, instead of nitrogen-based (i.e. where possible, quantities are now expressed in mmol C, instead of mmol N)
- A silicon cycle was added
- Variable nutrient ratios in detritus and fluff are now possible, as detritus and fluff are now put in separate variables for each element.
  This could cause problems with mass conservation in SELMA if detritus, sediment, and phytoplankton had not all the same nutrient ratio. 
- Calculation of chlorophyll is now exclusively based on phytoplankton C content and the Yc parameter (C:Chl ratio)
- An option to add buoyancy regulation in phytoplankton was added. Nutrient, temperature, and light thresholds can control vertical velocity.
- The maximum growth rate can now be scaled with temperature using the beta parameter, following Reynolds et al. 2001 (doi:10.1016/S0304-3800(01)00330-1)
- Light limitation now also follows Reynolds et al. 2001, using the alpha_light parameter
- Nutrient uptake parameter (alpha) can now differ per nutrient
- Zooplankton predation now also conserves mass if prey has a different nutrient ratio
- maxsed parameter has been removed
- Iron-bound phosphorus is included in the calculation of total-P
#### v1.0.1
- Compatible with FABM v1.0.4 (and it was added to FABM-plus, see below)
#### v1.1
- A DOM module was added, based on the DOMCAST model.
- Option to use 24h-averaged light for phytoplankton growth equations
- More phytoplankton temperature limitation options (WET/Lehman) and all handling of T-limitation occurs before nutrient and light limiation are checked. The "beta" parameter is no longer used by default, but only when tlim==3.
- An extra phytoplankton light limitation function (similar to the llim=1 option, but more gradual)
- Nutrient uptake constants are expressed in the concentration of each element. 
- Option to multiply phytoplankton lightlim and nutlim (rather than uses the minimum value of them)
- Update to some constants in the mineralisation equations
- Option to use anammox processes and sulphate-based mineralisation
- Nitrification rate can be specified
- Output of bottom diagnostics is varied over depth

## How to install?
Selmaprotbas can only be run together with the Framework for Aquatic Biogeochemical Models (FABM, "https://github.com/fabm-model/fabm"),
for example coupled to the General Ocean Turbulence Model (GOTM, "https://github.com/gotm-model/code"). 

To install Selmaprotbas, you first need to install a physical model which is also coupled to FABM. For example, when using GOTM, follow the
instructions on the GOTM homepage ("https://gotm.net/portfolio/software/").

### Method 1 (advised): 
When you run CMake, add -DFABM_INSTITUTES=selmaprotbas -DFABM_SELMAPROTBAS_BASE=<SELMAPROTBAS_DIR>, with <SELMAPROTBAS_DIR> the directory
in which the CMakeLists.txt file and the selmaprotbas folder on this GitHub page are placed. 

### Method 2:
Once the physical model is installed, go to the source code and copy the "selmaprotbas" folder on this GitHub page, with its contents,
in the "./code/extern/fabm/src/models" folder (next to the folder for the other included models, such as SELMA). Then, open the file
"./code/extern/fabm/src/CMakeLists.txt" and manually add "selmaprotbas" to the list of DEFAULT_INSTITUTES. 

After this is done, you need to run CMake again and compile the physical model. Now, Selmaprotbas should be ready for use.

### FABM-plus
Since 2025, this model is part of [FABM-plus](https://github.com/fabm-model/fabm-plus), and can be used together with a suite of physical models
that incorporate the FABM framework. If installing a physical model that comes with FABM-plus (e.g. [PyGETM](https://github.com/BoldingBruggeman/getm-rewrite)),
Selmaprotbas (and other biogeochemical models) will be directly included without having to do the installation steps described above.

## How to run the model?
You can run Selmaprotbas the same way you can run any other model coupled to a physical model by means of FABM. Here we will use GOTM as an example,
but the process is similar when using other physical models. 

In the gotm.yaml file, you need to switch the option "fabm/use" to "true". You also need to have a "fabm.yaml" file in the same directory as
gotm.yaml. Use the fabm.yaml file that we provide on this GitHub page. In this file you can change parameter and initialisation values,
or add/remove phytoplankton and zooplankton groups. If a parameter is not defined in the fabm.yaml file, it is set to its default (you
can therefore add more parameters than are in our provided file, or remove parameters that you plan to keep at the default anyways).

The only factor related to biogeochemistry you can't access in fabm.yaml file or the "fabm" section of gotm.yaml, is influx of biogeochemical
variables in the inflow. In the gotm.yaml file, you can add new subsections for each inflow (next to flow, temp, salt...) for e.g. selmaprotbas/po,
selmaprotbas/nn, etc. Follow the same structure as for the flow, temp, and salt subsections.  

If "fabm/use" is set to "true" in the gotm.yaml file, when you run GOTM, the model described in the fabm.yaml file is run coupled to it (in our
case, Selmaprotbas). In the GOTM output file, you will be able to access the results of the combined GOTM-FABM-Selmaprotbas run. 
