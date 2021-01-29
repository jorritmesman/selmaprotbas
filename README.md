# selmaprotbas
Adaptation to the SELMA model, coupled to the Framework for Aquatic Biogeochemical Models (FABM)

## How to install?
Selmaprotbas can only be run together with the Framework for Aquatic Biogeochemical Models (FABM, "https://github.com/fabm-model/fabm"),
coupled to the General Ocean Turbulence Model (GOTM, "https://github.com/gotm-model/code"). 

To install Selmaprotbas, you first need to install GOTM, following the instructions on the GOTM homepage ("https://gotm.net/portfolio/software/").
Once GOTM is installed, go to the source code and place the "selmaprotbas" folder on this GitHub page, with its contents, in the 
"./code/extern/fabm/src/models" folder (next to the folder for the other included models, such as SELMA). Then, open the file
"./code/extern/fabm/src/CMakeLists.txt" and manually add "selmaprotbas" to the list of DEFAULT_INSTITUTES. 

After this is done, you need to run CMake again and compile GOTM, following the installation instructions on the GOTM homepage. Now,
Selmaprotbas should be ready for use.

## How to run the model?
You can run Selmaprotbas the same way you can run any other model coupled to GOTM by means of FABM. In the gotm.yaml file, you need to switch the
option "fabm/use" to "true". You also need to have a "fabm.yaml" file in the same directory as gotm.yaml. Use the fabm.yaml file that we provide 
on this GitHub page. In this file you can change parameter and initialisation values, or add/remove phytoplankton and zooplankton groups. If a
parameter is not defined in the fabm.yaml file, it is set to its default (you can therefore add more parameters than are in our provided file,
or remove parameters that you plan to keep at the default anyways).

The only factor related to biogeochemistry you can't access in fabm.yaml file or the "fabm" section of gotm.yaml, is influx of biogeochemical
variables in the inflow. In the gotm.yaml file, you can add new subsections for each inflow (next to flow, temp, salt...) for e.g. selmaprotbas/po,
selmaprotbas/nn, etc. Follow the same structure as for the flow, temp, and salt subsections.  

If "fabm/use" is set to "true" in the gotm.yaml file, when you run GOTM, the model described in the fabm.yaml file is run coupled to it (in our
case Selmaprotbas). In the GOTM output file, you will be able to access the results of the combined GOTM-FABM-Selmaprotbas run. 
