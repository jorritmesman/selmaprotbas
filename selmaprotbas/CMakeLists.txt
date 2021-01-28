add_library(fabm_models_selmaprotbas OBJECT
            selmaprotbas_model_library.F90
            selmaprotbas.F90
            phytoplankton.F90
            zooplankton.F90
           )

add_dependencies(fabm_models_selmaprotbas fabm_base)

#-----------------------------------------------------------------------
# Copyright (C) 2017 - Karsten Bolding and Jorn Bruggeman (BB)         !
#-----------------------------------------------------------------------
