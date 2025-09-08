#include "fabm_driver.h"

module selmaprotbas_dom

   use fabm_types

   implicit none

   private
   ! This selmaprotbas module is built on the DOMCAST code base
   ! Original author(s) of DOMCAST: Members of the PROGNOS project
   ! Implemented into selmaprotbas by Jorrit Mesman

   type, extends(type_base_model), public :: type_selmaprotbas_dom
      ! Variable identifiers
	  type (type_state_variable_id)        :: id_dom_a, id_dom_b ! a is labile, b is semi-labile
	  type (type_state_variable_id)        :: id_o2, id_nn, id_aa, id_dd_c
	  type (type_dependency_id)            :: id_temp, id_swf
	  type (type_diagnostic_variable_id)   :: id_extinc, id_rate11, id_rate12, id_rate21, id_rate22, id_rate3a, id_rate3b, id_rate4a, id_rate4b, id_rate5, id_light !reaction rates
	  
	  ! Model parameters
	  real(rk) :: theta, km_o2, km_no3, k_inh_o2
	  real(rk) :: k_om1, k_om2, oc_dom, qy_dom, f_par, e_par, k_leach, q_leach
	  real(rk) :: k_floc, mole_per_weight, n_use_factor, ext_coef_a, ext_coef_b
	  logical  :: diagnostics
   contains
      procedure :: initialize
	  procedure :: do
      ! Reference model procedures here.
   end type

contains

   subroutine initialize(self, configunit)
      class (type_selmaprotbas_dom), intent(inout), target :: self
      integer,                       intent(in)            :: configunit
	  
	  ! Local variables
	  real(rk), parameter :: d_per_s = 1.0_rk / 86400.0_rk
	  real(rk), parameter :: y_per_s = d_per_s / 365.25_rk
 
      ! Register model parameters and variables here.
	  
	  ! Register model parameters
	  call self%get_parameter(self%km_o2, 'km_o2', 'mmol/m3', 'respiration', default=1.23e-2_rk)
	  call self%get_parameter(self%km_no3, 'km_no3', 'mmol/m3', 'denitrification', default=0.01_rk)
      call self%get_parameter(self%k_inh_o2, 'k_inh_o2', 'mmol/m3', 'inhibitation of denitrification by O2', default=0.33_rk)
	  call self%get_parameter(self%k_om1, 'k_om1', '1/yr', 'labile OM degradation rate', default=1.0_rk, scale_factor=y_per_s)
      call self%get_parameter(self%k_om2, 'k_om2', '1/yr', 'semi-labile OM2 degradation rate', default=0.1_rk, scale_factor=y_per_s)
	  call self%get_parameter(self%theta, 'theta',  '-', 'Temperature adjustment coefficient', default=1.047_rk)
	  call self%get_parameter(self%oc_dom, 'oc_dom', 'm2/mgDOM', 'Optical cross section of DOM', default=0.01_rk)
	  call self%get_parameter(self%qy_dom, 'qy_dom', 'mgDOM / mol PAR', 'Quantum yield', default=0.1_rk)
	  call self%get_parameter(self%f_par, 'f_par', '-', 'PAR fraction', default=0.45_rk)
      call self%get_parameter(self%e_par, 'e_par', 'J/mol', 'Energy of PAR photons', default=240800._rk)
	  call self%get_parameter(self%k_floc, 'k_floc', '1/d', 'flocculation rate', default=0.0006_rk, scale_factor=d_per_s)
	  call self%get_parameter(self%mole_per_weight, 'mole_per_weight', 'molC/gDOM', 'mol C per g DOM', default=0.0416_rk) ! Default assumes 50% C/DW weight ratio and 12.01 g/mole molar mass
	  call self%get_parameter(self%n_use_factor, 'n_use_factor', 'mmol N / mg DOM', 'Factor of bacterial mineralisation as N flux', default=0.8_rk)
	  call self%get_parameter(self%k_leach, 'k_leach', '1/d', 'Leaching constant, at 20 degrees', default=0.001_rk, scale_factor=d_per_s)
	  call self%get_parameter(self%q_leach, 'q_leach', '1/K', 'temperature dependence of leaching', default=0.02_rk)
	  call self%get_parameter(self%ext_coef_a, 'ext_coef_a', 'm2/mgDOM', 'Linear coefficient DOM light influence', default=0.1_rk)
	  call self%get_parameter(self%ext_coef_b, 'ext_coef_b', '-', 'Exponential coefficient DOM light influence', default=1.22_rk)
	  call self%get_parameter(self%diagnostics, 'diagnostics', '-', 'toggle diagnostic output', default=.false.)
	  
	  ! Register state variables
      call self%register_state_variable(self%id_dom_a, 'dom_a', 'mg/m3', 'DOM - labile', initial_value=0.0_rk, minimum=0.0_rk)
      call self%register_state_variable(self%id_dom_b, 'dom_b', 'mg/m3', 'DOM - semi-labile', initial_value=0.0_rk, minimum=0.0_rk)
	  
	  ! Register dependencies on external state variables
	  call self%register_state_dependency(self%id_o2, 'o2', 'mmol O2/m3', 'oxygen')
	  call self%register_state_dependency(self%id_nn,'nn','mmol N/m3', 'nitrate')
	  call self%register_state_dependency(self%id_aa,'aa','mmol N/m3', 'ammonium')
	  call self%register_state_dependency(self%id_dd_c,'dd_c','mmol C/m3', 'carbon detritus')
	  
	  ! Register light extinction diagnostic and light aggregation
	  call self%register_diagnostic_variable(self%id_extinc, 'extinc', '1/m', 'shading by CDOM')
	  call self%add_to_aggregate_variable(standard_variables%attenuation_coefficient_of_shortwave_flux, self%id_extinc)
	  
	  ! register diagnostic variable -> rates included only for debugging purposes
      if(self%diagnostics) then
		  call self%register_diagnostic_variable(self%id_rate11, 'rate11', 'mgDOM/m3/s', 'rate11 - dom_a_o2')
		  call self%register_diagnostic_variable(self%id_rate12, 'rate12', 'mgDOM/m3/s', 'rate12 - dom_b_o2')
		  call self%register_diagnostic_variable(self%id_rate21, 'rate21', 'mgDOM/m3/s', 'rate21 - dom_a_N')
		  call self%register_diagnostic_variable(self%id_rate22, 'rate22', 'mgDOM/m3/s', 'rate22 - dom_b_N')
		  call self%register_diagnostic_variable(self%id_rate3a, 'rate3a', 'mgDOM/m3/s', 'rate3a - photo_ox_a')
		  call self%register_diagnostic_variable(self%id_rate3b, 'rate3b', 'mgDOM/m3/s', 'rate3b - photo_ox_b')
		  call self%register_diagnostic_variable(self%id_rate4a, 'rate4a', 'mgDOM/m3/s', 'rate4a - dom_a_floc')
		  call self%register_diagnostic_variable(self%id_rate4b, 'rate4b', 'mgDOM/m3/s', 'rate4b - dom_b_floc')
		  call self%register_diagnostic_variable(self%id_rate5,  'rate5',  'mgDOM/m3/s', 'rate5 - leaching')
		  call self%register_diagnostic_variable(self%id_light,  'light',  'W/m2', 'light')
	  endif
	  
	  ! Register environmental dependencies
	  call self%register_dependency(self%id_temp, standard_variables%temperature)
	  call self%register_dependency(self%id_swf,  standard_variables%downwelling_shortwave_flux)
	  
	  call self%add_to_aggregate_variable(standard_variables%total_carbon, self%id_dom_a, scale_factor = self%mole_per_weight)
	  call self%add_to_aggregate_variable(standard_variables%total_carbon, self%id_dom_b, scale_factor = self%mole_per_weight)
	  
   end subroutine initialize

   subroutine do(self,_ARGUMENTS_DO_)
!
! !INPUT PARAMETERS:
   class (type_selmaprotbas_dom),intent(in) :: self
   _DECLARE_ARGUMENTS_DO_
!
! !LOCAL VARIABLES:
   real(rk)                   :: dom_a, dom_b 
   real(rk)                   :: temp, o2, nn, dd_c, swfz
   real(rk)                   :: temp_adj, rate_o2, rate_no3, R11, R12, R21, R22, R3a, R3b, R4a, R4b, r_leach
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops (if any)
   _LOOP_BEGIN_

   ! Retrieve current (local) state variable values
   _GET_(self%id_dom_a, dom_a) ! Labile DOM
   _GET_(self%id_dom_b, dom_b) ! Semi-labile DOM
   
   ! Retrieve current environmental conditions
   _GET_(self%id_temp, temp)    ! temperature
   _GET_(self%id_swf, swfz)             ! local shortwave radiation flux
   
   !retrieve domcast variables
   _GET_(self%id_o2, o2)
   _GET_(self%id_nn, nn)
   _GET_(self%id_dd_c, dd_c)
   
   temp_adj = self%theta ** (temp - 20.0_rk) ! temperature adjustment factor
   
   ! The unit of all rates is mgDOM/m3/s
   ! DOC mineralization by bacteria in the water column, controlled by O2
   rate_o2 = o2 / (self%km_o2 + o2) *  temp_adj
   R11 = self%k_om1 * rate_o2 * dom_a
   R12 = self%k_om2 * rate_o2 * dom_b  

   ! DOC mineralization by bacteria in the water column, controlled by NO3
   rate_no3 = nn / (self%km_no3 + nn) * self%k_inh_o2 / (self%k_inh_o2 + o2) * temp_adj
   R21 = self%k_om1 * rate_no3 * dom_a
   R22 = self%k_om2 * rate_no3 * dom_b
   
   ! Photo-oxidation and photo-mineralization
   R3a = self%oc_dom * self%qy_dom * self%f_par / self%e_par * swfz * dom_a
   R3b = self%oc_dom * self%qy_dom * self%f_par / self%e_par * swfz * dom_b
   ! Orig DOMCAST: R3 = self%oc_dom * self%qy_dom * self%f_par / self%e_par * swfz
   
   ! Flocculation
   R4a = self%k_floc * dom_a
   R4b = self%k_floc * dom_b ! DOMCAST had a "* swfz", but then the units would be wrong
   
   ! Source: leaching from detritus (POM). Only affects dom_a
   ! Similar to WET. Defaults of k_leach and q_leach are based on the WET defaults
   r_leach = dd_c * self%k_leach * exp(self%q_leach * (temp - 20.0_rk))
   
   ! Set light extinction
   _SET_DIAGNOSTIC_(self%id_extinc, self%ext_coef_a*dom_b**self%ext_coef_b)
   
   ! All processes degrade DOM pools
   _SET_ODE_(self%id_dom_a, -(R11 + R21 + R3a + R4a) + r_leach / self%mole_per_weight)
   _SET_ODE_(self%id_dom_b, -(R12 + R22 + R3b + R4b))
   
   ! O2 is consumed for DOM bacteria mineralization of both DOM pools
   _SET_ODE_(self%id_o2,  -(R11 + R12))

   ! Nitrate is consumed during DOM bacteria mineralization of both DOM pools
   _SET_ODE_(self%id_nn,  -(R21 + R22) * self%n_use_factor)
   
   ! Ammonium is created during bacteria mineralization
   _SET_ODE_(self%id_aa,  (R21 + R22) * self%n_use_factor)
   
   ! Carbon detritus is produced during flocculation
   ! Uses mole_per_weight ratio to convert from mg/m3 DOM to mmolC/m3 detritus
   _SET_ODE_(self%id_dd_c, (R4a + R4b) * self%mole_per_weight - r_leach)
   
   ! Export diagnostic variables -> rates included only for debugging purposes
   if(self%diagnostics) then
	   _SET_DIAGNOSTIC_(self%id_rate11, R11)
	   _SET_DIAGNOSTIC_(self%id_rate12, R12)
	   _SET_DIAGNOSTIC_(self%id_rate21, R21)
	   _SET_DIAGNOSTIC_(self%id_rate22, R22)
	   _SET_DIAGNOSTIC_(self%id_rate3a, R3a)
	   _SET_DIAGNOSTIC_(self%id_rate3b, R3b)
	   _SET_DIAGNOSTIC_(self%id_rate4a, R4a)
	   _SET_DIAGNOSTIC_(self%id_rate4b, R4b)
	   _SET_DIAGNOSTIC_(self%id_rate5, r_leach / self%mole_per_weight)
	   _SET_DIAGNOSTIC_(self%id_light, swfz)
   endif
   
   _LOOP_END_
   end subroutine do

end module
