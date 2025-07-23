#include "fabm_driver.h"

!-----------------------------------------------------------------------
!BOP
!
!                  MSI-ERGOM (Baltic Sea Ecosystem Model)                       
!                                                                            
!  Coded by Gennadi Lessin (PML) based upon code provided by Thomas Neumann
!  (IOW) and updates the older version of ERGOM included in GOTM-BIO.
!  This older version is still provided in FABM under the name gotm/ergom for           
!  historical/educational purposes.                                           
!                                                                            
!  A detailed description of the original model is given                     
!  in T. Neumann, W. Fennel and C. Kremp, 2002. Experimental simulations     
!  with an ecosystem model of the Baltic Sea: a nutrient load reduction      
!  experiment. Global Biogeochemical Cycles, 16 (3).                         
!  http://dx.doi.org/10.1029/2001GB001450.                                   
!                                                                            
!  The present version adds oxygen-dependent phosphorus dynamics     
!  between sediment and water and the effect of bio-resuspension, as         
!  described in T. Neumann and G. Schernewski, 2008. Eutrophication in the   
!  Baltic Sea and shifts in nitrogen fixation analyzed with a 3D ecosystem   
!  model, J. Mar. Sys., 74 (1.2), pp. 592-602. 
!
!  Revision history:
!  September 2015, by Dennis Trolle (AU):
!  Implemented a switch for choosing between fresh and marine (default) environments.
!  If "fresh" is selected, oxygen debt (negative O2) due to H2S production is disabled.
!  Added a sediment burial process, and a range of additional diagnostic variables to output, 
!  incl. chlorophyll a, oxygen and nutrients in mass concentration units. 
!  Updated yaml input file with new entries (e.g., sediment burial rate, and phytoplankton
!  carbon:chla ratios for individual phyto groups)
!  May 2016, by Dennis Trolle (AU):
!  Added option to switch on or off n-fixation by cyanobacteria
!  Added settling of diatoms to bottom sediments, where diatoms are converted to fluff once settled
!
! !INTERFACE:
   MODULE selmaprotbas_phytoplankton
!
! !DESCRIPTION:
!
! !USE:
   use fabm_types
   use fabm_expressions

   implicit none

   private
!
! !PUBLIC_DERIVED_TYPES:
  type,extends(type_base_model),public :: type_selmaprotbas_phytoplankton
      type (type_state_variable_id) :: id_c
      type (type_state_variable_id) :: id_aa,id_nn,id_po,id_o2,id_dd_c,id_dd_p,id_dd_n,id_dd_si,id_dic,id_si
      type (type_bottom_state_variable_id) :: id_fl_c,id_fl_p,id_fl_n,id_fl_si
      type (type_dependency_id) :: id_par, id_parmean
      type (type_dependency_id) :: id_temp
      type (type_horizontal_dependency_id) :: id_taub
      type (type_diagnostic_variable_id) :: id_chla
      type (type_diagnostic_variable_id) :: id_GPP
      type (type_diagnostic_variable_id) :: id_NPP

      real(rk) :: alpha_light, imin
      real(rk) :: alpha, alpha_n, alpha_p, alpha_si
      logical  :: nitrogen_fixation, use_24h_light, mult_llim_nutlim
      logical  :: buoyancy_regulation, buoy_temperature, buoy_nutrient
      real(rk) :: par_limit1, par_limit2, par_limit3, vert_vel1, vert_vel2, vert_vel3, vert_vel4
	    real(rk) :: buoy_temp_limit, vert_vel_temp
      real(rk) :: buoy_nutrient_limit, vert_vel_nutrient
      real(rk) :: rfr, rfn, rfs
      real(rk) :: r0
      real(rk) :: tll, beta, temp_opt, temp_sigma, temp_min
      integer :: tlim, llim
      real(rk) :: nb
      real(rk) :: deltao
      real(rk) :: Yc
      real(rk) :: sedrate
      real(rk) :: tau_crit
   contains
      procedure :: initialize
      procedure :: do
      procedure :: do_bottom
      procedure :: get_vertical_movement
   end type
!EOP
!-----------------------------------------------------------------------

   CONTAINS

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialise the selma/phytoplankton model
!
! !INTERFACE:
   subroutine initialize(self,configunit)
!
! !DESCRIPTION:
!   Here, the selma namelist is read and the variables exported by the model are registered with FABM
!
! !INPUT PARAMETERS:
   class(type_selmaprotbas_phytoplankton),intent(inout),target :: self
   integer,               intent(in)           :: configunit

   real(rk) :: c0, wz, kc
   real(rk),parameter :: secs_per_day = 86400._rk

   call self%get_parameter(c0,         'c0',    'mmol C/m3',   'background concentration',            default=0._rk)
   call self%get_parameter(self%par_limit1,   'par_limit1',    'W/m2',   'below this light intensity, buoyancy-regulating plankton change vertical velocity from vert_vel1 to vert_vel2', default=0._rk)
   call self%get_parameter(self%par_limit2,   'par_limit2',    'W/m2',   'below this light intensity, buoyancy-regulating plankton change vertical velocity to vert_vel3', default=self%par_limit1)
   call self%get_parameter(self%par_limit3,   'par_limit3',    'W/m2',   'below this light intensity, buoyancy-regulating plankton change vertical velocity to vert_vel4', default=self%par_limit2)
   call self%get_parameter(self%vert_vel1,   'vert_vel1',    'm/d',   'rate of vertical movement above par_limit1, if buoyancy_regulation = true', default=0._rk, scale_factor=1.0_rk/secs_per_day)
   call self%get_parameter(self%vert_vel2,   'vert_vel2',    'm/d',   'rate of vertical movement below par_limit1 and above par_limit2, if buoyancy_regulation = true', default=0._rk, scale_factor=1.0_rk/secs_per_day)
   call self%get_parameter(self%vert_vel3,   'vert_vel3',    'm/d',   'rate of vertical movement below par_limit2 and above par_limit3, if buoyancy_regulation = true', default=self%vert_vel2*secs_per_day, scale_factor=1.0_rk/secs_per_day)
   call self%get_parameter(self%vert_vel4,   'vert_vel4',    'm/d',   'rate of vertical movement below par_limit3, if buoyancy_regulation = true', default=self%vert_vel3*secs_per_day, scale_factor=1.0_rk/secs_per_day)
   call self%get_parameter(self%buoy_temp_limit,   'buoy_temp_limit',    'degrees C',   'if temperature falls below this temperature, vertical velocity is set to vert_vel_temp, if buoyancy_regulation = true', default=0._rk)
   call self%get_parameter(self%vert_vel_temp,   'vert_vel_temp',    'm/d',   'vertical velocity at temperatures below buoy_temp_limit, if buoyancy_regulation = true', default=self%vert_vel4*secs_per_day, scale_factor=1.0_rk/secs_per_day)
   call self%get_parameter(self%buoy_nutrient_limit,   'buoy_nutrient_limit',    '-',   'if nutrient limitation falls below this value (0-1), vertical velocity is set to vert_vel_nutrient, if buoyancy_regulation = true', default=0._rk)
   call self%get_parameter(self%vert_vel_nutrient,   'vert_vel_nutrient',    'm/d',   'vertical velocity at nutrient limitation below buoy_nutrient_limit, if buoyancy_regulation = true', default=self%vert_vel4*secs_per_day, scale_factor=1.0_rk/secs_per_day)
   call self%get_parameter(self%rfr,   'rfr',   'mol P/mol C', 'phosphorus : carbon ratio in phytoplankton',         default=1.0_rk/106.0_rk)
   call self%get_parameter(self%rfn,   'rfn',   'mol N/mol C', 'nitrogen : carbon ratio in phytoplankton',             default=16.0_rk/106.0_rk)
   call self%get_parameter(self%rfs,   'rfs',   'mol Si/mol C', 'silica : carbon ratio in phytoplankton',            default=0.000_rk)
   call self%get_parameter(self%alpha, 'alpha', 'mmol C/m3',   'half-saturation for nutrient uptake', default=1.65625_rk)
   call self%get_parameter(self%alpha_n, 'alpha_n', 'mmol C/m3',   'half-saturation for nitrogen uptake', default=self%alpha)
   call self%get_parameter(self%alpha_p, 'alpha_p', 'mmol C/m3',   'half-saturation for phosphorus uptake', default=self%alpha)
   call self%get_parameter(self%alpha_si, 'alpha_si', 'mmol C/m3',   'half-saturation for silica uptake', default=self%alpha)
   call self%get_parameter(self%r0,    'r0',    '1/d',         'maximum growth rate at 20 degrees C',  default=1.3_rk, scale_factor=1.0_rk/secs_per_day)
   call self%get_parameter(self%nitrogen_fixation,    'nitrogen_fixation', '', 'whether nitrogen fixation is used to acquire nitrogen', default=.false.)
   call self%get_parameter(self%buoyancy_regulation,    'buoyancy_regulation', '', 'whether cells can regulate vertical movement', default=.false.)
   call self%get_parameter(self%buoy_temperature,    'buoy_temperature', '', 'whether temperature can regulate buoyancy, if buoyancy_regulation is true', default=.false.)
   call self%get_parameter(self%buoy_nutrient,    'buoy_nutrient', '', 'whether nutrient limitation can regulate buoyancy, if buoyancy_regulation is true', default=.false.)
   call self%get_parameter(self%tlim,  'tlim',  '', 'temperature limitation of growth (0: none, 1: flagellate-style, 2: cyanobacteria-style, 3:PROTECH-style, 4:optimal-stdev, 5: Lehman-equation)', default=0)
   select case (self%tlim)
   case (1)
      call self%get_parameter(self%tll, 'tll', 'degrees C^2', 'half-saturation temperature, squared', default=100.0_rk)
   case (2)
      call self%get_parameter(self%tll, 'tll', 'degrees C', 'lower temperature limit', default=13.5_rk)
   case (3)
      call self%get_parameter(self%beta, 'beta', '', 'temperature growth correction factor', default=3.7_rk)
   case (4)
      call self%get_parameter(self%temp_opt, 'temp_opt', 'degrees C', 'optimum temperature of growth rate', default=20.0_rk)
	  call self%get_parameter(self%temp_sigma, 'temp_sigma', 'degrees C', 'growth rate temperature constant; sigma in Gaussian curve', default=11.6_rk)
   case (5)
      call self%get_parameter(self%temp_opt, 'temp_opt', 'degrees C', 'optimum temperature of growth rate', default=20.0_rk)
	  call self%get_parameter(self%temp_min, 'temp_min', 'degrees C', 'temperature below optimum where growth is 10% of maximum', default=0.0_rk)
   end select
   call self%get_parameter(self%llim,  'llim',  '', 'light limitation of growth (1: Reynolds, 2: Selma)', default=1)
   select case (self%llim)
   case (1)
      call self%get_parameter(self%alpha_light, 'alpha_light', 'd-1 [W/m2]-1', 'the slope of light-dependent growth', default=0.1_rk,scale_factor=1.0_rk/secs_per_day)
   case (2)
      call self%get_parameter(self%imin, 'imin', 'W/m2', 'minimal optimal light radiation', default=50._rk)
   end select
   call self%get_parameter(self%mult_llim_nutlim, 'mult_llim_nutlim', '-', 'multiply light and nutrient limitation', default=.false.)
   call self%get_parameter(self%nb,      'nb',      '1/d', 'excretion rate', default=0.01_rk, scale_factor=1.0_rk/secs_per_day)
   call self%get_parameter(self%deltao,  'deltao',  '1/d', 'mortality rate', default=0.02_rk, scale_factor=1.0_rk/secs_per_day)
   call self%get_parameter(self%Yc,      'Yc',      'mmol C/mg Chl a', 'carbon : chlorophyll a ratio', default=6.25_rk)
   call self%get_parameter(wz,           'wz',      'm/d',  'vertical velocity (positive: upwards/floating, negative: downwards/sinking)', default=0.0_rk, scale_factor=1.0_rk/secs_per_day)
   call self%get_parameter(kc,           'kc',      'm2/mmol C', 'specific light attenuation')
   !call self%get_parameter(self%sll,'sll','PSU', 'lower salinity limit', default=1.0_rk)
   !call self%get_parameter(self%sul,'sul','PSU', 'upper salinity limit', default=10.0_rk)
   call self%get_parameter(self%sedrate, 'sedrate', 'm/d', 'sedimentation rate', default=0.0_rk, scale_factor=1.0_rk/secs_per_day)
   call self%get_parameter(self%use_24h_light, 'use_24h_light', '-', 'use 24h averaged light for growth', default=.false.)

   call self%register_state_variable(self%id_c, 'c', 'mmol C/m3', 'concentration', minimum=0.0_rk, background_value=c0, vertical_movement=wz)
   call self%register_state_dependency(self%id_aa, 'aa', 'mmol N/m3', 'ammonium')
   call self%register_state_dependency(self%id_nn, 'nn', 'mmol N/m3', 'nitrate')
   call self%register_state_dependency(self%id_o2, 'o2', 'mmol O2/m3','oxygen')
   call self%register_state_dependency(self%id_po, 'po', 'mmol P/m3', 'phosphate')
   call self%register_state_dependency(self%id_si, 'si', 'mmol Si/m3', 'silica') ! there need to be the coupling in fabm.yaml with the selmaprotbas.F90 script
   call self%register_state_dependency(self%id_dd_c, 'dd_c', 'mmol C/m3', 'carbon detritus')
   call self%register_state_dependency(self%id_dd_p, 'dd_p', 'mmol P/m3', 'phosphorus detritus')
   call self%register_state_dependency(self%id_dd_n, 'dd_n', 'mmol N/m3', 'nitrogen detritus')
   call self%register_state_dependency(self%id_dd_si, 'dd_si', 'mmol Si/m3', 'silica detritus')
   call self%register_dependency(self%id_temp, standard_variables%temperature)
   call self%register_dependency(self%id_par,  standard_variables%downwelling_photosynthetic_radiative_flux)
   if (self%use_24h_light) then
	  call self%register_dependency(self%id_parmean, temporal_mean(self%id_par,period=86400._rk, resolution=3600._rk))
   end if
   
   if (self%sedrate>0.0_rk) then
      call self%get_parameter(self%tau_crit,'tau_crit','N/m2', 'critical shear stress', default=0.07_rk)
      call self%register_state_dependency(self%id_fl_c, 'fl_c', 'mmol C/m2', 'fluff')
      call self%register_state_dependency(self%id_fl_p, 'fl_p', 'mmol P/m2', 'phosphorus fluff')
      call self%register_state_dependency(self%id_fl_n, 'fl_n', 'mmol N/m2', 'nitrogen fluff')
      call self%register_state_dependency(self%id_fl_si, 'fl_si', 'mmol Si/m2', 'silica fluff')
      call self%register_dependency(self%id_taub, standard_variables%bottom_stress)
   end if
   call self%register_state_dependency(self%id_dic,standard_variables%mole_concentration_of_dissolved_inorganic_carbon, required=.false.)

   call self%add_to_aggregate_variable(standard_variables%total_nitrogen,   self%id_c, self%rfn)
   call self%add_to_aggregate_variable(standard_variables%total_phosphorus, self%id_c, self%rfr)
   call self%add_to_aggregate_variable(standard_variables%total_carbon,     self%id_c)
   call self%add_to_aggregate_variable(type_bulk_standard_variable(name='total_silica',units="mmol/m^3",aggregate_variable=.true.),self%id_c,scale_factor=self%rfs)
   call self%add_to_aggregate_variable(standard_variables%attenuation_coefficient_of_photosynthetic_radiative_flux, self%id_c, kc, include_background=.true.)

   call self%register_diagnostic_variable(self%id_chla, 'chla','mg chl a/m3', 'chlorophyll concentration')
   call self%register_diagnostic_variable(self%id_GPP,  'GPP', 'mmol/m3/d',   'gross primary production')
   call self%register_diagnostic_variable(self%id_NPP,  'NPP', 'mmol/m3/d',   'net primary production')

!  we create an aggregate variable for chlA
   call self%add_to_aggregate_variable(type_bulk_standard_variable(name='total_chlorophyll',units="mg/m^3",aggregate_variable=.true.),self%id_chla,scale_factor=1._rk)
   
   end subroutine initialize
!EOC


!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Right hand sides of selma model
!
! !INTERFACE:
   subroutine do(self,_ARGUMENTS_DO_)
!
! !INPUT PARAMETERS:
   class(type_selmaprotbas_phytoplankton), INTENT(IN) :: self
  _DECLARE_ARGUMENTS_DO_
!
! !LOCAL VARIABLES:
   real(rk) :: par, temp
   real(rk) :: r0_temp
   real(rk) :: c, cg
   real(rk) :: aa, nn, po, o2, si
   real(rk) :: lightlim, iopt
   real(rk) :: ntemp, ptemp, tempq, sitemp
   real(rk) :: nlim, plim, r, silim
   real(rk) :: lpn, lpd
   real(rk),parameter :: epsilon = 0.00000000001_rk
   real(rk),parameter :: secs_per_day = 86400._rk
   real(rk),parameter :: si_minimal = 0.00000000001_rk

   ! Enter spatial_loops (if any)
   _LOOP_BEGIN_

      _GET_(self%id_c,c)                  ! own concentration
      _GET_WITH_BACKGROUND_(self%id_c,cg) ! own concentration including background

      _GET_(self%id_aa,aa) ! ammonium
      _GET_(self%id_nn,nn) ! nitrate
      _GET_(self%id_po,po) ! phosphate
      _GET_(self%id_si,si) ! silica
      _GET_(self%id_o2,o2) ! oxygen

      if (self%use_24h_light) then
		_GET_(self%id_parmean, par)   ! 24-hour averaged local photosynthetically active radiation (W/m2)
	  else
		_GET_(self%id_par, par)   ! local photosynthetically active radiation (W/m2)
	  end if
	  
      _GET_(self%id_temp,temp) ! temperature (degrees Celsius)
	  
	  ! BOSP
	  ! Temperature-modification of growth rate, based on tlim parameter
	  ! 0: none, 1: flagellate-style, 2: cyanobacteria-style, 3:PROTECH-style, 4:optimal-stdev, 5: Lehman-equation
	  if (self%tlim == 0) then
		 r0_temp = self%r0
	  elseif (self%tlim == 1) then
         ! Flagellate-style: Type III [sigmoidal] functional response ranging between 1 and 2
		 ! Caveat: this function makes r0 the maximum growth rate at 0 degrees C and increases it by up to a factor 2
         tempq = max(temp,0.0_rk)**2
         r0_temp = self%r0 * (1.0_rk + tempq / (self%tll + tempq))
      elseif (self%tlim == 2) then
         ! Cyanobacteria-style: 0 at infinitely low temp, 0.5 at tll, 1.0 at infinitely high temp
         r0_temp = self%r0  / (1.0_rk + exp(self%tll - temp))
	  elseif (self%tlim == 3) then
	     ! Exponential change with temperature. 'beta' can be based on cell size aspects. Source: PROTBAS model
		 r0_temp = 10**(log10(self%r0) + self%beta * (1000._rk / 293._rk - 1000._rk / (273._rk + temp)))
	  elseif (self%tlim == 4) then
	     ! Optimal temperature with gaussian distribution. self%r0 is always the growth at 20 degrees. Source: WET (Water Ecosystems Tool) model.
		 r0_temp = self%r0 * exp(-0.5_rk / self%temp_sigma**2 * ((temp - self%temp_opt)**2 - (20.0_rk - self%temp_opt)**2))
	  elseif (self%tlim == 5) then
	     ! Similar to tlim == 4 (WET option), but in this equation, r0 is the growth at Topt (so not necessarily 20 degrees) and use of "t_min" (temp where growth is 10% of optimal) rather than "temp_sigma"
		 ! Source: Lehman et al. 1975 (doi:10.4319/lo.1975.20.3.0343)
		 r0_temp = self%r0 * exp(-2.3_rk * ((self%temp_opt - temp) / (self%temp_opt - self%temp_min))**2)
      end if
	  
	  ! Light limitation
	  if (self%llim == 1) then
		  if (par > r0_temp/self%alpha_light) then
			  lightlim = 1.0_rk
		  else
		     lightlim = par * self%alpha_light / r0_temp
		  end if
     elseif (self%llim == 2) then
        ! As in the original SELMA code
		  iopt = max(0.5_rk * par, self%imin)
        lightlim = par / iopt * exp(1.0_rk - par / iopt)
     end if
	  
	  ! EOSP
      
      ! Nitrogen limitation (type III functional response)
      if (self%nitrogen_fixation) then
         nlim = 1.0_rk
      else
         ntemp = (nn + aa)**2
         nlim = ntemp / (self%alpha_n * self%alpha_n * self%rfn * self%rfn + ntemp) ! MiMe eq. for IN
      end if

      ! Phosphorus limitation (type III functional response)
      ptemp = po**2
      plim = ptemp / (self%alpha_p * self%alpha_p * self%rfr * self%rfr + ptemp)
      
      ! Silica limitation
      ! A small fraction (si_minimal) was added to prevent division by zero in case of a Si concentration of zero and rfs=0.
      ! Note: if si is zero, there is now strong si limitation of phytoplankton groups with rfs>0
      
      sitemp = si**2
      silim = (sitemp + si_minimal) / (self%alpha_si * self%alpha_si * self%rfs * self%rfs + sitemp + si_minimal)
      
      ! Calculation actual growth rate
      if (self%mult_llim_nutlim) then
	     r = lightlim * min(nlim, plim, silim) * r0_temp
	  else
	     r = min(lightlim, nlim, plim, silim) * r0_temp
	  end if
      
      lpn = self%nb     ! excretion rate
      lpd = self%deltao ! mortality rate

      if (o2 <= 0.0_rk) then
         ! Anoxic: no growth or respiration, higher mortality
         r = 0.0_rk
         lpn = 0.0_rk
         lpd = 10._rk * self%deltao
      end if

      if (self%nitrogen_fixation) then
         ! Nitrogen acquired from dinitrogen gas (not tracked)
         _SET_ODE_(self%id_o2, + r * cg)
      else
         ! Nitrogen acquired from ammonium and nitrate pools (proportional to availability)
         _SET_ODE_(self%id_aa, - r * cg * self%rfn * aa/(nn + aa + epsilon))
         _SET_ODE_(self%id_nn, - r * cg * self%rfn * (1.0_rk - aa/(nn + aa + epsilon)))
         _SET_ODE_(self%id_o2,   r * cg * (nn * 1.302_rk + aa)/(nn + aa + epsilon))
      end if

      _SET_ODE_(self%id_o2, - lpn * c)
      _SET_ODE_(self%id_aa, + lpn * c * self%rfn)
      _SET_ODE_(self%id_po, self%rfr * (- r * cg + lpn * c))
      _SET_ODE_(self%id_si, self%rfs * (- r * cg + lpn * c))
      _SET_ODE_(self%id_c, r * cg - (lpn + lpd) * c)
      _SET_ODE_(self%id_dd_c, lpd * c)
      _SET_ODE_(self%id_dd_p, lpd * c * self%rfr)
      _SET_ODE_(self%id_dd_n, lpd * c * self%rfn)
      _SET_ODE_(self%id_dd_si, lpd * c * self%rfs)

      if (_AVAILABLE_(self%id_dic)) _SET_ODE_(self%id_dic, lpn * c - r * cg)

      _SET_DIAGNOSTIC_(self%id_chla, c/self%Yc) ! old comment before making script N-based: relation between carbon and nitrogen from Hecky et al 1993. The stoichiometry of carbon, nitrogen, and phosphorus in particulate matter of lakes and oceans. Limnology and Oceanography, 38: 709-724.
      _SET_DIAGNOSTIC_(self%id_GPP, secs_per_day * r * cg)
      _SET_DIAGNOSTIC_(self%id_NPP, secs_per_day *(r * cg - lpn * c))

   ! Leave spatial loops (if any)
   _LOOP_END_

   END subroutine do
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Right hand sides of benthic_predator model
!
! !INTERFACE:
   subroutine do_bottom(self,_ARGUMENTS_DO_BOTTOM_)
!
! !INPUT PARAMETERS:
   class (type_selmaprotbas_phytoplankton),intent(in) :: self
   _DECLARE_ARGUMENTS_DO_BOTTOM_
!
! !LOCAL VARIABLES:
   real(rk)                   :: c
   real(rk)                   :: taub
   real(rk)                   :: ll
!
!EOP
!-----------------------------------------------------------------------
!BOC
   if (self%sedrate == 0.0_rk) return

   ! Enter spatial loops over the horizontal domain (if any).
   _HORIZONTAL_LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
 !   if (self%fluff) then
   _GET_(self%id_c,c)
   _GET_HORIZONTAL_(self%id_taub,taub)

   ! Resuspension-sedimentation rate are computed as in GOTM-BIO
   ! Phytoplankton is assumed to become detritus/fluff as soon as it settles to bottom sediments,
   ! and can therefore not be resuspended from benthic layer

   ll=self%sedrate*max(0.0_rk, self%tau_crit-taub)/self%tau_crit

   ! Sediment resuspension, detritus settling, diatom settling, bio-resuspension, mineralization and burial
   _SET_BOTTOM_ODE_(self%id_fl_c,+ ll * c)
   _SET_BOTTOM_ODE_(self%id_fl_p,+ ll * c * self%rfr)
   _SET_BOTTOM_ODE_(self%id_fl_n,+ ll * c * self%rfn)
   _SET_BOTTOM_ODE_(self%id_fl_si,+ ll * c * self%rfs)
   _SET_BOTTOM_EXCHANGE_(self%id_c, -ll * c)

   ! Leave spatial loops over the horizontal domain (if any).
   _HORIZONTAL_LOOP_END_

   end subroutine do_bottom
!EOC

!BOSP (Beginning of SelmaProtbas code)
   subroutine get_vertical_movement(self,_ARGUMENTS_GET_VERTICAL_MOVEMENT_)
      class (type_selmaprotbas_phytoplankton),intent(in) :: self
      
      _DECLARE_ARGUMENTS_GET_VERTICAL_MOVEMENT_
      real(rk) :: par, temp
      real(rk) :: nlim, ntemp, plim, ptemp, silim, sitemp, nutlim
      real(rk) :: aa, nn, po, si
      real(rk),parameter :: si_minimal = 0.00000000001_rk
      
      if (.not.self%buoyancy_regulation) return
      
      ! spatial loop
	  
      _LOOP_BEGIN_
	     ! determine vertical velocity based on light level
         _GET_(self%id_par,par)
         if (par >= self%par_limit1) then
            _SET_VERTICAL_MOVEMENT_(self%id_c, self%vert_vel1)
         else if (par >= self%par_limit2) then
            _SET_VERTICAL_MOVEMENT_(self%id_c, self%vert_vel2)
         else if (par >= self%par_limit3) then
		    _SET_VERTICAL_MOVEMENT_(self%id_c, self%vert_vel3)
         else
		    _SET_VERTICAL_MOVEMENT_(self%id_c, self%vert_vel4)
         end if
         
		 ! overwrite vertical velocity if nutrient depletion can regulate buoyancy
       if(self%buoy_nutrient) then
          _GET_(self%id_nn,nn)
          _GET_(self%id_aa,aa)
          _GET_(self%id_po,po)
          _GET_(self%id_si,si)
		    
          !! recalculate nutrient limitation (copy code from earlier in this script)
          ! Nitrogen limitation (type III functional response)
          if (self%nitrogen_fixation) then
             nlim = 1.0_rk
          else
             ntemp = (nn + aa)**2
             nlim = ntemp / (self%alpha_n * self%alpha_n * self%rfn * self%rfn + ntemp) ! MiMe eq. for IN
          end if
          
          ! Phosphorus limitation (type III functional response)
          ptemp = po**2
          plim = ptemp / (self%alpha_p * self%alpha_p * self%rfr * self%rfr + ptemp)
         
          ! Silica limitation
          ! A small fraction (si_minimal) was added to prevent division by zero in case of a Si concentration of zero and rfs=0.
          ! Note: if si is zero, there is now strong si limitation of phytoplankton groups with rfs>0
         
          sitemp = si**2
          silim = (sitemp + si_minimal) / (self%alpha_si * self%alpha_si * self%rfs * self%rfs + sitemp + si_minimal)
         
          ! Calculation nutrient limitation
          nutlim = min(nlim, plim, silim)
          
          if(nutlim < self%buoy_nutrient_limit) then
             _SET_VERTICAL_MOVEMENT_(self%id_c, self%vert_vel_nutrient)
          end if
       end if
		 
		 ! overwrite vertical velocity if temperature can regulate buoyancy
		 if(self%buoy_temperature) then
		    _GET_(self%id_temp,temp)
          if(temp < self%buoy_temp_limit) then
             _SET_VERTICAL_MOVEMENT_(self%id_c, self%vert_vel_temp)
          end if
       end if
		 
      _LOOP_END_
   end subroutine get_vertical_movement
!EOSP
!-----------------------------------------------------------------------

  END MODULE selmaprotbas_phytoplankton

!-----------------------------------------------------------------------
