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
! !MODULE: selmaprotbas
!
! !INTERFACE:
   MODULE selmaprotbas_zooplankton
!
! !DESCRIPTION:
!
! !USE:
   use fabm_types
   use fabm_particle
   use fabm_expressions
   use fabm_builtin_models

   implicit none

   private
!
! !PUBLIC_DERIVED_TYPES:
  type,extends(type_particle_model),public :: type_selmaprotbas_zooplankton
      ! Variable identifiers
      type (type_model_id),      allocatable,dimension(:) :: id_prey
	  type (type_state_variable_id),allocatable,dimension(:) :: id_preyc
      type (type_dependency_id), allocatable,dimension(:) :: id_preyn,id_preyp,id_preys
	  type (type_state_variable_id)             :: id_c
      type (type_state_variable_id)             :: id_aa, id_po, id_dd_c, id_dd_p, id_dd_n, id_dd_si, id_dic, id_o2, id_si,id_dom_a
      type (type_dependency_id)                 :: id_temp

      ! Model parameters
      integer :: nprey
      real(rk), allocatable :: pref(:)
      real(rk) :: nue,sigma_b
      real(rk) :: iv,graz,toptz,zcl1
      real(rk) :: rfr,rfn,rfs
	  logical  :: couple_dom
	  real(rk) :: f_diss, mole_c_per_weight_dom
   contains
      procedure :: initialize
      procedure :: do
   end type
!EOP
!-----------------------------------------------------------------------

   CONTAINS

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialise the selma/zooplankton model
!
! !INTERFACE:
   subroutine initialize(self,configunit)
!
! !DESCRIPTION:
!   Here, the selma namelist is read and the variables exported by the model are registered with FABM
!
! !INPUT PARAMETERS:
   class(type_selmaprotbas_zooplankton),intent(inout),target :: self
   integer,               intent(in)           :: configunit
   real(rk),parameter :: secs_per_day = 86400._rk
   real(rk) :: c0
   integer :: iprey
   character(len=8) :: strprey

   call self%get_parameter(c0,       'c0',  'mmol C/m3',   'background concentration',    default=0._rk)
   call self%get_parameter(self%rfr, 'rfr', 'mol P/mol C', 'phosphorus : carbon ratio', default=1.0_rk/106.0_rk)
   call self%get_parameter(self%rfn, 'rfn', 'mol N/mol C', 'nitrogen : carbon ratio',     default=16.0_rk/106.0_rk)
   call self%get_parameter(self%rfs, 'rfs', 'mol Si/mol C', 'silicon : carbon ratio',     default=0.000_rk)
   call self%register_state_variable(self%id_c,'c','mmol C/m3', 'concentration', minimum=0.0_rk, background_value=c0)

   call self%get_parameter(self%nprey, 'nprey', '', 'number of prey', default=1) 
   allocate(self%id_prey(self%nprey))
   allocate(self%id_preyc(self%nprey))
   allocate(self%id_preyn(self%nprey))
   allocate(self%id_preyp(self%nprey))
   allocate(self%id_preys(self%nprey))
   allocate(self%pref(self%nprey))
   do iprey=1,self%nprey
      write (strprey,'(i0)') iprey
      call self%register_model_dependency(self%id_prey(iprey), 'prey' // trim(strprey))
	  call self%register_state_dependency(self%id_preyc(iprey),'preyc'//trim(strprey), 'mmol C/m3', 'prey '//trim(strprey))
	  call self%register_dependency(self%id_preyn(iprey),'preyn'//trim(strprey), 'mmol N/m3', 'prey '//trim(strprey))
	  call self%register_dependency(self%id_preyp(iprey),'preyp'//trim(strprey), 'mmol P/m3', 'prey '//trim(strprey))
	  call self%register_dependency(self%id_preys(iprey),'preys'//trim(strprey), 'mmol Si/m3', 'prey '//trim(strprey))
	  call self%request_coupling_to_model(self%id_preyc(iprey), self%id_prey(iprey),standard_variables%total_carbon)
	  call self%request_coupling_to_model(self%id_preyn(iprey), self%id_prey(iprey),standard_variables%total_nitrogen)
	  call self%request_coupling_to_model(self%id_preyp(iprey), self%id_prey(iprey),standard_variables%total_phosphorus)
	  call self%request_coupling_to_model(self%id_preys(iprey), self%id_prey(iprey),standard_variables%total_silicate)
	  
	  call self%get_parameter(self%pref(iprey), 'pref'//trim(strprey), '-', 'preference for prey '//trim(strprey), default=1.0_rk)
   end do

   call self%get_parameter(self%nue,     'nue',     'm3/d/mmol C',    'respiration rate',                default=0.001509_rk, scale_factor=1.0_rk/secs_per_day)
   call self%get_parameter(self%sigma_b, 'sigma_b', 'm3/d/mmol C',    'mortality rate',                  default=0.004528_rk, scale_factor=1.0_rk/secs_per_day)
   call self%get_parameter(self%iv,      'iv',      '1/(mmol C/m3)2', 'Ivlev constant, quadratic',       default=0.27341_rk)
   call self%get_parameter(self%graz,    'graz',    '1/d',            'grazing rate',                    default=0.5_rk, scale_factor=1.0_rk/secs_per_day)
   call self%get_parameter(self%toptz,   'toptz',   'deg C',          'optimal temperature for grazing', default=20._rk)
   call self%get_parameter(self%zcl1,    'zcl1',    '-',              'closure parameter',               default=50._rk)
   call self%get_parameter(self%couple_dom, 'couple_dom', '-',        'couple zooplankton to dom module',default=.false.)

   ! Register state variables
   call self%register_state_dependency(self%id_aa, 'aa', 'mmol N/m3', 'ammonium')
   call self%register_state_dependency(self%id_po, 'po', 'mmol P/m3', 'phosphate')
   call self%register_state_dependency(self%id_si, 'si', 'mmol Si/m3', 'silicon')
   call self%register_state_dependency(self%id_dd_c, 'dd_c', 'mmol C/m3', 'carbon detritus')
   call self%register_state_dependency(self%id_dd_p, 'dd_p', 'mmol P/m3', 'phosphorus detritus')
   call self%register_state_dependency(self%id_dd_n, 'dd_n', 'mmol N/m3', 'nitrogen detritus')
   call self%register_state_dependency(self%id_dd_si, 'dd_si', 'mmol Si/m3', 'silicon detritus')
   call self%register_state_dependency(self%id_o2, 'o2', 'mmol O2/m3','oxygen')
   
   ! DOM snippet
   call self%get_parameter(self%f_diss,'f_diss','-', 'fraction of zooplankton biomass that is dissolved, only used if DOM is coupled', default=0.2_rk, minimum=0.0_rk, maximum=1.0_rk)
   if (self%couple_dom) then
	  call self%register_state_dependency(self%id_dom_a, 'dom_a', 'mg/m3', 'DOM - labile')
	  call self%get_parameter(self%mole_c_per_weight_dom,'mole_c_per_weight_dom','molC/gDOM', 'mol C per g DOM', default=0.0416_rk) ! Default assumes 50% C/DW weight ratio and 12.01 g/mole molar mass
   else
      self%f_diss = 0.0_rk   ! if no DOM is provided, all biomass should go to detritus
   end if
   ! End DOM snippet

   ! Contribute to total nitrogen, phosphorus, carbon
   call self%add_to_aggregate_variable(standard_variables%total_nitrogen,   self%id_c, self%rfn)
   call self%add_to_aggregate_variable(standard_variables%total_phosphorus, self%id_c, self%rfr)
   call self%add_to_aggregate_variable(standard_variables%total_carbon,     self%id_c)
   call self%add_to_aggregate_variable(standard_variables%total_silicate, self%id_c, self%rfs)

   ! Register link to external DIC pool
   call self%register_state_dependency(self%id_dic,standard_variables%mole_concentration_of_dissolved_inorganic_carbon, required=.false.)

   ! Register environmental dependencies
   call self%register_dependency(self%id_temp, standard_variables%temperature)

   END subroutine initialize
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
   class(type_selmaprotbas_zooplankton), INTENT(IN) :: self
  _DECLARE_ARGUMENTS_DO_
!
! !LOCAL VARIABLES:
    real(rk) :: c, cg
    real(rk) :: prey,preyp,preyn,preys
    real(rk) :: temp, o2
    real(rk) :: food, zz0, food_eps, gg, lzd, lzn, graz_z, tlim
    real(rk) :: growth_red_p, growth_red_n, growth_red_si, growth_red_net
	real(rk) :: prey_rfn,prey_rfr,prey_rfs
    integer  :: iprey
    real(rk),parameter :: epsilon = 0.00000000001_rk
!EOP
!-----------------------------------------------------------------------
!BOC
    ! Enter spatial_loops (if any)
    _LOOP_BEGIN_

      ! Retrieve current (local) state variable values
      _GET_(self%id_c,c)                  ! own concentration
      _GET_WITH_BACKGROUND_(self%id_c,cg) ! own concentration including background

      food = 0
      do iprey=1,self%nprey
         _GET_(self%id_preyc(iprey),prey)
         _GET_(self%id_preyp(iprey),preyp)
		 _GET_(self%id_preyn(iprey),preyn)
		 _GET_(self%id_preys(iprey),preys)
		 
		 ! Calculate nutrient ratios of the prey
	     prey_rfr = preyp / max(prey, epsilon)
	     prey_rfn = preyn / max(prey, epsilon)
	     prey_rfs = preys / max(prey, epsilon)
		 
         ! BOSP
         ! Calculate if any of the ratios of the prey is lower than that of the zooplankton, which would result in reduced growth  
         ! P:C ratio
         if(prey_rfr < self%rfr) then
            growth_red_p = (self%rfr - prey_rfr)/self%rfr
         else
            growth_red_p = 0
         end if
         
         ! N:C ratio
         if(prey_rfn < self%rfn) then
            growth_red_n = (self%rfn - prey_rfn)/self%rfn
         else
            growth_red_n = 0
         end if
         
         ! Si:C ratio
         if(prey_rfs < self%rfs) then
            growth_red_si = (self%rfs - prey_rfs)/self%rfs
         else
            growth_red_si = 0
         end if
         
         growth_red_net = max(growth_red_p, growth_red_n, growth_red_si) ! Reduction of growth rate (fraction)
         
         ! EOSP
         
         food = food + self%pref(iprey) * prey * (1.000_rk - growth_red_net)
      end do

      ! Local environment
      _GET_(self%id_o2,o2)     ! oxygen (mmol O2/m3)
      _GET_(self%id_temp,temp) ! temperature (degrees Celsius)

      zz0 = self%zcl1*c*c

      food_eps = max(food, epsilon) ! Be sure food is positive
      gg = self%graz * (1.0_rk - exp(self%iv * food * food * (-1.0_rk))) !Grazing rate

      lzd = self%sigma_b ! Zooplankton mortality rate
      lzn = self%nue     ! Zooplankton respiration rate

      if (o2 <= 0.0_rk) then
         ! Anoxic conditions
         gg = 0.0_rk                ! No grazing
         lzd = 10._rk *self%sigma_b ! Higher mortality
         lzn = 0.0_rk               ! No respiration
      end if

      ! Zooplankton grazing depends on food availability and temperature
      tlim = (1.0_rk + 2.7183_rk/self%toptz/self%toptz * max(temp,0.0_rk)**2 * exp(1.0_rk - temp * 2.0_rk / self%toptz))
      graz_z = gg * cg/food_eps * tlim

      _SET_ODE_(self%id_o2, - lzn * zz0)
      _SET_ODE_(self%id_aa, self%rfn * lzn * zz0)
      _SET_ODE_(self%id_po, self%rfr * lzn * zz0)
	  _SET_ODE_(self%id_si, self%rfs * lzn * zz0)
      do iprey=1,self%nprey
         _GET_(self%id_preyc(iprey),prey)
         _GET_(self%id_preyp(iprey),preyp)
		 _GET_(self%id_preyn(iprey),preyn)
		 _GET_(self%id_preys(iprey),preys)
		 
		 prey_rfr = preyp / max(prey, epsilon)
	     prey_rfn = preyn / max(prey, epsilon)
	     prey_rfs = preys / max(prey, epsilon)
		 
         ! BOSP
         ! Calculate if any of the ratios of the prey is lower than that of the zooplankton, which would result in reduced growth  
         ! P:C ratio
         if(prey_rfr < self%rfr) then
            growth_red_p = (self%rfr - prey_rfr)/self%rfr
         else
            growth_red_p = 0.0_rk
         end if
         
         ! N:C ratio
         if(prey_rfn < self%rfn) then
            growth_red_n = (self%rfn - prey_rfn)/self%rfn
         else
            growth_red_n = 0.0_rk
         end if
         
         ! Si:C ratio
         if(prey_rfs < self%rfs) then
            growth_red_si = (self%rfs - prey_rfs)/self%rfs
         else
            growth_red_si = 0.0_rk
         end if
         
         growth_red_net = max(growth_red_p, growth_red_n, growth_red_si) ! Reduction of growth rate (fraction)
         
         ! Nutrients released as detritus from consumption of phytoplankton with different nutrient ratios
         ! If nutrient ratios of prey are lower than zooplankton's, prey is consumed at the same rate, but zooplankton grows less (since plankton nutrient ratios need to stay constant in the selmaprotbas model)
         ! Non-limiting nutrients of the phytoplankton that is consumed, but not used for growth, are released as detritus, according to the prey's nutrient ratios (1st term)
         ! If nutrient ratios of the prey are higher than zooplankton's, the excess nutrients in the amount of prey consumed, are also released as detritus, according to the difference of the prey's and zooplankton's nutrient ratios (2nd term)
         _SET_ODE_(self%id_dd_c, growth_red_net * graz_z * self%pref(iprey) * prey) ! carbon 
         _SET_ODE_(self%id_dd_p, (growth_red_net - growth_red_p) * self%pref(iprey) * graz_z * prey * prey_rfr + (1.000_rk - growth_red_net) * max(prey_rfr - self%rfr, 0.0_rk) * graz_z * prey * self%pref(iprey)) ! phosphorus
         _SET_ODE_(self%id_dd_n, (growth_red_net - growth_red_n) * self%pref(iprey) * graz_z * prey * prey_rfn + (1.000_rk - growth_red_net) * max(prey_rfn - self%rfn, 0.0_rk) * graz_z * prey * self%pref(iprey)) ! nitrogen
         _SET_ODE_(self%id_dd_si, (growth_red_net - growth_red_si) * self%pref(iprey) * graz_z * prey * prey_rfs + (1.000_rk - growth_red_net) * max(prey_rfs - self%rfs, 0.0_rk) * graz_z * prey * self%pref(iprey)) ! silicon
         
         ! EOSP
         
         _SET_ODE_(self%id_preyc(iprey), - graz_z * self%pref(iprey) * prey)
      end do
      _SET_ODE_(self%id_c, graz_z * food - (lzn + lzd) * zz0)
      _SET_ODE_(self%id_dd_c, + lzd * zz0 * (1.0_rk - self%f_diss))
      _SET_ODE_(self%id_dd_p, + lzd * zz0 * self%rfr)
      _SET_ODE_(self%id_dd_n, + lzd * zz0 * self%rfn)
      _SET_ODE_(self%id_dd_si, + lzd * zz0 * self%rfs)
	  
	  if (self%couple_dom) then
	     _SET_ODE_(self%id_dom_a, lzd * zz0 * self%f_diss / self%mole_c_per_weight_dom)
	  end if
	  
      if (_AVAILABLE_(self%id_dic)) _SET_ODE_(self%id_dic, lzn * zz0)

   ! Leave spatial loops (if any)
   _LOOP_END_

   END subroutine do
!EOC

!-----------------------------------------------------------------------

  END MODULE selmaprotbas_zooplankton

!-----------------------------------------------------------------------
