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

   implicit none

   private
!
! !PUBLIC_DERIVED_TYPES:
  type,extends(type_base_model),public :: type_selmaprotbas_zooplankton
      ! Variable identifiers
      type (type_state_variable_id)             :: id_c
      type (type_state_variable_id),allocatable :: id_prey(:)
      type (type_state_variable_id)             :: id_aa, id_po, id_dd_c, id_dd_p, id_dd_n, id_dd_si, id_dic, id_o2, id_si
      type (type_dependency_id)                 :: id_temp

      ! Model parameters
      integer :: nprey
      real(rk), allocatable :: pref(:)
      real(rk), allocatable :: prey_rfr(:)
      real(rk), allocatable :: prey_rfn(:)
      real(rk), allocatable :: prey_rfs(:)
      real(rk) :: nue,sigma_b
      real(rk) :: iv,graz,toptz,zcl1
      real(rk) :: rfr,rfn, rfs
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
   call self%get_parameter(self%rfs, 'rfs', 'mol Si/mol C', 'silica : carbon ratio',     default=0.000_rk)
   call self%register_state_variable(self%id_c,'c','mmol C/m3', 'concentration', minimum=0.0_rk, background_value=c0)

   call self%get_parameter(self%nprey, 'nprey', '', 'number of prey', default=1) 
   allocate(self%id_prey(self%nprey))
   allocate(self%pref(self%nprey))
   allocate(self%prey_rfr(self%nprey))
   allocate(self%prey_rfn(self%nprey))
   allocate(self%prey_rfs(self%nprey))
   do iprey=1,self%nprey
      write (strprey,'(i0)') iprey
      call self%register_state_dependency(self%id_prey(iprey),'prey'//trim(strprey), 'mmol C/m3', 'prey '//trim(strprey))
      call self%get_parameter(self%pref(iprey), 'pref'//trim(strprey), '-', 'preference for prey '//trim(strprey), default=1.0_rk)
      call self%get_parameter(self%prey_rfr(iprey), 'prey_rfr'//trim(strprey), '-', 'P:C ratio for prey '//trim(strprey), default=1.0_rk/106.0_rk)
      call self%get_parameter(self%prey_rfn(iprey), 'prey_rfn'//trim(strprey), '-', 'N:C ratio for prey '//trim(strprey), default=16.0_rk/106.0_rk)
      call self%get_parameter(self%prey_rfs(iprey), 'prey_rfs'//trim(strprey), '-', 'Si:C ratio for prey '//trim(strprey), default=0.000_rk)
   end do

   call self%get_parameter(self%nue,     'nue',     'm3/d/mmol C',    'respiration rate',                default=0.001509_rk, scale_factor=1.0_rk/secs_per_day)
   call self%get_parameter(self%sigma_b, 'sigma_b', 'm3/d/mmol C',    'mortality rate',                  default=0.004528_rk, scale_factor=1.0_rk/secs_per_day)
   call self%get_parameter(self%iv,      'iv',      '1/(mmol C/m3)2', 'Ivlev constant, quadratic',       default=0.27341_rk)
   call self%get_parameter(self%graz,    'graz',    '1/d',            'grazing rate',                    default=0.5_rk, scale_factor=1.0_rk/secs_per_day)
   call self%get_parameter(self%toptz,   'toptz',   'deg C',          'optimal temperature for grazing', default=20._rk)
   call self%get_parameter(self%zcl1,    'zcl1',    '-',              'closure parameter',               default=50._rk)

   ! Register state variables
   call self%register_state_dependency(self%id_aa, 'aa', 'mmol N/m3', 'ammonium')
   call self%register_state_dependency(self%id_po, 'po', 'mmol P/m3', 'phosphate')
   call self%register_state_dependency(self%id_si, 'si', 'mmol Si/m3', 'silica')
   call self%register_state_dependency(self%id_dd_c, 'dd_c', 'mmol C/m3', 'carbon detritus')
   call self%register_state_dependency(self%id_dd_p, 'dd_p', 'mmol P/m3', 'phosphorus detritus')
   call self%register_state_dependency(self%id_dd_n, 'dd_n', 'mmol N/m3', 'nitrogen detritus')
   call self%register_state_dependency(self%id_dd_si, 'dd_si', 'mmol Si/m3', 'silica detritus')
   call self%register_state_dependency(self%id_o2, 'o2', 'mmol O2/m3','oxygen')

   ! Contribute to total nitrogen, phosphorus, carbon
   call self%add_to_aggregate_variable(standard_variables%total_nitrogen,   self%id_c, self%rfn)
   call self%add_to_aggregate_variable(standard_variables%total_phosphorus, self%id_c, self%rfr)
   call self%add_to_aggregate_variable(standard_variables%total_carbon,     self%id_c)
   call self%add_to_aggregate_variable(type_bulk_standard_variable(name='total_silica',units="mmol/m^3",aggregate_variable=.true.),self%id_c,scale_factor=self%rfs)
   

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
    real(rk) :: prey
    real(rk) :: temp, o2
    real(rk) :: food, zz0, food_eps, gg, lzd, lzn, graz_z, tlim
    real(rk) :: growth_red_p, growth_red_n, growth_red_si, growth_red_net
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
         _GET_(self%id_prey(iprey),prey)
         
         ! BOSP
         ! Calculate if any of the ratios of the prey is lower than that of the zooplankton, which would result in reduced growth  
         ! P:C ratio
         if(self%prey_rfr(iprey) < self%rfr) then
            growth_red_p = (self%prey_rfr(iprey) - self%rfr)/self%rfr
         else
            growth_red_p = 0
         end if
         
         ! N:C ratio
         if(self%prey_rfn(iprey) < self%rfn) then
            growth_red_n = (self%prey_rfn(iprey) - self%rfn)/self%rfn
         else
            growth_red_n = 0
         end if
         
         ! Si:C ratio
         if(self%prey_rfs(iprey) < self%rfs) then
            growth_red_si = (self%prey_rfs(iprey) - self%rfs)/self%rfs
         else
            growth_red_si = 0
         end if
         
         growth_red_net = min(growth_red_p, growth_red_n, growth_red_si) ! Reduction of growth rate (fraction)
         
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
         _GET_(self%id_prey(iprey),prey)
         
         ! BOSP
         ! Calculate if any of the ratios of the prey is lower than that of the zooplankton, which would result in reduced growth  
         ! P:C ratio
         if(self%prey_rfr(iprey) < self%rfr) then
            growth_red_p = (self%prey_rfr(iprey) - self%rfr)/self%rfr
         else
            growth_red_p = 0.0_rk
         end if
         
         ! N:C ratio
         if(self%prey_rfn(iprey) < self%rfn) then
            growth_red_n = (self%prey_rfn(iprey) - self%rfn)/self%rfn
         else
            growth_red_n = 0.0_rk
         end if
         
         ! Si:C ratio
         if(self%prey_rfs(iprey) < self%rfs) then
            growth_red_si = (self%prey_rfs(iprey) - self%rfs)/self%rfs
         else
            growth_red_si = 0.0_rk
         end if
         
         growth_red_net = min(growth_red_p, growth_red_n, growth_red_si) ! Reduction of growth rate (fraction)
         
         ! Nutrients released as detritus from consumption of phytoplankton with different nutrient ratios
         ! If nutrient ratios of prey are lower than zooplankton's, prey is consumed at the same rate, but zooplankton grows less (since plankton nutrient ratios need to stay constant in the selmaprotbas model)
         ! Non-limiting nutrients of the phytoplankton that is consumed, but not used for growth, are released as detritus, according to the prey's nutrient ratios (1st term)
         ! If nutrient ratios of the prey are higher than zooplankton's, the excess nutrients in the amount of prey consumed, are also released as detritus, according to the difference of the prey's and zooplankton's nutrient ratios (2nd term)
         _SET_ODE_(self%id_dd_c, growth_red_net * graz_z * self%pref(iprey) * prey) ! carbon 
         _SET_ODE_(self%id_dd_p, (growth_red_net - growth_red_p) * self%pref(iprey) * graz_z * prey * self%prey_rfr(iprey) + (1.000_rk - growth_red_net) * max(self%prey_rfr(iprey) - self%rfr, 0.0_rk) * graz_z * prey * self%pref(iprey)) ! phosphorus
         _SET_ODE_(self%id_dd_n, (growth_red_net - growth_red_n) * self%pref(iprey) * graz_z * prey * self%prey_rfn(iprey) + (1.000_rk - growth_red_net) * max(self%prey_rfn(iprey) - self%rfn, 0.0_rk) * graz_z * prey * self%pref(iprey)) ! nitrogen
         _SET_ODE_(self%id_dd_si, (growth_red_net - growth_red_si) * self%pref(iprey) * graz_z * prey * self%prey_rfs(iprey) + (1.000_rk - growth_red_net) * max(self%prey_rfs(iprey) - self%rfs, 0.0_rk) * graz_z * prey * self%pref(iprey)) ! silica
         
         ! EOSP
         
         _SET_ODE_(self%id_prey(iprey), - graz_z * self%pref(iprey) * prey)
      end do
      _SET_ODE_(self%id_c, graz_z * food - (lzn + lzd) * zz0)
      _SET_ODE_(self%id_dd_c, + lzd * zz0)
      _SET_ODE_(self%id_dd_p, + lzd * zz0 * self%rfr)
      _SET_ODE_(self%id_dd_n, + lzd * zz0 * self%rfn)
      _SET_ODE_(self%id_dd_si, + lzd * zz0 * self%rfs)

      if (_AVAILABLE_(self%id_dic)) _SET_ODE_(self%id_dic, lzn * zz0)

   ! Leave spatial loops (if any)
   _LOOP_END_

   END subroutine do
!EOC

!-----------------------------------------------------------------------

  END MODULE selmaprotbas_zooplankton

!-----------------------------------------------------------------------
