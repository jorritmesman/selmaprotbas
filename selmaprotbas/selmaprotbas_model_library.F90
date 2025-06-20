module selmaprotbas_model_library

   use fabm_types, only: type_base_model_factory,type_base_model

   use selmaprotbas
   use selmaprotbas_phytoplankton
   use selmaprotbas_zooplankton
   use selmaprotbas_dom

   implicit none

   private

   type,extends(type_base_model_factory) :: type_factory
      contains
!KB      procedure :: initialize
      procedure :: create
   end type

   type (type_factory),save,target,public :: selmaprotbas_model_factory

contains

   subroutine create(self,name,model)
      class (type_factory),intent(in) :: self
      character(*),        intent(in) :: name
      class (type_base_model),pointer :: model

      select case (name)
         case ('selmaprotbas'); allocate(type_selmaprotbas::model)
         case ('phytoplankton'); allocate(type_selmaprotbas_phytoplankton::model)
         case ('zooplankton'); allocate(type_selmaprotbas_zooplankton::model)
		 case ('dom'); allocate(type_selmaprotbas_dom::model)
         case default
            call self%type_base_model_factory%create(name,model)
      end select

   end subroutine

end module
