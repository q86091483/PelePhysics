module react_type_module

  use network, only: nspec

  use amrex_fort_module, only : amrex_real
  implicit none

  type :: reaction_stat_t
     real(amrex_real) :: cost_value
     logical         :: reactions_succesful = .false.
  end type reaction_stat_t

  type :: react_t
     real(amrex_real) :: rhoedot_ext
     real(amrex_real), allocatable :: rhoY(:)
     real(amrex_real), allocatable :: rhoYdot_ext(:)
     real(amrex_real) :: rho
     real(amrex_real) :: T
     real(amrex_real) :: e
     integer :: i, j, k
  end type react_t
  
  interface build
     module procedure react_build
  end interface build

  interface destroy
     module procedure react_destroy
  end interface destroy


contains

  subroutine react_build(r)
    type(react_t), intent(inout) :: r
    allocate(r%rhoY(nspec))
    allocate(r%rhoYdot_ext(nspec))
  end subroutine react_build

  subroutine react_destroy(r)
    type(react_t), intent(inout) :: r
    deallocate(r%rhoY)
    deallocate(r%rhoYdot_ext)
  end subroutine react_destroy
  
end module react_type_module

