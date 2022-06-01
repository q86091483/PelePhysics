# AMReX
DIM        = 3
PRECISION  = DOUBLE
PROFILE    = FALSE
VERBOSE    = FALSE
DEBUG      = FALSE

# Compiler
COMP       = gnu
FCOMP      = gfortran
USE_MPI    = TRUE
USE_OMP    = FALSE
USE_CUDA   = TRUE
USE_HIP    = FALSE

# PelePhysics
FUEGO_GAS  = FALSE
TINY_PROFILE = TRUE

# define the location of the PELE_PHYSICS top directory
PELE_PHYSICS_HOME    ?= ../../..

# this flag activates the subcycling mode in the D/Cvode routines
DEFINES  += -DMOD_REACTOR

# Activates use of SUNDIALS
ifeq ($(USE_CUDA), TRUE)
  PELE_USE_KLU = FALSE
else
  ifeq ($(USE_HIP), TRUE)
    PELE_USE_KLU = FALSE
  else
    PELE_USE_KLU = FALSE
  endif
endif

ifeq ($(FUEGO_GAS), TRUE)
  Eos_Model       = Fuego
  Chemistry_Model = grimech30-noArN
  Transport_Model = Constant
else
  Eos_Model       = Manifold
  Chemistry_Model = Manifold3d
  Transport_Model = Manifold
endif

Bpack   := ./Make.package
Blocs   := .

include $(PELE_PHYSICS_HOME)/Testing/Exec/Make.PelePhysics
