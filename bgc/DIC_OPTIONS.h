#ifndef DIC_OPTIONS_H
#define DIC_OPTIONS_H
#include "PACKAGES_CONFIG.h"
#include "CPP_OPTIONS.h"

#ifdef ALLOW_DIC
C     Package-specific Options & Macros go here

#define DIC_BIOTIC
#define ALLOW_FE
#define MINFE
#undef SEDIMENT_FE_ELROD
#undef SEDIMENT_FE_DALE
#undef HYDROTHERMAL_FE

C Prognostic Iron Binding Ligands
#undef ALLOW_VARIABLE_LIGANDS

#define ALLOW_O2
#undef READ_PAR
#define AD_SAFE
#define DIC_NO_NEG
#undef ALLOW_DIC_COST
C these all need to be defined for coupling to atmospheric model
#undef USE_QSW
#undef USE_QSW_UNDERICE
#undef USE_PLOAD

C use surface salinity forcing (scaled by mean surf value) for DIC & ALK forcing
C could also set PTRACERS_EvPrRn(n)=0., may need PTRACERS_ref in data.ptracers
#define ALLOW_OLD_VIRTUALFLUX

C put back bugs related to Water-Vapour in carbonate chemistry & air-sea fluxes
#undef WATERVAP_BUG

C dissolution only below saturation horizon following method by Karsten Friis
#undef CAR_DISS

C Make sure carbon dissociation constants are on the same pH scale
#define CARBONCHEM_TOTALPHSCALE

C Include self-shading effect by phytoplankton
#undef LIGHT_CHL

#endif /* ALLOW_DIC */
#endif /* DIC_OPTIONS_H */

CEH3 ;;; Local Variables: ***
CEH3 ;;; mode:fortran ***
CEH3 ;;; End: ***
