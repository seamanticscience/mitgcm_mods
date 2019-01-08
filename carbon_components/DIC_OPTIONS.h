C $Header: /u/gcmpack/MITgcm/pkg/dic/DIC_OPTIONS.h,v 1.14 2016/01/11 21:46:55 jmc Exp $
C $Name:  $

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
C could also set PTRACERS_EvPrRn(n)=0., 
C may need PTRACERS_ref(1,n)=??? in data.ptracers
#define AD_SAFE
#define DIC_NO_NEG
#undef ALLOW_DIC_COST
C these all need to be defined for coupling to
C atmospheric model
#undef USE_QSW
#undef USE_QSW_UNDERICE
#undef USE_PLOAD

C use surface salinity forcing (scaled by mean surf value) for DIC & ALK forcing
#undef ALLOW_OLD_VIRTUALFLUX

C put back bugs related to Water-Vapour in carbonate chemistry & air-sea fluxes
#undef WATERVAP_BUG

C define to include carbon component calculations
#define DIC_CARBON_COMPONENTS
#define DIC_COMPONENT_DIAGS
#define DIC_INSITU_COMPONENTS

C define to include C13 calculations
#undef ALLOW_C13

C dissolution only below saturation horizon following method by Karsten Friis
#undef CAR_DISS

C Include self-shading effect by phytoplankton
#undef LIGHT_CHL

#endif /* ALLOW_DIC */
#endif /* DIC_OPTIONS_H */

CEH3 ;;; Local Variables: ***
CEH3 ;;; mode:fortran ***
CEH3 ;;; End: ***
