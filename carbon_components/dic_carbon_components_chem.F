C $Header: /u/gcmpack/MITgcm/pkg/dic/carbon_chem.F,v 1.23 2011/10/07 21:36:39 dfer Exp $
C $Name:  $

#include "DIC_OPTIONS.h"

C--  File carbon_chem.F:
C--   Contents
C--   o CALC_PCO2_APPROX
C--   o CARBON_COEFFS
C--   o CARBON_COEFFS_PRESSURE_DEP

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
#ifdef DIC_CARBON_COMPONENTS
C !ROUTINE: CALC_CSAT_APPROX

C !INTERFACE: ==========================================================
       SUBROUTINE CALC_CSAT_APPROX(
     I                       t,s,pco2eq,pt,sit,ta,
     I                       k1local,k2local,
     I                       k1plocal,k2plocal,k3plocal,
     I                       kslocal,kblocal,kwlocal,
     I                       ksilocal,kflocal,
     I                       k0local, fugflocal, fflocal,
     I                       btlocal,stlocal,ftlocal,
     U                       pHlocal,csatloc,co3local,
     I                       i,j,k,bi,bj,myIter,myThid )

C !DESCRIPTION:
C     *==========================================================*
C     | SUBROUTINE CALC_CSAT_APPROX                              |
C     *==========================================================*
C      CALCULATE Csat, EQUILIBRIUM DIC FOR GIVEN pCO2, T, Alk etc
C !USES: ========================================================*
      IMPLICIT NONE

C     == GLobal variables ==
#include "SIZE.h"
#include "DYNVARS.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "FFIELDS.h"
#include "DIC_VARS.h"

C     == Routine arguments ==
C       diclocal = total inorganic carbon (mol/m^3)
C             where 1 T = 1 metric ton = 1000 kg
C       ta  = total alkalinity (eq/m^3)
C       pt  = inorganic phosphate (mol/^3)
C       sit = inorganic silicate (mol/^3)
C       t   = temperature (degrees C)
C       s   = salinity (PSU)
        _RL  t, s, pt, sit, ta
        _RL  csatloc, pco2eq, pHlocal
        _RL  fflocal, btlocal, stlocal, ftlocal
        _RL  k1local, k2local
        _RL  k1plocal, k2plocal, k3plocal
        _RL  kslocal, kblocal, kwlocal, ksilocal, kflocal
        _RL  k0local, fugflocal
        _RL  co3local
        INTEGER i,j,k,bi,bj,myIter
        INTEGER myThid
CEOP

C     == Local variables ==
        _RL  phguess
        _RL  cag
        _RL  bohg
        _RL  hguess
        _RL  denom
        _RL  stuff, stuff2, stuff3
        _RL  hnew
        _RL  h3po4g, h2po4g, hpo4g, po4g
        _RL  siooh3g

c ---------------------------------------------------------------------
C Change units from the input of mol/m^3 -> mol/kg:
c (1 mol/m^3)  x (1 m^3/1024.5 kg)
c where the ocean mean surface density is 1024.5 kg/m^3
c Note: mol/kg are actually what the body of this routine uses
c for calculations.  Units are reconverted back to mol/m^3 at the
c end of this routine.
c To convert input in mol/m^3 -> mol/kg

        pt=pt*permil
        sit=sit*permil
        ta=ta*permil
c ---------------------------------------------------------------------
c set first guess and brackets for [H+] solvers
c first guess (for newton-raphson)
C             WRITE(standardMessageUnit,*) 
C     &      'CALC_CSAT_APPROX: ph guess'
        phguess = phlocal
C     - new approx method
C     - make estimate of htotal (hydrogen ion conc) using
C       appromate estimate of CA, carbonate alkalinity
        hguess = 10.0 _d 0**(-phguess)
        
C     - first estimate borate contribution using guess for [H+]
        bohg = (btlocal*kblocal)/(hguess+kblocal)

C     - first estimate of contribution from phosphate
C     based on Dickson and Goyet
        denom  = (hguess*hguess*hguess)
     &           + (k1plocal*hguess*hguess)
     &           + (k1plocal*k2plocal*hguess)
     &           + (k1plocal*k2plocal*k3plocal)
        h3po4g = (pt*hguess*hguess*hguess) / denom
        h2po4g = (pt*k1plocal*hguess*hguess) / denom
        hpo4g  = (pt*k1plocal*k2plocal*hguess) / denom
        po4g   = (pt*k1plocal*k2plocal*k3plocal) / denom

C     - estimate contribution from silicate
C     based on Dickson and Goyet
        siooh3g = (sit*ksilocal) / (ksilocal + hguess)

C     - now estimate carbonate alkalinity
        cag = ta - bohg - (kwlocal/hguess) + hguess
     &           - hpo4g - 2.0 _d 0*po4g + h3po4g
     &           - siooh3g

C     - estimate hydrogen ion conc
C        stuff = (k1local*fugflocal*k0local*pco2eq)/cag
        stuff = (k1local*fflocal*pco2eq)/cag  
        stuff2 = stuff*(stuff + 8.0 _d 0*k2local)  
        stuff3 = sqrt(stuff2)  
        hnew = 0.5 _d 0*(stuff + stuff3)  
        hguess = hnew  

C evaluate csat, equilibrium DIC concentration
        csatloc = (pco2eq*fflocal)
     &   *(1.0 _d 0 + (k1local/hnew) + ((k1local*k2local)/(hnew*hnew)))

C Output pt, sit, ta and csat in mol/m^3
         pt=pt/permil 
         sit=sit/permil 
         ta=ta/permil 
         csatloc=csatloc/permil  

C calc final pH
         phlocal = -log10(hnew)  

        RETURN
        END
        
C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

CBOP
C !ROUTINE: CARBON_COEFFS

C !INTERFACE: ==========================================================
      SUBROUTINE CARBON_COEFFS_SURF(
     I                   bi,bj,iMin,iMax,jMin,jMax,myThid)

C !DESCRIPTION:
C     *==========================================================*
C     | SUBROUTINE CARBON_COEFFS_SURF                                 |
C     | determine coefficients for surface carbon chemistry      |
C     | adapted from OCMIP2:  SUBROUTINE CO2CALC                 |
C     | mick follows, oct 1999                                   |
C     | minor changes to tidy, swd aug 2002                      |
C     *==========================================================*
C INPUT
C       diclocal = total inorganic carbon (mol/m^3)
C             where 1 T = 1 metric ton = 1000 kg
C       ta   = total alkalinity (eq/m^3)
C       pt   = inorganic phosphate (mol/^3)
C       sit  = inorganic silicate (mol/^3)
C       theta= temperature (degrees C)
C       salt = salinity (PSU)
C OUTPUT
C IMPORTANT: Some words about units - (JCO, 4/4/1999)
C     - Models carry tracers in mol/m^3 (on a per volume basis)
C     - Conversely, this routine, which was written by observationalists
C       (C. Sabine and R. Key), passes input arguments in umol/kg
C       (i.e., on a per mass basis)
C     - I have changed things slightly so that input arguments are in mol/m^3,
C     - Thus, all input concentrations (diclocal, ta, pt, and st) should be
C       given in mol/m^3  output arguments "co2star" and "dco2star"
C       are likewise be in mol/m^3.
C
C Apr 2011: fix vapour bug (following Bennington)
C--------------------------------------------------------------------------

C !USES: ===============================================================
        IMPLICIT NONE
C     == GLobal variables ==
#include "SIZE.h"
#include "DYNVARS.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "FFIELDS.h"
#include "DIC_VARS.h"

        INTEGER bi,bj,iMin,iMax,jMin,jMax,ksurf
        INTEGER myThid
CEOP

C LOCAL VARIABLES
        _RL  t
        _RL  s
        _RL  tk
        _RL  tk100
        _RL  tk1002
        _RL  dlogtk
        _RL  sqrtis
        _RL  sqrts
        _RL  s15
        _RL  scl
        _RL  s2
        _RL  invtk
        _RL  is
        _RL  is2
c add Bennington
        _RL  P1atm
        _RL  Rgas
        _RL  RT
        _RL  delta
        _RL  B1
        _RL  B
        _RL  tmpa1
        _RL  tmpa2
        _RL  tmpa3
        _RL  logKspc
        INTEGER i
        INTEGER j
        
        ksurf=1
C.....................................................................
C OCMIP note:
C Calculate all constants needed to convert between various measured
C carbon species. References for each equation are noted in the code.
C Once calculated, the constants are
C stored and passed in the common block "const". The original version
C of this code was based on the code by dickson in Version 2 of
C  Handbook of Methods C for the Analysis of the Various Parameters of
C the Carbon Dioxide System in Seawater , DOE, 1994 (SOP No. 3, p25-26).
C....................................................................

        do i=imin,imax
         do j=jmin,jmax
C           IF (debugMode) CALL DEBUG_MSG(
C     &         'INSIDE I/J LOOP',myThid)
          if (hFacC(i,j,ksurf,bi,bj).ne.0. _d 0) then
           t = theta(i,j,ksurf,bi,bj)
           s = salt(i,j,ksurf,bi,bj)
C terms used more than once
           tk = 273.15 _d 0 + t
           tk100 = tk/100. _d 0
           tk1002=tk100*tk100
           invtk=1.0 _d 0/tk
           dlogtk=log(tk)
           is=19.924 _d 0*s/(1000. _d 0-1.005 _d 0*s)
           is2=is*is
           sqrtis=sqrt(is)
           s2=s*s
           sqrts=sqrt(s)
           s15=s**1.5 _d 0
           scl=s/1.80655 _d 0
C -----------------------------------------------------------------------
C added by Val Bennington Nov 2010
C Fugacity Factor needed for non-ideality in ocean
C ff used for atmospheric correction for water vapor and pressure
C Weiss (1974) Marine Chemistry
           P1atm = 1.01325 _d 0 ! bars
           Rgas = 83.1451 _d 0 ! bar*cm3/(mol*K)
           RT = Rgas*tk
           delta = (57.7 _d 0 - 0.118 _d 0*tk)
           B1 = -1636.75 _d 0 + 12.0408 _d 0*tk - 0.0327957 _d 0*tk*tk
           B = B1 + 3.16528 _d 0*tk*tk*tk*(0.00001 _d 0)
           fugf3d(i,j,ksurf,bi,bj) = exp( (B+2. _d 0*delta) * P1atm / RT)
C------------------------------------------------------------------------
C f = k0(1-pH2O)*correction term for non-ideality
C Weiss & Price (1980, Mar. Chem., 8, 347-359  Eq 13 with table 6 values)
           ff3d(i,j,ksurf,bi,bj) = exp(-162.8301 _d 0 + 218.2968 _d 0/tk100
     &          + 90.9241 _d 0*log(tk100) - 1.47696 _d 0*tk1002 
     &          + s * (.025695 _d 0 - .025225 _d 0*tk100 
     &          + 0.0049867 _d 0*tk1002))
C------------------------------------------------------------------------
C K0 from Weiss 1974
           ak03d(i,j,ksurf,bi,bj) = exp(93.4517 _d 0/tk100 - 60.2409 _d 0
     &        + 23.3585 _d 0 * log(tk100) 
     &        + s * (0.023517 _d 0 - 0.023656 _d 0*tk100 
     &        + 0.0047036 _d 0*tk1002))
C------------------------------------------------------------------------
C k1 = [H][HCO3]/[H2CO3]
C k2 = [H][CO3]/[HCO3]
C Millero p.664 (1995) using Mehrbach et al. data on seawater scale
           ak13d(i,j,ksurf,bi,bj)=10.**(-1. _d 0*(3670.7 _d 0*invtk -
     &          62.008 _d 0 + 9.7944 _d 0*dlogtk -
     &          0.0118 _d 0 * s + 0.000116 _d 0*s2))
           ak23d(i,j,ksurf,bi,bj)=10.**(-1. _d 0*(1394.7 _d 0*invtk
     &          + 4.777 _d 0 - 0.0184 _d 0*s + 0.000118 _d 0*s2))
C------------------------------------------------------------------------
C kb = [H][BO2]/[HBO2]
C Millero p.669 (1995) using data from dickson (1990)
           akb3d(i,j,ksurf,bi,bj)=exp((-8966.90 _d 0- 2890.53 _d 0*sqrts
     &          - 77.942 _d 0*s + 1.728 _d 0*s15 - 0.0996 _d 0*s2)*invtk +
     &          (148.0248 _d 0 + 137.1942 _d 0*sqrts + 1.62142 _d 0*s) +
     &          (-24.4344 _d 0 - 25.085 _d 0*sqrts - 0.2474 _d 0*s) *
     &          dlogtk + 0.053105 _d 0*sqrts*tk)
C------------------------------------------------------------------------
C k1p = [H][H2PO4]/[H3PO4]
C DOE(1994) eq 7.2.20 with footnote using data from Millero (1974)
           ak1p3d(i,j,ksurf,bi,bj) = exp(-4576.752 _d 0*invtk 
     &          + 115.525 _d 0 - 18.453 _d 0*dlogtk +
     &          (-106.736 _d 0*invtk + 0.69171 _d 0)*sqrts +
     &          (-0.65643 _d 0*invtk - 0.01844 _d 0)*s)
C------------------------------------------------------------------------
C k2p = [H][HPO4]/[H2PO4]
C DOE(1994) eq 7.2.23 with footnote using data from Millero (1974))
           ak2p3d(i,j,ksurf,bi,bj) = exp(-8814.715 _d 0*invtk 
     &          + 172.0883 _d 0 - 27.927 _d 0*dlogtk +
     &          (-160.340 _d 0*invtk + 1.3566 _d 0) * sqrts +
     &          (0.37335 _d 0*invtk - 0.05778 _d 0) * s)
C------------------------------------------------------------------------
C k3p = [H][PO4]/[HPO4]
C DOE(1994) eq 7.2.26 with footnote using data from Millero (1974)
           ak3p3d(i,j,ksurf,bi,bj) = exp(-3070.75 _d 0*invtk 
     &          - 18.141 _d 0 + (17.27039 _d 0*invtk + 2.81197 _d 0)
     &          * sqrts + (-44.99486 _d 0*invtk - 0.09984 _d 0) * s)
C------------------------------------------------------------------------
C ksi = [H][SiO(OH)3]/[Si(OH)4]
C Millero p.671 (1995) using data from Yao and Millero (1995)
           aksi3d(i,j,ksurf,bi,bj) = exp(-8904.2 _d 0*invtk 
     &          + 117.385 _d 0 - 19.334 _d 0*dlogtk +
     &          (-458.79 _d 0*invtk + 3.5913 _d 0) * sqrtis +
     &          (188.74 _d 0*invtk - 1.5998 _d 0) * is +
     &          (-12.1652 _d 0*invtk + 0.07871 _d 0) * is2 +
     &          log(1.0 _d 0-0.001005 _d 0*s))
C------------------------------------------------------------------------
C kw = [H][OH]
C Millero p.670 (1995) using composite data
           akw3d(i,j,ksurf,bi,bj) = exp(-13847.26 _d 0*invtk 
     &          + 148.9652 _d 0 - 23.6521 _d 0*dlogtk +
     &          (118.67 _d 0*invtk - 5.977 _d 0 + 1.0495 _d 0 * dlogtk)
     &          * sqrts - 0.01615 _d 0 * s)
C------------------------------------------------------------------------
C ks = [H][SO4]/[HSO4]
C dickson (1990, J. chem. Thermodynamics 22, 113)
           aks3d(i,j,ksurf,bi,bj)=exp(-4276.1 _d 0*invtk + 141.328 _d 0 -
     &          23.093 _d 0*dlogtk +
     &   (-13856. _d 0*invtk + 324.57 _d 0 - 47.986 _d 0*dlogtk)*sqrtis+
     &   (35474. _d 0*invtk - 771.54 _d 0 + 114.723 _d 0*dlogtk)*is -
     &          2698. _d 0*invtk*is**1.5 _d 0 + 1776. _d 0*invtk*is2 +
     &          log(1.0 _d 0 - 0.001005 _d 0*s))
C------------------------------------------------------------------------
C kf = [H][F]/[HF]
C dickson and Riley (1979) -- change pH scale to total
           akf3d(i,j,ksurf,bi,bj)=exp(1590.2 _d 0*invtk - 12.641 _d 0 
     &            + 1.525 _d 0*sqrtis + log(1.0 _d 0 - 0.001005 _d 0*s) 
     &            + log(1.0 _d 0 + (0.1400 _d 0/96.062 _d 0)*(scl)
     &            /aks3d(i,j,ksurf,bi,bj)))
C------------------------------------------------------------------------
C Calculate concentrations for borate, sulfate, and fluoride
C Uppstrom (1974)
           bt3d(i,j,ksurf,bi,bj) = 0.000232 _d 0 * scl/10.811 _d 0
C Morris & Riley (1966)
           st3d(i,j,ksurf,bi,bj) = 0.14 _d 0 * scl/96.062 _d 0
C Riley (1965)
           ft3d(i,j,ksurf,bi,bj) = 0.000067 _d 0 * scl/18.9984 _d 0
C------------------------------------------------------------------------
C solubility product for calcite
C
c Following Takahashi (1982) GEOSECS handbook
C NOT SURE THIS IS WORKING???
C Ingle et al. (1973)
c          Ksp_T_Calc = ( -34.452 - 39.866*(s**0.333333)
c    &                  + 110.21*log(s) - 7.5752d-6 * (tk**2.0)
c    &                  ) * 1.0d-7
c with pressure dependence Culberson and Pytkowicz (1968)
c          xvalue  =  (36-0.20*t)*(pressc-1.0)/(83.143*tk)
c          Ksp_TP_Calc(i,j,bi,bj) = Ksp_T_Calc*exp(xvalue)
c
c
C Following Mucci (1983) - from Zeebe/Wolf-Gladrow equic.m
         tmpa1 = - 171.9065 _d 0 - (0.077993 _d 0*tk) + (2839.319 _d 0/tk)
     &            + (71.595 _d 0*log10(tk))
         tmpa2 = +(-0.77712 _d 0 + (0.0028426 _d 0*tk) + (178.34 _d 0/tk) )*sqrts
         tmpa3 = -(0.07711 _d 0*s) + (0.0041249 _d 0*s15)
         logKspc = tmpa1 + tmpa2 + tmpa3
         Ksp_TP_Calc3d(i,j,ksurf,bi,bj) = 10.0**logKspc
c        write(6,*)i,j,k,tmpa1,tmpa2,tmpa3,logkspc,Ksp_T_Calc
c with pressure dependence Culberson and Pytkowicz (1968)
c        xvalue  =  (36.0-0.20*t)*(pressc-1.0)/(83.143*tk)
c        Ksp_TP_Calc(i,j,bi,bj) = Ksp_T_Calc*exp(xvalue)

c alternative pressure depdendence
c following Millero (1995) but using info from Appendix A11 of
c Zeebe and Wolf-Gladrow (2001) book
c          dv = -48.6 - 0.5304*t
c          dk = -11.76d-3 - 0.3692*t
c          xvalue = - (dv/(bigR*tk))*pressc
c    &               + (0.5*dk/(bigR*tk))*pressc*pressc
c          Ksp_TP_Calc(i,j,bi,bj) = Ksp_T_Calc*exp(xvalue)

c alternative pressure dependence from Ingle (1975)
C Do it at the surface, ie 1 bar
C           pressc = 1. _d 0
C           zdum   = (pressc*10. _d 0 - 10. _d 0)/10. _d 0
C           xvalue = ( (48.8 _d 0 - 0.53 _d 0*t)*zdum
C     &                 + (-0.00588 _d 0 + 0.0001845 _d 0*t)*zdum*zdum)
C     &            / (188.93 _d 0*(t + 273.15 _d 0))
C
C           Ksp_TP_Calc3d(i,j,k,bi,bj) = Ksp_T_tmp*10**(xvalue)
         else
c add Bennington           
            fugf3d(i,j,ksurf,bi,bj)=0. _d 0
            ff3d(i,j,ksurf,bi,bj)=0. _d 0
            ak03d(i,j,ksurf,bi,bj)= 0. _d 0
            ak13d(i,j,ksurf,bi,bj)= 0. _d 0
            ak23d(i,j,ksurf,bi,bj)= 0. _d 0
            akb3d(i,j,ksurf,bi,bj)= 0. _d 0
            ak1p3d(i,j,ksurf,bi,bj) = 0. _d 0
            ak2p3d(i,j,ksurf,bi,bj) = 0. _d 0
            ak3p3d(i,j,ksurf,bi,bj) = 0. _d 0
            aksi3d(i,j,ksurf,bi,bj) = 0. _d 0
            akw3d(i,j,ksurf,bi,bj) = 0. _d 0
            aks3d(i,j,ksurf,bi,bj)= 0. _d 0
            akf3d(i,j,ksurf,bi,bj)= 0. _d 0
            bt3d(i,j,ksurf,bi,bj) = 0. _d 0
            st3d(i,j,ksurf,bi,bj) = 0. _d 0
            ft3d(i,j,ksurf,bi,bj) = 0. _d 0
            Ksp_TP_Calc3d(i,j,ksurf,bi,bj) = 0. _d 0
         endif
         end do
        end do

        RETURN
        END

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

CBOP
C !ROUTINE: CARBON_COEFFS_DEPTH

C !INTERFACE: ==========================================================
      SUBROUTINE CARBON_COEFFS_DEPTH(
     I                   bi,bj,iMin,iMax,jMin,jMax,
     I                   kMin,kMax,myThid)

C !DESCRIPTION:
C     *==========================================================*
C     | SUBROUTINE CARBON_COEFFS_DEPTH                           |
C     | determine coefficients for surface carbon chemistry      |
C     | adapted from OCMIP2:  SUBROUTINE CO2CALC                 |
C     | mick follows, oct 1999                                   |
C     | minor changes to tidy, swd aug 2002                      |
C     | MODIFIED FOR PRESSURE DEPENDENCE                         |
C     | Karsten Friis and Mick Follows 2004                      |
C     *==========================================================*
C INPUT
C       diclocal = total inorganic carbon (mol/m^3)
C             where 1 T = 1 metric ton = 1000 kg
C       ta  = total alkalinity (eq/m^3)
C       pt  = inorganic phosphate (mol/^3)
C       sit = inorganic silicate (mol/^3)
C       theta= temperature (degrees C)
C       salt = salinity (PSU)
C OUTPUT
C IMPORTANT: Some words about units - (JCO, 4/4/1999)
C     - Models carry tracers in mol/m^3 (on a per volume basis)
C     - Conversely, this routine, which was written by observationalists
C       (C. Sabine and R. Key), passes input arguments in umol/kg
C       (i.e., on a per mass basis)
C     - I have changed things slightly so that input arguments are in mol/m^3,
C     - Thus, all input concentrations (diclocal, ta, pt, and st) should be
C       given in mol/m^3  output arguments "co2star" and "dco2star"
C       are likewise be in mol/m^3.
C
C Apr 2011: fix vapour bug (following Bennington)
C
C NOW INCLUDES:
C PRESSURE DEPENDENCE of K1, K2, solubility product of calcite
C based on Takahashi, GEOSECS Atlantic Report, Vol. 1 (1981)
C
C--------------------------------------------------------------------------

C !USES: ===============================================================
        IMPLICIT NONE
C     == GLobal variables ==
#include "SIZE.h"
#include "DYNVARS.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "FFIELDS.h"
#include "DIC_VARS.h"
        INTEGER bi,bj,iMin,iMax,jMin,jMax,kMin,kMax
C K is depth index
        INTEGER myThid
CEOP

C LOCAL VARIABLES
        _RL  t
        _RL  s
        _RL  tk
        _RL  tk100
        _RL  tk1002
        _RL  dlogtk
        _RL  sqrtis
        _RL  sqrts
        _RL  s15
        _RL  scl
        _RL  s2
        _RL  invtk
        _RL  is
        _RL  is2
        INTEGER i
        INTEGER j
        INTEGER k
        INTEGER kk
        _RL  bdepth
        _RL  cdepth
        _RL  pressc
        _RL  Ksp_T_tmp
        _RL  xvalue
        _RL  zdum
        _RL  tmpa1
        _RL  tmpa2
        _RL  tmpa3
        _RL  logKspc
        _RL  dv
        _RL  dk
        _RL  pfactor
        _RL  bigR
c add Bennington
        _RL  P1atm
        _RL  Rgas
        _RL  RT
        _RL  delta
        _RL  B1
        _RL  B
C.....................................................................
C OCMIP note:
C Calculate all constants needed to convert between various measured
C carbon species. References for each equation are noted in the code.
C Once calculated, the constants are
C stored and passed in the common block "const". The original version
C of this code was based on the code by dickson in Version 2 of
C  Handbook of Methods C for the Analysis of the Various Parameters of
C the Carbon Dioxide System in Seawater , DOE, 1994 (SOP No. 3, p25-26).
C....................................................................

c determine pressure (bar) from depth
c 1 BAR at z=0m (atmos pressure)
c use UPPER surface of cell so top layer pressure = 0 bar
c for surface exchange coeffs

       do i=imin,imax
        do j=jmin,jmax
         do k = kMin,KMax
          bdepth = 0. _d 0
          cdepth = 0. _d 0
          pressc = 1. _d 0
          do kk = kMin,k
            cdepth = bdepth + 0.5 _d 0*drF(kk)
            bdepth = bdepth + drF(kk)
            pressc = 1. _d 0 + 0.1 _d 0*cdepth
          end do
          
          if (hFacC(i,j,k,bi,bj).ne.0. _d 0) then
           t = theta(i,j,k,bi,bj)
           s = salt (i,j,k,bi,bj)
           
C terms used more than once
           tk = 273.15 _d 0 + t
           tk100 = tk/100. _d 0
           tk1002=tk100*tk100
           invtk=1.0 _d 0/tk
           dlogtk=log(tk)
           is=19.924 _d 0*s/(1000. _d 0-1.005 _d 0*s)
           is2=is*is
           sqrtis=sqrt(is)
           s2=s*s
           sqrts=sqrt(s)
           s15=s**1.5 _d 0
           scl=s/1.80655 _d 0
           
C -----------------------------------------------------------------------
C added by Val Bennington Nov 2010
C Fugacity Factor needed for non-ideality in ocean
C ff used for atmospheric correction for water vapor and pressure
C Weiss (1974) Marine Chemistry
C           P1atm = 1.01325 _d 0 ! bars
           P1atm = pressc
           Rgas = 83.1451 _d 0 ! bar*cm3/(mol*K)
           RT = Rgas*tk
           delta = (57.7 _d 0 - 0.118 _d 0*tk)
           B1 = -1636.75 _d 0 + 12.0408 _d 0*tk - 0.0327957 _d 0*tk*tk
           B = B1 + 3.16528 _d 0*tk*tk*tk*(0.00001 _d 0)
           fugf3d(i,j,k,bi,bj) = exp( (B+2. _d 0*delta) * P1atm / RT)
C------------------------------------------------------------------------
C f = k0(1-pH2O)*correction term for non-ideality
C Weiss & Price (1980, Mar. Chem., 8, 347-359  Eq 13 with table 6 values)
           ff3d(i,j,k,bi,bj) = exp(-162.8301 _d 0 + 218.2968 _d 0/tk100
     &          + 90.9241 _d 0*log(tk100) - 1.47696 _d 0*tk1002
     &          + s * (.025695 _d 0 - .025225 _d 0*tk100
     &          + 0.0049867 _d 0*tk1002))
C------------------------------------------------------------------------
C K0 from Weiss 1974
           ak03d(i,j,k,bi,bj) = exp(93.4517 _d 0/tk100 - 60.2409 _d 0
     &        + 23.3585 _d 0 * log(tk100)
     &        + s * (0.023517 _d 0 - 0.023656 _d 0*tk100
     &        + 0.0047036 _d 0*tk1002))
C------------------------------------------------------------------------
C k1 = [H][HCO3]/[H2CO3]
C k2 = [H][CO3]/[HCO3]
C Millero p.664 (1995) using Mehrbach et al. data on seawater scale
           ak13d(i,j,k,bi,bj)=10**(-1 _d 0*(3670.7 _d 0*invtk -
     &          62.008 _d 0 + 9.7944 _d 0*dlogtk -
     &          0.0118 _d 0 * s + 0.000116 _d 0*s2))
           ak23d(i,j,k,bi,bj)=10**(-1*(1394.7 _d 0*invtk 
     &          + 4.777 _d 0 - 0.0184 _d 0*s + 0.000118 _d 0*s2))
C NOW PRESSURE DEPENDENCE:
c Following Takahashi (1981) GEOSECS report - quoting Culberson and
c Pytkowicz (1968)
c pressc = pressure in bars
           ak13d(i,j,k,bi,bj) = ak13d(i,j,k,bi,bj)*
     &         exp( (24.2 _d 0-0.085 _d 0*t)*(pressc-1.0 _d 0)/(83.143 _d 0*tk) )
c FIRST GO FOR K2: According to GEOSECS (1982) report
c          ak2(i,j,bi,bj) = ak2(i,j,bi,bj)*
c    &             exp( (26.4-0.040*t)*(pressc-1.0)/(83.143*tk) )
c SECOND GO FOR K2: corrected coeff according to CO2sys documentation
c                   E. Lewis and D. Wallace (1998) ORNL/CDIAC-105
           ak23d(i,j,k,bi,bj) = ak23d(i,j,k,bi,bj)*
     &         exp( (16.4 _d 0-0.040 _d 0*t)*(pressc-1.0 _d 0)/(83.143 _d 0*tk) )
C------------------------------------------------------------------------
C kb = [H][BO2]/[HBO2]
C Millero p.669 (1995) using data from dickson (1990)
           akb3d(i,j,k,bi,bj)=exp((-8966.90 _d 0- 2890.53 _d 0*sqrts -
     &          77.942 _d 0*s + 1.728 _d 0*s15 - 0.0996 _d 0*s2)*invtk +
     &          (148.0248 _d 0 + 137.1942 _d 0*sqrts + 1.62142 _d 0*s) +
     &          (-24.4344 _d 0 - 25.085 _d 0*sqrts - 0.2474 _d 0*s) *
     &          dlogtk + 0.053105 _d 0*sqrts*tk)

C Mick and Karsten - Dec 04
C ADDING pressure dependence based on Millero (1995), p675
C with additional info from CO2sys documentation (E. Lewis and
C D. Wallace, 1998 - see endnotes for commentary on Millero, 95)
           bigR = 83.145 _d 0
           dv = -29.48 _d 0 + 0.1622 _d 0*t + 2.608 _d -3*t*t
           dk = -2.84 _d -3
           pfactor = - (dv/(bigR*tk))*pressc
     &               + (0.5 _d 0*dk/(bigR*tk))*pressc*pressc
           akb3d(i,j,k,bi,bj) = akb3d(i,j,k,bi,bj)*exp(pfactor)
C------------------------------------------------------------------------
C k1p = [H][H2PO4]/[H3PO4]
C DOE(1994) eq 7.2.20 with footnote using data from Millero (1974)
           ak1p3d(i,j,k,bi,bj) = exp(-4576.752 _d 0*invtk 
     &          + 115.525 _d 0 - 18.453 _d 0*dlogtk +
     &          (-106.736 _d 0*invtk + 0.69171 _d 0)*sqrts +
     &          (-0.65643 _d 0*invtk - 0.01844 _d 0)*s)
C------------------------------------------------------------------------
C k2p = [H][HPO4]/[H2PO4]
C DOE(1994) eq 7.2.23 with footnote using data from Millero (1974))
           ak2p3d(i,j,k,bi,bj) = exp(-8814.715 _d 0*invtk 
     &          + 172.0883 _d 0 - 27.927 _d 0*dlogtk +
     &          (-160.34 _d 00*invtk + 1.3566 _d 0) * sqrts +
     &          (0.37335 _d 0*invtk - 0.05778 _d 0) * s)
C------------------------------------------------------------------------
C k3p = [H][PO4]/[HPO4]
C DOE(1994) eq 7.2.26 with footnote using data from Millero (1974)
           ak3p3d(i,j,k,bi,bj) = exp(-3070.75 _d 0*invtk 
     &          - 18.141 _d 0 + (17.27039 _d 0*invtk + 2.81197 _d 0) *
     &          sqrts + (-44.99486 _d 0*invtk - 0.09984 _d 0) * s)
C------------------------------------------------------------------------
C ksi = [H][SiO(OH)3]/[Si(OH)4]
C Millero p.671 (1995) using data from Yao and Millero (1995)
           aksi3d(i,j,k,bi,bj) = exp(-8904.2 _d 0*invtk 
     &          + 117.385 _d 0 - 19.334 _d 0*dlogtk +
     &          (-458.79 _d 0*invtk + 3.5913 _d 0) * sqrtis +
     &          (188.74 _d 0*invtk - 1.5998 _d 0) * is +
     &          (-12.1652 _d 0*invtk + 0.07871 _d 0) * is2 +
     &          log(1.0 _d 0-0.001005 _d 0*s))
C------------------------------------------------------------------------
C kw = [H][OH]
C Millero p.670 (1995) using composite data
           akw3d(i,j,k,bi,bj) = exp(-13847.26 _d 0*invtk + 148.9652 _d 0 -
     &          23.6521 _d 0*dlogtk +
     &          (118.67 _d 0*invtk - 5.977 _d 0 + 1.0495 _d 0 * dlogtk) *
     &          sqrts - 0.01615 _d 0 * s)
C------------------------------------------------------------------------
C ks = [H][SO4]/[HSO4]
C dickson (1990, J. chem. Thermodynamics 22, 113)
           aks3d(i,j,k,bi,bj)=exp(-4276.1 _d 0*invtk + 141.328 _d 0 -
     &          23.093 _d 0*dlogtk +
     &          (-13856 _d 0*invtk + 324.57 _d 0 - 47.986 _d 0*dlogtk)*sqrtis +
     &          (35474 _d 0*invtk - 771.54 _d 0 + 114.723 _d 0*dlogtk)*is -
     &          2698 _d 0*invtk*is**1.5 _d 0 + 1776 _d 0*invtk*is2 +
     &          log(1.0 _d 0 - 0.001005 _d 0*s))
C------------------------------------------------------------------------
C kf = [H][F]/[HF]
C dickson and Riley (1979) -- change pH scale to total
        akf3d(i,j,k,bi,bj)=exp(1590.2 _d 0*invtk - 12.641 _d 0 
     &          + 1.525 _d 0*sqrtis + log(1.0 _d 0 - 0.001005 _d 0*s)
     &          + log(1.0 _d 0 + (0.1400 _d 0/96.062 _d 0)*(scl)
     &          /aks3d(i,j,k,bi,bj)))
C------------------------------------------------------------------------
C Calculate concentrations for borate, sulfate, and fluoride
C Uppstrom (1974)
           bt3d(i,j,k,bi,bj) = 0.000232 _d 0 * scl/10.811 _d 0
C Morris & Riley (1966)
           st3d(i,j,k,bi,bj) = 0.14 _d 0 * scl/96.062 _d 0
C Riley (1965)
           ft3d(i,j,k,bi,bj) = 0.000067 _d 0 * scl/18.9984 _d 0
C------------------------------------------------------------------------
C solubility product for calcite
C
c Following Takahashi (1982) GEOSECS handbook
C NOT SURE THIS IS WORKING???
C Ingle et al. (1973)
c          Ksp_T_Calc = ( -34.452 - 39.866*(s**0.333333)
c    &                  + 110.21*log(s) - 7.5752d-6 * (tk**2.0)
c    &                  ) * 1.0d-7
c with pressure dependence Culberson and Pytkowicz (1968)
c          xvalue  =  (36-0.20*t)*(pressc-1.0)/(83.143*tk)
c          Ksp_TP_Calc(i,j,bi,bj) = Ksp_T_Calc*exp(xvalue)
c
c
C Following Mucci (1983) - from Zeebe/Wolf-Gladrow equic.m
         tmpa1 = - 171.9065 _d 0 - (0.077993 _d 0*tk) + (2839.319 _d 0/tk)
     &            + (71.595 _d 0*log10(tk))
         tmpa2 = +(-0.77712 _d 0 + (0.0028426 _d 0*tk) + (178.34 _d 0/tk) )*sqrts
         tmpa3 = -(0.07711 _d 0*s) + (0.0041249 _d 0*s15)
         logKspc = tmpa1 + tmpa2 + tmpa3
         Ksp_T_tmp = 10.0**logKspc
c        write(6,*)i,j,k,tmpa1,tmpa2,tmpa3,logkspc,Ksp_T_Calc
c with pressure dependence Culberson and Pytkowicz (1968)
c        xvalue  =  (36.0-0.20*t)*(pressc-1.0)/(83.143*tk)
c        Ksp_TP_Calc(i,j,bi,bj) = Ksp_T_Calc*exp(xvalue)

c alternative pressure depdendence
c following Millero (1995) but using info from Appendix A11 of
c Zeebe and Wolf-Gladrow (2001) book
c          dv = -48.6 - 0.5304*t
c          dk = -11.76d-3 - 0.3692*t
c          xvalue = - (dv/(bigR*tk))*pressc
c    &               + (0.5*dk/(bigR*tk))*pressc*pressc
c          Ksp_TP_Calc(i,j,bi,bj) = Ksp_T_Calc*exp(xvalue)

c alternative pressure dependence from Ingle (1975)

           zdum   = (pressc*10. _d 0 - 10. _d 0)/10. _d 0
           xvalue = ( (48.8 _d 0 - 0.53 _d 0*t)*zdum
     &                 + (-0.00588 _d 0 + 0.0001845 _d 0*t)*zdum*zdum)
     &            / (188.93 _d 0*(t + 273.15 _d 0))

           Ksp_TP_Calc3d(i,j,k,bi,bj) = Ksp_T_tmp*10**(xvalue)

C------------------------------------------------------------------------
         else
C add Bennington
            fugf3d(i,j,k,bi,bj)=0. _d 0
            ff3d(i,j,k,bi,bj)=0. _d 0
            ak03d(i,j,k,bi,bj)= 0. _d 0
            ak13d(i,j,k,bi,bj)= 0. _d 0
            ak23d(i,j,k,bi,bj)= 0. _d 0
            akb3d(i,j,k,bi,bj)= 0. _d 0
            ak1p3d(i,j,k,bi,bj) = 0. _d 0
            ak2p3d(i,j,k,bi,bj) = 0. _d 0
            ak3p3d(i,j,k,bi,bj) = 0. _d 0
            aksi3d(i,j,k,bi,bj) = 0. _d 0
            akw3d(i,j,k,bi,bj) = 0. _d 0
            aks3d(i,j,k,bi,bj)= 0. _d 0
            akf3d(i,j,k,bi,bj)= 0. _d 0
            bt3d(i,j,k,bi,bj) = 0. _d 0
            st3d(i,j,k,bi,bj) = 0. _d 0
            ft3d(i,j,k,bi,bj) = 0. _d 0
            Ksp_TP_Calc3d(i,j,k,bi,bj) = 0. _d 0
         endif
         end do
        end do
       end do
      return
      end
#endif /* DIC_CARBON_COMPONENTS */
