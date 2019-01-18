#include "DIC_OPTIONS.h"

C--  File carbon_chem.F:
C--   Contents
C--   o CALC_CSAT
C--   o CARBON_COEFFS_SURF
C--   o CARBON_COEFFS_DEPTH

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
#ifdef DIC_CARBON_COMPONENTS
C !ROUTINE: CALC_CSAT

C !INTERFACE: ==========================================================
       SUBROUTINE CALC_CSAT(
     I                       t,s,pco2eq,pt,sit,ta,
     U                       pHlocal,csatloc,
     I                       i,j,k,bi,bj,myIter,myThid )

C !DESCRIPTION:
C     *==========================================================*
C     | SUBROUTINE CALC_CSAT                            |
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
C     &      'CALC_CSAT: ph guess'
        phguess = phlocal
C     - new approx method
C     - make estimate of htotal (hydrogen ion conc) using
C       appromate estimate of CA, carbonate alkalinity
        hguess = 10.0 _d 0**(-phguess)
        
C     - first estimate borate contribution using guess for [H+]
        bohg = (bt(i,j,k,bi,bj)*akb(i,j,k,bi,bj))
     &        /(hguess+akb(i,j,k,bi,bj))

C     - first estimate of contribution from phosphate
C     based on Dickson and Goyet
        denom = (hguess*hguess*hguess)
     &        + (ak1p(i,j,k,bi,bj)*hguess*hguess)
     &        + (ak1p(i,j,k,bi,bj)*ak2p(i,j,k,bi,bj)*hguess)
     &        + (ak1p(i,j,k,bi,bj)*ak2p(i,j,k,bi,bj)
     &          *ak3p(i,j,k,bi,bj))
        h3po4g = (pt*hguess*hguess*hguess) / denom
        h2po4g = (pt*ak1p(i,j,k,bi,bj)*hguess*hguess) / denom
        hpo4g  = (pt*ak1p(i,j,k,bi,bj)*ak2p(i,j,k,bi,bj)
     &              *hguess) / denom
        po4g   = (pt*ak1p(i,j,k,bi,bj)*ak2p(i,j,k,bi,bj)
     &              *ak3p(i,j,k,bi,bj)) / denom

C     - estimate contribution from silicate
C     based on Dickson and Goyet
        siooh3g = (sit*aksi(i,j,k,bi,bj)) 
     &              / (aksi(i,j,k,bi,bj) + hguess)

C     - now estimate carbonate alkalinity
        cag = ta - bohg - (akw(i,j,k,bi,bj)/hguess) + hguess
     &           - hpo4g - 2.0 _d 0*po4g + h3po4g
     &           - siooh3g

C     - estimate hydrogen ion conc
C        stuff = (k1(i,j,k,bi,bj)*fugf(i,j,k,bi,bj)*k0(i,j,k,bi,bj)*pco2eq)/cag
        stuff = (ak1(i,j,k,bi,bj)*ff(i,j,k,bi,bj)*pco2eq)/cag  
        stuff2 = stuff*(stuff + 8.0 _d 0*ak2(i,j,k,bi,bj))  
        stuff3 = sqrt(stuff2)  
        hnew = 0.5 _d 0*(stuff + stuff3)  
        hguess = hnew  

C evaluate csat, equilibrium DIC concentration
        csatloc = (pco2eq*ff(i,j,k,bi,bj))
     &   *(1.0 _d 0 + (ak1(i,j,k,bi,bj)/hnew) 
     &   +((ak1(i,j,k,bi,bj)*ak2(i,j,k,bi,bj))/(hnew*hnew)))

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
     I                   bi,bj,iMin,iMax,jMin,jMax,
     I                   k,myThid)

C     Determine coefficients for surface carbon chemistry,
C     loaded into common block

      IMPLICIT NONE

C     MITgcm GLobal variables
#include "SIZE.h"
#include "DYNVARS.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "FFIELDS.h"
#include "DIC_VARS.h"

      INTEGER bi,bj,iMin,iMax,jMin,jMax,myThid
      INTEGER  i, j, k
      _RL t
      _RL s
      _RL t_k
      _RL t_k_o_100
      _RL t_k_o_100_2
      _RL dlog_t_k
      _RL sqrtis
      _RL sqrts
      _RL inv_t_k
      _RL ion_st
      _RL is_2
      _RL scl
      _RL s35
      _RL s_2
      _RL s_15
      _RL log_kgfw2kgsw
      _RL B1
      _RL B
      _RL delta
      _RL P1atm
      _RL RT
      _RL Rgas
C conversions for different pH scales
      _RL total2free
      _RL free2sw
      _RL total2sw
     
        DO i=imin,imax
         DO j=jmin,jmax
          IF (hFacC(i,j,k,bi,bj).gt.0. _d 0) THEN
           t = theta(i,j,k,bi,bj)
           s = salt (i,j,k,bi,bj)
C terms used more than once for:
C temperature
           t_k         = 273.15 _d 0 + t
           t_k_o_100   = t_k/100. _d 0
           t_k_o_100_2 = t_k_o_100*t_k_o_100
           inv_t_k=1.0 _d 0/t_k
           dlog_t_k=LOG(t_k)
C ionic strength (converted to kgsw)
           ion_st=19.924 _d 0*s/(1000. _d 0-1.005 _d 0*s)
           is_2=ion_st*ion_st
           sqrtis=sqrt(ion_st)
C salinity
           s_2  = s*s
           sqrts=sqrt(s)
           s_15 = s*sqrts
           scl  =s/1.80655 _d 0
           s35  =s/35. _d 0

           log_kgfw2kgsw = LOG(1. _d 0 - 0.001005 _d 0*s)
           
C These coefficients largely come from Munhoven (2013) SolveSAPHE Package
C -----------------------------------------------------------------------
C       Calculate the total borate concentration in mol/kg-SW
C       given the salinity of a sample
C       Ref: Uppström (1974), cited by  Dickson et al. (2007, chapter 5, p 10)
C            Millero (1982) cited in Millero (1995)
C       pH scale  : N/A
           bt(i,j,k,bi,bj) = 0.000232 _d 0* scl/10.811 _d 0
           
C -----------------------------------------------------------------------
C       Calculate the total fluoride concentration in mol/kg-SW
C       given the salinity of a sample
C       References: Riley (1965)
C       pH scale  : N/A
           ft(i,j,k,bi,bj) = 0.000067 _d 0 * scl/18.9984 _d 0
           
C -----------------------------------------------------------------------
C       Calculate the total sulfate concentration in mol/kg-SW
C       given the salinity of a sample
C       References: Morris & Riley (1966) quoted in Handbook (2007)
C       pH scale  : N/A
           st(i,j,k,bi,bj) = 0.14 _d 0 * scl/96.062 _d 0

C -----------------------------------------------------------------------
C       Calculate the total calcium concentration in mol/kg-SW
C       given the salinity of a sample
C       References: Culkin (1965) (???)
C       pH scale  : N/A
           cat(i,j,k,bi,bj) = 0.010282 _d 0*(s35)

C -----------------------------------------------------------------------
C       Calculate K0 in (mol/kg-SW)/atmosphere
C       References: Weiss (1979) [(mol/kg-SW)/atm]
C       pH scale  : N/A
C       Note      : currently no pressure correction
           ak0(i,j,k,bi,bj)  = EXP( 93.4517 _d 0/t_k_o_100 - 60.2409 _d 0
     &                 + 23.3585 _d 0*LOG(t_k_o_100)
     &                 + s * (0.023517 _d 0 - 0.023656 _d 0*t_k_o_100
     &                 + 0.0047036 _d 0*t_k_o_100_2))

C------------------------------------------------------------------------
C       Calculate f = k0(1-pH2O)*correction term for non-ideality
C       References: Weiss & Price (1980, Mar. Chem., 8, 347-359
C                   Eq 13 with table 6 values)
C       pH scale  : N/A
           ff(i,j,k,bi,bj) = exp(-162.8301 _d 0 + 218.2968 _d 0/t_k_o_100
     &          + 90.9241 _d 0*log(t_k_o_100) - 1.47696 _d 0*t_k_o_100_2
     &          + s * (.025695 _d 0 - .025225 _d 0*t_k_o_100
     &          + 0.0049867 _d 0*t_k_o_100_2))
     
C------------------------------------------------------------------------
C       Calculate Fugacity Factor needed for non-ideality in ocean
C       References: Weiss (1974) Marine Chemistry
C       pH scale  : N/A
           P1atm = 1.01325 _d 0 ! bars
           Rgas = 83.1451 _d 0  ! bar*cm3/(mol*K)
           RT = Rgas*t_k
           delta = (57.7 _d 0 - 0.118 _d 0*t_k)
           B1 = -1636.75 _d 0 + 12.0408 _d 0*t_k - 0.0327957 _d 0*t_k*t_k
           B  = B1 + 3.16528 _d 0*t_k*t_k*t_k*(0.00001 _d 0)

C   "x2" term often neglected (assumed=1) in applications of Weiss's (1974) eq.9
C    x2 = 1 - x1 = 1 - xCO2 (it is very close to 1, but not quite)
           fugf(i,j,k,bi,bj) = exp( (B+2. _d 0*delta) * P1atm / RT)
           
      IF ( selectK1K2const.EQ.1 ) THEN
C -----------------------------------------------------------------------
C       Calculate first dissociation constant of carbonic acid
C       in mol/kg-SW on the Seawater pH-scale.
C       References: Millero (1995, eq 35 -- pK1(MEHR));
C                   Mehrbach et al. (1973) data
C       pH scale:   Seawater
           ak1(i,j,k,bi,bj)=10.**(-1. _d 0*(3670.7 _d 0*inv_t_k
     &          - 62.008 _d 0 + 9.7944 _d 0*dlog_t_k
     &          - 0.0118 _d 0 * s + 0.000116 _d 0*s_2))

C       Calculate second dissociation constant of carbonic acid
C       in mol/kg-SW on the Seawater pH-scale.
C       References: Millero (1995, eq 36 -- pK2(MEHR))
C                   Mehrbach et al. (1973) data
C       pH scale:   Seawater
           ak2(i,j,k,bi,bj)=10.**(-1. _d 0*(1394.7 _d 0*inv_t_k
     &           + 4.777 _d 0 - 0.0184 _d 0*s + 0.000118 _d 0*s_2))

      ELSEIF ( selectK1K2const.EQ.2 ) THEN
C -----------------------------------------------------------------------
C       Calculate first dissociation constant of carbonic acid
C       in mol/kg-SW on the Total pH-scale.
C       References: Roy et al. (1993) -- also Handbook (1994)
C       pH scale  : Total
C       Note      : converted here from mol/kg-H2O to mol/kg-SW
           ak1(i,j,k,bi,bj)  =  EXP(-2307.1255 _d 0*inv_t_k + 2.83655 _d 0
     &              - 1.5529413 _d 0*dlog_t_k
     &              + (-4.0484 _d 0*inv_t_k - 0.20760841)*sqrts
     &              + 0.08468345*s
     &              - 0.00654208*s_15
     &              + log_kgfw2kgsw)

C       Calculate second dissociation constant of carbonic acid
C       in mol/kg-SW on the Total pH-scale.
C       References: Roy et al. (1993) -- also Handbook (1994)
C       pH scale  : Total
C       Note      : converted here from mol/kg-H2O to mol/kg-SW
           ak2(i,j,k,bi,bj) = EXP(-3351.6106 _d 0*inv_t_k - 9.226508 _d 0
     &              - 0.2005743 _d 0*dlog_t_k
     &              + ( -23.9722 _d 0*inv_t_k - 0.106901773 _d 0)*sqrts
     &              + 0.1130822*s - 0.00846934 _d 0*s_15
     &              + log_kgfw2kgsw)

      ELSEIF ( selectK1K2const.EQ.3 ) THEN
C -----------------------------------------------------------------------
C       Calculate first dissociation constant of carbonic acid
C       in mol/kg-SW on the Seawater pH-scale.
C       References: Millero (1995, eq 50 -- ln K1(COM))
C       pH scale:   Seawater
           ak1(i,j,k,bi,bj)  = EXP(2.18867 _d 0 - 2275.0360 _d 0*inv_t_k
     &              - 1.468591 _d 0*dlog_t_k
     &              + ( -0.138681 _d 0 - 9.33291 _d 0*inv_t_k)*sqrts
     &              + 0.0726483 _d 0*s - 0.00574938 _d 0*s_15)

C       Calculate second dissociation constant of carbonic acid
C       in mol/kg-SW on the Seawater pH-scale.
C       References: Millero (1995, eq 51 -- ln K2(COM))
C       pH scale:   Seawater
           ak2(i,j,k,bi,bj)  = EXP(-0.84226 _d 0 - 3741.1288 _d 0*inv_t_k
     &              -  1.437139 _d 0*dlog_t_k
     &              + (-0.128417 _d 0 - 24.41239 _d 0*inv_t_k)*sqrts
     &              +  0.1195308 _d 0*s - 0.00912840 _d 0*s_15)

      ELSEIF ( selectK1K2const.EQ.4 ) THEN
C -----------------------------------------------------------------------
C       Calculate first dissociation constant of carbonic acid
C       in mol/kg-SW on the Total pH-scale.
C       Suitable when 2 < T < 35 and 19 < S < 43
C       References: Luecker et al. (2000) -- also Handbook (2007)
C       pH scale:   Total
           ak1(i,j,k,bi,bj)  = 10. _d 0**( 61.2172 _d 0
     &                 - 3633.86 _d 0*inv_t_k - 9.67770 _d 0*dlog_t_k
     &                 + s*(0.011555 - s*0.0001152 _d 0))

C       Calculate second dissociation constant of carbonic acid
C       in mol/kg-SW on the Total pH-scale.
C       Suitable when 2 < T < 35 and 19 < S < 43
C       References: Luecker et al. (2000) -- also Handbook (2007)
C       pH scale:   Total
           ak2(i,j,k,bi,bj)  = 10. _d 0**(-25.9290 _d 0
     &                 - 471.78 _d 0*inv_t_k + 3.16967 _d 0*dlog_t_k
     &                 + s*(0.01781 _d 0 - s*0.0001122 _d 0))

      ELSEIF ( selectK1K2const.EQ.5 ) THEN
C -----------------------------------------------------------------------
C       Calculate first dissociation constant of carbonic acid
C       in mol/kg-SW on the Seawater pH-scale.
C       Suitable when 0 < T < 50 and 1 < S < 50
C       References: Millero (2010, Mar. Fresh Wat. Res.)
C                   Millero (1979) pressure correction
C       pH scale:   Seawater
           ak1(i,j,k,bi,bj) = 10.0 _d 0**(-1*(6320.813 _d 0*inv_t_k
     &      + 19.568224 _d 0*dlog_t_k -126.34048 _d 0
     &      + 13.4038 _d 0*sqrts + 0.03206 _d 0*s - (5.242 _d -5)*s_2
     &      + (-530.659 _d 0*sqrts - 5.8210 _d 0*s)*inv_t_k
     &      -2.0664 _d 0*sqrts*dlog_t_k))

C       Calculate second dissociation constant of carbonic acid
C       in mol/kg-SW on the Seawater pH-scale.
C       Suitable when 0 < T < 50 and 1 < S < 50
C       References: Millero (2010, Mar. Fresh Wat. Res.)
C                   Millero (1979) pressure correction
C       pH scale:   Seawater
           ak2(i,j,k,bi,bj) = 10.0 _d 0**(-1*(5143.692 _d 0*inv_t_k
     &     + 14.613358 _d 0*dlog_t_k -90.18333 _d 0
     &     + 21.3728 _d 0*sqrts + 0.1218 _d 0*s - (3.688 _d -4)*s_2
     &     + (-788.289 _d 0*sqrts - 19.189 _d 0*s)*inv_t_k
     &     -3.374 _d 0*sqrts*dlog_t_k))

      ELSEIF ( selectK1K2const.EQ.6 ) THEN
C -----------------------------------------------------------------------
C       Calculate first dissociation constant of carbonic acid
C       in mol/kg-SW on the Seawater pH-scale.
C       Suitable when 0 < T < 50 and 1 < S < 50
C       References: Waters, Millero, Woosley (Mar. Chem., 165, 66-67, 2014)
C                   Millero (1979) pressure correction
C       pH scale:   Seawater
           ak1(i,j,k,bi,bj) = 10.0 _d 0**(-1*(6320.813 _d 0*inv_t_k
     &     + 19.568224 _d 0*dlog_t_k -126.34048 _d 0
     &     + 13.409160 _d 0*sqrts + 0.031646 _d 0*s - (5.1895 _d -5)*s_2
     &     + (-531.3642 _d 0*sqrts - 5.713 _d 0*s)*inv_t_k
     &     -2.0669166 _d 0*sqrts*dlog_t_k))

C       Calculate second dissociation constant of carbonic acid
C       in mol/kg-SW on the Seawater pH-scale.
C       Suitable when 0 < T < 50 and 1 < S < 50
C       References: Waters, Millero, Woosley (Mar. Chem., 165, 66-67, 2014)
C                   Millero (1979) pressure correction
C       pH scale:   Seawater
           ak2(i,j,k,bi,bj) = 10.0 _d 0**(-1*
     &     ( 5143.692 _d 0*inv_t_k + 14.613358 _d 0*dlog_t_k
     &      - 90.18333 _d 0 + 21.225890 _d 0*sqrts + 0.12450870 _d 0*s
     &      - (3.7243 _d -4)*s_2
     &      + (-779.3444 _d 0*sqrts - 19.91739 _d 0*s)*inv_t_k
     &      - 3.3534679 _d 0*sqrts*dlog_t_k ) )

      ENDIF /* selectK1K2const */
      
C -----------------------------------------------------------------------
C       Calculate boric acid dissociation constant KB
C       in mol/kg-SW on the total pH-scale.
C       References: Dickson (1990, eq. 23) -- also Handbook (2007, eq. 37)
C       pH scale  : total
           akb(i,j,k,bi,bj)  = EXP(( -8966.90 _d 0 - 2890.53 _d 0*sqrts
     &      -77.942 _d 0*s + 1.728 _d 0*s_15 - 0.0996 _d 0*s_2 )*inv_t_k
     &      + (148.0248 _d 0 + 137.1942 _d 0*sqrts + 1.62142 _d 0*s)
     &      + (-24.4344 _d 0 - 25.085 _d 0*sqrts - 0.2474 _d 0*s)* 
     &      dlog_t_k + 0.053105 _d 0*sqrts*t_k )

C -----------------------------------------------------------------------
C       Calculate the first dissociation constant
C       of phosphoric acid (H3PO4) in seawater
C       References: Yao and Millero (1995)
C       pH scale  : Seawater
           ak1p(i,j,k,bi,bj) = EXP(115.54 _d 0
     &              - 4576.752 _d 0*inv_t_k - 18.453 _d 0*dlog_t_k
     &              + ( 0.69171 _d 0 -  106.736 _d 0*inv_t_k)*sqrts
     &              + (-0.01844 _d 0 -  0.65643 _d 0*inv_t_k)*s)

C -----------------------------------------------------------------------
C       Calculate the second dissociation constant
C       of phosphoric acid (H3PO4) in seawater
C       References: Yao and Millero (1995)
C       pH scale  : Seawater
           ak2p(i,j,k,bi,bj) = EXP( 172.1033 _d 0
     &              - 8814.715 _d 0*inv_t_k
     &              -   27.927 _d 0*dlog_t_k
     &              + (  1.3566 _d 0 -  160.340 _d 0*inv_t_k)*sqrts
     &              + (-0.05778 _d 0 +  0.37335 _d 0*inv_t_k)*s)

C -----------------------------------------------------------------------
C       Calculate the third dissociation constant
C       of phosphoric acid (H3PO4) in seawater
C       References: Yao and Millero (1995)
C       pH scale  : Seawater
           ak3p(i,j,k,bi,bj) = EXP(-18.126 _d 0  -  3070.75 _d 0*inv_t_k
     &                + ( 2.81197 _d 0 + 17.27039 _d 0*inv_t_k)*sqrts
     &                + (-0.09984 _d 0 - 44.99486 _d 0*inv_t_k)*s)

C -----------------------------------------------------------------------
C       Calculate the first dissociation constant
C       of silicic acid (H4SiO4) in seawater
C       References: Yao and Millero (1995) cited by Millero (1995)
C       pH scale  : Seawater (according to Dickson et al, 2007)
C       Note      : converted here from mol/kg-H2O to mol/kg-sw
           aksi(i,j,k,bi,bj) = EXP(
     &           117.40 _d 0 - 8904.2 _d 0*inv_t_k
     &         - 19.334 _d 0 * dlog_t_k
     &         + ( 3.5913 _d 0 -  458.79 _d 0*inv_t_k) * sqrtis
     &         + (-1.5998 _d 0 +  188.74 _d 0*inv_t_k) * ion_st
     &         + (0.07871 _d 0 - 12.1652 _d 0*inv_t_k) * ion_st*ion_st
     &         + log_kgfw2kgsw)

C -----------------------------------------------------------------------
C       Calculate the dissociation constant
C       of ammonium in sea-water [mol/kg-SW]
C       References: Yao and Millero (1995)
C       pH scale  : Seawater
C           akn(i,j,k,bi,bj)  = EXP(-0.25444 _d 0 -  6285.33 _d 0*inv_t_k
C     &                + 0.0001635 _d 0 * t_k
C     &                + ( 0.46532 _d 0 - 123.7184 _d 0*inv_t_k)*sqrts
C     &                + (-0.01992 _d 0 +  3.17556 _d 0*inv_t_k)*s)

C -----------------------------------------------------------------------
C       Calculate the dissociation constant of hydrogen sulfide in sea-water
C       References: Millero et al. (1988) (cited by Millero (1995)
C       pH scale  : - Seawater (according to Yao and Millero, 1995,
C                               p. 82: "refitted if necessary")
C                   - Total (according to Lewis and Wallace, 1998)
C       Note      : we stick to Seawater here for the time being
C       Note      : the fits from Millero (1995) and Yao and Millero (1995)
C                   derive from Millero et al. (1998), with all the coefficients
C                   multiplied by -ln(10)
C           akhs(i,j,k,bi,bj) = EXP( 225.838 _d 0 - 13275.3 _d 0*inv_t_k
C     &               - 34.6435 _d 0 * dlog_t_k
C     &               +  0.3449 _d 0*sqrts -  0.0274 _d 0*s)

C -----------------------------------------------------------------------
C       Calculate the dissociation constant of hydrogen sulfate (bisulfate)
C       References: Dickson (1990) -- also Handbook (2007)
C       pH scale  : free
C       Note      : converted here from mol/kg-H2O to mol/kg-SW
           aks(i,j,k,bi,bj) = EXP(141.328 _d 0
     &                -   4276.1 _d 0*inv_t_k -  23.093 _d 0*dlog_t_k
     &                + ( 324.57 _d 0 - 13856. _d 0*inv_t_k
     &                -   47.986 _d 0*dlog_t_k) * sqrtis
     &                + (-771.54 _d 0 + 35474. _d 0*inv_t_k
     &                +  114.723 _d 0*dlog_t_k) * ion_st
     &                - 2698. _d 0*inv_t_k*ion_st**1.5 _d 0
     &                + 1776. _d 0*inv_t_k*ion_st*ion_st
     &                + log_kgfw2kgsw)
     
C -----------------------------------------------------------------------
C       Calculate the dissociation constant \beta_{HF} [(mol/kg-SW)^{-1}]
C       in (mol/kg-SW)^{-1}, where
C         \beta_{HF} = \frac{ [HF] }{ [H^{+}] [F^{-}] }
C       References: Dickson and Riley (1979)
C       pH scale  : free
C       Note      : converted here from mol/kg-H2O to mol/kg-SW
           akf(i,j,k,bi,bj) = EXP(1590.2 _d 0*inv_t_k - 12.641 _d 0
     &                 + 1.525 _d 0*sqrtis
     &                 + log_kgfw2kgsw)           
           
C -----------------------------------------------------------------------
C       Calculate the water dissociation constant Kw in (mol/kg-SW)^2
C       References: Millero (1995) for value at pressc = 0
C       pH scale  : Seawater
           akw(i,j,k,bi,bj) =  EXP(148.9802 _d 0 - 13847.26 _d 0*inv_t_k
     &              -   23.6521 _d 0*dlog_t_k
     &              + (-5.977 _d 0 + 118.67 _d 0*inv_t_k
     &              + 1.0495 _d 0*dlog_t_k)*sqrts
     &              -   0.01615 _d 0*s )

C -----------------------------------------------------------------------
C pH scale conversion factors and conversions. Everything on total pH scale
           total2free = 1. _d 0/
     &                 (1. _d 0 + st(i,j,k,bi,bj)/aks(i,j,k,bi,bj))

           free2sw = 1. _d 0
     &                 + (st(i,j,k,bi,bj)/ aks(i,j,k,bi,bj))
     &                 + (ft(i,j,k,bi,bj)/(akf(i,j,k,bi,bj)
     &                   *total2free))

           total2sw = total2free * free2sw

           aphscale(i,j,k,bi,bj) = 1. _d 0 
     &                    + st(i,j,k,bi,bj)/aks(i,j,k,bi,bj)

      IF ( selectK1K2const.NE.2 .AND. selectK1K2const.NE.4 ) THEN
C Convert to the total pH scale
           ak1(i,j,k,bi,bj)  = ak1(i,j,k,bi,bj)/total2sw
           ak2(i,j,k,bi,bj)  = ak2(i,j,k,bi,bj)/total2sw
      ENDIF
           ak1p(i,j,k,bi,bj) = ak1p(i,j,k,bi,bj)/total2sw
           ak2p(i,j,k,bi,bj) = ak2p(i,j,k,bi,bj)/total2sw
           ak3p(i,j,k,bi,bj) = ak3p(i,j,k,bi,bj)/total2sw
           aksi(i,j,k,bi,bj) = aksi(i,j,k,bi,bj)/total2sw
C           akn (i,j,k,bi,bj) = akn (i,j,k,bi,bj)/total2sw
C           akhs(i,j,k,bi,bj) = akhs(i,j,k,bi,bj)/total2sw
           aks (i,j,k,bi,bj) = aks (i,j,k,bi,bj)/total2free
           akf (i,j,k,bi,bj) = akf (i,j,k,bi,bj)/total2free
           akw (i,j,k,bi,bj) = akw (i,j,k,bi,bj)/total2free
           
C -----------------------------------------------------------------------
C       Calculate the stoichiometric solubility product
C       of calcite in seawater
C       References: Mucci (1983)
C       pH scale  : N/A
C       Units     : (mol/kg-SW)^2

           Ksp_TP_Calc(i,j,k,bi,bj) = 10. _d 0**(-171.9065 _d 0
     &             - 0.077993 _d 0*t_k
     &             + 2839.319 _d 0*inv_t_k + 71.595 _d 0*dlog_t_k
     &             + ( -0.77712 _d 0 + 0.0028426 _d 0*t_k
     &             + 178.34 _d 0*inv_t_k)*sqrts
     &             - 0.07711 _d 0*s + 0.0041249 _d 0*s_15)

C -----------------------------------------------------------------------
C       Calculate the stoichiometric solubility product
C       of aragonite in seawater
C       References: Mucci (1983)
C       pH scale  : N/A
C       Units     : (mol/kg-SW)^2

           Ksp_TP_Arag(i,j,k,bi,bj) = 10. _d 0**(-171.945 _d 0
     &             - 0.077993 _d 0*t_k
     &             + 2903.293 _d 0*inv_t_k + 71.595 _d 0*dlog_t_k
     &             + ( -0.068393 _d 0 + 0.0017276 _d 0*t_k
     &             + 88.135 _d 0*inv_t_k)*sqrts
     &             - 0.10018 _d 0*s + 0.0059415 _d 0*s_15)
         ELSE
           bt(i,j,k,bi,bj)  = 0. _d 0
           st(i,j,k,bi,bj)  = 0. _d 0
           ft(i,j,k,bi,bj)  = 0. _d 0
           cat(i,j,k,bi,bj) = 0. _d 0
           fugf(i,j,k,bi,bj)= 0. _d 0
           ff(i,j,k,bi,bj)  = 0. _d 0
           ak0(i,j,k,bi,bj) = 0. _d 0
           ak1(i,j,k,bi,bj) = 0. _d 0
           ak2(i,j,k,bi,bj) = 0. _d 0
           akb(i,j,k,bi,bj) = 0. _d 0
           ak1p(i,j,k,bi,bj)= 0. _d 0
           ak2p(i,j,k,bi,bj)= 0. _d 0
           ak3p(i,j,k,bi,bj)= 0. _d 0
           aksi(i,j,k,bi,bj)= 0. _d 0
           akw(i,j,k,bi,bj) = 0. _d 0
           aks(i,j,k,bi,bj) = 0. _d 0
           akf(i,j,k,bi,bj) = 0. _d 0
           Ksp_TP_Calc(i,j,k,bi,bj) = 0. _d 0
           Ksp_TP_Arag(i,j,k,bi,bj) = 0. _d 0
           aphscale(i,j,k,bi,bj)    = 0. _d 0
         ENDIF
         END DO
        END DO
      RETURN
      END
C     END SUBROUTINE CARBON_COEFFS_SURF
C -----------------------------------------------------------------------           
           
CBOP
C !ROUTINE: CARBON_COEFFS_DEPTH

C !INTERFACE: ==========================================================
      SUBROUTINE CARBON_COEFFS_DEPTH(
     I                   bi,bj,iMin,iMax,jMin,jMax,
     I                   k,myThid)

C     Add depth dependence to carbon chemistry coefficients loaded into
C     common block. Corrections occur on the Seawater pH scale and are
C     converted back to the total scale
      IMPLICIT NONE

C     MITgcm GLobal variables
#include "SIZE.h"
#include "DYNVARS.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "FFIELDS.h"
#include "DIC_VARS.h"

      INTEGER bi,bj,iMin,iMax,jMin,jMax,k
C K is depth index
        INTEGER myThid
CEOP

C LOCAL VARIABLES
      INTEGER i, j, kLev
      _RL  bdepth
      _RL  cdepth
      _RL  pressc
      _RL t
      _RL s
      _RL zds
      _RL t_k
      _RL t_k_o_100
      _RL t_k_o_100_2
      _RL dlog_t_k
      _RL sqrtis
      _RL sqrts
      _RL s_2
      _RL s_15
      _RL inv_t_k
      _RL ion_st
      _RL is_2
      _RL scl
      _RL zrt
      _RL B1
      _RL B
      _RL delta
      _RL Pzatm
      _RL zdvi
      _RL zdki
      _RL pfac
C pH scale converstions
      _RL total2free_surf
      _RL free2sw_surf
      _RL total2sw_surf
      _RL total2free
      _RL free2sw
      _RL total2sw

c determine pressure (bar) from depth
c 1 BAR at z=0m (atmos pressure)
c use UPPER surface of cell so top layer pressure = 0 bar
c for surface exchange coeffs

cmick..............................
c        write(6,*)'k ',k

        bdepth = 0. _d 0
        cdepth = 0. _d 0
        pressc = 1. _d 0
        DO kLev = 1,k
            cdepth = bdepth + 0.5 _d 0*drF(kLev)
            bdepth = bdepth + drF(kLev)
            pressc = 1. _d 0 + 0.1 _d 0*cdepth
        END DO
cmick...................................................
c        write(6,*)'depth,pressc ',cdepth,pressc
cmick....................................................

        DO i=imin,imax
         DO j=jmin,jmax
          IF (hFacC(i,j,k,bi,bj).gt.0. _d 0) THEN
           t = theta(i,j,k,bi,bj)
           s = salt (i,j,k,bi,bj)
C terms used more than once for:
C temperature
           t_k = 273.15 _d 0 + t
           zrt= 83.14472 _d 0 * t_k
           t_k_o_100 = t_k/100. _d 0
           t_k_o_100_2=t_k_o_100*t_k_o_100
           inv_t_k=1.0 _d 0/t_k
           dlog_t_k=log(t_k)
C ionic strength
           ion_st=19.924 _d 0*s/(1000. _d 0-1.005 _d 0*s)
           is_2=ion_st*ion_st
           sqrtis=sqrt(ion_st)
C salinity
           s_2=s*s
           sqrts=sqrt(s)
           s_15=s*sqrts
           scl=s/1.80655 _d 0
           zds = s - 34.8 _d 0

           total2free_surf = 1. _d 0/
     &                 (1. _d 0 + st(i,j,k,bi,bj)/aks(i,j,k,bi,bj))

           free2sw_surf = 1. _d 0
     &                 + st(i,j,k,bi,bj)/ aks(i,j,k,bi,bj)
     &                 + ft(i,j,k,bi,bj)/(akf(i,j,k,bi,bj)
     &                  *total2free_surf)

           total2sw_surf = total2free_surf * free2sw_surf

C------------------------------------------------------------------------
C       Recalculate Fugacity Factor needed for non-ideality in ocean
C       with pressure dependence.
C       Reference : Weiss (1974) Marine Chemistry
C       pH scale  : N/A
           Pzatm = 1.01325 _d 0+pressc ! bars
           delta = (57.7 _d 0 - 0.118 _d 0*t_k)
           B1 = -1636.75 _d 0 + 12.0408 _d 0*t_k - 0.0327957 _d 0*t_k*t_k
           B  = B1 + 3.16528 _d 0*t_k*t_k*t_k*(0.00001 _d 0)
C   "x2" term often neglected (assumed=1) in applications of Weiss's (1974) eq.9
C    x2 = 1 - x1 = 1 - xCO2 (it is very close to 1, but not quite)
           fugf(i,j,k,bi,bj) = exp( (B+2. _d 0*delta) * Pzatm / zrt)

C -----------------------------------------------------------------------
C       Apply pressure dependence to the dissociation constant of hydrogen
C       sulfate (bisulfate).  Ref: Millero (1995) for pressure correction
           zdvi         =  -18.03 _d 0 + t*(0.0466 _d 0 + t*0.316 _d -3)
           zdki         = ( -4.53 _d 0 + t*0.0900 _d 0)*1. _d -3
           pfac = (-zdvi + zdki*pressc/2. _d 0)*pressc/zrt
           aks(i,j,k,bi,bj) = total2free_surf*aks(i,j,k,bi,bj) 
     &                                       * exp(pfac)

           total2free = 1. _d 0/
     &                 (1. _d 0 + st(i,j,k,bi,bj)/aks(i,j,k,bi,bj))

C free2sw has an additional component from fluoride
           free2sw    = 1. _d 0
     &                 + st(i,j,k,bi,bj)/ aks(i,j,k,bi,bj)

           aks(i,j,k,bi,bj) = aks(i,j,k,bi,bj)/total2free

           aphscale(i,j,k,bi,bj) = 1. _d 0 
     &                    + st(i,j,k,bi,bj)/aks(i,j,k,bi,bj)
C -----------------------------------------------------------------------
C       Apply pressure dependence to dissociation constant for hydrogen fluoride
C       References: Millero (1995) for pressure correction
           zdvi   =   -9.78 _d 0 + t*(-0.0090 _d 0 - t*0.942 _d -3)
           zdki   = ( -3.91 _d 0 + t*0.054 _d 0)*1. _d -3
           pfac   = (-zdvi + zdki*pressc/2. _d 0)*pressc/zrt
           akf(i,j,k,bi,bj) = total2free_surf*akf(i,j,k,bi,bj) 
     &                                       * exp(pfac)

           free2sw = free2sw
     &               + ft(i,j,k,bi,bj)/akf(i,j,k,bi,bj)
           total2sw = total2free * free2sw

           akf(i,j,k,bi,bj) = akf(i,j,k,bi,bj)/total2free

C -----------------------------------------------------------------------
C       Apply pressure dependence to 1rst dissociation constant of carbonic acid
C       References: Millero (1982) pressure correction
           zdvi =  -25.50 _d 0 -0.151 _d 0*zds + 0.1271 _d 0*t
           zdki = ( -3.08 _d 0 -0.578 _d 0*zds + 0.0877 _d 0*t)*1. _d -3
           pfac = (-zdvi + zdki*pressc/2. _d 0)*pressc/zrt
           ak1(i,j,k,bi,bj) = (total2free_surf*ak1(i,j,k,bi,bj)
     &                      * exp(pfac))/total2sw

C -----------------------------------------------------------------------
C       Apply pressure dependence to 2nd dissociation constant of carbonic acid
C       References: Millero (1979) pressure correction
           zdvi = -15.82 _d 0 + 0.321 _d 0*zds - 0.0219 _d 0*t
           zdki = ( 1.13 _d 0 - 0.314 _d 0*zds - 0.1475 _d 0*t)*1. _d -3
           pfac = (-zdvi + zdki*pressc/2. _d 0)*pressc/zrt
           ak2(i,j,k,bi,bj) = (total2free_surf*ak2(i,j,k,bi,bj)
     &                      * exp(pfac))/total2sw

C -----------------------------------------------------------------------
C       Apply pressure dependence to the boric acid dissociation constant KB
C       References: Millero (1979) pressure correction
           zdvi      = -29.48 _d 0 + 0.295 _d 0*zds + 0.1622 _d 0*t
     &                 - 0.002608 _d 0*t*t
           zdki      = (-2.84 _d 0 + 0.354 _d 0*zds)*1. _d -3
           pfac =  (-zdvi + zdki*pressc/2. _d 0)*pressc/zrt
           akb(i,j,k,bi,bj) = (total2free_surf*akb(i,j,k,bi,bj)
     &                      * exp(pfac))/total2sw

C -----------------------------------------------------------------------
C       Apply pressure dependence to water dissociation constant Kw in
C       (mol/kg-SW)^2. Ref.: Millero (pers. comm. 1996) for pressure correction
           zdvi    =  -20.02 _d 0 + 0.1119 _d 0*t - 0.1409 _d -2*t*t
           zdki    = ( -5.13 _d 0 + 0.0794 _d 0*t)*1. _d -3
           pfac = (-zdvi + zdki*pressc/2. _d 0)*pressc/zrt
           akw(i,j,k,bi,bj) = (total2free_surf*akw(i,j,k,bi,bj)
     &                      * exp(pfac))/total2sw

C -----------------------------------------------------------------------
C       Apply pressure dependence to the first dissociation constant
C       of phosphoric acid (H3PO4) in seawater
C       References: Millero (1995) for pressure correction
           zdvi      =  -14.51 _d 0 + 0.1211 _d 0*t - 0.321 _d -3*t*t
           zdki      = ( -2.67 _d 0 + 0.0427 _d 0*t)*1. _d -3
           pfac = (-zdvi + zdki*pressc/2. _d 0)*pressc/zrt
           ak1p(i,j,k,bi,bj) = (total2free_surf*ak1p(i,j,k,bi,bj)
     &                       * exp(pfac))/total2sw

C -----------------------------------------------------------------------
C       Apply pressure dependence to the second dissociation constant
C       of phosphoric acid (H3PO4) in seawater
C       References: Millero (1995) for pressure correction
           zdvi       =  -23.12 _d 0 + 0.1758 _d 0*t -2.647 _d -3*t*t
           zdki       = ( -5.15 _d 0 +   0.09 _d 0*t)*1. _d -3
           pfac = (-zdvi + zdki*pressc/2. _d 0)*pressc/zrt
           ak2p(i,j,k,bi,bj) = (total2free_surf*ak2p(i,j,k,bi,bj)
     &                       * exp(pfac))/total2sw

C -----------------------------------------------------------------------
C       Apply pressure dependence to the third dissociation constant
C       of phosphoric acid (H3PO4) in seawater
C       References: Millero (1995) for pressure correction
           zdvi      =  -26.57 _d 0 + 0.2020 _d 0*t -3.042 _d -3*t*t
           zdki      = ( -4.08 _d 0 + 0.0714 _d 0*t)*1. _d -3
           pfac = (-zdvi + zdki*pressc/2. _d 0)*pressc/zrt
           ak3p(i,j,k,bi,bj) = (total2free_surf*ak3p(i,j,k,bi,bj)
     &                       * exp(pfac))/total2sw

C -----------------------------------------------------------------------
C       Apply pressure dependence to the first dissociation constant
C       of silicic acid (H4SiO4) in seawater
C       References: Millero (1979) pressure correction. Note: Pressure
C        correction estimated to be the same as borate (Millero, 1995)
           zdvi      = -29.48 _d 0 + 0.295 _d 0*zds + 0.1622 _d 0*t
     &                 - 0.002608 _d 0*t*t
           zdki      = (-2.84 _d 0 + 0.354 _d 0*zds)*1. _d -3
           pfac =  (-zdvi + zdki*pressc/2. _d 0)*pressc/zrt
           aksi(i,j,k,bi,bj) = (total2free_surf*aksi(i,j,k,bi,bj)
     &                       * exp(pfac))/total2sw

C -----------------------------------------------------------------------
C       Apply pressure dependence to the dissociation constant of hydrogen
C       sulfide in sea-water
C       References: Millero (1995) for pressure correction
C           zdvi         =  -14.80 _d 0 + t*(0.0020 _d 0 - t*0.400 _d -3)
C           zdki         = (  2.89 _d 0 + t*0.054 _d 0)*1. _d -3
C           pfac  = (-zdvi + zdki*pressc/2. _d 0)*pressc/zrt
C           akhs(i,j,k,bi,bj) = (total2free_surf*akhs(i,j,k,bi,bj)
C     &                       * exp(pfac))/total2sw

C -----------------------------------------------------------------------
C       Apply pressure dependence to the dissociation constant
C       of ammonium in sea-water [mol/kg-SW]
C       References: Millero (1995) for pressure correction
C           zdvi         =  -26.43 _d 0 + t*(0.0889 _d 0 - t*0.905 _d -3)
C           zdki         = ( -5.03 _d 0 + t*0.0814 _d 0)*1. _d -3
C           pfac  = (-zdvi + zdki*pressc/2. _d 0)*pressc/zrt
C           akn(i,j,k,bi,bj) = (total2free_surf*akn(i,j,k,bi,bj)
C     &                      * exp(pfac))/total2sw

C -----------------------------------------------------------------------
C       Apply pressure dependence to the stoichiometric solubility product
C       of calcite in seawater
C       References: Millero (1995) for pressure correction
           zdvi      =  -48.76 _d 0 + 0.5304 _d 0*t
           zdki      = (-11.76 _d 0 + 0.3692 _d 0*t)*1. _d -3
           pfac = (-zdvi + zdki*pressc/2. _d 0)*pressc/zrt
           Ksp_TP_Calc(i,j,k,bi,bj) = Ksp_TP_Calc(i,j,k,bi,bj) 
     &                                      * exp(pfac)

C -----------------------------------------------------------------------
C       Apply pressure dependence to the stoichiometric solubility product
C       of aragonite in seawater
C       References: Millero (1979) for pressure correction
           zdvi      =  -48.76 _d 0 + 0.5304 _d 0*t  + 2.8 _d 0
           zdki      = (-11.76 _d 0 + 0.3692 _d 0*t)*1. _d -3
           pfac = (-zdvi + zdki*pressc/2. _d 0)*pressc/zrt
           Ksp_TP_Arag(i,j,k,bi,bj) = Ksp_TP_Arag(i,j,k,bi,bj)
     &                                      * exp(pfac)
         ELSE
           bt(i,j,k,bi,bj)  = 0. _d 0
           st(i,j,k,bi,bj)  = 0. _d 0
           ft(i,j,k,bi,bj)  = 0. _d 0
           cat(i,j,k,bi,bj) = 0. _d 0
           fugf(i,j,k,bi,bj)= 0. _d 0
           ff(i,j,k,bi,bj)  = 0. _d 0
           ak0(i,j,k,bi,bj) = 0. _d 0
           ak1(i,j,k,bi,bj) = 0. _d 0
           ak2(i,j,k,bi,bj) = 0. _d 0
           akb(i,j,k,bi,bj) = 0. _d 0
           ak1p(i,j,k,bi,bj)= 0. _d 0
           ak2p(i,j,k,bi,bj)= 0. _d 0
           ak3p(i,j,k,bi,bj)= 0. _d 0
           aksi(i,j,k,bi,bj)= 0. _d 0
           akw(i,j,k,bi,bj) = 0. _d 0
           aks(i,j,k,bi,bj) = 0. _d 0
           akf(i,j,k,bi,bj) = 0. _d 0
           Ksp_TP_Calc(i,j,k,bi,bj) = 0. _d 0
           Ksp_TP_Arag(i,j,k,bi,bj) = 0. _d 0
           aphscale(i,j,k,bi,bj)    = 0. _d 0
         ENDIF
       ENDDO
      ENDDO
      RETURN
      END
C     END SUBROUTINE CARBON_COEFFS_DEPTH
C      =========================================================================
#endif /* DIC_CARBON_COMPONENTS */
