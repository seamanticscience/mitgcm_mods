C $Header: /u/gcmpack/MITgcm/pkg/dic/calcite_saturation.F,v 1.11 2011/10/07 21:36:39 dfer Exp $
C $Name:  $

#include "DIC_OPTIONS.h"

CBOP
C !ROUTINE: CAR_FLUX

C !INTERFACE: ==========================================================
      SUBROUTINE CALCITE_SATURATION(PTR_DIC, PTR_ALK, PTR_PO4,
     I           bi,bj,imin,imax,jmin,jmax,
     I           myIter,myTime,myThid)

C !DESCRIPTION:
C  Calculate carbonate fluxes

C !USES: ===============================================================
      IMPLICIT NONE
#include "SIZE.h"
#include "DYNVARS.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "DIC_VARS.h"

C !INPUT PARAMETERS: ===================================================
C  myThid               :: thread number
C  myIter               :: current timestep
C  myTime               :: current time
C  bioac                :: biological productivity
       _RL  PTR_DIC(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
       _RL  PTR_ALK(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
       _RL  PTR_PO4(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)

      INTEGER myIter
      _RL myTime
      INTEGER myThid
       INTEGER imin, imax, jmin, jmax, bi, bj

C !OUTPUT PARAMETERS: ===================================================

#ifdef ALLOW_PTRACERS
#ifdef DIC_BIOTIC

C !LOCAL VARIABLES: ====================================================
C  i,j,k                  :: loop indices
c  ko                     :: loop-within-loop index
c depth_u, depth_l        :: depths of upper and lower interfaces
c flux_u, flux_l          :: flux through upper and lower interfaces
c zbase                   :: depth of bottom of current productive layer
       INTEGER I,J,k
#ifndef DIC_CARBON_COMPONENTS       
       _RL carbonate
       _RL calcium
#endif
       _RL silicaTEST
       _RL po4local
       _RL diclocal
       _RL alklocal
       _RL pCO2local
       _RL pHlocal
       INTEGER CO3ITER
       INTEGER CO3ITERmax
CEOP


cmick...................................................
       write(6,*)'myIter ',myIter,'  CALLED CALCITEcd_SATURATION'
c      write(6,*)'WARNING calcite_sat needs 3d silica & H0 set=7.9'
c       write(6,*)'        - & Fixed first guess of deep pH to 7.9'
cmick....................................................

c determine carbonate ion concentration through full domain
c determine calcite saturation state
       DO k=1,Nr

#ifdef DIC_CARBON_COMPONENTS
C Fill 3d values using coeffs from dic_carbon_components_chem
C Calculate the coefficients and then apply separate pressure 
C   corrections.
#ifdef ALLOW_DEBUG
          IF (debugMode) CALL DEBUG_CALL(
     &            'CARBON_COEFFS_SURF',myThid)
#endif
           IF (k .GT. 1 ) THEN
             CALL CARBON_COEFFS_SURF(
     I              bi,bj,iMin,iMax,jMin,jMax,
     I              k,myThid)
           ENDIF
       
#ifdef ALLOW_DEBUG
          IF (debugMode) CALL DEBUG_CALL(
     &           'CARBON_COEFFS_DEPTH',myThid)
#endif
           IF (k .GT. 1 ) THEN
             CALL CARBON_COEFFS_DEPTH(
     I                   bi,bj,iMin,iMax,jMin,jMax,
     I                   k,myThid)
           ENDIF
#else
        CALL CARBON_COEFFS_PRESSURE_DEP(
     I                       theta,salt,
     I                       bi,bj,iMin,iMax,jMin,jMax,
     I                       k,myThid)
#endif

        DO j=jmin,jmax
         DO i=imin,imax

           if ( hFacC(i,j,k,bi,bj) .gt. 0. _d 0 ) then
#ifndef DIC_CARBON_COMPONENTS
             calcium = 1.028 _d -2*salt(i,j,k,bi,bj)/35. _d 0
#endif

c 30 micromol = 0.03 mol m-3
             silicaTEST = 0.03 _d 0
             po4local = PTR_PO4(i,j,k)
             diclocal = PTR_DIC(i,j,k)
             alklocal = PTR_ALK(i,j,k)
             pHlocal = 7.9 _d 0

CMICK - TEMPORARY!!!!!
CMICK silica = fixed
CMICK silica = fixed
C
CMICK -DEC 04
CMICK- NOW ITERATE pH SOLVER AT DEPTH ONLY
CMICK  TO ENSURE ACCURATE ESTIMATE OF CO3 AT DEPTH
CMICK - NOTE Si STILL USING A UNIFORM DUMMY VALUE
             CO3itermax = 10
CMICK - SO NOW WE ITERATE, UPDATING THE ESTIMATE OF pH and CO3--
CMICK - SINCE WE CALL THIS FOR DEEP OCEAN INFREQUENTLY (MONTHLY?)
CMIKC - CAN AFFORD TO MAKE SEVERAL ITERATIONS...
             DO CO3iter = 1, CO3itermax
#ifdef DIC_CARBON_COMPONENTS
            CALL CALC_PCO2_APPROX(
     I        theta(i,j,k,bi,bj),salt(i,j,k,bi,bj),
     I        diclocal, po4local,silicaTEST,alklocal,
     I        ak1(i,j,k,bi,bj),ak2(i,j,k,bi,bj),
     I        ak1p(i,j,k,bi,bj),ak2p(i,j,k,bi,bj),ak3p(i,j,k,bi,bj),
     I        aks(i,j,k,bi,bj),akb(i,j,k,bi,bj),akw(i,j,k,bi,bj),
     I        aksi(i,j,k,bi,bj),akf(i,j,k,bi,bj),
     I        ak0(i,j,k,bi,bj), fugf(i,j,k,bi,bj),
     I        ff(i,j,k,bi,bj),
     I        bt(i,j,k,bi,bj),st(i,j,k,bi,bj),ft(i,j,k,bi,bj),
     U        pH_glob(i,j,k,bi,bj),pCO2local,
     U        co3(i,j,k,bi,bj),hco3(i,j,k,bi,bj),
     U        co2aq(i,j,k,bi,bj),
     I        i,j,k,bi,bj,myIter,myThid )
#else
            CALL CALC_PCO2_APPROX(
     I        surftemp(i,j),surfsalt(i,j),
     I        surfdic(i,j), surfphos(i,j),
     I        surfsi(i,j),surfalk(i,j),
     I        ak1(i,j,bi,bj),ak2(i,j,bi,bj),
     I        ak1p(i,j,bi,bj),ak2p(i,j,bi,bj),ak3p(i,j,bi,bj),
     I        aks(i,j,bi,bj),akb(i,j,bi,bj),akw(i,j,bi,bj),
     I        aksi(i,j,bi,bj),akf(i,j,bi,bj),
     I        ak0(i,j,bi,bj), fugf(i,j,bi,bj),
     I        ff(i,j,bi,bj),
     I        bt(i,j,bi,bj),st(i,j,bi,bj),ft(i,j,bi,bj),
     U        pH(i,j,bi,bj),pCO2(i,j,bi,bj),
     U        co3(i,j,bi,bj),hco3(i,j,bi,bj),co2aq(i,j,bi,bj),
     I        i,j,k,bi,bj,myIter,myThid )
#endif
c........................................................
c               if(i .eq. 76 .and. j .eq. 36  .and. k .eq. 15) then
c                 write(6,*)'Iteration, pH = ',CO3iter,pHlocal
c               endif
c........................................................
             END DO

#ifdef DIC_CARBON_COMPONENTS
             omegaC(i,j,k,bi,bj) = cat(i,j,k,bi,bj) 
     &                           * co3(i,j,k,bi,bj)*permil /
     &                          Ksp_TP_Calc(i,j,k,bi,bj)
#else             
             omegaC(i,j,k,bi,bj) = calcium * carbonate*permil /
     &                   Ksp_TP_Calc(i,j,bi,bj)
#endif


cmick...................................................
c            if(omegaC(i,j,k,bi,bj) .eq. 0.) then
c             if(i .eq. 76 .and. j .eq. 36  .and. k .eq. 15) then
c               write(6,*)'i,j,k,KS,CO3,pHCa,T,S,hfacc,omega',
c     &                 i,j,k,
c     &                 Ksp_TP_Calc(i,j,bi,bj),
c     &                 carbonate,calcium,pHlocal,
c     &                 theta(i,j,k,bi,bj),salt(i,j,k,bi,bj),
c     &                 hfacc(i,j,k,bi,bj),omegaC(i,j,k,bi,bj)
c              write(6,*)'Ksp_TP_Calc',
c     &                 Ksp_TP_Calc(i,j,bi,bj)
c               write(6,*)'dic, alk, po4 ',
c     &                 diclocal, alklocal,po4local
c               write(6,*)'k1, k2, k1p, k2p, k3p ',
c     &                 ak1(i,j,bi,bj),ak2(i,j,bi,bj),
c     &                ak1p(i,j,bi,bj),ak2p(i,j,bi,bj),ak3p(i,j,bi,bj)
c               write(6,*)'ks, kb, kw, ksi ',
c     &                aks(i,j,bi,bj),akb(i,j,bi,bj),akw(i,j,bi,bj),
c     &                aksi(i,j,bi,bj)
c               write(6,*)'akf, ff, bt, st, ft ',
c     &                akf(i,j,bi,bj),ff(i,j,bi,bj),
c     &                bt(i,j,bi,bj),st(i,j,bi,bj),ft(i,j,bi,bj)
c             end if
cmick....................................................
           else
             omegaC(i,j,k,bi,bj) = 0. _d 0
           endif

         ENDDO
        ENDDO

       ENDDO
c
#endif
#endif
       RETURN
       END
