C $Header: /u/gcmpack/MITgcm/pkg/dic/dic_surfforcing_init.F,v 1.33 2011/10/07 21:36:39 dfer Exp $
C $Name:  $

#include "DIC_OPTIONS.h"
#include "PTRACERS_OPTIONS.h"

CBOP
C !ROUTINE: DIC_SURFFORCING_INIT

C !INTERFACE: ==========================================================
      SUBROUTINE DIC_SURFFORCING_INIT(
     I          myThid)

C !DESCRIPTION:
C  Calculate first guess of pH

C !USES: ===============================================================
      IMPLICIT NONE
#include "SIZE.h"
#include "DYNVARS.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "FFIELDS.h"
#include "DIC_VARS.h"
#include "PTRACERS_SIZE.h"
#include "PTRACERS_PARAMS.h"
#include "PTRACERS_FIELDS.h"
#include "DIC_LOAD.h"

C !INPUT PARAMETERS: ===================================================
C  myThid               :: thread number
      INTEGER  myThid

#ifdef ALLOW_DIC

C !LOCAL VARIABLES: ====================================================
       INTEGER i,j,k,kLev,pco2_init_iter
       INTEGER intimeP, intime0, intime1
       _RL aWght, bWght, co3dummy
C Number of iterations for pCO2 solvers...
C Solubility relation coefficients
C local variables for carbon chem
      INTEGER iMin,iMax,jMin,jMax, bi, bj
      _RL surfalk(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL surfphos(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL surfsi(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL surftemp(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL surfsalt(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL surfdic(1-OLx:sNx+OLx,1-OLy:sNy+OLy) 
      INTEGER iprt,jprt
      LOGICAL pH_isLoaded
CEOP

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      kLev=1

      CALL DIC_ATMOS(0, startTime, nIter0, myThid )
      
#ifdef DIC_CARBON_COMPONENTS
        IF ( dic_pCO2_ref .EQ. UNSET_RL ) THEN
C capture initial pCO2 or set from specified reference
         DO bj = myByLo(myThid), myByHi(myThid)
          DO bi = myBxLo(myThid), myBxHi(myThid)
           DO j=1-Oly,sNy+Oly
            DO i=1-Olx,sNx+Olx
              dic_pCO2_init(i,j,bi,bj)=AtmospCO2(i,j,bi,bj)
            ENDDO
           ENDDO
          ENDDO
         ENDDO
        ELSE
         DO bj = myByLo(myThid), myByHi(myThid)
          DO bi = myBxLo(myThid), myBxHi(myThid)
           DO j=1-Oly,sNy+Oly
            DO i=1-Olx,sNx+Olx
              dic_pCO2_init(i,j,bi,bj)=dic_pCO2_ref
            ENDDO
           ENDDO
          ENDDO
         ENDDO         
        ENDIF
#endif /* DIC_CARBON_COMPONENTS */

#ifdef ALLOW_C13
C Fill total_atmos_c13, atmos_13pco2
		IF ( nIter0.GT.PTRACERS_Iter0 .OR.
	 &       (nIter0.EQ.PTRACERS_Iter0 .AND. pickupSuff.NE.' ')
	 &      ) THEN
            CALL DIC_READ_13C_PICKUP( nIter0, myThid )
        ELSE
C Initialize from dic_p13co2 supplied in data.dic
         DO bj = myByLo(myThid), myByHi(myThid)
          DO bi = myBxLo(myThid), myBxHi(myThid)
           DO j=1-Oly,sNy+Oly
            DO i=1-Olx,sNx+Olx
              total_atmos_c13(i,j,bi,bj)=dic_p13co2*1.77 _d 20
              atmos_13pco2   (i,j,bi,bj)=dic_p13co2
            ENDDO
           ENDDO
          ENDDO
         ENDDO
        ENDIF
#endif /* ALLOW_C13 */

ccccccccccccccccccccccccccccccccccccccccc
      IF ( periodicExternalForcing ) THEN
#ifndef DIC_CARBON_COMPONENTS
c read in silica field
         CALL LEF_ZERO( silica0,myThid )
         CALL LEF_ZERO( silica1,myThid )

C--   Now calculate whether it is time to update the forcing arrays
       CALL GET_PERIODIC_INTERVAL(
     O                   intimeP, intime0, intime1, bWght, aWght,
     I                   externForcingCycle, externForcingPeriod,
     I                   deltaTclock, startTime, myThid )

       _BARRIER
       _BEGIN_MASTER(myThid)
        WRITE(standardMessageUnit,'(A,I10,A,2(2I5,A))')
     &   ' DIC_SURFFORCING_INIT, it=', nIter0,
     &   ' : Reading new data, i0,i1=', intime0, intime1
       _END_MASTER(myThid)


       IF ( DIC_silicaFile .NE. ' '  ) THEN
         CALL READ_REC_XY_RS( DIC_silicaFile,silica0,intime0,
     &        nIter0,myThid )
         CALL READ_REC_XY_RS( DIC_silicaFile,silica1,intime1,
     &        nIter0,myThid )
       ENDIF

#ifdef ALLOW_OFFLINE
       IF ( useOffLine ) THEN
         CALL OFFLINE_FIELDS_LOAD( startTime, nIter0, myThid )
       ENDIF
#endif

       _EXCH_XY_RS(silica0, myThid )
       _EXCH_XY_RS(silica1, myThid )

       IF ( DIC_silicaFile .NE. ' '  ) THEN
        DO bj = myByLo(myThid), myByHi(myThid)
         DO bi = myBxLo(myThid), myBxHi(myThid)
          DO j=1-Oly,sNy+Oly
           DO i=1-Olx,sNx+Olx
             SILICA(i,j,bi,bj)= bWght*silica0(i,j,bi,bj)
     &                        + aWght*silica1(i,j,bi,bj)
           ENDDO
          ENDDO
         ENDDO
        ENDDO
       ENDIF
#endif /* DIC_CARBON_COMPONENTS */
c end periodicExternalForcing
      ENDIF

C =================================================================

      jMin=1
      jMax=sNy
      iMin=1
      iMax=sNx

      DO bj=myByLo(myThid),myByHi(myThid)
       DO bi=myBxLo(myThid),myBxHi(myThid)
         DO j=1-OLy,sNy+OLy
          DO i=1-Olx,sNx+OLx
            pH(i,j,bi,bj) = 8. _d 0
          ENDDO
         ENDDO
       ENDDO
      ENDDO

      DO bj=myByLo(myThid),myByHi(myThid)
       DO bi=myBxLo(myThid),myBxHi(myThid)
        DO j=1-OLy,sNy+OLy
         DO i=1-OLx,sNx+OLx
          ak0(i,j,bi,bj)=0. _d 0
          ak1(i,j,bi,bj)=0. _d 0
          ak2(i,j,bi,bj)=0. _d 0
          akw(i,j,bi,bj)=0. _d 0
          akb(i,j,bi,bj)=0. _d 0
          akf(i,j,bi,bj)=0. _d 0
          ak1p(i,j,bi,bj)=0. _d 0
          ak2p(i,j,bi,bj)=0. _d 0
          ak3p(i,j,bi,bj)=0. _d 0
          aksi(i,j,bi,bj)=0. _d 0
          fugf(i,j,bi,bj)=0. _d 0
          ff(i,j,bi,bj)=0. _d 0
          ft(i,j,bi,bj)=0. _d 0
          st(i,j,bi,bj)=0. _d 0
          bt(i,j,bi,bj)=0. _d 0
         ENDDO
        ENDDO
       ENDDO
      ENDDO

      pH_isLoaded = .FALSE.
      IF ( nIter0.GT.PTRACERS_Iter0 .OR.
     &    (nIter0.EQ.PTRACERS_Iter0 .AND. pickupSuff.NE.' ')
     &   ) THEN
C       Read pH from a pickup file if needed
        CALL DIC_READ_PICKUP(
     O                        pH_isLoaded,
     I                        nIter0, myThid )
      ENDIF

      DO bj=myByLo(myThid),myByHi(myThid)
       DO bi=myBxLo(myThid),myBxHi(myThid)

C determine inorganic carbon chem coefficients
        DO j=jMin,jMax
         DO i=iMin,iMax

#ifdef DIC_BIOTIC
#ifdef DIC_BOUNDS
             surfalk(i,j) = max(0.4 _d 0,
     &                          min(10. _d 0,PTRACER(i,j,klev,bi,bj,2)))
     &                          * maskC(i,j,kLev,bi,bj)
             surfphos(i,j)  = max(1.0 _d -11,
     &                          min(1._d -1, PTRACER(i,j,klev,bi,bj,3)))
     &                          * maskC(i,j,kLev,bi,bj)
#else
             surfalk(i,j) = PTRACER(i,j,klev,bi,bj,2)
     &                          * maskC(i,j,kLev,bi,bj)
             surfphos(i,j)= PTRACER(i,j,klev,bi,bj,3)
     &                          * maskC(i,j,kLev,bi,bj)
#endif
#else
             surfalk(i,j) = 2.366595 _d 0 *salt(i,j,kLev,bi,bj)/35. _d 0
     &                          * maskC(i,j,kLev,bi,bj)
             surfphos(i,j)= 5.1225 _d -4 * maskC(i,j,kLev,bi,bj)
#endif
#ifdef DIC_CARBON_COMPONENTS            
C FOR 3d NON-INTERACTIVE Si
             surfsi(i,j)   = SILICA(i,j,kLev,bi,bj)
     &                          * maskC(i,j,kLev,bi,bj)
#else
C FOR 2d NON-INTERACTIVE Si
             surfsi(i,j)   = SILICA(i,j,bi,bj)
     &                          * maskC(i,j,kLev,bi,bj)
#endif

#ifdef DIC_BOUNDS
            surftemp(i,j) = max(-4. _d 0,
     &                          min(50. _d 0, theta(i,j,kLev,bi,bj)))
            surfsalt(i,j) = max(4. _d 0,
     &                          min(50. _d 0, salt(i,j,kLev,bi,bj)))
            surfdic(i,j)  = max(0.4 _d 0,
     &                         min(10. _d 0, PTRACER(i,j,klev,bi,bj,1)))
#else
            surftemp(i,j) = theta(i,j,kLev,bi,bj)
            surfsalt(i,j) = salt(i,j,kLev,bi,bj)
            surfdic(i,j)  = PTRACER(i,j,kLev,bi,bj,1)
     &                    * maskC(i,j,kLev,bi,bj)
#endif
         ENDDO
        ENDDO

        CALL CARBON_COEFFS(
     I                       surftemp,surfsalt,
     I                       bi,bj,iMin,iMax,jMin,jMax,myThid)

C====================================================================

        IF ( .NOT.pH_isLoaded ) THEN
C set guess of pH for first step here

          WRITE(standardMessageUnit,*) 'QQ: pCO2 approximation method'
c first approximation
C$TAF LOOP = parallel
          DO j=jMin,jMax
C$TAF LOOP = parallel
           DO i=iMin,iMax
            IF ( maskC(i,j,kLev,bi,bj) .GT. 0. _d 0) THEN
C$TAF init dic_surf = static, 10
             DO pco2_init_iter=1,10
C$TAF STORE pH(i,j,bi,bj)                              = dic_surf
C$TAF STORE surfalk(i,j), surfphos(i,j), surfsi(i,j)   = dic_surf
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
     I        i,j,kLev,bi,bj,pco2_init_iter,myThid )
             ENDDO
            ENDIF
           ENDDO
          ENDDO
C          iprt = MIN(20,sNx)
C          jprt = MIN(20,sNy)
C          WRITE(standardMessageUnit,*) 'QQ first guess pH',
C     &        pH(iprt,jprt,bi,bj)
C     &        theta(iprt,jprt,1,bi,bj), salt(iprt,jprt,1,bi,bj),
C     &        surfdic(iprt,jprt), surfphos(iprt,jprt),
C     &        surfsi(iprt,jprt),surfalk(iprt,jprt)
          CALL PLOT_FIELD_XYRL(pH, 
     &         'First Guess pH (dic_surfforcing_init)', 1, myThid)
        ENDIF

C     end bi,bj loops
       ENDDO
      ENDDO

#endif /* ALLOW_DIC */
      RETURN
      END