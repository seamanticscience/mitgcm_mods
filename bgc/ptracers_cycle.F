C $Header: /u/gcmpack/MITgcm/pkg/ptracers/ptracers_cycle.F,v 1.7 2012/03/08 17:13:31 jmc Exp $
C $Name:  $

#include "PTRACERS_OPTIONS.h"

CBOP
C !ROUTINE: PTRACERS_CYCLE

C !INTERFACE: ==========================================================
      SUBROUTINE PTRACERS_CYCLE( bi,bj,myTime,myIter,myThid )

C !DESCRIPTION:
C     Cycle time-stepping arrays

C !USES: ===============================================================
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "PTRACERS_SIZE.h"
#include "PTRACERS_PARAMS.h"
#include "PTRACERS_START.h"
#include "PTRACERS_FIELDS.h"
#include "GAD.h"

C !INPUT PARAMETERS: ===================================================
C  bi,bj                :: tile indices
C  myTime               :: model time
C  myIter               :: time-step number
C  myThid               :: thread number
      INTEGER i,j,k,bi,bj      
      _RL myTime
      INTEGER myIter
      INTEGER myThid

C !OUTPUT PARAMETERS: ==================================================
C  none

#ifdef ALLOW_PTRACERS
#ifdef ALLOW_DIAGNOSTICS
C     !FUNCTIONS:
      LOGICAL  DIAGNOSTICS_IS_ON
      EXTERNAL DIAGNOSTICS_IS_ON
      CHARACTER*4 GAD_DIAG_SUFX
      EXTERNAL    GAD_DIAG_SUFX
#endif /* ALLOW_DIAGNOSTICS */

C !LOCAL VARIABLES: ====================================================
C  iTracer       :: loop indices
      INTEGER iTracer
      INTEGER GAD_TR
      CHARACTER*8 diagName
      CHARACTER*4 diagSufx
      LOGICAL diag_tend
      _RL tmpMk(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
CEOP

C Loop over passive tracers
      DO iTracer=1,PTRACERS_numInUse
       IF ( PTRACERS_StepFwd(iTracer) ) THEN
#ifdef ALLOW_DIAGNOSTICS
C JML at this point, pTracer array contains previous tracer concentration, while
C    gPtr contains the update tracer concentration due to explicit AND implicit advection
C BUT NOT Convective adjustment if ivdc_kappa>0 (i.e. not using KPP).
        GAD_TR = GAD_TR1 + iTracer - 1

        IF ( useDiagnostics ) THEN
          diagSufx = GAD_DIAG_SUFX( GAD_TR, myThid )
          diagName = 'g'//diagSufx
          diag_tend = DIAGNOSTICS_IS_ON( diagName, myThid )
        ENDIF
        
        IF ( DIAGNOSTICS_IS_ON(diagName,myThid) ) THEN
           DO k=1,Nr
            DO j = 1,sNy
             DO i = 1,sNx
               tmpMk(i,j,k,bi,bj) = (gPtr(i,j,k,bi,bj,iTracer)
     &                              -pTracer(i,j,k,bi,bj,iTracer))          
     &                              /dTtracerLev(k)
             ENDDO
            ENDDO
           ENDDO
           CALL DIAGNOSTICS_FILL(tmpMk,diagName,0,Nr,1,bi,bj,myThid)
        ENDIF
#endif /* ALLOW_DIAGNOSTICS */

C Cycle time-stepping arrays
        CALL CYCLE_TRACER(
     I             bi,bj,
     U             pTracer(1-OLx,1-OLy,1,1,1,iTracer),
     U             gPtr(1-OLx,1-OLy,1,1,1,iTracer),
     I             myTime,myIter,myThid )

C End of tracer loop
       ENDIF
      ENDDO

#endif /* ALLOW_PTRACERS */

      RETURN
      END
