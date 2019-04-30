C $Header: /u/gcmpack/MITgcm/pkg/salt_plume/salt_plume_forcing_surf.F,v 1.4 2010/07/06 00:14:29 dimitri Exp $
C $Name:  $

#include "SALT_PLUME_OPTIONS.h"

CBOP
C     !ROUTINE: SALT_PLUME_FORCING_SURF_PTR
C     !INTERFACE:
      SUBROUTINE SALT_PLUME_FORCING_SURF_PTR(
     I                            bi, bj, iMin, iMax, jMin, jMax,
     I                            myTime,myIter,myThid )

C     !DESCRIPTION: \bv
C     *==========================================================*
C     | S/R SALT_PLUME_FORCING_SURF
C     | o saltPlume is the amount of salt rejected by ice while freezing;
C     |   Some of this will also be rejected tracers and so this is 
C     |   subtracted from surfaceForcingPTR and will be redistributed
C     |   to multiple vertical levels later on as per Duffy et al. (GRL 1999)
C     *==========================================================*
C     \ev

C     !USES:
      IMPLICIT NONE
C     == Global variables ==
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "FFIELDS.h"
#include "DYNVARS.h"
#include "SALT_PLUME.h"
#include "PTRACERS_SIZE.h"
#include "PTRACERS_PARAMS.h"
#include "PTRACERS_FIELDS.h"
#include "PTRACERS_START.h"

C     !INPUT PARAMETERS:
C     bi,bj                :: tile indices
C     myTime               :: model time
C     myIter               :: time-step number
C     myThid               :: thread number
      INTEGER bi, bj, iMin, iMax, jMin, jMax, iTrc
      _RL myTime
      _RL plumeflux(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      INTEGER myIter
      INTEGER myThid
#ifdef ALLOW_DIAGNOSTICS
      CHARACTER*8 diagName
      LOGICAL  DIAGNOSTICS_IS_ON
      EXTERNAL DIAGNOSTICS_IS_ON
#endif
CEOP

#if defined (ALLOW_SALT_PLUME) && (defined SALT_PLUME_DIC)

C     !LOCAL VARIABLES:
C     i,j                  :: loop indices
C     ks                   :: surface level index
      INTEGER i, j, ks
      IF ( usingPCoords ) THEN
         ks = Nr
      ELSE
         ks = 1
      ENDIF

      DO iTrc=1,PTRACERS_numInUse
       IF (convertFW2Salt .EQ. -1.) THEN
C-    use local surface tracer field to calculate forcing term:
        DO j = jMin, jMax
         DO i = iMin, iMax
          IF (maskC(i,j,ks,bi,bj).NE.0.0 _d 0) THEN
           plumeflux(i,j)=saltPlumeFlux(i,j,bi,bj) * mass2rUnit
     &     *( pTracer(i,j,ks,bi,bj,iTrc)/SALT(i,j,ks,bi,bj))

           surfaceForcingPTr(i,j,bi,bj,iTrc) =
     &         surfaceForcingPTr(i,j,bi,bj,iTrc)
     &       - plumeflux(i,j)
          ELSE
           plumeflux(i,j)=0. _d 0
          ENDIF
         ENDDO
        ENDDO
       ELSE
C-    use uniform tracer value to calculate forcing term:
        DO j = jMin, jMax
         DO i = iMin, iMax
          IF (maskC(i,j,ks,bi,bj).NE.0.0 _d 0) THEN
           plumeflux(i,j)=saltPlumeFlux(i,j,bi,bj) * mass2rUnit
     &       *( pTracer(i,j,ks,bi,bj,iTrc)/convertFW2Salt )
     
           surfaceForcingPTr(i,j,bi,bj,iTrc) =
     &         surfaceForcingPTr(i,j,bi,bj,iTrc)
     &       - plumeflux(i,j)
          ELSE
           plumeflux(i,j)=0. _d 0
          ENDIF
         ENDDO
        ENDDO 
C-    end local-surface-tracer / uniform-value distinction
       ENDIF
       
#ifdef ALLOW_DIAGNOSTICS
       diagName = '        '
       WRITE(diagName,'(A6,A2)') 'SPflxT',PTRACERS_ioLabel(iTrc)
       IF ( DIAGNOSTICS_IS_ON( diagName , myThid ) ) THEN
        CALL DIAGNOSTICS_FILL( plumeflux,diagName,1,1,2,bi,bj,myThid )
       ENDIF
#endif /* ALLOW_DIAGNOSTICS */       
      ENDDO
#endif /* ALLOW_SALT_PLUME */

      RETURN
      END
