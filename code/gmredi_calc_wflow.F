C $Header: /u/gcmpack/MITgcm/pkg/gmredi/gmredi_calc_wflow.F,v 1.5 2008/05/30 02:50:16 gforget Exp $
C $Name:  $

#include "GMREDI_OPTIONS.h"

CBOP
C     !ROUTINE: GMREDI_CALC_WFLOW
C     !INTERFACE:
      SUBROUTINE GMREDI_CALC_WFLOW(
     U                       wFld, rTrans,
     I                       k, bi, bj, myThid )
C     !DESCRIPTION:
C     Add GM-bolus flow to Eulerian vertical transport.

C     !USES:
      IMPLICIT NONE

C     == GLobal variables ==
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "GMREDI.h"

C     !INPUT/OUTPUT PARAMETERS:
C     == Routine arguments ==
C     wFld   :: vertical volume transport (updated)
C     rTrans :: vertical volume transport (updated)
C     k      :: level index
C     bi,bj  :: tile indices
C     myThid :: thread number
      _RL wFld  (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL wEddy (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL rTrans(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL eTrans(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      INTEGER k, bi, bj, myThid

#ifdef ALLOW_GMREDI
#ifdef GM_BOLUS_ADVEC
C     !LOCAL VARIABLES:
C     == Local variables ==
C     i, j :: loop indices
      INTEGER i, j
      _RL delPsi, tmpFac0, tmpFac1
CEOP
      tmpFac0=0. _d 0
      tmpFac1=1. _d 0

#ifdef ALLOW_DEBUG
      IF (debugMode) CALL DEBUG_ENTER('GMREDI_CALC_WFLOW',myThid)
#endif
         
      IF (GM_AdvForm .AND. .NOT.GM_AdvSeparate
     & .AND. .NOT.GM_InMomAsStress) THEN

       DO j=1-Oly,sNy+Oly-1
        DO i=1-Olx,sNx+Olx-1
         delPsi = ( dyG(i+1,j,bi,bj)*GM_PsiX(i+1,j,k,bi,bj)
     &             -dyG( i ,j,bi,bj)*GM_PsiX( i ,j,k,bi,bj)
     &             +dxG(i,j+1,bi,bj)*GM_PsiY(i,j+1,k,bi,bj)
     &             -dxG(i, j ,bi,bj)*GM_PsiY(i, j ,k,bi,bj)
     &            )*maskC(i,j,k,bi,bj)
         wEddy(i,j)  = delPsi*recip_rA(i,j,bi,bj)
         wFld(i,j)   = wFld(i,j) + wEddy(i,j)
         eTrans(i,j) = delPsi
         rTrans(i,j) = rTrans(i,j) + delPsi
        ENDDO
       ENDDO
       
       IF (useDiagnostics) THEN
        IF (k .EQ. 2) THEN
C Hack to fill surface, even though this routine isnt called when k=1
         CALL DIAGNOSTICS_SCALE_FILL(wFld  ,tmpFac0,tmpFac1,
     &                                'GM_W_RES',1,1,2,bi,bj,myThid)
         CALL DIAGNOSTICS_SCALE_FILL(wEddy ,tmpFac0,tmpFac1,
     &                                'GM_W_EDD',1,1,2,bi,bj,myThid)
         CALL DIAGNOSTICS_SCALE_FILL(rTrans,tmpFac0,tmpFac1,
     &                                'GM_WR_TR',1,1,2,bi,bj,myThid)
         CALL DIAGNOSTICS_SCALE_FILL(eTrans,tmpFac0,tmpFac1,
     &                                'GM_WE_TR',1,1,2,bi,bj,myThid)
        ENDIF

        CALL DIAGNOSTICS_FILL(wFld,  'GM_W_RES',k,1,2,bi,bj,myThid)
        CALL DIAGNOSTICS_FILL(wEddy, 'GM_W_EDD',k,1,2,bi,bj,myThid) 
        CALL DIAGNOSTICS_FILL(rTrans,'GM_WR_TR',k,1,2,bi,bj,myThid)
        CALL DIAGNOSTICS_FILL(eTrans,'GM_WE_TR',k,1,2,bi,bj,myThid)
       ENDIF 
      ENDIF     
#endif /* GM_BOLUS_ADVEC */
#endif /* ALLOW_GMREDI */

      RETURN
      END
