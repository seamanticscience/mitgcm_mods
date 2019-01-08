#include "GMREDI_OPTIONS.h"

CBOP
C     !ROUTINE: GMREDI_RESIDUAL_FLOW
C     !INTERFACE:
      SUBROUTINE GMREDI_RESIDUAL_FLOW(
     U                  uFld, vFld, wFld,
     I                  bi, bj, myIter, myThid )
C     !DESCRIPTION:
C     Add GM-bolus velocity to Eulerian velocity to get Residual Mean velocity.

C     !USES:
      IMPLICIT NONE

C     == GLobal variables ==
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "GMREDI.h"
#ifdef ALLOW_EDDYPSI
# include "DYNVARS.h"
# include "FFIELDS.h"
#endif

C     !INPUT/OUTPUT PARAMETERS:
C     == Routine arguments ==
C     uFld   :: zonal      velocity (updated)
C     vFld   :: meridional velocity (updated)
C     wFld   :: vertical volume transport (updated)
C     bi,bj  :: tile indices
C     myIter :: my Iteration number
C     myThid :: my Thread Id number
      INTEGER bi, bj, myIter, myThid
      _RL uFld  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL vFld  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL wFld  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL uEddy (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL vEddy (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL wEddy (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL uETrans(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL vETrans(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL wETrans(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RS xA    (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RS yA    (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RS maskUp(1-OLx:sNx+OLx,1-OLy:sNy+OLy)      
      _RL uTrans(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL vTrans(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL wTrans(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      
#ifdef ALLOW_GMREDI
#ifdef GM_BOLUS_ADVEC

C     !LOCAL VARIABLES:
C     == Local variables ==
C     i, j, k :: loop indices
      INTEGER i, j, k
      INTEGER kp1
      _RL maskp1
      _RL delPsi
#ifdef ALLOW_EDDYPSI
      _RL ustar, vstar
#endif
CEOP

      IF ( GM_AdvForm .AND. .NOT.GM_AdvSeparate
     &     .AND. .NOT.GM_InMomAsStress ) THEN

       DO k=1,Nr
        kp1 = MIN(k+1,Nr)
        maskp1 = 1.
        IF (k.GE.Nr) maskp1 = 0.

C--   Calculate tracer cell face open areas
        DO j=1-OLy,sNy+OLy
         DO i=1-OLx,sNx+OLx
          xA(i,j) = _dyG(i,j,bi,bj)*deepFacC(k)
     &           *drF(k)*_hFacW(i,j,k,bi,bj)
          yA(i,j) = _dxG(i,j,bi,bj)*deepFacC(k)
     &           *drF(k)*_hFacS(i,j,k,bi,bj)
         ENDDO
        ENDDO
        
        DO j=1-OLy,sNy+OLy
         DO i=1-OLx,sNx+OLx
           delPsi = GM_PsiX(i,j,kp1,bi,bj)*maskp1
     &            - GM_PsiX(i,j, k, bi,bj)
           uEddy(i,j,k)= delPsi*recip_drF(k)*_recip_hFacW(i,j,k,bi,bj)
           uETrans(i,j,k)= dyG(i,j,bi,bj)*delPsi*maskW(i,j,k,bi,bj)

           uFld(i,j,k)  = uFld(i,j,k) + uEddy(i,j,k)
C--   Calculate "volume transports" through tracer cell faces.
C     anelastic: scaled by rhoFacC (~ mass transport)
           uTrans(i,j,k)= uFld(i,j,k)*xA(i,j)*rhoFacC(k)
         ENDDO
        ENDDO
        DO j=1-OLy,sNy+OLy
         DO i=1-OLx,sNx+OLx
           delPsi = GM_PsiY(i,j,kp1,bi,bj)*maskp1
     &            - GM_PsiY(i,j, k, bi,bj)
           vEddy(i,j,k)= delPsi*recip_drF(k)*_recip_hFacS(i,j,k,bi,bj)
           vETrans(i,j,k)= dxG(i,j,bi,bj)*delPsi*maskS(i,j,k,bi,bj)
         
           vFld(i,j,k)  = vFld(i,j,k) + vEddy(i,j,k)
C--   Calculate "volume transports" through tracer cell faces.
C     anelastic: scaled by rhoFacC (~ mass transport)
           vTrans(i,j,k)= vFld(i,j,k)*yA(i,j)*rhoFacC(k)
         ENDDO
        ENDDO
        DO j=1-OLy,sNy+OLy-1
         DO i=1-OLx,sNx+OLx-1
           delPsi = ( dyG(i+1,j,bi,bj)*GM_PsiX(i+1,j,k,bi,bj)
     &               -dyG( i ,j,bi,bj)*GM_PsiX( i ,j,k,bi,bj)
     &               +dxG(i,j+1,bi,bj)*GM_PsiY(i,j+1,k,bi,bj)
     &               -dxG(i, j ,bi,bj)*GM_PsiY(i, j ,k,bi,bj)
     &              )*maskC(i,j,k,bi,bj)
           wEddy(i,j,k)= delPsi*recip_rA(i,j,bi,bj)
           wETrans(i,j,k) = delPsi

           wFld(i,j,k)  = wFld(i,j,k) + wEddy(i,j,k)
           
C--   Calculate "volume transports" through tracer cell faces.
C     anelastic: scaled by rhoFacC (~ mass transport)           
           IF (k.EQ.1) THEN
C-      Surface interface :
              maskUp(i,j) = 0. _d 0
              wTrans(i,j,k) = 0. _d 0
           ELSE
              maskUp(i,j)  = maskC(i,j,k-1,bi,bj)*maskC(i,j,k,bi,bj)
              wTrans(i,j,k)= wFld(i,j,k)*rA(i,j,bi,bj)*maskUp(i,j)
     &                              *deepFac2F(k)*rhoFacF(k)
           ENDIF
         ENDDO
        ENDDO

       ENDDO

#ifdef ALLOW_DIAGNOSTICS
C JML Residual velocity diagnostics
        IF (useDiagnostics) THEN
C Eddy bolus velocity diagnostics
          CALL DIAGNOSTICS_FILL(uEddy,  'GM_U_EDD',0,Nr,1,bi,bj,myThid)
          CALL DIAGNOSTICS_FILL(vEddy,  'GM_V_EDD',0,Nr,1,bi,bj,myThid)
          CALL DIAGNOSTICS_FILL(wEddy,  'GM_W_EDD',0,Nr,1,bi,bj,myThid)
C Eddy bolus transport diagnostics
          CALL DIAGNOSTICS_FILL(uETrans,'GM_UE_TR',0,Nr,1,bi,bj,myThid)
          CALL DIAGNOSTICS_FILL(vETrans,'GM_VE_TR',0,Nr,1,bi,bj,myThid)
          CALL DIAGNOSTICS_FILL(wETrans,'GM_WE_TR',0,Nr,1,bi,bj,myThid)
C Residual velocity diagnostics          
          CALL DIAGNOSTICS_FILL(uFld,   'GM_U_RES',0,Nr,1,bi,bj,myThid)
          CALL DIAGNOSTICS_FILL(vFld,   'GM_V_RES',0,Nr,1,bi,bj,myThid)
          CALL DIAGNOSTICS_FILL(wFld,   'GM_W_RES',0,Nr,1,bi,bj,myThid)
C Residual transport diagnostics
          CALL DIAGNOSTICS_FILL(uTrans, 'GM_UR_TR',0,Nr,1,bi,bj,myThid)
          CALL DIAGNOSTICS_FILL(vTrans, 'GM_VR_TR',0,Nr,1,bi,bj,myThid)
          CALL DIAGNOSTICS_FILL(wTrans, 'GM_WR_TR',0,Nr,1,bi,bj,myThid)
        ENDIF
#endif /* ALLOW_DIAGNOSTICS */      
#ifdef ALLOW_EDDYPSI
      ELSEIF( GM_AdvForm .AND. .NOT.GM_AdvSeparate
     &        .AND. GM_InMomAsStress ) THEN

C     Calculate the mean velocity from the residual and bolus
       DO k=1,Nr
        kp1 = MIN(k+1,Nr)
        maskp1 = 1.
        IF (k.GE.Nr) maskp1 = 0.

        DO j=1-OLy,sNy+OLy
         DO i=1-OLx,sNx+OLx
          delPsi = GM_PsiX(i,j,kp1,bi,bj)*maskp1
     &           - GM_PsiX(i,j, k, bi,bj)
          uEddy(i,j,k) = delPsi*recip_drF(k)*_recip_hFacW(i,j,k,bi,bj)
          uETrans(i,j,k)= dyG(i,j,bi,bj)*delPsi*hFacW(i,j,k,bi,bj)

          uEulerMean(i,j,k,bi,bj) = uVel(i,j,k,bi,bj) - uEddy(i,j,k)
         ENDDO
        ENDDO
        DO j=1-OLy,sNy+OLy
         DO i=1-OLx,sNx+OLx
          delPsi = GM_PsiY(i,j,kp1,bi,bj)*maskp1
     &           - GM_PsiY(i,j, k, bi,bj)
          vEddy(i,j,k) = delPsi*recip_drF(k)*_recip_hFacS(i,j,k,bi,bj)
          vETrans(i,j,k)= dxG(i,j,bi,bj)*delPsi*hFacS(i,j,k,bi,bj)

          vEulerMean(i,j,k,bi,bj) = vVel(i,j,k,bi,bj) - vEddy(i,j,k)
         ENDDO
        ENDDO

       ENDDO

#ifdef ALLOW_DIAGNOSTICS
       IF ( useDiagnostics ) THEN
        CALL DIAGNOSTICS_FILL(uEulerMean,'U_EulerM',0,Nr,1,bi,bj,myThid)
        CALL DIAGNOSTICS_FILL(vEulerMean,'V_EulerM',0,Nr,1,bi,bj,myThid)
        CALL DIAGNOSTICS_FILL(uEddy,  'GM_U_EDD'   ,0,Nr,1,bi,bj,myThid)
        CALL DIAGNOSTICS_FILL(vEddy , 'GM_V_EDD'   ,0,Nr,1,bi,bj,myThid)
        CALL DIAGNOSTICS_FILL(uETrans,'GM_UE_TR'   ,0,Nr,1,bi,bj,myThid)
        CALL DIAGNOSTICS_FILL(vETrans,'GM_VE_TR'   ,0,Nr,1,bi,bj,myThid)
       ENDIF
#endif /* ALLOW_DIAGNOSTICS */
#endif /* ALLOW_EDDYPSI */

      ENDIF
#endif /* GM_BOLUS_ADVEC */
#endif /* ALLOW_GMREDI */

      RETURN
      END
