C $Header: /u/gcmpack/MITgcm/pkg/gmredi/gmredi_calc_tensor.F,v 1.40 2011/07/13 22:59:53 jmc Exp $
C $Name:  $

#include "GMREDI_OPTIONS.h"
#ifdef ALLOW_KPP
# include "KPP_OPTIONS.h"
#endif

CBOP
C     !ROUTINE: GMREDI_CALC_TENSOR
C     !INTERFACE:
      SUBROUTINE GMREDI_CALC_TENSOR(
     I             iMin, iMax, jMin, jMax,
     I             sigmaX, sigmaY, sigmaR,
     I             bi, bj, myTime, myIter, myThid )

C     !DESCRIPTION: \bv
C     *==========================================================*
C     | SUBROUTINE GMREDI_CALC_TENSOR
C     | o Calculate tensor elements for GM/Redi tensor.
C     *==========================================================*
C     *==========================================================*
C     \ev

C     !USES:
      IMPLICIT NONE

C     == Global variables ==
#include "SIZE.h"
#include "GRID.h"
#include "DYNVARS.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GMREDI.h"
#include "GMREDI_TAVE.h"
#ifdef ALLOW_KPP
# include "KPP.h"
#endif

#ifdef ALLOW_AUTODIFF_TAMC
#include "tamc.h"
#include "tamc_keys.h"
#endif /* ALLOW_AUTODIFF_TAMC */

C     !INPUT/OUTPUT PARAMETERS:
C     == Routine arguments ==
C     bi, bj    :: tile indices
C     myTime    :: Current time in simulation
C     myIter    :: Current iteration number in simulation
C     myThid    :: My Thread Id. number
C
      INTEGER iMin,iMax,jMin,jMax
      _RL sigmaX(1-Olx:sNx+Olx,1-Oly:sNy+Oly,Nr)
      _RL sigmaY(1-Olx:sNx+Olx,1-Oly:sNy+Oly,Nr)
      _RL sigmaR(1-Olx:sNx+Olx,1-Oly:sNy+Oly,Nr)
      INTEGER bi, bj
      _RL     myTime
      INTEGER myIter
      INTEGER myThid
CEOP

#ifdef ALLOW_GMREDI

C     !LOCAL VARIABLES:
C     == Local variables ==
      INTEGER i,j,k,ks
      _RL SlopeX(1-Olx:sNx+Olx,1-Oly:sNy+Oly)
      _RL SlopeY(1-Olx:sNx+Olx,1-Oly:sNy+Oly)
      _RL dSigmaDx(1-Olx:sNx+Olx,1-Oly:sNy+Oly)
      _RL dSigmaDy(1-Olx:sNx+Olx,1-Oly:sNy+Oly)
      _RL dSigmaDr(1-Olx:sNx+Olx,1-Oly:sNy+Oly)
      _RL dSigmaDrref(1-Olx:sNx+Olx,1-Oly:sNy+Oly)
      _RL SlopeSqr(1-Olx:sNx+Olx,1-Oly:sNy+Oly)
      _RL taperFct(1-Olx:sNx+Olx,1-Oly:sNy+Oly)
      _RL ldd97_LrhoC(1-Olx:sNx+Olx,1-Oly:sNy+Oly)
      _RL ldd97_LrhoW(1-Olx:sNx+Olx,1-Oly:sNy+Oly)
      _RL ldd97_LrhoS(1-Olx:sNx+Olx,1-Oly:sNy+Oly)
      _RL Cspd, LrhoInf, LrhoSup, fCoriLoc
      _RL Kgm_tmp, isopycK, bolus_K

      INTEGER kLow_W (1-Olx:sNx+Olx,1-Oly:sNy+Oly)
      INTEGER kLow_S (1-Olx:sNx+Olx,1-Oly:sNy+Oly)
      _RL locMixLayer(1-Olx:sNx+Olx,1-Oly:sNy+Oly)
      _RL baseSlope  (1-Olx:sNx+Olx,1-Oly:sNy+Oly)
      _RL hTransLay  (1-Olx:sNx+Olx,1-Oly:sNy+Oly)
      _RL recipLambda(1-Olx:sNx+Olx,1-Oly:sNy+Oly)
      INTEGER  km1
#if ( defined (GM_NON_UNITY_DIAGONAL) || defined (GM_EXTRA_DIAGONAL) )
      INTEGER kp1
      _RL maskp1
#endif

#ifdef GM_VISBECK_VARIABLE_K
#ifdef OLD_VISBECK_CALC
      _RL Ssq(1-Olx:sNx+Olx,1-Oly:sNy+Oly)
#else
      _RL dSigmaH, dSigmaR, dSigmaRref
      _RL Sloc, M2loc, N2ref, N2rat
      _RL N2tap, GMtap
#endif
      _RL recipMaxSlope
      _RL deltaH, integrDepth
      _RL N2loc, SNloc
#endif /* GM_VISBECK_VARIABLE_K */

#ifdef ALLOW_DIAGNOSTICS
      LOGICAL  doDiagRediFlx
      LOGICAL  DIAGNOSTICS_IS_ON
      EXTERNAL DIAGNOSTICS_IS_ON
#if ( defined (GM_NON_UNITY_DIAGONAL) || defined (GM_EXTRA_DIAGONAL) )
      _RL dTdz
      _RL tmp1k(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
#endif
#endif /* ALLOW_DIAGNOSTICS */

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

#ifdef ALLOW_AUTODIFF_TAMC
          act1 = bi - myBxLo(myThid)
          max1 = myBxHi(myThid) - myBxLo(myThid) + 1
          act2 = bj - myByLo(myThid)
          max2 = myByHi(myThid) - myByLo(myThid) + 1
          act3 = myThid - 1
          max3 = nTx*nTy
          act4 = ikey_dynamics - 1
          igmkey = (act1 + 1) + act2*max1
     &                      + act3*max1*max2
     &                      + act4*max1*max2*max3
#endif /* ALLOW_AUTODIFF_TAMC */

#ifdef ALLOW_DIAGNOSTICS
      doDiagRediFlx = .FALSE.
      IF ( useDiagnostics ) THEN
        doDiagRediFlx = DIAGNOSTICS_IS_ON('GM_KuzTz', myThid )
        doDiagRediFlx = doDiagRediFlx .OR.
     &                  DIAGNOSTICS_IS_ON('GM_KvzTz', myThid )
      ENDIF
#endif

#ifdef GM_VISBECK_VARIABLE_K
      recipMaxSlope = 0. _d 0
      IF ( GM_Visbeck_maxSlope.GT.0. _d 0 ) THEN
        recipMaxSlope = 1. _d 0 / GM_Visbeck_maxSlope
      ENDIF
      DO k=1,Nr
       DO j=1-Oly,sNy+Oly
        DO i=1-Olx,sNx+Olx
          GM_diffK_3d(i,j,k,bi,bj) = 0. _d 0
        ENDDO     
       ENDDO
      ENDDO
#endif

C--   set ldd97_Lrho (for tapering scheme ldd97):
      IF ( GM_taper_scheme.EQ.'ldd97' .OR.
     &     GM_taper_scheme.EQ.'fm07' ) THEN
       Cspd = 2. _d 0
       LrhoInf = 15. _d 3
       LrhoSup = 100. _d 3
C-     Tracer point location (center):
       DO j=1-Oly,sNy+Oly
        DO i=1-Olx,sNx+Olx
         IF (fCori(i,j,bi,bj).NE.0.) THEN
           ldd97_LrhoC(i,j) = Cspd/ABS(fCori(i,j,bi,bj))
         ELSE
           ldd97_LrhoC(i,j) = LrhoSup
         ENDIF
         ldd97_LrhoC(i,j) = MAX(LrhoInf,MIN(ldd97_LrhoC(i,j),LrhoSup))
        ENDDO
       ENDDO
C-     U point location (West):
       DO j=1-Oly,sNy+Oly
        kLow_W(1-Olx,j) = 0
        ldd97_LrhoW(1-Olx,j) = LrhoSup
        DO i=1-Olx+1,sNx+Olx
         kLow_W(i,j) = MIN(kLowC(i-1,j,bi,bj),kLowC(i,j,bi,bj))
         fCoriLoc = op5*(fCori(i-1,j,bi,bj)+fCori(i,j,bi,bj))
         IF (fCoriLoc.NE.0.) THEN
           ldd97_LrhoW(i,j) = Cspd/ABS(fCoriLoc)
         ELSE
           ldd97_LrhoW(i,j) = LrhoSup
         ENDIF
         ldd97_LrhoW(i,j) = MAX(LrhoInf,MIN(ldd97_LrhoW(i,j),LrhoSup))
        ENDDO
       ENDDO
C-     V point location (South):
       DO i=1-Olx+1,sNx+Olx
         kLow_S(i,1-Oly) = 0
         ldd97_LrhoS(i,1-Oly) = LrhoSup
       ENDDO
       DO j=1-Oly+1,sNy+Oly
        DO i=1-Olx,sNx+Olx
         kLow_S(i,j) = MIN(kLowC(i,j-1,bi,bj),kLowC(i,j,bi,bj))
         fCoriLoc = op5*(fCori(i,j-1,bi,bj)+fCori(i,j,bi,bj))
         IF (fCoriLoc.NE.0.) THEN
           ldd97_LrhoS(i,j) = Cspd/ABS(fCoriLoc)
         ELSE
           ldd97_LrhoS(i,j) = LrhoSup
         ENDIF
         ldd97_LrhoS(i,j) = MAX(LrhoInf,MIN(ldd97_LrhoS(i,j),LrhoSup))
        ENDDO
       ENDDO
      ELSE
C-    Just initialize to zero (not use anyway)
       DO j=1-Oly,sNy+Oly
        DO i=1-Olx,sNx+Olx
          ldd97_LrhoC(i,j) = 0. _d 0
          ldd97_LrhoW(i,j) = 0. _d 0
          ldd97_LrhoS(i,j) = 0. _d 0
        ENDDO
       ENDDO
      ENDIF

#ifdef GM_BOLUS_ADVEC
      DO k=1,Nr
       DO j=1-Oly,sNy+Oly
        DO i=1-Olx,sNx+Olx
         GM_PsiX(i,j,k,bi,bj)  = 0. _d 0
         GM_PsiY(i,j,k,bi,bj)  = 0. _d 0
        ENDDO
       ENDDO
      ENDDO
#endif /* GM_BOLUS_ADVEC */
#ifdef ALLOW_AUTODIFF_TAMC
      DO k=1,Nr
       DO j=1-Oly,sNy+Oly
        DO i=1-Olx,sNx+Olx
         Kwx(i,j,k,bi,bj)  = 0. _d 0
         Kwy(i,j,k,bi,bj)  = 0. _d 0
         Kwz(i,j,k,bi,bj)  = 0. _d 0
# ifdef GM_NON_UNITY_DIAGONAL
         Kux(i,j,k,bi,bj)  = 0. _d 0
         Kvy(i,j,k,bi,bj)  = 0. _d 0
# endif
# ifdef GM_EXTRA_DIAGONAL
         Kuz(i,j,k,bi,bj)  = 0. _d 0
         Kvz(i,j,k,bi,bj)  = 0. _d 0
# endif
        ENDDO
       ENDDO
      ENDDO
#endif /* ALLOW_AUTODIFF_TAMC */

C--   Initialise Mixed Layer related array:
      DO j=1-Oly,sNy+Oly
       DO i=1-Olx,sNx+Olx
         hTransLay(i,j) = R_low(i,j,bi,bj)
         baseSlope(i,j) =  0. _d 0
         recipLambda(i,j) = 0. _d 0
         locMixLayer(i,j) = 0. _d 0
       ENDDO
      ENDDO
#ifdef ALLOW_KPP
      IF ( useKPP ) THEN
       DO j=1-Oly,sNy+Oly
        DO i=1-Olx,sNx+Olx
         locMixLayer(i,j) = KPPhbl(i,j,bi,bj)
        ENDDO
       ENDDO
      ELSE
#else
      IF ( .TRUE. ) THEN
#endif
       DO j=1-Oly,sNy+Oly
        DO i=1-Olx,sNx+Olx
         locMixLayer(i,j) = hMixLayer(i,j,bi,bj)
        ENDDO
       ENDDO
      ENDIF

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
C-- 1rst loop on k : compute Tensor Coeff. at W points.
C JML have to fill in the surface level for HMM2011, use Smax instead of Sloc
#ifdef GM_VISBECK_VARIABLE_K
#ifndef OLD_VISBECK_CALC
       ks = 1
#ifdef ALLOW_AUTODIFF_TAMC
       kkey = (igmkey-1)*Nr + 1
       DO j=1-Oly,sNy+Oly
        DO i=1-Olx,sNx+Olx
         dSigmaDx(i,j)     = 0. _d 0
         dSigmaDy(i,j)     = 0. _d 0
         dSigmaDr(i,j)     = 0. _d 0
        ENDDO
       ENDDO
#endif /* ALLOW_AUTODIFF_TAMC */

       DO j=1-Oly+1,sNy+Oly-1
        DO i=1-Olx+1,sNx+Olx-1
C      Gradient of Sigma at rVel points, probably shouldnt use k+1
         dSigmaDx(i,j)=op25*( sigmaX(i+1,j,ks+1)+sigmaX(i,j,ks+1)
     &                       +sigmaX(i+1,j, ks )+sigmaX(i,j, ks )
     &                      )*maskC(i,j,ks,bi,bj)
         dSigmaDy(i,j)=op25*( sigmaY(i,j+1,ks+1)+sigmaY(i,j,ks+1)
     &                       +sigmaY(i,j+1, ks )+sigmaY(i,j, ks )
     &                      )*maskC(i,j,ks,bi,bj)
c        dSigmaDr(i,j)=sigmaR(i,j,k)
        ENDDO
       ENDDO
       
        DO j=1-Oly,sNy+Oly
         DO i=1-Olx,sNx+Olx
          dSigmaDr(i,j) = MIN( sigmaR(i,j,ks), 0. _d 0 )
         ENDDO
        ENDDO
       
        DO j=1-Oly+1,sNy+Oly-1
         DO i=1-Olx+1,sNx+Olx-1
           IF ( maskC(i,j,ks,bi,bj).NE.0. ) THEN
C--      compute: ( M^2 * S )^1/2   (= S*N since S=M^2/N^2 )
C        a 5 points average gives a more "homogeneous" formulation
C        (same stencil and same weights as for dSigmaH calculation)
           dSigmaR = ( dSigmaDr(i,j)*4. _d 0
     &               + dSigmaDr(i-1,j)
     &               + dSigmaDr(i+1,j)
     &               + dSigmaDr(i,j-1)
     &               + dSigmaDr(i,j+1)
     &               )/( 4. _d 0
     &                 + maskC(i-1,j,ks,bi,bj)
     &                 + maskC(i+1,j,ks,bi,bj)
     &                 + maskC(i,j-1,ks,bi,bj)
     &                 + maskC(i,j+1,ks,bi,bj)
     &                 )
     
           IF ( GM_useFM2006N2 ) THEN
C  N^2 vertical dependence on the constant diffusivity
C calculate reference stratification at 2nd level for N2ref
                dSigmaRref = ( dSigmaDrref(i,j)*4. _d 0
     &               + dSigmaDrref(i-1,j)
     &               + dSigmaDrref(i+1,j)
     &               + dSigmaDrref(i,j-1)
     &               + dSigmaDrref(i,j+1)
     &               )/( 4. _d 0
     &                 + maskC(i-1,j,ks,bi,bj)
     &                 + maskC(i+1,j,ks,bi,bj)
     &                 + maskC(i,j-1,ks,bi,bj)
     &                 + maskC(i,j+1,ks,bi,bj)
     &                 )            
            
                N2loc = -gravity*recip_rhoConst*dSigmaR
                N2ref = -gravity*recip_rhoConst*dSigmaRref

                N2rat = N2loc/N2ref
C  JML use the gkw91 tapering scheme to limit N2rat and avoid singularities
                N2tap = MIN (1. _d 0, 1. _d 0/(N2rat*N2rat))
                GM_diffK_3d(i,j,ks,bi,bj) = GM_Visbeck_maxVal_K
     &                     * N2tap * N2rat

           ELSE IF ( GM_useHMM2011 ) THEN
C--     f/sqrt(Ri) = M^2/N^2 * N where M^2 and N^2 are horizontal & vertical 
C          gradient of buoyancy, but basically, dont use the deltaH 
C          dimensionless fraction to calculate depth average.

C       Calculate terms for mean Richardson number which is used
C       in the "variable K" parameterisaton:
                dSigmaH = dSigmaDx(i,j)*dSigmaDx(i,j)
     &             + dSigmaDy(i,j)*dSigmaDy(i,j)
     
             IF ( dSigmaH .GT. 0. _d 0 ) THEN
                 dSigmaH = SQRT( dSigmaH )
C-       compute slope, limited by GM_Visbeck_maxSlope:
              IF ( -dSigmaR.GT.dSigmaH*recipMaxSlope ) THEN
                  Sloc = dSigmaH / ( -dSigmaR )
              ELSE
                  Sloc = GM_Visbeck_maxSlope
              ENDIF

C              M2loc = gravity*recip_rhoConst*dSigmaH
              N2loc = -gravity*recip_rhoConst*dSigmaR

              IF ( N2loc.GT.0. _d 0 ) THEN
                  SNloc = Sloc*SQRT(N2loc)
              ELSE
                  SNloc = 0. _d 0
              ENDIF
             ELSE
                  SNloc = 0. _d 0
             ENDIF
             
            GM_diffK_3d(i,j,ks,bi,bj) = GM_Visbeck_alpha
     &              *GM_Visbeck_length*GM_Visbeck_length*SNloc
C endif for useFM2006 or HMM2011 calculation
           ENDIF
C endif for maskC          
          ENDIF
         ENDDO
        ENDDO
#endif /* ndef OLD_VISBECK_CALC */
#endif /* GM_VISBECK_VARIABLE_K */      

C JML Now loop normally through the depth levels
      DO k=Nr,2,-1

#ifdef ALLOW_AUTODIFF_TAMC
       kkey = (igmkey-1)*Nr + k
       DO j=1-Oly,sNy+Oly
        DO i=1-Olx,sNx+Olx
         SlopeX(i,j)       = 0. _d 0
         SlopeY(i,j)       = 0. _d 0
         dSigmaDx(i,j)     = 0. _d 0
         dSigmaDy(i,j)     = 0. _d 0
         dSigmaDr(i,j)     = 0. _d 0
         SlopeSqr(i,j)     = 0. _d 0
         taperFct(i,j)     = 0. _d 0
        ENDDO
       ENDDO
#endif /* ALLOW_AUTODIFF_TAMC */

       DO j=1-Oly+1,sNy+Oly-1
        DO i=1-Olx+1,sNx+Olx-1
C      Gradient of Sigma at rVel points
         dSigmaDx(i,j)=op25*( sigmaX(i+1,j,k-1)+sigmaX(i,j,k-1)
     &                       +sigmaX(i+1,j, k )+sigmaX(i,j, k )
     &                      )*maskC(i,j,k,bi,bj)
         dSigmaDy(i,j)=op25*( sigmaY(i,j+1,k-1)+sigmaY(i,j,k-1)
     &                       +sigmaY(i,j+1, k )+sigmaY(i,j, k )
     &                      )*maskC(i,j,k,bi,bj)
c        dSigmaDr(i,j)=sigmaR(i,j,k)
        ENDDO
       ENDDO

#ifdef GM_VISBECK_VARIABLE_K
#ifndef OLD_VISBECK_CALC
C JML Add additional option for 3d richardson number parameterisation
C     as in Hoffman and Morales Maqueda (2011) GRL, 10.1029/2010GL045972
C     or the vertical variation in GMK as in Ferreira and Marshall (2006)
      IF ( .NOT. ( GM_useHMM2011 .OR. GM_useFM2006N2 ) ) THEN
       IF ( GM_Visbeck_alpha.GT.0. .AND.
     &      -rC(k-1).LT.GM_Visbeck_depth ) THEN

        DO j=1-Oly,sNy+Oly
         DO i=1-Olx,sNx+Olx
          dSigmaDr(i,j) = MIN( sigmaR(i,j,k), 0. _d 0 )
         ENDDO
        ENDDO

C--     Depth average of f/sqrt(Ri) = M^2/N^2 * N
C       M^2 and N^2 are horizontal & vertical gradient of buoyancy.

C       Calculate terms for mean Richardson number which is used
C       in the "variable K" parameterisaton:
C       compute depth average from surface down to the bottom or
C       GM_Visbeck_depth, whatever is the shallower.

        DO j=1-Oly+1,sNy+Oly-1
         DO i=1-Olx+1,sNx+Olx-1
          IF ( maskC(i,j,k,bi,bj).NE.0. ) THEN
           integrDepth = -rC( kLowC(i,j,bi,bj) )
C-      in 2 steps to avoid mix of RS & RL type in min fct. arguments
           integrDepth = MIN( integrDepth, GM_Visbeck_depth )
C-      to recover "old-visbeck" form with Visbeck_minDepth = Visbeck_depth
           integrDepth = MAX( integrDepth, GM_Visbeck_minDepth )
C       Distance between level center above and the integration depth
           deltaH = integrDepth + rC(k-1)
C       If negative then we are below the integration level
C       (cannot be the case with 2 conditions on maskC & -rC(k-1))
C       If positive we limit this to the distance from center above
           deltaH = MIN( deltaH, drC(k) )
C       Now we convert deltaH to a non-dimensional fraction
           deltaH = deltaH/( integrDepth+rC(1) )

C--      compute: ( M^2 * S )^1/2   (= S*N since S=M^2/N^2 )
C        a 5 points average gives a more "homogeneous" formulation
C        (same stencil and same weights as for dSigmaH calculation)
           dSigmaR = ( dSigmaDr(i,j)*4. _d 0
     &               + dSigmaDr(i-1,j)
     &               + dSigmaDr(i+1,j)
     &               + dSigmaDr(i,j-1)
     &               + dSigmaDr(i,j+1)
     &               )/( 4. _d 0
     &                 + maskC(i-1,j,k,bi,bj)
     &                 + maskC(i+1,j,k,bi,bj)
     &                 + maskC(i,j-1,k,bi,bj)
     &                 + maskC(i,j+1,k,bi,bj)
     &                 )
           dSigmaH = dSigmaDx(i,j)*dSigmaDx(i,j)
     &             + dSigmaDy(i,j)*dSigmaDy(i,j)
           IF ( dSigmaH .GT. 0. _d 0 ) THEN
             dSigmaH = SQRT( dSigmaH )
C-       compute slope, limited by GM_Visbeck_maxSlope:
             IF ( -dSigmaR.GT.dSigmaH*recipMaxSlope ) THEN
              Sloc = dSigmaH / ( -dSigmaR )
             ELSE
              Sloc = GM_Visbeck_maxSlope
             ENDIF
             M2loc = gravity*recip_rhoConst*dSigmaH
c            SNloc = SQRT( Sloc*M2loc )
             N2loc = -gravity*recip_rhoConst*dSigmaR
c            N2loc = -gravity*recip_rhoConst*dSigmaDr(i,j)
             IF ( N2loc.GT.0. _d 0 ) THEN
               SNloc = Sloc*SQRT(N2loc)
             ELSE
               SNloc = 0. _d 0
             ENDIF
           ELSE
             SNloc = 0. _d 0
           ENDIF
           GM_diffK_3d(i,j,1,bi,bj) = GM_diffK_3d(i,j,1,bi,bj)
     &       +deltaH*GM_Visbeck_alpha
     &              *GM_Visbeck_length*GM_Visbeck_length*SNloc
          ENDIF
         ENDDO
        ENDDO
       ENDIF
C JML Add additional option for 3d richardson number parameterisation
C     as in Hoffman and Morales Maqueda (2011) GRL, 10.1029/2010GL045972
      ELSE
        DO j=1-Oly,sNy+Oly
         DO i=1-Olx,sNx+Olx
          dSigmaDr(i,j)    = MIN( sigmaR(i,j,k), 0. _d 0 )
          dSigmaDrref(i,j) = MIN( sigmaR(i,j,2), 0. _d 0 )
         ENDDO
        ENDDO

        DO j=1-Oly+1,sNy+Oly-1
         DO i=1-Olx+1,sNx+Olx-1
          IF ( maskC(i,j,k,bi,bj).NE.0. ) THEN
C--      compute: ( M^2 * S )^1/2   (= S*N since S=M^2/N^2 )
C        a 5 points average gives a more "homogeneous" formulation
C        (same stencil and same weights as for dSigmaH calculation)
           dSigmaR = ( dSigmaDr(i,j)*4. _d 0
     &               + dSigmaDr(i-1,j)
     &               + dSigmaDr(i+1,j)
     &               + dSigmaDr(i,j-1)
     &               + dSigmaDr(i,j+1)
     &               )/( 4. _d 0
     &                 + maskC(i-1,j,k,bi,bj)
     &                 + maskC(i+1,j,k,bi,bj)
     &                 + maskC(i,j-1,k,bi,bj)
     &                 + maskC(i,j+1,k,bi,bj)
     &                 )
     
           IF ( GM_useFM2006N2 ) THEN
C  N^2 vertical dependence on the constant diffusivity
C calculate reference stratification at 2nd level for N2ref
                dSigmaRref = ( dSigmaDrref(i,j)*4. _d 0
     &               + dSigmaDrref(i-1,j)
     &               + dSigmaDrref(i+1,j)
     &               + dSigmaDrref(i,j-1)
     &               + dSigmaDrref(i,j+1)
     &               )/( 4. _d 0
     &                 + maskC(i-1,j,2,bi,bj)
     &                 + maskC(i+1,j,2,bi,bj)
     &                 + maskC(i,j-1,2,bi,bj)
     &                 + maskC(i,j+1,2,bi,bj)
     &                 )            
            
                N2loc = -gravity*recip_rhoConst*dSigmaR
                N2ref = -gravity*recip_rhoConst*dSigmaRref

                N2rat = N2loc/N2ref
C  JML use the gkw91 tapering scheme to limit N2rat and avoid singularities
                N2tap = MIN (1. _d 0, 1. _d 0/(N2rat*N2rat))
                GM_diffK_3d(i,j,k,bi,bj) = GM_Visbeck_maxVal_K
     &                     * N2tap * N2rat

           ELSE IF ( GM_useHMM2011 ) THEN
C--     f/sqrt(Ri) = M^2/N^2 * N where M^2 and N^2 are horizontal & vertical 
C          gradient of buoyancy, but basically, dont use the deltaH 
C          dimensionless fraction to calculate depth average.

C       Calculate terms for mean Richardson number which is used
C       in the "variable K" parameterisaton:
                dSigmaH = dSigmaDx(i,j)*dSigmaDx(i,j)
     &             + dSigmaDy(i,j)*dSigmaDy(i,j)
     
             IF ( dSigmaH .GT. 0. _d 0 ) THEN
                 dSigmaH = SQRT( dSigmaH )
C-       compute slope, limited by GM_Visbeck_maxSlope:
              IF ( -dSigmaR.GT.dSigmaH*recipMaxSlope ) THEN
                  Sloc = dSigmaH / ( -dSigmaR )
              ELSE
                  Sloc = GM_Visbeck_maxSlope
              ENDIF

C              M2loc = gravity*recip_rhoConst*dSigmaH
              N2loc = -gravity*recip_rhoConst*dSigmaR

              IF ( N2loc.GT.0. _d 0 ) THEN
                  SNloc = Sloc*SQRT(N2loc)
              ELSE
                  SNloc = 0. _d 0
              ENDIF
             ELSE
                  SNloc = 0. _d 0
             ENDIF
             
            GM_diffK_3d(i,j,k,bi,bj) = GM_Visbeck_alpha
     &              *GM_Visbeck_length*GM_Visbeck_length*SNloc
C endif for useFM2006 or HMM2011 calculation
           ENDIF
C endif for maskC          
          ENDIF
         ENDDO
        ENDDO
C endif for JML variable parameterisation statement
      ENDIF  
#endif /* ndef OLD_VISBECK_CALC */
#endif /* GM_VISBECK_VARIABLE_K */
       DO j=1-Oly,sNy+Oly
        DO i=1-Olx,sNx+Olx
         dSigmaDr(i,j)=sigmaR(i,j,k)
        ENDDO
       ENDDO

#ifdef ALLOW_AUTODIFF_TAMC
CADJ STORE dSigmaDx(:,:)       = comlev1_bibj_k, key=kkey, byte=isbyte
CADJ STORE dSigmaDy(:,:)       = comlev1_bibj_k, key=kkey, byte=isbyte
CADJ STORE dSigmaDr(:,:)       = comlev1_bibj_k, key=kkey, byte=isbyte
CADJ STORE baseSlope(:,:)      = comlev1_bibj_k, key=kkey, byte=isbyte
CADJ STORE hTransLay(:,:)      = comlev1_bibj_k, key=kkey, byte=isbyte
CADJ STORE recipLambda(:,:)    = comlev1_bibj_k, key=kkey, byte=isbyte
#endif /* ALLOW_AUTODIFF_TAMC */

C     Calculate slopes for use in tensor, taper and/or clip
       CALL GMREDI_SLOPE_LIMIT(
     O             SlopeX, SlopeY,
     O             SlopeSqr, taperFct,
     U             hTransLay, baseSlope, recipLambda,
     U             dSigmaDr,
     I             dSigmaDx, dSigmaDy,
     I             ldd97_LrhoC, locMixLayer, rF,
     I             kLowC(1-Olx,1-Oly,bi,bj),
     I             k, bi, bj, myTime, myIter, myThid )

       DO j=1-Oly+1,sNy+Oly-1
        DO i=1-Olx+1,sNx+Olx-1
C      Mask Iso-neutral slopes
         SlopeX(i,j)=SlopeX(i,j)*maskC(i,j,k,bi,bj)
         SlopeY(i,j)=SlopeY(i,j)*maskC(i,j,k,bi,bj)
         SlopeSqr(i,j)=SlopeSqr(i,j)*maskC(i,j,k,bi,bj)
        ENDDO
       ENDDO

#ifdef ALLOW_AUTODIFF_TAMC
CADJ STORE SlopeX(:,:)       = comlev1_bibj_k, key=kkey, byte=isbyte
CADJ STORE SlopeY(:,:)       = comlev1_bibj_k, key=kkey, byte=isbyte
CADJ STORE SlopeSqr(:,:)     = comlev1_bibj_k, key=kkey, byte=isbyte
CADJ STORE dSigmaDr(:,:)     = comlev1_bibj_k, key=kkey, byte=isbyte
CADJ STORE taperFct(:,:)     = comlev1_bibj_k, key=kkey, byte=isbyte
#endif /* ALLOW_AUTODIFF_TAMC */

C      Components of Redi/GM tensor
       DO j=1-Oly+1,sNy+Oly-1
        DO i=1-Olx+1,sNx+Olx-1
          Kwx(i,j,k,bi,bj)= SlopeX(i,j)*taperFct(i,j)
          Kwy(i,j,k,bi,bj)= SlopeY(i,j)*taperFct(i,j)
          Kwz(i,j,k,bi,bj)= SlopeSqr(i,j)*taperFct(i,j)
        ENDDO
       ENDDO

#ifdef GM_VISBECK_VARIABLE_K
#ifdef OLD_VISBECK_CALC
       DO j=1-Oly+1,sNy+Oly-1
        DO i=1-Olx+1,sNx+Olx-1

C- note (jmc) : moved here since only used in VISBECK_VARIABLE_K
C           but do not know if *taperFct (or **2 ?) is necessary
        Ssq(i,j)=SlopeSqr(i,j)*taperFct(i,j)

C--     Depth average of M^2/N^2 * N

C       Calculate terms for mean Richardson number
C       which is used in the "variable K" parameterisaton.
C       Distance between interface above layer and the integration depth
        deltaH=abs(GM_Visbeck_depth)-abs(rF(k))
C       If positive we limit this to the layer thickness
        integrDepth = drF(k)
        deltaH=min(deltaH,integrDepth)
C       If negative then we are below the integration level
        deltaH=max(deltaH, 0. _d 0)
C       Now we convert deltaH to a non-dimensional fraction
        deltaH=deltaH/GM_Visbeck_depth

        IF ( Ssq(i,j).NE.0. .AND. dSigmaDr(i,j).NE.0. ) THEN
         N2loc = -gravity*recip_rhoConst*dSigmaDr(i,j)
         SNloc = SQRT(Ssq(i,j)*N2loc )
         GM_diffK_3d(i,j,1,bi,bj) = GM_diffK_3d(i,j,1,bi,bj)
     &       +deltaH*GM_Visbeck_alpha
     &              *GM_Visbeck_length*GM_Visbeck_length*SNloc
        ENDIF

        ENDDO
       ENDDO
#endif /* OLD_VISBECK_CALC */
#endif /* GM_VISBECK_VARIABLE_K */

C-- end 1rst loop on vertical level index k
      ENDDO


#ifdef GM_VISBECK_VARIABLE_K
#ifdef ALLOW_AUTODIFF_TAMC
CADJ STORE GM_diffK_3d(:,:,bi,bj) = comlev1_bibj, key=igmkey, byte=isbyte
#endif
C-     Limit range that KapGM can take
          IF ( GM_useHMM2011 .OR. GM_useFM2006N2 ) THEN
C JML Add additional option for 3d GMK parameterisations
             DO k=1,Nr
               DO j=1-Oly+1,sNy+Oly-1
                 DO i=1-Olx+1,sNx+Olx-1
C Taper GM_diffK_3d from either method using gkw91 scheme
                  GMtap = MIN (1. _d 0, 
     &       GM_Visbeck_maxVal_K*GM_Visbeck_maxVal_K/
     &       GM_diffK_3d(i,j,k,bi,bj)*GM_diffK_3d(i,j,k,bi,bj))
                  
                  GM_diffK_3d(i,j,k,bi,bj) = GM_diffK_3d(i,j,k,bi,bj)
     &       * GMtap 

                  GM_diffK_3d(i,j,k,bi,bj)=
     &       MIN( MAX( GM_diffK_3d(i,j,k,bi,bj), GM_Visbeck_minVal_K ),
     &       GM_Visbeck_maxVal_K )
                 ENDDO
               ENDDO
             ENDDO
          ELSE
           IF ( GM_Visbeck_alpha.GT.0. ) THEN
               DO j=1-Oly+1,sNy+Oly-1
                 DO i=1-Olx+1,sNx+Olx-1
             GM_diffK_3d(i,j,1,bi,bj)=
     &       MIN( MAX( GM_diffK_3d(i,j,1,bi,bj), GM_Visbeck_minVal_K ),
     &            GM_Visbeck_maxVal_K )
                 ENDDO
               ENDDO
           ENDIF
          ENDIF
cph( NEW
#ifdef ALLOW_AUTODIFF_TAMC
CADJ STORE GM_diffK_3d(:,:,bi,bj) = comlev1_bibj, key=igmkey, byte=isbyte
#endif
cph)
#endif /* GM_VISBECK_VARIABLE_K */

C-    express the Tensor in term of Diffusivity (= m**2 / s )
      DO k=1,Nr
#ifdef ALLOW_AUTODIFF_TAMC
       kkey = (igmkey-1)*Nr + k
# if (defined (GM_NON_UNITY_DIAGONAL) || \
      defined (GM_VISBECK_VARIABLE_K))
CADJ STORE Kwx(:,:,k,bi,bj) = comlev1_bibj_k, key=kkey, byte=isbyte
CADJ STORE Kwy(:,:,k,bi,bj) = comlev1_bibj_k, key=kkey, byte=isbyte
CADJ STORE Kwz(:,:,k,bi,bj) = comlev1_bibj_k, key=kkey, byte=isbyte
# endif
#endif
       km1 = MAX(k-1,1)
       isopycK = GM_isopycK
     &         *(GM_isoFac1d(km1)+GM_isoFac1d(k))*op5
       bolus_K = GM_background_K
     &         *(GM_bolFac1d(km1)+GM_bolFac1d(k))*op5
       DO j=1-Oly+1,sNy+Oly-1
        DO i=1-Olx+1,sNx+Olx-1
#ifdef ALLOW_KAPREDI_CONTROL
         Kgm_tmp = kapredi(i,j,k,bi,bj)
#else
         Kgm_tmp = isopycK*GM_isoFac2d(i,j,bi,bj)
#endif
#ifdef ALLOW_KAPGM_CONTROL
     &           + GM_skewflx*kapgm(i,j,k,bi,bj)
#else
     &           + GM_skewflx*bolus_K*GM_bolFac2d(i,j,bi,bj)
#endif
#ifdef GM_VISBECK_VARIABLE_K
      IF ( GM_useHMM2011 .OR. GM_useFM2006N2 ) THEN
C JML Add additional option for 3d GMK parameterisations
         Kgm_tmp=Kgm_tmp + GM_diffK_3d(i,j,k,bi,bj)*(1. _d 0 + GM_skewflx)
      ELSE
         Kgm_tmp=Kgm_tmp + GM_diffK_3d(i,j,1,bi,bj)*(1. _d 0 + GM_skewflx)
      ENDIF
#endif
         Kwx(i,j,k,bi,bj)= Kgm_tmp*Kwx(i,j,k,bi,bj)
         Kwy(i,j,k,bi,bj)= Kgm_tmp*Kwy(i,j,k,bi,bj)
#ifdef ALLOW_KAPREDI_CONTROL
         Kgm_tmp = kapredi(i,j,k,bi,bj)
#else
         Kgm_tmp = isopycK*GM_isoFac2d(i,j,bi,bj)
#endif
#ifdef GM_VISBECK_VARIABLE_K
      IF ( GM_useHMM2011 .OR. GM_useFM2006N2 ) THEN
C JML Add additional option for 3d GMK parameterisations
         Kgm_tmp=Kgm_tmp + GM_diffK_3d(i,j,k,bi,bj)
      ELSE
         Kgm_tmp=Kgm_tmp + GM_diffK_3d(i,j,1,bi,bj)
      ENDIF
#endif
         Kwz(i,j,k,bi,bj)  = Kgm_tmp*Kwz(i,j,k,bi,bj)
        ENDDO
       ENDDO
      ENDDO

#ifdef ALLOW_DIAGNOSTICS
      IF ( useDiagnostics .AND. GM_taper_scheme.EQ.'fm07' ) THEN
       CALL DIAGNOSTICS_FILL( hTransLay, 'GM_hTrsL', 0,1,2,bi,bj,myThid)
       CALL DIAGNOSTICS_FILL( baseSlope, 'GM_baseS', 0,1,2,bi,bj,myThid)
       CALL DIAGNOSTICS_FILL(recipLambda,'GM_rLamb', 0,1,2,bi,bj,myThid)
      ENDIF
#endif /* ALLOW_DIAGNOSTICS */

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
C--   Calculate Stream-Functions used in Advective Form:

#ifdef GM_BOLUS_ADVEC
      IF (GM_AdvForm) THEN
#ifdef GM_BOLUS_BVP
       IF (GM_UseBVP) THEN
        CALL GMREDI_CALC_PSI_BVP(
     I             bi, bj, iMin, iMax, jMin, jMax,
     I             sigmaX, sigmaY, sigmaR,
     I             myThid )
       ELSE
#endif
        CALL GMREDI_CALC_PSI_B(
     I              bi, bj, iMin, iMax, jMin, jMax,
     I              sigmaX, sigmaY, sigmaR,
     I              ldd97_LrhoW, ldd97_LrhoS,
     I              myThid )
#ifdef GM_BOLUS_BVP
       ENDIF
#endif
      ENDIF
#endif

#ifndef GM_EXCLUDE_SUBMESO
      IF ( GM_useSubMeso .AND. GM_AdvForm ) THEN
        CALL SUBMESO_CALC_PSI(
     I              bi, bj, iMin, iMax, jMin, jMax,
     I              sigmaX, sigmaY, sigmaR,
     I              locMixLayer,
     I              myIter, myThid )
      ENDIF
#endif /* ndef GM_EXCLUDE_SUBMESO */

#if ( defined (GM_NON_UNITY_DIAGONAL) || defined (GM_EXTRA_DIAGONAL) )
C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
C-- 2nd  k loop : compute Tensor Coeff. at U point

#ifdef ALLOW_KPP
      IF ( useKPP ) THEN
       DO j=1-Oly,sNy+Oly
        DO i=2-Olx,sNx+Olx
         locMixLayer(i,j) = ( KPPhbl(i-1,j,bi,bj)
     &                      + KPPhbl( i ,j,bi,bj) )*op5
        ENDDO
       ENDDO
      ELSE
#else
      IF ( .TRUE. ) THEN
#endif
       DO j=1-Oly,sNy+Oly
        DO i=2-Olx,sNx+Olx
         locMixLayer(i,j) = ( hMixLayer(i-1,j,bi,bj)
     &                      + hMixLayer( i ,j,bi,bj) )*op5
        ENDDO
       ENDDO
      ENDIF
      DO j=1-Oly,sNy+Oly
       DO i=1-Olx,sNx+Olx
         hTransLay(i,j) =  0.
         baseSlope(i,j) =  0.
         recipLambda(i,j)= 0.
       ENDDO
       DO i=2-Olx,sNx+Olx
         hTransLay(i,j) = MAX( R_low(i-1,j,bi,bj), R_low(i,j,bi,bj) )
       ENDDO
      ENDDO

      DO k=Nr,1,-1
       kp1 = MIN(Nr,k+1)
       maskp1 = 1. _d 0
       IF (k.GE.Nr) maskp1 = 0. _d 0
#ifdef ALLOW_AUTODIFF_TAMC
       kkey = (igmkey-1)*Nr + k
#endif

C     Gradient of Sigma at U points
       DO j=1-Oly+1,sNy+Oly-1
        DO i=1-Olx+1,sNx+Olx-1
         dSigmaDx(i,j)=sigmaX(i,j,k)
     &                       *_maskW(i,j,k,bi,bj)
         dSigmaDy(i,j)=op25*( sigmaY(i-1,j+1,k)+sigmaY(i,j+1,k)
     &                       +sigmaY(i-1, j ,k)+sigmaY(i, j ,k)
     &                      )*_maskW(i,j,k,bi,bj)
         dSigmaDr(i,j)=op25*( sigmaR(i-1,j, k )+sigmaR(i,j, k )
     &                      +(sigmaR(i-1,j,kp1)+sigmaR(i,j,kp1))*maskp1
     &                      )*_maskW(i,j,k,bi,bj)
        ENDDO
       ENDDO

#ifdef ALLOW_AUTODIFF_TAMC
CADJ STORE SlopeSqr(:,:)       = comlev1_bibj_k, key=kkey, byte=isbyte
CADJ STORE dSigmaDx(:,:)       = comlev1_bibj_k, key=kkey, byte=isbyte
CADJ STORE dSigmaDy(:,:)       = comlev1_bibj_k, key=kkey, byte=isbyte
CADJ STORE dSigmaDr(:,:)       = comlev1_bibj_k, key=kkey, byte=isbyte
CADJ STORE locMixLayer(:,:)    = comlev1_bibj_k, key=kkey, byte=isbyte
CADJ STORE baseSlope(:,:)      = comlev1_bibj_k, key=kkey, byte=isbyte
CADJ STORE hTransLay(:,:)      = comlev1_bibj_k, key=kkey, byte=isbyte
CADJ STORE recipLambda(:,:)    = comlev1_bibj_k, key=kkey, byte=isbyte
#endif /* ALLOW_AUTODIFF_TAMC */

C     Calculate slopes for use in tensor, taper and/or clip
       CALL GMREDI_SLOPE_LIMIT(
     O             SlopeX, SlopeY,
     O             SlopeSqr, taperFct,
     U             hTransLay, baseSlope, recipLambda,
     U             dSigmaDr,
     I             dSigmaDx, dSigmaDy,
     I             ldd97_LrhoW, locMixLayer, rC,
     I             kLow_W,
     I             k, bi, bj, myTime, myIter, myThid )

#ifdef ALLOW_AUTODIFF_TAMC
CADJ STORE SlopeSqr(:,:)       = comlev1_bibj_k, key=kkey, byte=isbyte
CADJ STORE taperFct(:,:)       = comlev1_bibj_k, key=kkey, byte=isbyte
#endif /* ALLOW_AUTODIFF_TAMC */

#ifdef GM_NON_UNITY_DIAGONAL
c      IF ( GM_nonUnitDiag ) THEN
        DO j=1-Oly+1,sNy+Oly-1
         DO i=1-Olx+1,sNx+Olx-1
          Kux(i,j,k,bi,bj) =
#ifdef ALLOW_KAPREDI_CONTROL
     &     ( kapredi(i,j,k,bi,bj)
#else
     &     ( GM_isopycK*GM_isoFac1d(k)
     &        *op5*(GM_isoFac2d(i-1,j,bi,bj)+GM_isoFac2d(i,j,bi,bj))
#endif
     &     )*taperFct(i,j)
#ifdef GM_VISBECK_VARIABLE_K
      IF ( GM_useHMM2011 .OR. GM_useFM2006N2 ) THEN
C JML Add additional option for 3d GMK parameterisations
            Kux(i,j,k,bi,bj) = Kux(i,j,k,bi,bj) + (
     &        op5*(GM_diffK_3d(i,j,k,bi,bj)+GM_diffK_3d(i-1,j,k,bi,bj))
     &        )*taperFct(i,j)
      ELSE
            Kux(i,j,k,bi,bj) = Kux(i,j,k,bi,bj) + (
     &        op5*(GM_diffK_3d(i,j,1,bi,bj)+GM_diffK_3d(i-1,j,1,bi,bj))
     &        )*taperFct(i,j)
      ENDIF
#endif
         ENDDO
        ENDDO
#ifdef ALLOW_AUTODIFF_TAMC
# ifdef GM_EXCLUDE_CLIPPING
CADJ STORE Kux(:,:,k,bi,bj)  = comlev1_bibj_k, key=kkey, byte=isbyte
# endif
#endif
        DO j=1-Oly+1,sNy+Oly-1
         DO i=1-Olx+1,sNx+Olx-1
          Kux(i,j,k,bi,bj) = MAX( Kux(i,j,k,bi,bj), GM_Kmin_horiz )
         ENDDO
        ENDDO
c      ENDIF
#endif /* GM_NON_UNITY_DIAGONAL */

#ifdef GM_EXTRA_DIAGONAL

#ifdef ALLOW_AUTODIFF_TAMC
CADJ STORE SlopeX(:,:)       = comlev1_bibj_k, key=kkey, byte=isbyte
CADJ STORE taperFct(:,:)     = comlev1_bibj_k, key=kkey, byte=isbyte
#endif /* ALLOW_AUTODIFF_TAMC */
       IF ( GM_ExtraDiag ) THEN
        DO j=1-Oly+1,sNy+Oly-1
         DO i=1-Olx+1,sNx+Olx-1
          Kuz(i,j,k,bi,bj) =
#ifdef ALLOW_KAPREDI_CONTROL
     &     ( kapredi(i,j,k,bi,bj)
#else
     &     ( GM_isopycK*GM_isoFac1d(k)
     &        *op5*(GM_isoFac2d(i-1,j,bi,bj)+GM_isoFac2d(i,j,bi,bj))
#endif
#ifdef ALLOW_KAPGM_CONTROL
     &     - GM_skewflx*kapgm(i,j,k,bi,bj)
#else
     &     - GM_skewflx*GM_background_K*GM_bolFac1d(k)
     &        *op5*(GM_bolFac2d(i-1,j,bi,bj)+GM_bolFac2d(i,j,bi,bj))
#endif
     &     )*SlopeX(i,j)*taperFct(i,j)
#ifdef GM_VISBECK_VARIABLE_K
      IF ( GM_useHMM2011 .OR. GM_useFM2006N2 ) THEN
C JML Add additional option for 3d GMK parameterisations
            Kuz(i,j,k,bi,bj) = Kuz(i,j,k,bi,bj) + (
     &        op5*(GM_diffK_3d(i,j,k,bi,bj)+GM_diffK_3d(i-1,j,k,bi,bj))*GM_advect
     &        )*SlopeX(i,j)*taperFct(i,j)
      ELSE
            Kuz(i,j,k,bi,bj) = Kuz(i,j,k,bi,bj) + (
     &        op5*(GM_diffK_3d(i,j,1,bi,bj)+GM_diffK_3d(i-1,j,1,bi,bj))*GM_advect
     &        )*SlopeX(i,j)*taperFct(i,j)
      ENDIF
#endif
         ENDDO
        ENDDO
       ENDIF
#endif /* GM_EXTRA_DIAGONAL */

#ifdef ALLOW_DIAGNOSTICS
       IF (doDiagRediFlx) THEN
        km1 = MAX(k-1,1)
        DO j=1,sNy
         DO i=1,sNx+1
C         store in tmp1k Kuz_Redi
#ifdef ALLOW_KAPREDI_CONTROL
          tmp1k(i,j) = ( kapredi(i,j,k,bi,bj)
#else
          tmp1k(i,j) = ( GM_isopycK*GM_isoFac1d(k)
     &        *op5*(GM_isoFac2d(i-1,j,bi,bj)+GM_isoFac2d(i,j,bi,bj))
#endif
     &                 )*SlopeX(i,j)*taperFct(i,j)
#ifdef GM_VISBECK_VARIABLE_K
      IF ( GM_useHMM2011 .OR. GM_useFM2006N2 ) THEN
C JML Add additional option for 3d GMK parameterisations
            tmp1k(i,j) = tmp1k(i,j) + (
     &        op5*(GM_diffK_3d(i,j,k,bi,bj)+GM_diffK_3d(i-1,j,k,bi,bj))
     &        )*SlopeX(i,j)*taperFct(i,j)
      ELSE
            tmp1k(i,j) = tmp1k(i,j) + (
     &        op5*(GM_diffK_3d(i,j,1,bi,bj)+GM_diffK_3d(i-1,j,1,bi,bj))
     &        )*SlopeX(i,j)*taperFct(i,j)
      ENDIF
#endif
C#endif  --necessary?
         ENDDO
        ENDDO
        DO j=1,sNy
         DO i=1,sNx+1
C-        Vertical gradients interpolated to U points
          dTdz = (
     &     +recip_drC(k)*
     &       ( maskC(i-1,j,k,bi,bj)*
     &           (theta(i-1,j,km1,bi,bj)-theta(i-1,j,k,bi,bj))
     &        +maskC( i ,j,k,bi,bj)*
     &           (theta( i ,j,km1,bi,bj)-theta( i ,j,k,bi,bj))
     &       )
     &     +recip_drC(kp1)*
     &       ( maskC(i-1,j,kp1,bi,bj)*
     &           (theta(i-1,j,k,bi,bj)-theta(i-1,j,kp1,bi,bj))
     &        +maskC( i ,j,kp1,bi,bj)*
     &           (theta( i ,j,k,bi,bj)-theta( i ,j,kp1,bi,bj))
     &       )      ) * 0.25 _d 0
           tmp1k(i,j) = dyG(i,j,bi,bj)*drF(k)
     &                * _hFacW(i,j,k,bi,bj)
     &                * tmp1k(i,j) * dTdz
         ENDDO
        ENDDO
        CALL DIAGNOSTICS_FILL(tmp1k, 'GM_KuzTz', k,1,2,bi,bj,myThid)
       ENDIF
#endif /* ALLOW_DIAGNOSTICS */

C-- end 2nd  loop on vertical level index k
      ENDDO

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
C-- 3rd  k loop : compute Tensor Coeff. at V point

#ifdef ALLOW_KPP
      IF ( useKPP ) THEN
       DO j=2-Oly,sNy+Oly
        DO i=1-Olx,sNx+Olx
         locMixLayer(i,j) = ( KPPhbl(i,j-1,bi,bj)
     &                      + KPPhbl(i, j ,bi,bj) )*op5
        ENDDO
       ENDDO
      ELSE
#else
      IF ( .TRUE. ) THEN
#endif
       DO j=2-Oly,sNy+Oly
        DO i=1-Olx,sNx+Olx
         locMixLayer(i,j) = ( hMixLayer(i,j-1,bi,bj)
     &                      + hMixLayer(i, j ,bi,bj) )*op5
        ENDDO
       ENDDO
      ENDIF
      DO j=1-Oly,sNy+Oly
       DO i=1-Olx,sNx+Olx
         hTransLay(i,j) =  0.
         baseSlope(i,j) =  0.
         recipLambda(i,j)= 0.
       ENDDO
      ENDDO
      DO j=2-Oly,sNy+Oly
       DO i=1-Olx,sNx+Olx
         hTransLay(i,j) = MAX( R_low(i,j-1,bi,bj), R_low(i,j,bi,bj) )
       ENDDO
      ENDDO

C     Gradient of Sigma at V points
      DO k=Nr,1,-1
       kp1 = MIN(Nr,k+1)
       maskp1 = 1. _d 0
       IF (k.GE.Nr) maskp1 = 0. _d 0
#ifdef ALLOW_AUTODIFF_TAMC
       kkey = (igmkey-1)*Nr + k
#endif

       DO j=1-Oly+1,sNy+Oly-1
        DO i=1-Olx+1,sNx+Olx-1
         dSigmaDx(i,j)=op25*( sigmaX(i, j ,k) +sigmaX(i+1, j ,k)
     &                       +sigmaX(i,j-1,k) +sigmaX(i+1,j-1,k)
     &                      )*_maskS(i,j,k,bi,bj)
         dSigmaDy(i,j)=sigmaY(i,j,k)
     &                       *_maskS(i,j,k,bi,bj)
         dSigmaDr(i,j)=op25*( sigmaR(i,j-1, k )+sigmaR(i,j, k )
     &                      +(sigmaR(i,j-1,kp1)+sigmaR(i,j,kp1))*maskp1
     &                      )*_maskS(i,j,k,bi,bj)
        ENDDO
       ENDDO

#ifdef ALLOW_AUTODIFF_TAMC
CADJ STORE dSigmaDx(:,:)       = comlev1_bibj_k, key=kkey, byte=isbyte
CADJ STORE dSigmaDy(:,:)       = comlev1_bibj_k, key=kkey, byte=isbyte
CADJ STORE dSigmaDr(:,:)       = comlev1_bibj_k, key=kkey, byte=isbyte
CADJ STORE baseSlope(:,:)      = comlev1_bibj_k, key=kkey, byte=isbyte
CADJ STORE hTransLay(:,:)      = comlev1_bibj_k, key=kkey, byte=isbyte
CADJ STORE recipLambda(:,:)    = comlev1_bibj_k, key=kkey, byte=isbyte
#endif /* ALLOW_AUTODIFF_TAMC */

C     Calculate slopes for use in tensor, taper and/or clip
       CALL GMREDI_SLOPE_LIMIT(
     O             SlopeX, SlopeY,
     O             SlopeSqr, taperFct,
     U             hTransLay, baseSlope, recipLambda,
     U             dSigmaDr,
     I             dSigmaDx, dSigmaDy,
     I             ldd97_LrhoS, locMixLayer, rC,
     I             kLow_S,
     I             k, bi, bj, myTime, myIter, myThid )

cph(
#ifdef ALLOW_AUTODIFF_TAMC
cph(
CADJ STORE taperfct(:,:)       = comlev1_bibj_k, key=kkey, byte=isbyte
cph)
#endif /* ALLOW_AUTODIFF_TAMC */
cph)

#ifdef GM_NON_UNITY_DIAGONAL
c      IF ( GM_nonUnitDiag ) THEN
        DO j=1-Oly+1,sNy+Oly-1
         DO i=1-Olx+1,sNx+Olx-1
          Kvy(i,j,k,bi,bj) =
#ifdef ALLOW_KAPREDI_CONTROL
     &     ( kapredi(i,j,k,bi,bj)
#else
     &     ( GM_isopycK*GM_isoFac1d(k)
     &        *op5*(GM_isoFac2d(i,j-1,bi,bj)+GM_isoFac2d(i,j,bi,bj))
#endif
     &     )*taperFct(i,j)
#ifdef GM_VISBECK_VARIABLE_K
      IF ( GM_useHMM2011 .OR. GM_useFM2006N2 ) THEN
C JML Add additional option for 3d GMK parameterisations
            Kvy(i,j,k,bi,bj) = Kvy(i,j,k,bi,bj) + (
     &        op5*(GM_diffK_3d(i,j,k,bi,bj)+GM_diffK_3d(i,j-1,k,bi,bj))
     &        )*taperFct(i,j)
      ELSE
            Kvy(i,j,k,bi,bj) = Kvy(i,j,k,bi,bj) + (
     &        op5*(GM_diffK_3d(i,j,1,bi,bj)+GM_diffK_3d(i,j-1,1,bi,bj))
     &        )*taperFct(i,j)
      ENDIF
#endif

         ENDDO
        ENDDO
#ifdef ALLOW_AUTODIFF_TAMC
# ifdef GM_EXCLUDE_CLIPPING
CADJ STORE Kvy(:,:,k,bi,bj)  = comlev1_bibj_k, key=kkey, byte=isbyte
# endif
#endif
        DO j=1-Oly+1,sNy+Oly-1
         DO i=1-Olx+1,sNx+Olx-1
          Kvy(i,j,k,bi,bj) = MAX( Kvy(i,j,k,bi,bj), GM_Kmin_horiz )
         ENDDO
        ENDDO
c      ENDIF
#endif /* GM_NON_UNITY_DIAGONAL */

#ifdef GM_EXTRA_DIAGONAL

#ifdef ALLOW_AUTODIFF_TAMC
CADJ STORE SlopeY(:,:)       = comlev1_bibj_k, key=kkey, byte=isbyte
CADJ STORE taperFct(:,:)     = comlev1_bibj_k, key=kkey, byte=isbyte
#endif /* ALLOW_AUTODIFF_TAMC */
       IF ( GM_ExtraDiag ) THEN
        DO j=1-Oly+1,sNy+Oly-1
         DO i=1-Olx+1,sNx+Olx-1
          Kvz(i,j,k,bi,bj) =
#ifdef ALLOW_KAPREDI_CONTROL
     &     ( kapredi(i,j,k,bi,bj)
#else
     &     ( GM_isopycK*GM_isoFac1d(k)
     &        *op5*(GM_isoFac2d(i,j-1,bi,bj)+GM_isoFac2d(i,j,bi,bj))
#endif
#ifdef ALLOW_KAPGM_CONTROL
     &     - GM_skewflx*kapgm(i,j,k,bi,bj)
#else
     &     - GM_skewflx*GM_background_K*GM_bolFac1d(k)
     &        *op5*(GM_bolFac2d(i,j-1,bi,bj)+GM_bolFac2d(i,j,bi,bj))
#endif
     &     )*SlopeY(i,j)*taperFct(i,j)
#ifdef GM_VISBECK_VARIABLE_K
      IF ( GM_useHMM2011 .OR. GM_useFM2006N2 ) THEN
C JML Add additional option for 3d GMK parameterisations
            Kvz(i,j,k,bi,bj) = Kvz(i,j,k,bi,bj) + (
     &        op5*(GM_diffK_3d(i,j,k,bi,bj)+GM_diffK_3d(i,j-1,k,bi,bj))*GM_advect
     &        )*SlopeY(i,j)*taperFct(i,j)
      ELSE
            Kvz(i,j,k,bi,bj) = Kvz(i,j,k,bi,bj) + (
     &        op5*(GM_diffK_3d(i,j,1,bi,bj)+GM_diffK_3d(i,j-1,1,bi,bj))*GM_advect
     &        )*SlopeY(i,j)*taperFct(i,j)
      ENDIF
#endif
         ENDDO
        ENDDO
       ENDIF
#endif /* GM_EXTRA_DIAGONAL */

#ifdef ALLOW_DIAGNOSTICS
       IF (doDiagRediFlx) THEN
        km1 = MAX(k-1,1)
        DO j=1,sNy+1
         DO i=1,sNx
C         store in tmp1k Kvz_Redi
#ifdef ALLOW_KAPREDI_CONTROL
          tmp1k(i,j) = ( kapredi(i,j,k,bi,bj)
#else
          tmp1k(i,j) = ( GM_isopycK*GM_isoFac1d(k)
     &        *op5*(GM_isoFac2d(i,j-1,bi,bj)+GM_isoFac2d(i,j,bi,bj))
#endif
     &        )*SlopeY(i,j)*taperFct(i,j)
#ifdef GM_VISBECK_VARIABLE_K
      IF ( GM_useHMM2011 .OR. GM_useFM2006N2 ) THEN
C JML Add additional option for 3d GMK parameterisations
            tmp1k(i,j) = tmp1k(i,j) + (
     &        op5*(GM_diffK_3d(i,j,k,bi,bj)+GM_diffK_3d(i,j-1,k,bi,bj))
     &        )*SlopeY(i,j)*taperFct(i,j)
      ELSE
            tmp1k(i,j) = tmp1k(i,j) + (
     &        op5*(GM_diffK_3d(i,j,1,bi,bj)+GM_diffK_3d(i,j-1,1,bi,bj))
     &        )*SlopeY(i,j)*taperFct(i,j)
      ENDIF
#endif
         ENDDO
        ENDDO
        DO j=1,sNy+1
         DO i=1,sNx
C-        Vertical gradients interpolated to U points
          dTdz = (
     &     +recip_drC(k)*
     &       ( maskC(i,j-1,k,bi,bj)*
     &           (theta(i,j-1,km1,bi,bj)-theta(i,j-1,k,bi,bj))
     &        +maskC(i, j ,k,bi,bj)*
     &           (theta(i, j ,km1,bi,bj)-theta(i, j ,k,bi,bj))
     &       )
     &     +recip_drC(kp1)*
     &       ( maskC(i,j-1,kp1,bi,bj)*
     &           (theta(i,j-1,k,bi,bj)-theta(i,j-1,kp1,bi,bj))
     &        +maskC(i, j ,kp1,bi,bj)*
     &           (theta(i, j ,k,bi,bj)-theta(i, j ,kp1,bi,bj))
     &       )      ) * 0.25 _d 0
           tmp1k(i,j) = dxG(i,j,bi,bj)*drF(k)
     &                * _hFacS(i,j,k,bi,bj)
     &                * tmp1k(i,j) * dTdz
         ENDDO
        ENDDO
        CALL DIAGNOSTICS_FILL(tmp1k, 'GM_KvzTz', k,1,2,bi,bj,myThid)
       ENDIF
#endif /* ALLOW_DIAGNOSTICS */

C-- end 3rd  loop on vertical level index k
      ENDDO

#endif /* GM_NON_UNITY_DIAGONAL || GM_EXTRA_DIAGONAL */

#ifdef ALLOW_TIMEAVE
C--   Time-average
      IF ( taveFreq.GT.0. ) THEN

         CALL TIMEAVE_CUMULATE( GM_Kwx_T, Kwx, Nr,
     &                          deltaTclock, bi, bj, myThid )
         CALL TIMEAVE_CUMULATE( GM_Kwy_T, Kwy, Nr,
     &                          deltaTclock, bi, bj, myThid )
         CALL TIMEAVE_CUMULATE( GM_Kwz_T, Kwz, Nr,
     &                          deltaTclock, bi, bj, myThid )
#ifdef GM_VISBECK_VARIABLE_K
       IF ( GM_Visbeck_alpha.NE.0. .or. GM_useFM2006N2 ) THEN
C JML Add additional option for 3d richardson number parameterisations
         CALL TIMEAVE_CUMULATE( Visbeck_K_T, GM_diffK_3d, Nr,
     &                          deltaTclock, bi, bj, myThid )
       ENDIF
#endif
#ifdef GM_BOLUS_ADVEC
       IF ( GM_AdvForm ) THEN
         CALL TIMEAVE_CUMULATE( GM_PsiXtave, GM_PsiX, Nr,
     &                          deltaTclock, bi, bj, myThid )
         CALL TIMEAVE_CUMULATE( GM_PsiYtave, GM_PsiY, Nr,
     &                          deltaTclock, bi, bj, myThid )
       ENDIF
#endif
       GM_timeAve(bi,bj) = GM_timeAve(bi,bj)+deltaTclock

      ENDIF
#endif /* ALLOW_TIMEAVE */

#ifdef ALLOW_DIAGNOSTICS
      IF ( useDiagnostics ) THEN
        CALL GMREDI_DIAGNOSTICS_FILL(bi,bj,myThid)
      ENDIF
#endif /* ALLOW_DIAGNOSTICS */

#endif /* ALLOW_GMREDI */

      RETURN
      END

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

CBOP
C     !ROUTINE: GMREDI_CALC_TENSOR_DUMMY
C     !INTERFACE:
      SUBROUTINE GMREDI_CALC_TENSOR_DUMMY(
     I             iMin, iMax, jMin, jMax,
     I             sigmaX, sigmaY, sigmaR,
     I             bi, bj, myTime, myIter, myThid )

C     !DESCRIPTION: \bv
C     *==========================================================*
C     | SUBROUTINE GMREDI_CALC_TENSOR_DUMMY
C     | o Calculate tensor elements for GM/Redi tensor.
C     *==========================================================*
C     \ev

C     !USES:
      IMPLICIT NONE

C     == Global variables ==
#include "SIZE.h"
#include "EEPARAMS.h"
#include "GMREDI.h"

C     !INPUT/OUTPUT PARAMETERS:
      _RL sigmaX(1-Olx:sNx+Olx,1-Oly:sNy+Oly,Nr)
      _RL sigmaY(1-Olx:sNx+Olx,1-Oly:sNy+Oly,Nr)
      _RL sigmaR(1-Olx:sNx+Olx,1-Oly:sNy+Oly,Nr)
      INTEGER iMin,iMax,jMin,jMax
      INTEGER bi, bj
      _RL     myTime
      INTEGER myIter
      INTEGER myThid
CEOP

#ifdef ALLOW_GMREDI
C     !LOCAL VARIABLES:
      INTEGER i, j, k

      DO k=1,Nr
       DO j=1-Oly+1,sNy+Oly-1
        DO i=1-Olx+1,sNx+Olx-1
         Kwx(i,j,k,bi,bj) = 0.0
         Kwy(i,j,k,bi,bj) = 0.0
         Kwz(i,j,k,bi,bj) = 0.0
        ENDDO
       ENDDO
      ENDDO
#endif /* ALLOW_GMREDI */

      RETURN
      END
