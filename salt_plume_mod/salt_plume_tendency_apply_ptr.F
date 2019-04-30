C $Header: /u/gcmpack/MITgcm/pkg/salt_plume/salt_plume_tendency_apply_s.F,v 1.7 2012/01/30 11:19:02 mlosch Exp $
C $Name:  $

#include "SALT_PLUME_OPTIONS.h"

CBOP 0
C     !ROUTINE: SALT_PLUME_TENDENCY_APPLY_PTR
C     !INTERFACE:
      SUBROUTINE SALT_PLUME_TENDENCY_APPLY_PTR(
     &     iMin, iMax, jMin, jMax,
     &     bi,bj,kLev,iTrc,myTime,myThid)

C     !DESCRIPTION:
C     Add salt_plume tendency terms to PTRACER tendency.
C     Routine works for one level at a time.
C     SaltPlume is the amount of salt rejected by ice while freezing;
C     it is here redistributed to multiple vertical levels as per
C     Duffy et al. (GRL 1999).

C     !INPUT PARAMETERS:
      IMPLICIT NONE
#include "SIZE.h"
#include "GRID.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "DYNVARS.h"
#include "SALT_PLUME.h"
#include "PTRACERS_SIZE.h"
#include "PTRACERS_PARAMS.h"
#include "PTRACERS_FIELDS.h"

C     !INPUT PARAMETERS:
      integer iMin, iMax, jMin, jMax, kLev, iTrc, bi, bj, myThid
      _RL myTime
#ifdef ALLOW_DIAGNOSTICS
      CHARACTER*8 diagName
      LOGICAL  DIAGNOSTICS_IS_ON
      EXTERNAL DIAGNOSTICS_IS_ON
#endif

CEOP

#if defined (ALLOW_SALT_PLUME) && (defined SALT_PLUME_DIC)

C     !LOCAL VARIABLES:
      integer i, j
      _RL minusone
      parameter(minusone = -1.)
      _RL plumefrac(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL plumetend(1-OLx:sNx+OLx,1-OLy:sNy+OLy)

#ifdef TARGET_NEC_SX
      integer imt
      parameter( imt=(sNx+2*OLx)*(sNy+2*OLy) )
      _RL plumekb2D(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
#else
      integer two2
      parameter(two2 = 2)
      _RL plumekb(two2), SPdepth(two2)
#endif
#ifdef TARGET_NEC_SX
C     The vector version computes plumekb2D at each grid point, but this
C     is still faster than non-vector code.
      IF ( kLev .LT. Nr ) THEN
       DO j=1-Oly,sNy+Oly
        DO i=1-Olx,sNx+Olx
         plumekb2D(i,j)=abs(rF(kLev))
        ENDDO
       ENDDO
       CALL SALT_PLUME_FRAC(
     I      imt,minusone,SaltPlumeDepth(1-Olx,1-Oly,bi,bj),
     U      plumekb2D,
     I      myTime, 1, myThid )
       DO j=1-Oly,sNy+Oly
        DO i=1-Olx,sNx+Olx
         plumefrac(I,J) = plumekb2D(i,j)*maskC(i,j,klev,bi,bj)
         plumekb2D(i,j) = abs(rF(kLev+1))
        ENDDO
       ENDDO
       CALL SALT_PLUME_FRAC(
     I      imt,minusone,SaltPlumeDepth(1-Olx,1-Oly,bi,bj),
     U      plumekb2D,
     I      myTime, 1, myThid )
       DO j=1-Oly,sNy+Oly
        DO i=1-Olx,sNx+Olx
C     Penetrating saltplume fraction:
         IF ( (maskC(i,j,kLev,bi,bj).NE.0.0 _d 0) .AND.
     &        ( SaltPlumeDepth(i,j,bi,bj) .GT. -rF(kLev)) ) THEN
          plumefrac(i,j)     = plumefrac(i,j)
     &         -plumekb2D(i,j)*maskC(i,j,klev,bi,bj)
          IF (convertFW2Salt .EQ. -1.) THEN
C-    use local surface tracer field to calculate forcing term:
           plumetend(i,j)=saltPlumeFlux(i,j,bi,bj)*plumefrac(I,J)
     &        *mass2rUnit
     &        *( pTracer(i,j,kLev,bi,bj,iTrc) / SALT(i,j,kLev,bi,bj))
     &        *recip_drF(kLev)*_recip_hFacC(i,j,kLev,bi,bj)
            
            gPtr(i,j,kLev,bi,bj,iTrc) = gPtr(i,j,kLev,bi,bj,iTrc)
     &        +  plumetend(i,j)
          ELSE
C-    use uniform tracer value to calculate forcing term:
           plumetend(i,j)=saltPlumeFlux(i,j,bi,bj)*plumefrac(I,J)
     &        *mass2rUnit
     &        *( pTracer(i,j,kLev,bi,bj,iTrc) / convertFW2Salt)
     &        *recip_drF(kLev)*_recip_hFacC(i,j,kLev,bi,bj) 
     
            gPtr(i,j,kLev,bi,bj,iTrc) = gPtr(i,j,kLev,bi,bj,iTrc)
     &        +  plumetend(i,j)          
          ENDIF
         ELSE
          plumefrac(i,j)=0.0 _d 0    
          plumetend(i,j)=0.0 _d 0
         ENDIF   
        ENDDO
       ENDDO
      ENDIF
#else
      DO j=jMin,jMax
       DO i=iMin,iMax
C Penetrating saltplume fraction:
        IF ( (maskC(i,j,kLev,bi,bj).NE.0.0 _d 0) .AND.
     &        ( SaltPlumeDepth(i,j,bi,bj) .GT. -rF(kLev)) ) THEN
         plumefrac(I,J) = 0. _d 0
         plumekb(1)=abs(rF(klev))
         plumekb(2)=abs(rF(klev+1))
         SPdepth(1)=SaltPlumeDepth(i,j,bi,bj)
         SPdepth(2)=SaltPlumeDepth(i,j,bi,bj)
         CALL SALT_PLUME_FRAC(
     I                   two2,minusone,SPdepth,
     U                   plumekb,
     I                   myTime, 1, myThid )
         plumefrac(I,J)=(plumekb(1)-plumekb(2))*maskC(i,j,klev,bi,bj)

          IF (convertFW2Salt .EQ. -1.) THEN
C-    use local surface tracer field to calculate forcing term:
           plumetend(i,j)=saltPlumeFlux(i,j,bi,bj)*plumefrac(I,J)
     &        *mass2rUnit
     &        *( pTracer(i,j,kLev,bi,bj,iTrc) / SALT(i,j,kLev,bi,bj))
     &        *recip_drF(kLev)*_recip_hFacC(i,j,kLev,bi,bj)
            
            gPtr(i,j,kLev,bi,bj,iTrc) = gPtr(i,j,kLev,bi,bj,iTrc)
     &        +  plumetend(i,j) 
          ELSE
C-    use uniform tracer value to calculate forcing term:
           plumetend(i,j)=saltPlumeFlux(i,j,bi,bj)*plumefrac(I,J)
     &        *mass2rUnit
     &        *( pTracer(i,j,kLev,bi,bj,iTrc) / convertFW2Salt)
     &        *recip_drF(kLev)*_recip_hFacC(i,j,kLev,bi,bj) 
     
            gPtr(i,j,kLev,bi,bj,iTrc) = gPtr(i,j,kLev,bi,bj,iTrc)
     &        +  plumetend(i,j)     
          ENDIF
        ELSE
           plumefrac(i,j)=0.0 _d 0
           plumetend(i,j)=0.0 _d 0
        ENDIF
       ENDDO
      ENDDO
#endif /* TARGET_NEC_SX */

#ifdef ALLOW_DIAGNOSTICS
      diagName = '        '
       WRITE(diagName,'(A6,A2)') 'SPtndT',PTRACERS_ioLabel(iTrc)
      IF ( DIAGNOSTICS_IS_ON( diagName , myThid ) ) THEN
       CALL DIAGNOSTICS_FILL (
     &      plumetend,diagName,kLev,1,2,bi,bj,myThid )
      ENDIF
#endif /* ALLOW_DIAGNOSTICS */     
#endif /* ALLOW_SALT_PLUME */
      
      RETURN
      END
