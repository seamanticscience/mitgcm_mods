C $Header: /u/gcmpack/MITgcm/pkg/salt_plume/salt_plume_diagnostics_init.F,v 1.3 2010/01/12 21:31:15 jmc Exp $
C $Name:  $

#include "SALT_PLUME_OPTIONS.h"

CBOP 0
C     !ROUTINE: SALT_PLUME_DIAGNOSTICS_INIT
C     !INTERFACE:
      SUBROUTINE SALT_PLUME_DIAGNOSTICS_INIT( myThid )

C     !DESCRIPTION:
C     Initialize list of all available diagnostics

C     !USES:
      IMPLICIT NONE
#include "SIZE.h"      
#include "EEPARAMS.h"

#ifdef ALLOW_PTRACERS
#include "PTRACERS_SIZE.h"
#include "PTRACERS_PARAMS.h"
#endif

C     !INPUT/OUTPUT PARAMETERS:
C     myThid ::  my Thread Id number
      INTEGER myThid
CEOP

#ifdef ALLOW_SALT_PLUME
#ifdef ALLOW_DIAGNOSTICS

C     !LOCAL VARIABLES:
C     === Local variables ===
C     diagNum   :: diagnostics number in the (long) list of available diag.
C     diagName  :: local short name (8c) of a diagnostics
C     diagCode  :: local parser field with characteristics of the diagnostics
C              cf head of S/R DIAGNOSTICS_INIT_EARLY or DIAGNOSTICS_MAIN_INIT
C     diagUnits :: local string (16c): physical units of a diagnostic field
C     diagTitle :: local string (80c): description of field in diagnostic
      INTEGER       diagNum
      CHARACTER*8   diagName
      CHARACTER*16  diagCode
      CHARACTER*16  diagUnits
      CHARACTER*(80) diagTitle

#ifdef ALLOW_PTRACERS
      INTEGER      iTrc, ilnb, tracerId
      CHARACTER*7  trcUnits
      CHARACTER*9  flxUnits, wUnits
      CHARACTER*30 locName
      CHARACTER*4  diagSufx
C     Functions ::
      INTEGER     ILNBLNK
      EXTERNAL    ILNBLNK
#endif

      diagName  = 'PLUMEKB '
      diagTitle = 'fractional plume: [0-1] (unitless)'
      diagUnits = '                '
      diagCode  = 'SM      MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I       diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'SPHeff  '
      diagTitle = 'Salt Plume effective ice thickness (iceVol/oceArea)'
      diagUnits = 'm'
      diagCode  = 'SM      U1      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I       diagName, diagCode, diagUnits, diagTitle, 0, myThid )
     
      diagName  = 'SPdHeff '
      diagTitle = 'Salt Plume change in effective ice thickness'
      diagUnits = 'm/s'
      diagCode  = 'SM      U1      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I       diagName, diagCode, diagUnits, diagTitle, 0, myThid )
     
      diagName  = 'oceSPflx'
      diagTitle = 'net surface Salt flux rejected into the ocean during'
     &          //' freezing, (+=down),'
     &          //' >0 increases salinity'
      diagUnits = 'g/m^2/s         '
      diagCode  = 'SM      U1      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I       diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'oceSPDep'
      diagTitle = 'Salt plume depth based on density criterion (>0)'
      diagUnits = 'm               '
      diagCode  = 'SM      U1      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I       diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'oceSPtnd'
      diagTitle = 'salt tendency due to salt plume flux' 
     &          //' >0 increases salinity'
      diagUnits = 'g/m^2/s         '
      diagCode  = 'SM      MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I       diagName, diagCode, diagUnits, diagTitle, 0, myThid )

#ifdef ALLOW_PTRACERS
      DO iTrc=1,PTRACERS_numInUse
C--    Set default name & tracer Units:
       WRITE(locName,'(A,I4,A)') 'Tracer ',iTrc
       trcUnits  = '(kg/kg)'
C-     use name & units from data.ptracers :
       ilnb = ILNBLNK(PTRACERS_units(iTrc))
       IF ( ilnb.GE.1 ) THEN
         ilnb = LEN(trcUnits)
         trcUnits = PTRACERS_units(iTrc)(1:ilnb)
       ENDIF
       ilnb = ILNBLNK(PTRACERS_long_names(iTrc))
       IF ( ilnb.GE.1 ) THEN
         ilnb = MIN(LEN(locName),ilnb)
         WRITE(locName,'(A)') PTRACERS_long_names(iTrc)(1:ilnb)
       ELSE
         ilnb = ILNBLNK(PTRACERS_names(iTrc))
         IF ( ilnb.GE.1 ) THEN
           ilnb = MIN(LEN(locName),ilnb)
           WRITE(locName,'(A)') PTRACERS_names(iTrc)(1:ilnb)
         ENDIF
       ENDIF
       ilnb = MAX(ILNBLNK(locName),1)

       WRITE(diagName,'(A6,A2)') 'SPflxT',PTRACERS_ioLabel(iTrc)
       WRITE(diagTitle,'(2A)') 'Flux rejected during freezing, of ',
     &                 locName(1:ilnb)
       diagUnits = trcUnits//'/m^2/s'
      diagCode  = 'SM      U1      '
       CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I           diagName, diagCode, diagUnits, diagTitle, 0, myThid )          

       WRITE(diagName,'(A6,A2)') 'SPtndT',PTRACERS_ioLabel(iTrc)
      WRITE(diagTitle,'(2A)') 'Tendency due to salt plume flux, of ', 
     &                 locName(1:ilnb)
       diagUnits = trcUnits//'.m^3/s'
      diagCode  = 'SM      MR      '
       CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I           diagName, diagCode, diagUnits, diagTitle, 0, myThid )          
      ENDDO

#endif /* ALLOW_PTRACERS */
#endif /* ALLOW_DIAGNOSTICS */
#endif /* ALLOW_SALT_PLUME */

      RETURN
      END
