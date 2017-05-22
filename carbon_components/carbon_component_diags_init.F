C $Header: /u/gcmpack/MITgcm/pkg/dic/dic_diagnostics_init.F,v 1.1 2010/05/14 21:52:33 jahn Exp $
C $Name:  $

#include "GCHEM_OPTIONS.h"
#include "DIC_OPTIONS.h"
C !INTERFACE: ==========================================================
      SUBROUTINE CARBON_COMPONENT_DIAGS_INIT(myThid )

C !DESCRIPTION:
C define diagnostics for dic package
C experiment

C !USES: ===============================================================
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"

C !INPUT PARAMETERS: ===================================================
C  myThid               :: thread number
      INTEGER myThid
CEOP

#ifdef ALLOW_DIAGNOSTICS

C     !LOCAL VARIABLES:
      INTEGER       diagNum
      INTEGER       diagMate
      CHARACTER*8   diagName
      CHARACTER*16  diagCode
      CHARACTER*16  diagUnits
      CHARACTER*(80) diagTitle

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
C     Define diagnostics Names :
#if ((defined DIC_CARBON_COMPONENTS) && (defined DIC_COMPONENT_DIAGS))
      diagName  = 'CSATINIT'
      diagTitle = 'Saturated DIC at initial atm pCO2'
      diagUnits = 'mol C/m3      '
      diagCode  = 'SM P    MR    '
      diagMate  = diagNum + 2
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I       diagName, diagCode, diagUnits, diagTitle, diagMate, myThid)

      diagName  = 'CRESINIT'
      diagTitle = 'Residual DIC at initial atm pCO2'
      diagUnits = 'mol C/m3      '
      diagCode  = 'SM      MR    '
      diagMate  = diagNum
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I       diagName, diagCode, diagUnits, diagTitle, diagMate, myThid)

      diagName  = 'CSATSIT '
      diagTitle = 'Saturated DIC at atm pCO2 of last surface contact'
      diagUnits = 'mol C/m3      '
      diagCode  = 'SM P    MR    '
      diagMate  = diagNum + 2
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I       diagName, diagCode, diagUnits, diagTitle, diagMate, myThid)

      diagName  = 'CRESSIT '
      diagTitle = 'Residual DIC at atm pCO2 of last surface contact'
      diagUnits = 'mol C/m3      '
      diagCode  = 'SM      MR    '
      diagMate  = diagNum
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I       diagName, diagCode, diagUnits, diagTitle, diagMate, myThid)

      diagName  = 'CSAT    '
      diagTitle = 'Saturated DIC at current atm pCO2'
      diagUnits = 'mol C/m3      '
      diagCode  = 'SM P    MR    '
      diagMate  = diagNum + 2
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I       diagName, diagCode, diagUnits, diagTitle, diagMate, myThid)
     
      diagName  = 'CRES    '
      diagTitle = 'Residual DIC at current atm pCO2'
      diagUnits = 'mol C/m3      '
      diagCode  = 'SM      MR    '
      diagMate  = diagNum
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I       diagName, diagCode, diagUnits, diagTitle, diagMate, myThid)
     
      diagName  = 'CSOFT   '
      diagTitle = 'Soft-tissue carbon concentration'
      diagUnits = 'mol C/m3      '
      diagCode  = 'SM      MR    '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I       diagName, diagCode, diagUnits, diagTitle, 0, myThid )
     
      diagName  = 'CCARB   '
      diagTitle = 'Carbonate concentration'
      diagUnits = 'mol C/m3      '
      diagCode  = 'SM      MR    '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I       diagName, diagCode, diagUnits, diagTitle, 0, myThid )
     
      diagName  = 'PREG    '
      diagTitle = 'Regenerated Phosphate concentration'
      diagUnits = 'mol P/m3      '
      diagCode  = 'SM      MR    '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I       diagName, diagCode, diagUnits, diagTitle, 0, myThid )
     
      diagName  = 'AREG    '
      diagTitle = 'Regenerated Alkalinity concentration'
      diagUnits = 'mol eq/m3      '
      diagCode  = 'SM      MR    '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I       diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'SIPRE   '
      diagTitle = 'Preformed Silica concentration'
      diagUnits = 'mol Si/m3'
      diagCode  = 'SM      MR    '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I       diagName, diagCode, diagUnits, diagTitle, 0, myThid )
     
      diagName  = 'Silica  '
      diagTitle = 'Non-interactive 3d silica field'
      diagUnits = 'mol Si/m3 '
      diagCode  = 'SM P    MR    '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I       diagName, diagCode, diagUnits, diagTitle, 0, myThid )
     
       diagName  ='pH_GLOB '
      diagTitle = 'Global contemporary pH field'
      diagUnits = ' '
      diagCode  = 'SM P    MR    '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I       diagName, diagCode, diagUnits, diagTitle, 0, myThid )
     
      diagName  = 'pH_INIT  '
      diagTitle = 'Global pH field at initial pCO2'
      diagUnits = ' '
      diagCode  = 'SM P    MR    '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I       diagName, diagCode, diagUnits, diagTitle, 0, myThid )
     
      diagName  = 'pH_INSIT'
      diagTitle = 'Global pH field at in situ pCO2'
      diagUnits = ' '
      diagCode  = 'SM P    MR    '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I       diagName, diagCode, diagUnits, diagTitle, 0, myThid )     
#endif
C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

#endif /* ALLOW_DIAGNOSTICS */

      RETURN
      END
