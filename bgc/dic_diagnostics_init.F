#include "GCHEM_OPTIONS.h"
#include "DIC_OPTIONS.h"

C !INTERFACE: ==========================================================
      SUBROUTINE DIC_DIAGNOSTICS_INIT(myThid )

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
      CHARACTER*8   diagName
      CHARACTER*16  diagCode
      CHARACTER*16  diagUnits
      CHARACTER*(80) diagTitle

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
C     Define diagnostics Names :

      diagName  = 'DICBIOA '
      diagTitle = 'Biological Productivity (mol/m3/s)'
      diagUnits = 'mol/m3/sec      '
      diagCode  = 'SMRP    MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I       diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'DICCARB '
      diagTitle = 'Carbonate chg-biol prod and remin (mol eq/m3/s)'
      diagUnits = 'mol eq/m3/sec   '
      diagCode  = 'SMRP    MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I       diagName, diagCode, diagUnits, diagTitle, 0, myThid )
     
      diagName  = 'OMEGAC  '
      diagTitle = 'Carbonate saturation'
      diagUnits = 'mol eq/m3/s   '
      diagCode  = 'SM      MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I       diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'DICPFLUX'
      diagTitle = 'Phosphate flux'
      diagUnits = 'mol P/m3/s   '
      diagCode  = 'SM P    MR    '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I       diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'CEXPFLX '
      diagTitle = 'Particulate organic carbon export flux'
      diagUnits = 'mol C/m3/s   '
      diagCode  = 'SM P    MR    '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I       diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'CARBFLUX'
      diagTitle = 'Carbonate export flux'
      diagUnits = 'mol C/m3/s   '
      diagCode  = 'SM P    MR    '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I       diagName, diagCode, diagUnits, diagTitle, 0, myThid )
     
      diagName  = 'DICTFLX '
      diagTitle = 'Tendency of DIC due to air-sea exch (mol/m3/s)'
      diagUnits = 'mol/m3/s      '
      diagCode  = 'SM P    L1      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I       diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'DICCFLX '
      diagTitle = 'Flux of CO2 - air-sea exch (mol/m2/s)'
      diagUnits = 'mol/m2/s      '
      diagCode  = 'SM P    L1      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I       diagName, diagCode, diagUnits, diagTitle, 0, myThid )
     
      diagName  = 'DICATCAR'
      diagTitle = 'Atmospheric carbon inventory (mol C)'
      diagUnits = 'mol C    '
      diagCode  = 'SM P    L1      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I       diagName, diagCode, diagUnits, diagTitle, 0, myThid )
     
      diagName  = 'DICATCO2'
      diagTitle = 'Atmospheric pCO2 (atm)'
      diagUnits = 'atm      '
      diagCode  = 'SM P    L1      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I       diagName, diagCode, diagUnits, diagTitle, 0, myThid )
     
      diagName  = 'DICAFLX '
      diagTitle = 'Tendency of ALK due to air-sea exch (mol/m3/s)'
      diagUnits = 'mol/m3/s      '
      diagCode  = 'SM P    L1      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I       diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'DICRDOP '
      diagTitle = 'Remineralization of DOP'
      diagUnits = 'mol P/m3/s   '
      diagCode  = 'SM      MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I       diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'DICGDIC '
      diagTitle = 'Tendency of DIC concentration'
      diagUnits = 'mol C/m3/s   '
      diagCode  = 'SM      MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I       diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'DICGDOP '
      diagTitle = 'Tendency of DOP concentration'
      diagUnits = 'mol P/m3/s   '
      diagCode  = 'SM      MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I       diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'DICGPO4 '
      diagTitle = 'Tendency of PO4 concentration'
      diagUnits = 'mol P/m3/s   '
      diagCode  = 'SM      MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I       diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'DICGALK '
      diagTitle = 'Tendency of ALK concentration'
      diagUnits = 'mol eq/m3/s   '
      diagCode  = 'SM      MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I       diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'DICCFLXK'
      diagTitle = 'CO2 air sea flux coefficient (m/s)'
      diagUnits = 'm/s      '
      diagCode  = 'SM P    L1      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I       diagName, diagCode, diagUnits, diagTitle, 0, myThid )
     
      diagName  = 'DICFLXCA'
      diagTitle = 'CO2 Flux, Atmospheric conc'
      diagUnits = 'mol/m3      '
      diagCode  = 'SM P    L1      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I       diagName, diagCode, diagUnits, diagTitle, 0, myThid )
     
      diagName  = 'DICFLXCO'
      diagTitle = 'CO2 Flux, Ocean conc'
      diagUnits = 'mol/m3      '
      diagCode  = 'SM P    L1      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I       diagName, diagCode, diagUnits, diagTitle, 0, myThid )
	 
      diagName  = 'DICVCFLX'
      diagTitle = 'Flux of CO2 - Virtual Flux'
      diagUnits = 'mol/m3/s      '
      diagCode  = 'SM P    L1      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I       diagName, diagCode, diagUnits, diagTitle, 0, myThid )
     
      diagName  = 'DICPCO2 '
      diagTitle = 'Partial Pressure of CO2 (atm)'
      diagUnits = 'atm             '
      diagCode  = 'SM P    M1      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I       diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'DICPHAV '
      diagTitle = 'pH (dimensionless)'
      diagUnits = 'dimensionless   '
      diagCode  = 'SM P    M1      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I       diagName, diagCode, diagUnits, diagTitle, 0, myThid )

#ifdef ALLOW_O2
      diagName  = 'DICOFLX '
      diagTitle = 'Tendency of O2 due to air-sea exch (mol/m3/s)'
      diagUnits = 'mol/m3/s      '
      diagCode  = 'SM      L1      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I       diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'DICGO2  '
      diagTitle = 'Tendency of O2 concentration'
      diagUnits = 'mol O2/m3/s      '
      diagCode  = 'SM      MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I       diagName, diagCode, diagUnits, diagTitle, 0, myThid )
#endif

#ifdef ALLOW_FE     
      diagName  = 'DICFEFLX'
      diagTitle = 'Deposition of Aeolian Fe in dust (mol/m3/s)'
      diagUnits = 'mol/m3/s      '
      diagCode  = 'SM P    L1      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I       diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'DICGFE  '
      diagTitle = 'Tendency of FE concentration'
      diagUnits = 'mol FE/m3/s      '
      diagCode  = 'SM      MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I       diagName, diagCode, diagUnits, diagTitle, 0, myThid )

#ifdef SUBSURFACE_FE_SOURCES
C JML include iron source due to sediment input into bottom layer            
      diagName  = 'DICSEDFE'
      diagTitle = 'Input of sedimentary Fe (mol/m3/s)'
      diagUnits = 'mol/m3/s      '
      diagCode  = 'SM P    L1      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I       diagName, diagCode, diagUnits, diagTitle, 0, myThid )

C JML include iron source due to hydrothermal input into bottom layer                 
      diagName  = 'DICHTVFE'
      diagTitle = 'Input of Hydrothermal Vent Fe (mol/m3/s)'
      diagUnits = 'mol/m3/s      '
      diagCode  = 'SM P    L1      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I       diagName, diagCode, diagUnits, diagTitle, 0, myThid )
#endif /* SUBSURFACE_FE_SOURCES */

#ifdef ALLOW_VARIABLE_LIGANDS
      diagName  = 'DICGLIG '
      diagTitle = 'Tendency of Ligand concentration'
      diagUnits = 'mol Lig/m3/s      '
      diagCode  = 'SM      MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I       diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'LIGLAMDA'
      diagTitle = 'Timescale of ligand degradation'
      diagUnits = '1/s      '
      diagCode  = 'SM      MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I       diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'LIGINCR '
      diagTitle = 'Ligand production rate'
      diagUnits = 'mol Lig/m3/s      '
      diagCode  = 'SM      MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I       diagName, diagCode, diagUnits, diagTitle, 0, myThid )
     
      diagName  = 'LIGDECR '
      diagTitle = 'Ligand degradation rate'
      diagUnits = 'mol Lig/m3/s      '
      diagCode  = 'SM      MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I       diagName, diagCode, diagUnits, diagTitle, 0, myThid )
#endif /* ALLOW_VARIABLE_LIGANDS */

#endif /* ALLOW_FE */
C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

#endif /* ALLOW_DIAGNOSTICS */

      RETURN
      END
