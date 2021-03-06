#include "DIC_OPTIONS.h"

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP
C     !ROUTINE: DIC_INIT_FIXED

C     !INTERFACE:
      SUBROUTINE DIC_INIT_FIXED( myThid )

C     !DESCRIPTION:
C     Initialize fixed quantities

C     !USES:
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "DIC_VARS.h"
#include "DIC_ATMOS.h"
#include "DIC_LOAD.h"
#ifdef DIC_CARBON_COMPONENTS
#include "PTRACERS_SIZE.h"
#include "PTRACERS_FIELDS.h"
#include "PTRACERS_PARAMS.h"
#endif
C     !INPUT PARAMETERS:
C     myThid       :: my Thread Id number
      INTEGER myThid
CEOP

#ifdef ALLOW_DIC
      INTEGER k
      CHARACTER*(MAX_LEN_MBUF) msgBuf
      CHARACTER*(MAX_LEN_MBUF) fldName
      INTEGER iUnit, nTrac

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

      _BEGIN_MASTER(myThid)

C set up coefficients for DIC chemistry
C define Schmidt no. coefficients for CO2
      sca1 = 2073.1 _d 0
      sca2 = -125.62 _d 0
      sca3 =    3.6276 _d 0
      sca4 =   -0.043219 _d 0
C define Schmidt no. coefficients for O2
C based on Keeling et al [GBC, 12, 141, (1998)]
      sox1 = 1638.0 _d 0
      sox2 = -81.83 _d 0
      sox3 =    1.483 _d 0
      sox4 =   -0.008004 _d 0

C coefficients for determining saturation O2
      oA0=  2.00907 _d 0
      oA1=  3.22014 _d 0
      oA2=  4.05010 _d 0
      oA3=  4.94457 _d 0
      oA4= -2.56847 _d -1
      oA5=  3.88767 _d 0
      oB0= -6.24523 _d -3
      oB1= -7.37614 _d -3
      oB2= -1.03410 _d -2
      oB3= -8.17083 _d -3
      oC0= -4.88682 _d -7

C     Set other constant/flag

C                   dic_int2=number entries to read
C                   dic_int3=start timestep,
C                   dic_int4=timestep between file entries)
      IF ( dic_int1.EQ.2.OR.dic_int1.eq.4) then
        CALL MDSFINDUNIT( iUnit, mythid )
        OPEN(UNIT=iUnit,FILE='co2atmos.dat',STATUS='old')
        DO k=1,dic_int2
          READ(iUnit,*) co2atmos(k)
          WRITE(standardMessageUnit,*) 'co2atmos',co2atmos(k)
        ENDDO
        CLOSE(iUnit)
      ENDIF

#ifdef DIC_BIOTIC

#ifdef ALLOW_FE
#ifdef ALLOW_VARIABLE_LIGANDS
C Given the ratio of steady state bio production divided by the steady state
C  ligand concentration (i.e. real world estimates for gamma and lambda) then
C  find the reference timescale of ligand degredation, lambda_ref.
        IF (gamma_lig.NE.0. _d 0.AND.
     &     lambda_over_gamma.NE.0. _d 0) THEN
           lambda_ref=lambda_over_gamma*gamma_lig
        ELSE
           WRITE(msgBuf,'(A)')
     &      'DIC_INIT_FIXED: gamma_lig and gamma_over_lambda are zero'
           CALL PRINT_ERROR( msgBuf, myThid )
           STOP 'ABNORMAL END: S/R DIC_INIT_FIXED: 
     &       variable ligand error' 
        ENDIF
#endif /* ALLOW_VARIABLE_LIGANDS */
#endif /* ALLOW_FE */

#ifdef ALLOW_AUTODIFF
cph nlev=Nr used in TAF store directive
      nlev=Nr
#else
      k = 1
      DO WHILE ( k .LT. Nr .AND. rC(k) .GE. -zcrit )
        k = k+1
      ENDDO
      nlev = k
#endif
C--   Print out nlev values :
      iUnit = standardMessageUnit
      WRITE(msgBuf,'(A)') ' '
      CALL PRINT_MESSAGE(msgBuf,iUnit,SQUEEZE_RIGHT,myThid)
      WRITE(msgBuf,'(A)') '// ==================================='
      CALL PRINT_MESSAGE(msgBuf,iUnit,SQUEEZE_RIGHT,myThid)
      WRITE(msgBuf,'(A)') '//  DIC_INIT_FIXED parameters :'
      CALL PRINT_MESSAGE(msgBuf,iUnit,SQUEEZE_RIGHT,myThid)

       CALL WRITE_0D_I( nlev, INDEX_NONE, 'nlev =',
     &   ' /* Number of level over which Bio act is computed */')
#endif /* DIC_BIOTIC */

      _END_MASTER(myThid)

#ifdef ALLOW_MNC
      IF ( useMNC ) THEN
        CALL DIC_MNC_INIT( myThid )
      ENDIF
#endif /* ALLOW_MNC */

#ifdef DIC_CARBON_COMPONENTS
C Find the locations of tracers in the ptracers array
       DO nTrac=1,PTRACERS_numInUse
        fldName=PTRACERS_names(nTrac)
        
        SELECT CASE (fldName(1:3))
         CASE ("dic")
          iDIC=nTrac
          CALL WRITE_0D_I( iDIC, INDEX_NONE, 'iDIC =',
     &   ' /* Index of the DIC tracer */')          
         CASE ("alk")
          iALK=nTrac
         CALL WRITE_0D_I( iALK, INDEX_NONE, 'iALK =',
     &   ' /* Index of the ALK tracer */')          
         CASE ("po4")
          iPO4=nTrac   
         CALL WRITE_0D_I( iPO4, INDEX_NONE, 'iPO4 =',
     &   ' /* Index of the PO4 tracer */')          
         CASE ("dop")
          iDOP=nTrac    
         CALL WRITE_0D_I( iDOP, INDEX_NONE, 'iDOP =',
     &   ' /* Index of the DOP tracer */')          
         CASE ("o2")
          iDO2=nTrac     
         CALL WRITE_0D_I( iDO2, INDEX_NONE, 'iDO2 =',
     &   ' /* Index of the O2 tracer */')          
         CASE ("fe")
          iDFE=nTrac
         CALL WRITE_0D_I( iDFE, INDEX_NONE, 'iDFE =',
     &   ' /* Index of the FE tracer */')          
         CASE ("cpre")
          iCPRE=nTrac
         CALL WRITE_0D_I( iCPRE, INDEX_NONE, 'iCPRE =',
     &   ' /* Index of the CPRE tracer */')          
         CASE ("apre")
          iAPRE=nTrac   
         CALL WRITE_0D_I( iAPRE, INDEX_NONE, 'iAPRE =',
     &   ' /* Index of the APRE tracer */')          
         CASE ("ppre")
          iPPRE=nTrac 
         CALL WRITE_0D_I( iPPRE, INDEX_NONE, 'iPPRE =',
     &   ' /* Index of the PPRE tracer */')          
         CASE ("opre")
          iOPRE=nTrac
         CALL WRITE_0D_I( iOPRE, INDEX_NONE, 'iOPRE =',
     &   ' /* Index of the OPRE tracer */')          
         CASE ("wm_age")
          iAGE=nTrac  
         CALL WRITE_0D_I( iAGE, INDEX_NONE, 'iAGE =',
     &   ' /* Index of the AGE tracer */')          
         CASE ("atmpco2")
          iPCO2=nTrac
         CALL WRITE_0D_I( iPCO2, INDEX_NONE, 'iPCO2 =',
     &   ' /* Index of the PCO2 tracer */')          
         CASE ("csat")
          iCSAT=nTrac
         CALL WRITE_0D_I( iCSAT, INDEX_NONE, 'iCSAT =',
     &   ' /* Index of the CSAT tracer */')          
         CASE ("cdis")
          iCDIS=nTrac
         CALL WRITE_0D_I( iCDIS, INDEX_NONE, 'iCDIS =',
     &   ' /* Index of the CDIS tracer */')          
         END SELECT
        ENDDO 
#endif

#ifdef ALLOW_DIAGNOSTICS
      IF ( useDIAGNOSTICS ) THEN
        CALL DIC_DIAGNOSTICS_INIT( myThid )
        CALL DIC_COMPONENT_DIAGS_INIT( myThid )
      ENDIF
#endif /* ALLOW_DIAGNOSTICS */

#endif /* ALLOW_DIC */

      RETURN
      END
