#include "DIC_OPTIONS.h"

C     *==========================================================*
C     | DIC_VARS.h
C     | o Abiotic Carbon Variables
C     *==========================================================*

C     AtmospCO2   :: Atmospheric pCO2 (atm).
C     AtmosP      :: Atmospheric Pressure loaded from file (atm).
C     pH          :: surface ocean pH (acidity) for pCO2 
C                       calculations.
C     pCO2        :: surface ocean partial pressure of CO2 (atm).
C     FluxCO2     :: Air-sea flux of CO2 (mol/m2/s).
C     wind        :: Wind speed loaded from file for air-sea 
C                       flux calculations (m/s).
C     FIce        :: Fraction of sea ice cover loaded from file 
C                       (or set by thice/seaice)
C                       for air-sea flux calculations.
C     Silica      :: Surface ocean concentration of silicate for 
C                       pCO2 calculations. Read in from file (mol/m3).
C     Kwexch_Pre  :: Common part of piston velocity used for 
C                       for air-sea CO2 and O2 flux calculations.
       COMMON /CARBON_NEEDS/
     &              AtmospCO2, AtmosP, pH, pCO2,
     &              co3, hco3, co2aq, FluxCO2,
     &              wind, FIce, Silica, Kwexch_Pre
      _RL  AtmospCO2(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL  AtmosP(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL  pH(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL  pCO2(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL  co3(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL  hco3(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL  co2aq(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL  FluxCO2(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL  wind(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL  FIce(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL  Silica(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL  Kwexch_Pre(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)

C Store dissociation and carbon chemistry coefficients for 
C    pCO2 solvers (see carbon_chem.F). 
C ak0 =  [H2CO2]/pCO2        
C   [Weiss 1974]
C ak1 =  [H][HCO3]/[H2CO3]
C   [Millero p.664 (1995) using Mehrbach et al. data on seawater scale]
C ak2 =  [H][CO3]/[HCO3]
C   [Millero p.664 (1995) using Mehrbach et al. data on seawater scale]
C akw =  [H][OH]
C   [Millero p.670 (1995) using composite data]
C akb =  [H][BO2]/[HBO2]
C   [Millero p.669 (1995) using data from dickson (1990)]
C aks =  [H][SO4]/[HSO4]
C   [dickson (1990, J. chem. Thermodynamics 22, 113)]
C akf =  [H][F]/[HF]
C   [dickson and Riley (1979)]
C ak1p = [H][H2PO4]/[H3PO4]
C   [DOE(1994) eq 7.2.20 with footnote using data from Millero (1974)]
C ak2p = [H][HPO4]/[H2PO4]
C   [DOE(1994) eq 7.2.23 with footnote using data from Millero (1974)]
C ak3p = [H][PO4]/[HPO4]
C   [DOE(1994) eq 7.2.26 with footnote using data from Millero (1974)]
C aksi = [H][SiO(OH)3]/[Si(OH)4]
C   [Millero p.671 (1995) using data from Yao and Millero (1995)]
C ft  = estimated fluoride concentration 
C   [Riley (1965)]
C st  = estimated sulphate concentration 
C   [Morris & Riley (1966)]
C bt  = estimated borate concentration 
C   [Uppstrom (1974)]
C fugf correct for non-ideality in ocean 
C   [Weiss (1974) Marine Chemistry]
C ff used for water vapor and pressure correction
C   [Weiss & Price (1980, Mar. Chem., 8, 347-359; Eq 13 with table 6 values)]
C Ksp_TP_Calc solubility product for calcite
C   [Following Mucci (1983) with pressure dependence from Ingle (1975)]
       COMMON /CARBON_CHEM/
     &                     ak0,ak1,ak2,akw,akb,aks,akf,
     &                     ak1p,ak2p,ak3p,aksi, fugf,
     &                     ff,ft,st,bt, Ksp_TP_Calc
      _RL  ak0(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL  ak1(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL  ak2(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL  akw(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL  akb(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL  aks(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL  akf(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL  ak1p(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL  ak2p(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL  ak3p(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL  aksi(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL  ff(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
C Fugacity Factor added by Val Bennington Nov. 2010
      _RL  fugf(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL  ft(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL  st(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL  bt(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL  Ksp_TP_Calc(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)

C Store dissociation coefficients for O2 solubility
       COMMON /OXYGEN_CHEM/
     &              oA0,oA1,oA2,oA3,oA4,oA5,
     &              oB0,oB1,oB2,oB3,
     &              oC0
      _RL oA0,oA1,oA2,oA3,oA4,oA5
      _RL oB0,oB1,oB2,oB3
      _RL oC0

C Store surface mean tracer values for virtual FW fluxes
C permil : is conversion factor for mol/m3 to mol/kg
C          assumes uniform (surface) density
C Pa2Atm : for conversion of atmospheric pressure
C          when coming from atmospheric model
       COMMON /GLOBAL_SURF_MEAN/
     &                          gsm_alk,gsm_s,gsm_t,gsm_dic,
     &                          gsm_c14,permil,Pa2Atm
      _RL  gsm_alk
      _RL  gsm_s
      _RL  gsm_t
      _RL  gsm_DIC
      _RL  gsm_C14
      _RL  permil
      _RL  Pa2Atm

C Schmidt number coefficients for air sea gas exchange of
C   carbon and oxygen
      COMMON /DIC_SCHMIDT_NO/
     &                    sca1, sca2, sca3, sca4,
     &                    sox1, sox2, sox3, sox4
      _RL  sca1
      _RL  sca2
      _RL  sca3
      _RL  sca4
      _RL  sox1
      _RL  sox2
      _RL  sox3
      _RL  sox4

C--   COMMON /DIC_FILENAMES/
C  DIC_windFile    :: file name of wind speeds
C  DIC_atmospFile  :: file name of atmospheric pressure
C  DIC_iceFile     :: file name of seaice fraction
C  DIC_ironFile    :: file name of aeolian iron flux
C  DIC_silicaFile  :: file name of surface silica
C  DIC_parFile     :: file name of photosynthetically available radiation (PAR)
C  DIC_chlaFile    :: file name of chlorophyll climatology
C  DIC_forcingPeriod :: periodic forcing parameter specific for dic (seconds)
C  DIC_forcingCycle  :: periodic forcing parameter specific for dic (seconds)
C  dic_int*          :: Handle the atmospheric boundary condition for pCO2 
C  dic_int1:
C          0=use default pCO2 (278.d-6)
C          1=use constant value of dic_pCO2 read in from data.dic
C          2=read in from file named co2atmos.dat
C          3=interact with well mixed atmospheric box (can use dic_pCO2 as initial value)
C  dic_int2          :: number pCO2 entries to read from file
C  dic_int3          :: start timestep
C  dic_int4          :: timestep between file entries
C  dic_pCO2          :: atmospheric pCO2 to be read from data.dic
C JML include iron source due to hydrothermal input into bottom layer  
      COMMON /DIC_FILENAMES/
     &        DIC_windFile, DIC_atmospFile, DIC_iceFile,
     &        DIC_ironFile, DIC_silicaFile, DIC_parFile,
     &        DIC_chlaFile, DIC_hydroventFile,
     &        DIC_forcingPeriod, DIC_forcingCycle,
     &        dic_pCO2, dic_int1, dic_int2, dic_int3, dic_int4

      CHARACTER*(MAX_LEN_FNAM) DIC_windFile
      CHARACTER*(MAX_LEN_FNAM) DIC_atmospFile
      CHARACTER*(MAX_LEN_FNAM) DIC_iceFile
      CHARACTER*(MAX_LEN_FNAM) DIC_ironFile
      CHARACTER*(MAX_LEN_FNAM) DIC_hydroventFile
      CHARACTER*(MAX_LEN_FNAM) DIC_silicaFile
      CHARACTER*(MAX_LEN_FNAM) DIC_parFile
      CHARACTER*(MAX_LEN_FNAM) DIC_chlaFile
      _RL     DIC_forcingPeriod
      _RL     DIC_forcingCycle
      _RL     dic_pCO2
      INTEGER dic_int1
      INTEGER dic_int2
      INTEGER dic_int3
      INTEGER dic_int4

#ifdef DIC_BIOTIC
C     *==========================================================*
C     | o Biological Carbon Variables
C     *==========================================================*

C For averages of the output using TIMEAVE
C  BIOave                 :: biological productivity [mol P/m3/s]
C  CARave                 :: carbonate changes due to biological
C                             productivity and remineralization [mol C/m3/s]
C  SURave                 :: tendency of DIC due to air-sea exchange
C                            and virtual flux [mol C/m3/s]
C  SUROave                :: tendency of O2 due to air-sea exchange [mol O2/m3/s]
C  pCO2ave                :: surface ocean pCO2 [uatm]
C  pHave                  :: surface ocean pH
C  fluxCO2ave             :: Air-sea flux of CO2 [mol C/m2/s] 
C  omegaCave              :: Local saturation state with respect to calcite 
C  pfluxave               :: changes to PO4 due to flux and remineralization [mol P/m3/s] 
C  epfluxave              :: export flux of PO4 through each layer [mol P/m2/s] 
C  cfluxave               :: carbonate changes due to flux and remineralization [mol C/m3/s]
C  DIC_timeAve            :: period over which DIC averages are calculated [s]

C Values for biotic biogeochemistry 
C          (many set in data.dic, defaults to values in dic_readparms.F)
C  par                    :: photosynthetically active radiation  (light available for 
C                            phytoplankton growth) [W/m2]
C  alphaUniform, alpha    :: maximum rate of biological activity [mol P/m3/s]
C                            read in alphaUniform and filled in 2d array alpha
C  rainRatioUniform, rain_ratio :: inorganic/organic particulate carbon rain ratio (PIC/POC)
C                            read in rainRatioUniform and filled in 2d array rain_ratio
C  InputFe                :: aeolian deposition of TOTAL IRON in dust [mol/m2/s]
C  omegaC                 :: Local saturation state with respect to calcite 
C  CHL                    :: chlorophyll climatology data for self-shading effect [mg/m3]
C  Kpo4, KFE, lit0        :: half saturation constants for phosphate, iron and light 
C							 [mol P/m3, mol Fe/m3, W/m2]
C  DOPfraction            :: fraction of new production going to DOP  
C  zcrit, nlev            :: Minimum Depth (m) over which biological activity
C                            is computed. determines nlev as the indice of the
C                            first layer deeper than -zcrit
C  KRemin                 :: remineralization power law coeffient
C  KDOPremin              :: DOP remineralization rate [1/s]
C  zca                    :: scale depth for CaCO3 remineralization [m] 
C  R_op,R_cp,R_NP, R_FeP  :: stochiometric ratios of nutrients 
C                            (assumption of stoichometry of plankton and particulate 
C                             and dissolved organic matter)
C  O2crit                 :: critical oxygen level [mol/m3] 
C  alpfe                  :: solubility of aeolian fe [fraction]
C  KScav                  :: iron scavenging rate [1/s]
C  ligand_stab            :: ligand-free iron stability constant [m3/mol]
C  ligand_tot             :: uniform, invariant total free ligand conc [mol/m3]
C  freefemax              :: max soluble free iron [mol/m3]
C  fesedflux_pcm          :: ratio of sediment iron to sinking organic matter
C  FeIntSec               :: Sediment Fe flux, intersect value in:
C                                    Fe_flux = fesedflux_pcm*pflux + FeIntSec
C  parfrac                :: fraction of Qsw that is PAR
C  k0                     :: light attentuation coefficient of water [1/m]
C  kchl                   :: light attentuation coefficient of chlorophyll [m2/mg]
C  alphamax, alphamin     :: not used (legacy adjoint param)
C  calpha                 :: not used (legacy adjoint param)
C  crain_ratio            :: not used (legacy adjoint param)
C  cInputFe               :: not used (legacy adjoint param)
C  calpfe                 :: not used (legacy adjoint param)
C  feload                 :: not used (legacy adjoint param)
C  cfeload                :: not used (legacy adjoint param)
C  QSW_underice           :: is Qsw is masked by ice fraction?
C JML Added some extras for particle dependent Fe scavenging rates
C       and prognostic Ligands
      COMMON /BIOTIC_NEEDS/
     &     BIOave, CARave, SURave, SUROave, pCO2ave, pHave,
     &     fluxCO2ave, omegaCave, pfluxave, epfluxave, cfluxave,
     &     DIC_timeAve,
     &     alpha, rain_ratio, InputFe, omegaC, CHL,
     &     Kpo4, DOPfraction, zcrit, KRemin, O2crit,
     &     KDOPremin,zca,R_op,R_cp,R_NP, R_FeP, R_SIP,
     &     alpfe, ligand_stab, ligand_tot, KFE, freefemax, 
     &     KScav, KScav_yr, KScav_surf, KScav_poc_yr, KScav_poc,
     &     R_pop2poc, R_dust2fe, poc_wsp,  dust_wsp,
     &     scav_exp, scav_ratio, KScav_dust_yr, 
     &     KScav_dust, KScav_background, KScav_background_yr,
     &     sedfe_factor, sedfe_offset, fe_sed_depth_max,
     &     HydroInputHe3, solfe, R_FeHe3, fe_vent_depth_min,
     &     par, parfrac, k0, kchl, lit0,
     &     alphaUniform_yr, alphaUniform, rainRatioUniform,
     &     alphamax, alphamin,
     &     calpha, crain_ratio, cInputFe, calpfe, feload, cfeload,
     &     nlev, QSW_underice,
     &     selectFeSedSource, selectFeHTVSource
     
      INTEGER nlev
      INTEGER selectFeSedSource, selectFeHTVSource
      
C     For averages
      _RL BIOave(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL CARave(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL SURave(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL SUROave(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL pCO2ave(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL pHave(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL fluxCO2ave(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL OmegaCave(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL pfluxave(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL epfluxave(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL cfluxave(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL DIC_timeAve(nSx,nSy)

C     values for biogeochemistry
C JML Added some extras for particle dependent Fe scavenging rates,
C       prognostic Ligands, and Elrod (via Hajoon) sediment source
C   CHL           :: chlorophyll climatology [mg/m3]
C   sedfe_factor :: ratio of sediment iron to sinking organic matter
C   sedfe_offset      :: Sediment Fe flux, intersect value in:
C                    Fe_flux = sedfe_factor*pflux + sedfe_offset
      _RL par(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL alpha(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL rain_ratio(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL InputFe(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL HydroInputHe3(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL omegaC(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL CHL(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL Kpo4
      _RL DOPfraction
      _RL zcrit
      _RL KRemin
      _RL KDOPremin
      _RL zca
      _RL R_op
      _RL R_cp
      _RL R_NP
      _RL R_FeP
      _RL R_SIP
      _RL O2crit
      _RL alpfe
      _RL ligand_stab
      _RL ligand_tot
      _RL KFe
      _RL freefemax
C Scavenging values
      _RL KScav_yr
      _RL KScav
      _RL KScav_surf
      _RL KScav_poc_yr
      _RL KScav_poc
      _RL R_pop2poc
      _RL R_dust2fe
      _RL poc_wsp
      _RL dust_wsp
      _RL scav_exp
      _RL scav_ratio
      _RL KScav_dust_yr
      _RL KScav_dust
      _RL KScav_background_yr
      _RL KScav_background
      _RL sedfe_factor
      _RL sedfe_offset
      _RL R_FeHe3
      _RL solfe
      _RL fe_vent_depth_min
      _RL fe_sed_depth_max
C     values for light limited bio activity
C   k0      :: Light attentuation coefficient for water [1/m]
C   kchl    :: Light attentuation coefficient fct of chlorophyll [m2/mg]
      _RL k0, kchl, parfrac, lit0
      _RL alphaUniform_yr
      _RL alphaUniform
      _RL rainRatioUniform
      _RL alphamax, alphamin
      _RL calpha
      _RL crain_ratio
      _RL cInputFe(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL calpfe
      _RL cfeload
      _RL feload(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)

      LOGICAL QSW_underice

#ifdef ALLOW_VARIABLE_LIGANDS
      COMMON /PROGNOSTIC_LIGANDS/
     &     gamma_lig, lambda_over_gamma, lambda_ref,
     &     lig_q10, lig_tref
     
      _RL gamma_lig
      _RL lambda_over_gamma
      _RL lambda_ref
      _RL lig_q10
      _RL lig_tref
#endif
      
#endif /* DIC_BIOTIC */

CEH3 ;;; Local Variables: ***
CEH3 ;;; mode:fortran ***
CEH3 ;;; End: ***
