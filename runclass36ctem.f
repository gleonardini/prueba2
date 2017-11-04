      PROGRAM RUNCLASS36CTEM

C Esto es una prueba
C     REVISION HISTORY:
C
C     * JUL 2 2015
C     * JOE MELTON : Took many calculations out of this driver and into subroutines. Introduced
C                    modular structure and made CTEM vars into pointers. Harmonized CLASS v. 3.6.1
C                    code with CTEM code and into this driver.
C
C     * JAN 14 2014
C     * JOE MELTON : Harmonized the field capacity and wilting point calculations between CLASS and CTEM.
C                    took the code out of runclassctem and it is now fully done within CLASSB. Harmonized names too.
C
C     * JUN 2014
C     * RUDRA SHRESTHA : ADD IN WETLAND CODE
C
C     * JUL 2013
C     * JOE MELTON : REMOVED CTEM1 AND CTEM2 OPTIONS, REPLACED WITH CTEM_ON. INTRODUCE
C                                  MODULES. RESTRUCTURE OUTPUTS AND CTEM VARIABLE DECLARATIONS
C                                  OUTPUTS ARE NOW THE SAME FOR BOTH MOSAIC AND COMPOSITE MODES
C
C     * DEC 2012
C     * JOE MELTON : REMOVED GOTO STATEMENTS, CLEANED UP AND FIXED INCONSISTENCIES
C                                  IN HOW INPUT DATA READ IN. ALSO MADE LUC WORK FOR
C                                  BOTH COMPOSITE AND MOSAIC APPROACHES.
C
C     * OCT 2012
C     * YIRAN PENG AND JOE MELTON: BRING IN COMPETITION TO 3.6 AND MAKE IT
C                                  SO THE MODEL CAN START FROM BARE GROUND
C                                  OR FROM THE INI AND CTM FILES INPUTS
C
C     * SEP 2012
C     * JOE MELTON: COUPLED CLASS3.6 AND CTEM
C
C     * NOV 2011
C     * YIRAN PENG AND VIVEK ARORA: COUPLED CLASS3.5 AND CTEM
C
C     * SEPT 8, 2009
C     * RONG LI AND VIVEK ARORA: COUPLED CLASS3.4 AND CTEM
C
C=======================================================================

C     * DIMENSION STATEMENTS.

C     * FIRST SET OF DEFINITIONS:
C     * BACKGROUND VARIABLES, AND PROGNOSTIC AND DIAGNOSTIC
C     * VARIABLES NORMALLY PROVIDED BY AND/OR USED BY THE GCM.
C     * THE SUFFIX "ROT" REFERS TO VARIABLES EXISTING ON THE
C     * MOSAIC GRID ON THE CURRENT LATITUDE CIRCLE.  THE SUFFIX
C     * "GAT" REFERS TO THE SAME VARIABLES AFTER THEY HAVE UNDERGONE
C     * A "GATHER" OPERATION IN WHICH THE TWO MOSAIC DIMENSIONS
C     * ARE COLLAPSED INTO ONE.  THE SUFFIX "ROW" REFERS BOTH TO
C     * GRID-CONSTANT INPUT VARIABLES. AND TO GRID-AVERAGED
C     * DIAGNOSTIC VARIABLES.
C
C     * THE FIRST DIMENSION ELEMENT OF THE "ROT" VARIABLES
C     * REFERS TO THE NUMBER OF GRID CELLS ON THE CURRENT
C     * LATITUDE CIRCLE.  IN THIS STAND-ALONE VERSION, THIS
C     * NUMBER IS ARBITRARILY SET TO THREE, TO ALLOW UP TO THREE
C     * SIMULTANEOUS TESTS TO BE RUN.  THE SECOND DIMENSION
C     * ELEMENT OF THE "ROT" VARIABLES REFERS TO THE MAXIMUM
C     * NUMBER OF TILES IN THE MOSAIC.  IN THIS STAND-ALONE
C     * VERSION, THIS NUMBER IS SET TO EIGHT.  THE FIRST
C     * DIMENSION ELEMENT IN THE "GAT" VARIABLES IS GIVEN BY
C     * THE PRODUCT OF THE FIRST TWO DIMENSION ELEMENTS IN THE
C     * "ROT" VARIABLES.

C     The majority of CTEM parameters are stored in ctem_params.f90. We access them
c     through use statements for modules:
      use ctem_params,        only : initpftpars,nlat,nmos,ilg,nmon,
     1                               ican, ignd,icp1, icc, iccp1,
     2                               monthend, mmday,modelpft, l2max,
     3                                deltat, abszero, monthdays,seed,
     4                                crop, NBS

      use landuse_change,     only : initialize_luc, readin_luc

      use ctem_statevars,     only : vrot,vgat,c_switch,initrowvars,
     1                               class_out,resetclassmon,
     2                               resetclassyr,resetmidmonth,
     3                               resetmonthend_g,ctem_grd_mo,
     4                               resetyearend_g,ctem_grd_yr,
     5                               resetclassaccum,ctem_grd,
     6                               ctem_tile, ctem_tile_mo,
     7                               ctem_tile_yr,resetgridavg,
     8                               resetmonthend_m,resetyearend_g,
     9                               resetyearend_m

      use io_driver,          only : read_from_ctm, create_outfiles,
     1                               write_ctm_rs, class_monthly_aw,
     2                               ctem_annual_aw,ctem_monthly_aw,
     3                               close_outfiles,ctem_daily_aw

c
      implicit none
C
C     * INTEGER CONSTANTS.
C
      INTEGER IDISP,IZREF,ISLFD,IPCP,IWF,IPAI,IHGT,IALC,
     1        IALS,IALG,N,ITG,ITC,ITCG,isnoalb,igralb

      INTEGER NLTEST,NMTEST,NCOUNT,NDAY,
     1        IMONTH,NDMONTH,NT,
     2        IHOUR,IMIN,IDAY,IYEAR,NML,NMW,JLAT,
     3        NLANDCS,NLANDGS,NLANDC,NLANDG,NLANDI,I,J,K,L,M,
     4        NTLD
C
      INTEGER IHOUR2,IMIN2,IDAY2,IYEAR2
C
      INTEGER K1,K2,K3,K4,K5,K6,K7,K8,K9,K10,K11,ITA,ITCAN,ITD,
     1        ITAC,ITS,ITSCR,ITD2,ITD3,ITD4,ISTEPS,NFS,NDRY,NAL,NFT
      REAL(KIND = 8)TAHIST(200),TCHIST(200),TACHIST(200),TDHIST(200),
     1     TSHIST(200),TSCRHIST(200),TD2HIST(200),TD3HIST(200),
     2     TD4HIST(200),PAICAN(ILG)

      INTEGER*4 TODAY(3), NOW(3)
C
C     * LAND SURFACE PROGNOSTIC VARIABLES.
C
      REAL(KIND = 8),DIMENSION(NLAT,NMOS,IGND) ::
     1        TBARROT,   THLQROT,   THICROT
C
      REAL(KIND = 8),DIMENSION(NLAT,NMOS) ::
     1        TPNDROT,   ZPNDROT,   TBASROT,
     2        ALBSROT,   TSNOROT,   RHOSROT,
     3        SNOROT ,   TCANROT,   RCANROT,
     4        SCANROT,   GROROT ,   CMAIROT,
     5        TACROT ,   QACROT ,   WSNOROT,
     6        REFROT,    BCSNROT
C
      REAL(KIND = 8)    TSFSROT(NLAT,NMOS,4)
C
      REAL(KIND = 8),DIMENSION(ILG,IGND) ::
     1        TBARGAT, THLQGAT, THICGAT
C
      REAL(KIND = 8),DIMENSION(ILG) ::
     1        TPNDGAT,   ZPNDGAT,   TBASGAT,
     2        ALBSGAT,   TSNOGAT,   RHOSGAT,
     3        SNOGAT ,   TCANGAT,   RCANGAT,
     4        SCANGAT,   GROGAT ,   CMAIGAT,
     5        TACGAT ,   QACGAT ,   WSNOGAT,
     6        REFGAT,    BCSNGAT

C
      REAL(KIND = 8)    TSFSGAT(ILG,4)
C
C     * GATHER-SCATTER INDEX ARRAYS.
C
      INTEGER  ILMOS (ILG),JLMOS(ILG),IWMOS(ILG),JWMOS (ILG)
C
C     * CANOPY AND SOIL INFORMATION ARRAYS.
C     * (THE LAST DIMENSION OF MOST OF THESE ARRAYS IS GIVEN BY
C     * THE NUMBER OF SOIL LAYERS (IGND), THE NUMBER OF BROAD
C     * VEGETATION CATEGORIES (ICAN), OR ICAN+1.
C
      REAL(KIND = 8),DIMENSION(NLAT,NMOS,ICP1) ::
     1              FCANROT,  LNZ0ROT,
     2              ALVCROT,  ALICROT
C
      REAL(KIND = 8),DIMENSION(NLAT,NMOS,ICAN) ::
     1              PAMXROT,  PAMNROT,
     2              CMASROT,  ROOTROT,
     3              RSMNROT,  QA50ROT,
     4              VPDAROT,  VPDBROT,
     5              PSGAROT,  PSGBROT,
     6              PAIDROT,  HGTDROT,
     7              ACVDROT,  ACIDROT
C
      REAL(KIND = 8),DIMENSION(ILG,ICP1) ::
     1              FCANGAT,  LNZ0GAT,
     2              ALVCGAT,  ALICGAT
C
      REAL(KIND = 8),DIMENSION(ILG,ICAN) ::
     1              PAMXGAT,  PAMNGAT,
     2              CMASGAT,  ROOTGAT,
     3              RSMNGAT,  QA50GAT,
     4              VPDAGAT,  VPDBGAT,
     5              PSGAGAT,  PSGBGAT,
     6              PAIDGAT,  HGTDGAT,
     7              ACVDGAT,  ACIDGAT
C
      REAL(KIND = 8),DIMENSION(NLAT,NMOS,IGND) ::
     1        THPROT ,  THRROT ,  THMROT ,
     2        BIROT  ,  PSISROT,  GRKSROT,
     3        THRAROT,  HCPSROT,
     4        TCSROT ,  THFCROT,  PSIWROT,
     5        THLWROT,  DLZWROT,  ZBTWROT
C
      REAL(KIND = 8),DIMENSION(NLAT,NMOS) ::
     1        DRNROT ,   XSLPROT,   GRKFROT,
     2        WFSFROT,   WFCIROT,   ALGWROT,
     3        ALGDROT,   ASVDROT,   ASIDROT,
     4        AGVDROT,   AGIDROT,   ZSNLROT,
     5        ZPLGROT,   ZPLSROT,   ZSNOROT,
     6        ALGWVROT,  ALGWNROT,  ALGDVROT,
     7        ALGDNROT,  EMISROT

C
      REAL(KIND = 8),DIMENSION(NLAT,NMOS,NBS) ::
     1        SALBROT,   CSALROT
C
      REAL(KIND = 8),DIMENSION(NLAT,NBS) ::
     1        FSDBROL,   FSFBROL,   FSSBROL

C
      REAL(KIND = 8),DIMENSION(ILG,IGND) ::
     1        THPGAT ,  THRGAT ,  THMGAT ,
     2        BIGAT  ,  PSISGAT,  GRKSGAT,
     3        THRAGAT,  HCPSGAT,
     4        TCSGAT ,  THFCGAT,  PSIWGAT,
     5        THLWGAT,  DLZWGAT,  ZBTWGAT
C
      REAL(KIND = 8),DIMENSION(ILG) ::
     1        DRNGAT ,   XSLPGAT,   GRKFGAT,
     2        WFSFGAT,   WFCIGAT,   ALGWGAT,
     3        ALGDGAT,   ASVDGAT,   ASIDGAT,
     4        AGVDGAT,   AGIDGAT,   ZSNLGAT,
     5        ZPLGGAT,   ZPLSGAT,   ALGWVGAT,
     6        ALGWNGAT,  ALGDVGAT,  ALGDNGAT,
     7        EMISGAT
C
      REAL(KIND = 8)SANDROT(NLAT,NMOS,IGND), CLAYROT(NLAT,NMOS,IGND),
     1        ORGMROT(NLAT,NMOS,IGND), SOCIROT(NLAT,NMOS),
     2        SDEPROT(NLAT,NMOS),FAREROT(NLAT,NMOS)
C
      INTEGER MIDROT (NLAT,NMOS),ISNDROT(NLAT,NMOS,IGND),
     1        ISNDGAT( ILG,IGND),IGDRROT(NLAT,NMOS),
     2        IGDRGAT( ILG)
C
      REAL(KIND = 8),DIMENSION(ILG,NBS) ::
     1        FSDBGAT,   FSFBGAT,   FSSBGAT,
     2        SALBGAT,   CSALGAT
C
C     * ARRAYS ASSOCIATED WITH COMMON BLOCKS.
C
      REAL(KIND = 8) THPORG (3),THRORG (3),THMORG(3),BORG(3),
     1      PSISORG(3), GRKSORG(  3)
C
      REAL(KIND = 8)  CANEXT(ICAN), XLEAF (ICAN), ZORAT (ICAN),
     1      DELZ  (IGND), ZBOT  (IGND),
     2      GROWYR (  18,4,2)
C
C     * ATMOSPHERIC AND GRID-CONSTANT INPUT VARIABLES.
C
      REAL(KIND = 8),DIMENSION(NLAT) ::
     1      ZRFMROW,   ZRFHROW,   ZDMROW ,   ZDHROW ,
     2      ZBLDROW,   FSVHROW,   FSIHROW,   RADJROW,
     3      CSZROW ,   FDLROW ,   ULROW  ,   VLROW  ,
     4      TAROW  ,   QAROW  ,   PRESROW,   PREROW ,
     5      PADRROW,   VPDROW ,   TADPROW,   RHOAROW,
     6      RPCPROW,   TRPCROW,   SPCPROW,   TSPCROW,
     7      RHSIROW,   FCLOROW,   DLONROW,   UVROW  ,
     8      XDIFFUS,   GCROW  ,   Z0ORROW,   GGEOROW,
     9      RPREROW,   SPREROW,   VMODROW,   DLATROW,
     A      FSSROW,    PRENROW,   CLDTROW,   FSGROL,
     B      FLGROL,    GUSTROL,   DEPBROW
C
      REAL(KIND = 8),DIMENSION(NLAT) ::
     1	    FSSROW2,    FDLROW2,    PREROW2,    TAROW2,    
     2      QAROW2,     UVROW2,     PRESROW2
C
      REAL(KIND = 8),DIMENSION(ILG) ::
     1      ZRFMGAT,   ZRFHGAT,   ZDMGAT ,   ZDHGAT ,
     2      ZBLDGAT,   FSVHGAT,   FSIHGAT,   RADJGAT,
     3      CSZGAT ,   FDLGAT ,   ULGAT  ,   VLGAT  ,
     4      TAGAT  ,   QAGAT  ,   PRESGAT,   PREGAT ,
     5      PADRGAT,   VPDGAT ,   TADPGAT,   RHOAGAT,
     6      RPCPGAT,   TRPCGAT,   SPCPGAT,   TSPCGAT,
     7      RHSIGAT,   FCLOGAT,   DLONGAT,   Z0ORGAT,
     8      GGEOGAT,   VMODGAT,   FSGGAT,    FLGGAT,
     9      GUSTGAT,   DEPBGAT,   GTBS,      SFCUBS,
     +      SFCVBS,    USTARBS,   TCSNOW,    GSNOW

C
C     * LAND SURFACE DIAGNOSTIC VARIABLES.
C
      REAL(KIND = 8),DIMENSION(NLAT,NMOS) ::
     1      CDHROT ,   CDMROT ,   HFSROT ,   TFXROT ,
     2      QEVPROT,   QFSROT ,   QFXROT ,   PETROT ,
     3      GAROT  ,   EFROT  ,   GTROT  ,   QGROT  ,
     4      ALVSROT,   ALIRROT,   FSNOROT,   SFRHROT,
     5      SFCTROT,   SFCUROT,   SFCVROT,   SFCQROT,
     6      FSGVROT,   FSGSROT,   FSGGROT,   FLGVROT,
     7      FLGSROT,   FLGGROT,   HFSCROT,   HFSSROT,
     8      HFSGROT,   HEVCROT,   HEVSROT,   HEVGROT,
     9      HMFCROT,   HMFNROT,   HTCCROT,   HTCSROT,
     A      PCFCROT,   PCLCROT,   PCPNROT,   PCPGROT,
     B      QFGROT ,   QFNROT ,   QFCLROT,   QFCFROT,
     C      ROFROT ,   ROFOROT,   ROFSROT,   ROFBROT,
     D      TROFROT,   TROOROT,   TROSROT,   TROBROT,
     E      ROFCROT,   ROFNROT,   ROVGROT,   WTRCROT,
     F      WTRSROT,   WTRGROT,   DRROT  ,   WTABROT,
     G      ILMOROT,   UEROT  ,   HBLROT
C
      REAL(KIND = 8),DIMENSION(ILG) ::
     1      CDHGAT ,   CDMGAT ,   HFSGAT ,   TFXGAT ,
     2      QEVPGAT,   QFSGAT ,   QFXGAT ,   PETGAT ,
     3      GAGAT  ,   EFGAT  ,   GTGAT  ,   QGGAT  ,
     4      ALVSGAT,   ALIRGAT,   FSNOGAT,   SFRHGAT,
     5      SFCTGAT,   SFCUGAT,   SFCVGAT,   SFCQGAT,
     6      FSGVGAT,   FSGSGAT,   FSGGGAT,   FLGVGAT,
     7      FLGSGAT,   FLGGGAT,   HFSCGAT,   HFSSGAT,
     8      HFSGGAT,   HEVCGAT,   HEVSGAT,   HEVGGAT,
     9      HMFCGAT,   HMFNGAT,   HTCCGAT,   HTCSGAT,
     A      PCFCGAT,   PCLCGAT,   PCPNGAT,   PCPGGAT,
     B      QFGGAT ,   QFNGAT ,   QFCLGAT,   QFCFGAT,
     C      ROFGAT ,   ROFOGAT,   ROFSGAT,   ROFBGAT,
     D      TROFGAT,   TROOGAT,   TROSGAT,   TROBGAT,
     E      ROFCGAT,   ROFNGAT,   ROVGGAT,   WTRCGAT,
     F      WTRSGAT,   WTRGGAT,   DRGAT  ,   WTABGAT,
     G      ILMOGAT,   UEGAT  ,   HBLGAT ,   QLWOGAT,
     H      FTEMP  ,   FVAP   ,   RIB
C
      REAL(KIND = 8),DIMENSION(NLAT) ::
     1      CDHROW ,   CDMROW ,   HFSROW ,   TFXROW ,
     2      QEVPROW,   QFSROW ,   QFXROW ,   PETROW ,
     3      GAROW  ,   EFROW  ,   GTROW  ,   QGROW  ,
     4      ALVSROW,   ALIRROW,   FSNOROW,   SFRHROW,
     5      SFCTROW,   SFCUROW,   SFCVROW,   SFCQROW,
     6      FSGVROW,   FSGSROW,   FSGGROW,   FLGVROW,
     7      FLGSROW,   FLGGROW,   HFSCROW,   HFSSROW,
     8      HFSGROW,   HEVCROW,   HEVSROW,   HEVGROW,
     9      HMFCROW,   HMFNROW,   HTCCROW,   HTCSROW,
     A      PCFCROW,   PCLCROW,   PCPNROW,   PCPGROW,
     B      QFGROW ,   QFNROW ,   QFCLROW,   QFCFROW,
     C      ROFROW ,   ROFOROW,   ROFSROW,   ROFBROW,
     D      ROFCROW,   ROFNROW,   ROVGROW,   WTRCROW,
     E      WTRSROW,   WTRGROW,   DRROW  ,   WTABROW,
     F      ILMOROW,   UEROW  ,   HBLROW
C
      REAL(KIND = 8)HMFGROT(NLAT,NMOS,IGND),HTCROT (NLAT,NMOS,IGND),
     1        QFCROT (NLAT,NMOS,IGND), GFLXROT(NLAT,NMOS,IGND),
     2        HMFGGAT(ILG,IGND), HTCGAT (ILG,IGND),
     3        QFCGAT (ILG,IGND), GFLXGAT(ILG,IGND),
     4        HMFGROW(NLAT,IGND),HTCROW (NLAT,IGND),
     5        QFCROW (NLAT,IGND),GFLXROW(NLAT,IGND)
C
      INTEGER     ITCTROT(NLAT,NMOS,6,50),  ITCTGAT(ILG,6,50)
      INTEGER     ISUM(6)

C
C    *VARIABLES FOR MEP
C
      REAL(KIND = 8)  RN(ILG)
      REAL(KIND = 8)  H_MEP(ILG)
      REAL(KIND = 8)  LE_MEP(ILG)
      REAL(KIND = 8)  G_MEP(ILG)
C

C     * ARRAYS USED FOR OUTPUT AND DISPLAY PURPOSES.
C     * (THE SUFFIX "ACC" REFERS TO ACCUMULATOR ARRAYS USED IN
C     * CALCULATING TIME AVERAGES.)

      CHARACTER     TITLE1*4,     TITLE2*4,     TITLE3*4,
     1              TITLE4*4,     TITLE5*4,     TITLE6*4
      CHARACTER     NAME1*4,      NAME2*4,      NAME3*4,
     1              NAME4*4,      NAME5*4,      NAME6*4
      CHARACTER     PLACE1*4,     PLACE2*4,     PLACE3*4,
     1              PLACE4*4,     PLACE5*4,     PLACE6*4

      REAL(KIND = 8),DIMENSION(NLAT) ::
     1              PREACC ,   GTACC  ,   QEVPACC,
     2              HFSACC ,   ROFACC ,   SNOACC ,
     3              ALVSACC,   ALIRACC,   FSINACC,
     4              FLINACC,   TAACC  ,   UVACC  ,
     5              PRESACC,   QAACC  ,
     6              EVAPACC,   FLUTACC,   OVRACC ,
     7              HMFNACC,   WTBLACC,   WSNOACC,
     8              RHOSACC,   TSNOACC,   TCANACC,
     9              RCANACC,   SCANACC,   GROACC ,
     A              CANARE ,   SNOARE

      REAL(KIND = 8)  TBARACC(NLAT,IGND), THLQACC(NLAT,IGND),
     1              THICACC(NLAT,IGND), THALACC(NLAT,IGND)
C
!     * ARRAYS DEFINED TO PASS INFORMATION BETWEEN THE THREE MAJOR
C     * SUBSECTIONS OF CLASS ("CLASSA", "CLASST" AND "CLASSW").

      REAL(KIND = 8),DIMENSION(ILG,IGND) ::
     1        TBARC  ,     TBARG  ,     TBARCS ,
     2        TBARGS ,     THLIQC ,     THLIQG ,
     3        THICEC ,     THICEG ,     FROOT  ,
     4        HCPC   ,     HCPG   ,     FROOTS ,
     5        TCTOPC ,     TCBOTC ,
     6        TCTOPG ,     TCBOTG
C
      REAL(KIND = 8)  FC(ILG),FG(ILG),FCS(ILG),FGS(ILG),
     1      RBCOEF (ILG), ZSNOW  (ILG),
     2      FSVF   (ILG), FSVFS  (ILG),
     3      ALVSCN (ILG), ALIRCN (ILG), ALVSG  (ILG), ALIRG  (ILG),
     4      ALVSCS (ILG), ALIRCS (ILG), ALVSSN (ILG), ALIRSN (ILG),
     5      ALVSGC (ILG), ALIRGC (ILG), ALVSSC (ILG), ALIRSC (ILG),
     6      TRVSCN (ILG), TRIRCN (ILG), TRVSCS (ILG), TRIRCS (ILG),
     7      RC     (ILG), RCS    (ILG), FRAINC (ILG), FSNOWC (ILG),
     8      FRAICS (ILG), FSNOCS (ILG),
     9      CMASSC (ILG), CMASCS (ILG), DISP   (ILG), DISPS  (ILG),
     A      ZOMLNC (ILG), ZOELNC (ILG), ZOMLNG (ILG), ZOELNG (ILG),
     B      ZOMLCS (ILG), ZOELCS (ILG), ZOMLNS (ILG), ZOELNS (ILG),
     C      TRSNOWC (ILG), CHCAP  (ILG), CHCAPS (ILG),
     D      GZEROC (ILG), GZEROG (ILG), GZROCS (ILG), GZROGS (ILG),
     E      G12C   (ILG), G12G   (ILG), G12CS  (ILG), G12GS  (ILG),
     F      G23C   (ILG), G23G   (ILG), G23CS  (ILG), G23GS  (ILG),
     G      QFREZC (ILG), QFREZG (ILG), QMELTC (ILG), QMELTG (ILG),
     I      EVAPC  (ILG), EVAPCG (ILG), EVAPG  (ILG), EVAPCS (ILG),
     J      EVPCSG (ILG), EVAPGS (ILG), TCANO  (ILG), TCANS  (ILG),
     K      RAICAN (ILG), SNOCAN (ILG), RAICNS (ILG), SNOCNS (ILG),
     L      CWLCAP (ILG), CWFCAP (ILG), CWLCPS (ILG), CWFCPS (ILG),
     M      TSNOCS (ILG), TSNOGS (ILG), RHOSCS (ILG), RHOSGS (ILG),
     N      WSNOCS (ILG), WSNOGS (ILG),
     O      TPONDC (ILG), TPONDG (ILG), TPNDCS (ILG), TPNDGS (ILG),
     P      ZPLMCS (ILG), ZPLMGS (ILG), ZPLIMC (ILG), ZPLIMG (ILG)
C
      REAL(KIND = 8)  ALTG(ILG,NBS),ALSNO(ILG,NBS),TRSNOWG(ILG,NBS)

C
C     * DIAGNOSTIC ARRAYS USED FOR CHECKING ENERGY AND WATER
C     * BALANCES.
C
      REAL(KIND = 8) CTVSTP(ILG),CTSSTP(ILG),CT1STP(ILG),CT2STP(ILG),
     1     CT3STP(ILG),WTVSTP(ILG),WTSSTP(ILG),WTGSTP(ILG)
C
C     * CONSTANTS AND TEMPORARY VARIABLES.
C
      REAL(KIND = 8) DEGLON,DAY,DECL,HOUR,COSZ,CUMSNO,EVAPSUM,
     1     QSUMV,QSUMS,QSUM1,QSUM2,QSUM3,WSUMV,WSUMS,WSUMG,ALTOT,
     2     FSSTAR,FLSTAR,QH,QE,BEG,SNOMLT,ZSN,TCN,TSN,TPN,GTOUT,TAC,
     3     ALTOT_YR,TSURF,ALAVG,ALMAX,ACTLYR,FTAVG,FTMAX,FTABLE

C
C     * COMMON BLOCK PARAMETERS.
C
      REAL(KIND = 8) X1,X2,X3,X4,G,GAS,X5,X6,CPRES,GASV,X7,
     1    CPI,X8,CELZRO,X9,
     2    X10,X11,X12,X13,X14,X15,SIGMA,X16,DELTIM,DELT,TFREZ,
     3    RGAS,RGASV,GRAV,SBC,VKC,CT,VMIN,TCW,TCICE,TCSAND,TCCLAY,
     4    TCOM,TCDRYS,RHOSOL,RHOOM,HCPW,HCPICE,HCPSOL,HCPOM,HCPSND,
     5    HCPCLY,SPHW,SPHICE,SPHVEG,SPHAIR,RHOW,RHOICE,TCGLAC,CLHMLT,
     6    CLHVAP,PI,ZOLNG,ZOLNS,ZOLNI,ZORATG,ALVSI,ALIRI,ALVSO,ALIRO,
     7    ALBRCK,DELTA,CGRAV,CKARM,CPD,AS,ASX,CI,BS,BETA,FACTN,HMIN,
     8    ANGMAX,A,B
C
c================= CTEM array declaration ===============================\
c
c     Local variables for coupling CLASS and CTEM
c
      integer ictemmod

      integer strlen !strlen can go in arg reader subroutine.
      character*80   titlec1 !, titlec2, titlec3
      character*80   argbuff
      character*160  command
c
       integer   lopcount,  isumc,
     1           k1c,       k2c,     jhhstd,
     2           jhhendd,   jdstd,   jdendd,      jhhsty,
     3           jhhendy,   jdsty,   jdendy,
     4           month1,
     5           month2,      xday,  ctemloop,nummetcylyrs,
     6           ncyear,  co2yr,
     5           spinfast,   nol2pfts(4),
     6           popyr,
     7           metcylyrst, metcycendyr, climiyear, popcycleyr,
     8           cypopyr, lucyr, cylucyr, endyr,bigpftc(2),
     9           obswetyr, cywetldyr, trans_startyr, jmosty,
     +           obslghtyr
c
      real(KIND = 8), pointer ::  fsstar_g
      real(KIND = 8), pointer ::  flstar_g
      real(KIND = 8), pointer ::  qh_g
      real(KIND = 8), pointer ::  qe_g
      real(KIND = 8), pointer ::  snomlt_g
      real(KIND = 8), pointer ::  beg_g
      real(KIND = 8), pointer ::  gtout_g
      real(KIND = 8), pointer ::  tpn_g
      real(KIND = 8), pointer ::  altot_g
      real(KIND = 8), pointer ::  tcn_g
      real(KIND = 8), pointer ::  tsn_g
      real(KIND = 8), pointer ::  zsn_g

       real(KIND = 8) co2concin,popdin(nlat),setco2conc, sumfare,
     1           temp_var, barefrac,  todfrac(ilg,icc), barf(nlat)

      real(KIND = 8) grclarea(ilg), crop_temp_frac(ilg,2)
c
c     Competition related variables

       real(KIND = 8) fsinacc_gat(ilg), flutacc_gat(ilg),
     1      flinacc_gat(ilg),
     2      alswacc_gat(ilg), allwacc_gat(ilg), pregacc_gat(ilg),
     3      altot_gat,        fsstar_gat,       flstar_gat,
     $      netrad_gat(ilg),  preacc_gat(ilg)
c
!       These go into CTEM but that is about it...
       real(KIND = 8) tcurm(ilg),srpcuryr(ilg),dftcuryr(ilg),
     1      tmonth(12,ilg), anpcpcur(ilg),anpecur(ilg),
     2      gdd5cur(ilg),   surmncur(ilg),defmncur(ilg),
     3      srplscur(ilg),  defctcur(ilg)

c
!       These go into CTEM but that is about it...
       real(KIND = 8) lyglfmasgat(ilg,icc),   geremortgat(ilg,icc),
     1      intrmortgat(ilg,icc),     lambdagat(ilg,icc),
         !2            rnded_pft(icc),
     3            ccgat(ilg,icc),         mmgat(ilg,icc) !,
!     4            temparray(icc),                  temp

      real(KIND = 8)  xdiffusgat(ilg) ! the corresponding ROW is CLASS's XDIFFUS

!     For these below, the corresponding ROWs are defined by CLASS

      real(KIND = 8)  sdepgat(ilg),       orgmgat(ilg,ignd),
     1      sandgat(ilg,ignd),  claygat(ilg,ignd)

      ! CLASS Monthly Outputs:

      real(KIND = 8), pointer, dimension(:) :: ALVSACC_MO
      real(KIND = 8), pointer, dimension(:) :: ALIRACC_MO
      real(KIND = 8), pointer, dimension(:) :: FLUTACC_MO
      real(KIND = 8), pointer, dimension(:) :: FSINACC_MO
      real(KIND = 8), pointer, dimension(:) :: FLINACC_MO
      real(KIND = 8), pointer, dimension(:) :: HFSACC_MO
      real(KIND = 8), pointer, dimension(:) :: QEVPACC_MO
      real(KIND = 8), pointer, dimension(:) :: SNOACC_MO
      real(KIND = 8), pointer, dimension(:) :: WSNOACC_MO
      real(KIND = 8), pointer, dimension(:) :: ROFACC_MO
      real(KIND = 8), pointer, dimension(:) :: PREACC_MO
      real(KIND = 8), pointer, dimension(:) :: EVAPACC_MO
      real(KIND = 8), pointer, dimension(:) :: TAACC_MO

      real(KIND = 8), pointer :: FSSTAR_MO
      real(KIND = 8), pointer :: FLSTAR_MO
      real(KIND = 8), pointer :: QH_MO
      real(KIND = 8), pointer :: QE_MO

      real(KIND = 8), pointer, dimension(:,:) :: TBARACC_MO
      real(KIND = 8), pointer, dimension(:,:) :: THLQACC_MO
      real(KIND = 8), pointer, dimension(:,:) :: THICACC_MO

    ! CLASS yearly output for class grid-mean

      real(KIND = 8), pointer, dimension(:) :: ALVSACC_YR
      real(KIND = 8), pointer, dimension(:) :: ALIRACC_YR
      real(KIND = 8), pointer, dimension(:) :: FLUTACC_YR
      real(KIND = 8), pointer, dimension(:) :: FSINACC_YR
      real(KIND = 8), pointer, dimension(:) :: FLINACC_YR
      real(KIND = 8), pointer, dimension(:) :: HFSACC_YR
      real(KIND = 8), pointer, dimension(:) :: QEVPACC_YR
      real(KIND = 8), pointer, dimension(:) :: ROFACC_YR
      real(KIND = 8), pointer, dimension(:) :: PREACC_YR
      real(KIND = 8), pointer, dimension(:) :: EVAPACC_YR
      real(KIND = 8), pointer, dimension(:) :: TAACC_YR

      real(KIND = 8), pointer :: FSSTAR_YR
      real(KIND = 8), pointer :: FLSTAR_YR
      real(KIND = 8), pointer :: QH_YR
      real(KIND = 8), pointer :: QE_YR

      logical, pointer :: ctem_on
      logical, pointer :: parallelrun
      logical, pointer :: mosaic
      logical, pointer :: cyclemet
      logical, pointer :: dofire
      logical, pointer :: run_model
      logical, pointer :: met_rewound
      logical, pointer :: reach_eof
      logical, pointer :: compete
      logical, pointer :: start_bare
      logical, pointer :: rsfile
      logical, pointer :: lnduseon
      logical, pointer :: co2on
      logical, pointer :: popdon
      logical, pointer :: inibioclim
      logical, pointer :: start_from_rs
      logical, pointer :: dowetlands
      logical, pointer :: obswetf
      logical, pointer :: transient_run

      ! ROW vars:
      logical, pointer, dimension(:,:,:) :: pftexistrow
      integer, pointer, dimension(:,:,:) :: colddaysrow
      integer, pointer, dimension(:,:) :: icountrow
      integer, pointer, dimension(:,:,:) :: lfstatusrow
      integer, pointer, dimension(:,:,:) :: pandaysrow

      integer, pointer, dimension(:) :: stdalngrd

      real(KIND = 8), pointer, dimension(:,:) :: tcanrs
      real(KIND = 8), pointer, dimension(:,:) :: tsnors
      real(KIND = 8), pointer, dimension(:,:) :: tpndrs
      real(KIND = 8), pointer, dimension(:,:,:) :: csum
      real(KIND = 8), pointer, dimension(:,:,:) :: tbaraccrow_m
      real(KIND = 8), pointer, dimension(:,:) :: tcanoaccrow_m
      real(KIND = 8), pointer, dimension(:,:) :: uvaccrow_m
      real(KIND = 8), pointer, dimension(:,:) :: vvaccrow_m

      real(KIND = 8), pointer, dimension(:,:,:) :: ailcminrow         !
      real(KIND = 8), pointer, dimension(:,:,:) :: ailcmaxrow         !
      real(KIND = 8), pointer, dimension(:,:,:) :: dvdfcanrow         !
      real(KIND = 8), pointer, dimension(:,:,:) :: gleafmasrow        !
      real(KIND = 8), pointer, dimension(:,:,:) :: bleafmasrow        !
      real(KIND = 8), pointer, dimension(:,:,:) :: stemmassrow        !
      real(KIND = 8), pointer, dimension(:,:,:) :: rootmassrow        !
      real(KIND = 8), pointer, dimension(:,:,:) :: pstemmassrow       !
      real(KIND = 8), pointer, dimension(:,:,:) :: pgleafmassrow      !
      real(KIND = 8), pointer, dimension(:,:,:) :: fcancmxrow
      real(KIND = 8), pointer, dimension(:,:) :: gavglairow
      real(KIND = 8), pointer, dimension(:,:,:) :: zolncrow
      real(KIND = 8), pointer, dimension(:,:,:) :: ailcrow
      real(KIND = 8), pointer, dimension(:,:,:) :: ailcgrow
      real(KIND = 8), pointer, dimension(:,:,:) :: ailcgsrow
      real(KIND = 8), pointer, dimension(:,:,:) :: fcancsrow
      real(KIND = 8), pointer, dimension(:,:,:) :: fcancrow
      real(KIND = 8), pointer, dimension(:,:) :: co2concrow
      real(KIND = 8), pointer, dimension(:,:,:) :: co2i1cgrow
      real(KIND = 8), pointer, dimension(:,:,:) :: co2i1csrow
      real(KIND = 8), pointer, dimension(:,:,:) :: co2i2cgrow
      real(KIND = 8), pointer, dimension(:,:,:) :: co2i2csrow
      real(KIND = 8), pointer, dimension(:,:,:) :: ancsvegrow
      real(KIND = 8), pointer, dimension(:,:,:) :: ancgvegrow
      real(KIND = 8), pointer, dimension(:,:,:) :: rmlcsvegrow
      real(KIND = 8), pointer, dimension(:,:,:) :: rmlcgvegrow
      real(KIND = 8), pointer, dimension(:,:,:) :: slairow
      real(KIND = 8), pointer, dimension(:,:,:) :: ailcbrow
      real(KIND = 8), pointer, dimension(:,:) :: canresrow
      real(KIND = 8), pointer, dimension(:,:,:) :: flhrlossrow

      real(KIND = 8), pointer, dimension(:,:,:) :: grwtheffrow
      real(KIND = 8), pointer, dimension(:,:,:) :: lystmmasrow
      real(KIND = 8), pointer, dimension(:,:,:) :: lyrotmasrow
      real(KIND = 8), pointer, dimension(:,:,:) :: tymaxlairow
      real(KIND = 8), pointer, dimension(:,:) :: vgbiomasrow
      real(KIND = 8), pointer, dimension(:,:) :: gavgltmsrow
      real(KIND = 8), pointer, dimension(:,:) :: gavgscmsrow
      real(KIND = 8), pointer, dimension(:,:,:) :: stmhrlosrow
      real(KIND = 8), pointer, dimension(:,:,:,:) :: rmatcrow
      real(KIND = 8), pointer, dimension(:,:,:,:) :: rmatctemrow
      real(KIND = 8), pointer, dimension(:,:,:) :: litrmassrow
      real(KIND = 8), pointer, dimension(:,:,:) :: soilcmasrow
      real(KIND = 8), pointer, dimension(:,:,:) :: vgbiomas_vegrow

      real(KIND = 8), pointer, dimension(:,:,:) :: emit_co2row
      real(KIND = 8), pointer, dimension(:,:,:) :: emit_corow
      real(KIND = 8), pointer, dimension(:,:,:) :: emit_ch4row
      real(KIND = 8), pointer, dimension(:,:,:) :: emit_nmhcrow
      real(KIND = 8), pointer, dimension(:,:,:) :: emit_h2row
      real(KIND = 8), pointer, dimension(:,:,:) :: emit_noxrow
      real(KIND = 8), pointer, dimension(:,:,:) :: emit_n2orow
      real(KIND = 8), pointer, dimension(:,:,:) :: emit_pm25row
      real(KIND = 8), pointer, dimension(:,:,:) :: emit_tpmrow
      real(KIND = 8), pointer, dimension(:,:,:) :: emit_tcrow
      real(KIND = 8), pointer, dimension(:,:,:) :: emit_ocrow
      real(KIND = 8), pointer, dimension(:,:,:) :: emit_bcrow
      real(KIND = 8), pointer, dimension(:,:) :: burnfracrow
      real(KIND = 8), pointer, dimension(:,:,:) :: burnvegfrow
      real(KIND = 8), pointer, dimension(:,:) :: probfirerow
      real(KIND = 8), pointer, dimension(:,:) :: btermrow
      real(KIND = 8), pointer, dimension(:,:) :: ltermrow
      real(KIND = 8), pointer, dimension(:,:) :: mtermrow

      real(KIND = 8), pointer, dimension(:) :: extnprobgrd
      real(KIND = 8), pointer, dimension(:) :: prbfrhucgrd
      real(KIND = 8), pointer, dimension(:,:) :: mlightnggrd

      real(KIND = 8), pointer, dimension(:,:,:) :: bmasvegrow
      real(KIND = 8), pointer, dimension(:,:,:) :: cmasvegcrow
      real(KIND = 8), pointer, dimension(:,:,:) :: veghghtrow
      real(KIND = 8), pointer, dimension(:,:,:) :: rootdpthrow
      real(KIND = 8), pointer, dimension(:,:) :: rmlrow
      real(KIND = 8), pointer, dimension(:,:) :: rmsrow
      real(KIND = 8), pointer, dimension(:,:,:) :: tltrleafrow
      real(KIND = 8), pointer, dimension(:,:,:) :: tltrstemrow
      real(KIND = 8), pointer, dimension(:,:,:) :: tltrrootrow
      real(KIND = 8), pointer, dimension(:,:,:) :: leaflitrrow
      real(KIND = 8), pointer, dimension(:,:,:) :: roottemprow
      real(KIND = 8), pointer, dimension(:,:,:) :: afrleafrow
      real(KIND = 8), pointer, dimension(:,:,:) :: afrstemrow
      real(KIND = 8), pointer, dimension(:,:,:) :: afrrootrow
      real(KIND = 8), pointer, dimension(:,:,:) :: wtstatusrow
      real(KIND = 8), pointer, dimension(:,:,:) :: ltstatusrow
      real(KIND = 8), pointer, dimension(:,:) :: rmrrow

      REAL(KIND = 8), pointer, dimension(:) :: wetfracgrd
      REAL(KIND = 8), pointer, dimension(:,:) :: ch4wet1row
      REAL(KIND = 8), pointer, dimension(:,:) :: ch4wet2row
      REAL(KIND = 8), pointer, dimension(:,:) :: wetfdynrow
      REAL(KIND = 8), pointer, dimension(:,:) :: ch4dyn1row
      REAL(KIND = 8), pointer, dimension(:,:) :: ch4dyn2row
      REAL(KIND = 8), pointer, dimension(:,:) :: wetfrac_mon

      REAL(KIND = 8), pointer, dimension(:,:) :: lucemcomrow
      REAL(KIND = 8), pointer, dimension(:,:) :: lucltrinrow
      REAL(KIND = 8), pointer, dimension(:,:) :: lucsocinrow

      REAL(KIND = 8), pointer, dimension(:,:) :: npprow
      REAL(KIND = 8), pointer, dimension(:,:) :: neprow
      REAL(KIND = 8), pointer, dimension(:,:) :: nbprow
      REAL(KIND = 8), pointer, dimension(:,:) :: gpprow
      REAL(KIND = 8), pointer, dimension(:,:) :: hetroresrow
      REAL(KIND = 8), pointer, dimension(:,:) :: autoresrow
      REAL(KIND = 8), pointer, dimension(:,:) :: soilcresprow
      REAL(KIND = 8), pointer, dimension(:,:) :: rmrow
      REAL(KIND = 8), pointer, dimension(:,:) :: rgrow
      REAL(KIND = 8), pointer, dimension(:,:) :: litresrow
      REAL(KIND = 8), pointer, dimension(:,:) :: socresrow
      REAL(KIND = 8), pointer, dimension(:,:) :: dstcemlsrow
      REAL(KIND = 8), pointer, dimension(:,:) :: litrfallrow
      REAL(KIND = 8), pointer, dimension(:,:) :: humiftrsrow

      REAL(KIND = 8), pointer, dimension(:,:,:) :: gppvegrow
      REAL(KIND = 8), pointer, dimension(:,:,:) :: nepvegrow
      REAL(KIND = 8), pointer, dimension(:,:,:) :: nbpvegrow
      REAL(KIND = 8), pointer, dimension(:,:,:) :: nppvegrow
      REAL(KIND = 8), pointer, dimension(:,:,:) :: hetroresvegrow
      REAL(KIND = 8), pointer, dimension(:,:,:) :: autoresvegrow
      REAL(KIND = 8), pointer, dimension(:,:,:) :: litresvegrow
      REAL(KIND = 8), pointer, dimension(:,:,:) :: soilcresvegrow
      REAL(KIND = 8), pointer, dimension(:,:,:) :: rmlvegaccrow
      REAL(KIND = 8), pointer, dimension(:,:,:) :: rmsvegrow
      REAL(KIND = 8), pointer, dimension(:,:,:) :: rmrvegrow
      REAL(KIND = 8), pointer, dimension(:,:,:) :: rgvegrow

      REAL(KIND = 8), pointer, dimension(:,:,:) :: rothrlosrow
      REAL(KIND = 8), pointer, dimension(:,:,:) :: pfcancmxrow
      REAL(KIND = 8), pointer, dimension(:,:,:) :: nfcancmxrow
      REAL(KIND = 8), pointer, dimension(:,:,:) :: alvsctmrow
      REAL(KIND = 8), pointer, dimension(:,:,:) :: paicrow
      REAL(KIND = 8), pointer, dimension(:,:,:) :: slaicrow
      REAL(KIND = 8), pointer, dimension(:,:,:) :: alirctmrow
      REAL(KIND = 8), pointer, dimension(:,:) :: cfluxcgrow
      REAL(KIND = 8), pointer, dimension(:,:) :: cfluxcsrow
      REAL(KIND = 8), pointer, dimension(:,:) :: dstcemls3row
      REAL(KIND = 8), pointer, dimension(:,:,:) :: anvegrow
      REAL(KIND = 8), pointer, dimension(:,:,:) :: rmlvegrow

      ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>.
      ! GAT version:

      logical, pointer, dimension(:,:) :: pftexistgat
      integer, pointer, dimension(:,:) :: colddaysgat
      integer, pointer, dimension(:) :: icountgat
      integer, pointer, dimension(:,:) :: lfstatusgat
      integer, pointer, dimension(:,:) :: pandaysgat

      integer, pointer, dimension(:) :: stdalngat

      REAL(KIND = 8), pointer, dimension(:) :: lightng

      REAL(KIND = 8), pointer, dimension(:,:) :: ailcmingat         !
      REAL(KIND = 8), pointer, dimension(:,:) :: ailcmaxgat         !
      REAL(KIND = 8), pointer, dimension(:,:) :: dvdfcangat         !
      REAL(KIND = 8), pointer, dimension(:,:) :: gleafmasgat        !
      REAL(KIND = 8), pointer, dimension(:,:) :: bleafmasgat        !
      REAL(KIND = 8), pointer, dimension(:,:) :: stemmassgat        !
      REAL(KIND = 8), pointer, dimension(:,:) :: rootmassgat        !
      REAL(KIND = 8), pointer, dimension(:,:) :: pstemmassgat       !
      REAL(KIND = 8), pointer, dimension(:,:) :: pgleafmassgat      !
      REAL(KIND = 8), pointer, dimension(:,:) :: fcancmxgat
      REAL(KIND = 8), pointer, dimension(:) :: gavglaigat
      REAL(KIND = 8), pointer, dimension(:,:) :: zolncgat
      REAL(KIND = 8), pointer, dimension(:,:) :: ailcgat
      REAL(KIND = 8), pointer, dimension(:,:) :: ailcggat
      REAL(KIND = 8), pointer, dimension(:,:) :: ailcgsgat
      REAL(KIND = 8), pointer, dimension(:,:) :: fcancsgat
      REAL(KIND = 8), pointer, dimension(:,:) :: fcancgat
      REAL(KIND = 8), pointer, dimension(:) :: co2concgat
      REAL(KIND = 8), pointer, dimension(:,:) :: co2i1cggat
      REAL(KIND = 8), pointer, dimension(:,:) :: co2i1csgat
      REAL(KIND = 8), pointer, dimension(:,:) :: co2i2cggat
      REAL(KIND = 8), pointer, dimension(:,:) :: co2i2csgat
      REAL(KIND = 8), pointer, dimension(:,:) :: ancsveggat
      REAL(KIND = 8), pointer, dimension(:,:) :: ancgveggat
      REAL(KIND = 8), pointer, dimension(:,:) :: rmlcsveggat
      REAL(KIND = 8), pointer, dimension(:,:) :: rmlcgveggat
      REAL(KIND = 8), pointer, dimension(:,:) :: slaigat
      REAL(KIND = 8), pointer, dimension(:,:) :: ailcbgat
      REAL(KIND = 8), pointer, dimension(:) :: canresgat
      REAL(KIND = 8), pointer, dimension(:,:) :: flhrlossgat

      REAL(KIND = 8), pointer, dimension(:,:) :: grwtheffgat
      REAL(KIND = 8), pointer, dimension(:,:) :: lystmmasgat
      REAL(KIND = 8), pointer, dimension(:,:) :: lyrotmasgat
      REAL(KIND = 8), pointer, dimension(:,:) :: tymaxlaigat
      REAL(KIND = 8), pointer, dimension(:) :: vgbiomasgat
      REAL(KIND = 8), pointer, dimension(:) :: gavgltmsgat
      REAL(KIND = 8), pointer, dimension(:) :: gavgscmsgat
      REAL(KIND = 8), pointer, dimension(:,:) :: stmhrlosgat
      REAL(KIND = 8), pointer, dimension(:,:,:) :: rmatcgat
      REAL(KIND = 8), pointer, dimension(:,:,:) :: rmatctemgat
      REAL(KIND = 8), pointer, dimension(:,:) :: litrmassgat
      REAL(KIND = 8), pointer, dimension(:,:) :: soilcmasgat
      REAL(KIND = 8), pointer, dimension(:,:) :: vgbiomas_veggat

      REAL(KIND = 8), pointer, dimension(:,:) :: emit_co2gat
      REAL(KIND = 8), pointer, dimension(:,:) :: emit_cogat
      REAL(KIND = 8), pointer, dimension(:,:) :: emit_ch4gat
      REAL(KIND = 8), pointer, dimension(:,:) :: emit_nmhcgat
      REAL(KIND = 8), pointer, dimension(:,:) :: emit_h2gat
      REAL(KIND = 8), pointer, dimension(:,:) :: emit_noxgat
      REAL(KIND = 8), pointer, dimension(:,:) :: emit_n2ogat
      REAL(KIND = 8), pointer, dimension(:,:) :: emit_pm25gat
      REAL(KIND = 8), pointer, dimension(:,:) :: emit_tpmgat
      REAL(KIND = 8), pointer, dimension(:,:) :: emit_tcgat
      REAL(KIND = 8), pointer, dimension(:,:) :: emit_ocgat
      REAL(KIND = 8), pointer, dimension(:,:) :: emit_bcgat
      REAL(KIND = 8), pointer, dimension(:) :: burnfracgat
      REAL(KIND = 8), pointer, dimension(:,:) :: burnvegfgat
      REAL(KIND = 8), pointer, dimension(:) :: probfiregat
      REAL(KIND = 8), pointer, dimension(:) :: btermgat
      REAL(KIND = 8), pointer, dimension(:) :: ltermgat
      REAL(KIND = 8), pointer, dimension(:) :: mtermgat

      REAL(KIND = 8), pointer, dimension(:) :: extnprobgat
      REAL(KIND = 8), pointer, dimension(:) :: prbfrhucgat
      REAL(KIND = 8), pointer, dimension(:,:) :: mlightnggat

      REAL(KIND = 8), pointer, dimension(:,:) :: bmasveggat
      REAL(KIND = 8), pointer, dimension(:,:) :: cmasvegcgat
      REAL(KIND = 8), pointer, dimension(:,:) :: veghghtgat
      REAL(KIND = 8), pointer, dimension(:,:) :: rootdpthgat
      REAL(KIND = 8), pointer, dimension(:) :: rmlgat
      REAL(KIND = 8), pointer, dimension(:) :: rmsgat
      REAL(KIND = 8), pointer, dimension(:,:) :: tltrleafgat
      REAL(KIND = 8), pointer, dimension(:,:) :: tltrstemgat
      REAL(KIND = 8), pointer, dimension(:,:) :: tltrrootgat
      REAL(KIND = 8), pointer, dimension(:,:) :: leaflitrgat
      REAL(KIND = 8), pointer, dimension(:,:) :: roottempgat
      REAL(KIND = 8), pointer, dimension(:,:) :: afrleafgat
      REAL(KIND = 8), pointer, dimension(:,:) :: afrstemgat
      REAL(KIND = 8), pointer, dimension(:,:) :: afrrootgat
      REAL(KIND = 8), pointer, dimension(:,:) :: wtstatusgat
      REAL(KIND = 8), pointer, dimension(:,:) :: ltstatusgat
      REAL(KIND = 8), pointer, dimension(:) :: rmrgat

      REAL(KIND = 8), pointer, dimension(:,:) :: wetfrac_sgrd
      REAL(KIND = 8), pointer, dimension(:) :: ch4wet1gat
      REAL(KIND = 8), pointer, dimension(:) :: ch4wet2gat
      REAL(KIND = 8), pointer, dimension(:) :: wetfdyngat
      REAL(KIND = 8), pointer, dimension(:) :: ch4dyn1gat
      REAL(KIND = 8), pointer, dimension(:) :: ch4dyn2gat

      REAL(KIND = 8), pointer, dimension(:) :: lucemcomgat
      REAL(KIND = 8), pointer, dimension(:) :: lucltringat
      REAL(KIND = 8), pointer, dimension(:) :: lucsocingat

      REAL(KIND = 8), pointer, dimension(:) :: nppgat
      REAL(KIND = 8), pointer, dimension(:) :: nepgat
      REAL(KIND = 8), pointer, dimension(:) :: nbpgat
      REAL(KIND = 8), pointer, dimension(:) :: gppgat
      REAL(KIND = 8), pointer, dimension(:) :: hetroresgat
      REAL(KIND = 8), pointer, dimension(:) :: autoresgat
      REAL(KIND = 8), pointer, dimension(:) :: soilcrespgat
      REAL(KIND = 8), pointer, dimension(:) :: rmgat
      REAL(KIND = 8), pointer, dimension(:) :: rggat
      REAL(KIND = 8), pointer, dimension(:) :: litresgat
      REAL(KIND = 8), pointer, dimension(:) :: socresgat
      REAL(KIND = 8), pointer, dimension(:) :: dstcemlsgat
      REAL(KIND = 8), pointer, dimension(:) :: litrfallgat
      REAL(KIND = 8), pointer, dimension(:) :: humiftrsgat

      REAL(KIND = 8), pointer, dimension(:,:) :: gppveggat
      REAL(KIND = 8), pointer, dimension(:,:) :: nepveggat
      REAL(KIND = 8), pointer, dimension(:,:) :: nbpveggat
      REAL(KIND = 8), pointer, dimension(:,:) :: nppveggat
      REAL(KIND = 8), pointer, dimension(:,:) :: hetroresveggat
      REAL(KIND = 8), pointer, dimension(:,:) :: autoresveggat
      REAL(KIND = 8), pointer, dimension(:,:) :: litresveggat
      REAL(KIND = 8), pointer, dimension(:,:) :: soilcresveggat
      REAL(KIND = 8), pointer, dimension(:,:) :: rmlvegaccgat
      REAL(KIND = 8), pointer, dimension(:,:) :: rmsveggat
      REAL(KIND = 8), pointer, dimension(:,:) :: rmrveggat
      REAL(KIND = 8), pointer, dimension(:,:) :: rgveggat

      REAL(KIND = 8), pointer, dimension(:,:) :: rothrlosgat
      REAL(KIND = 8), pointer, dimension(:,:) :: pfcancmxgat
      REAL(KIND = 8), pointer, dimension(:,:) :: nfcancmxgat
      REAL(KIND = 8), pointer, dimension(:,:) :: alvsctmgat
      REAL(KIND = 8), pointer, dimension(:,:) :: paicgat
      REAL(KIND = 8), pointer, dimension(:,:) :: slaicgat
      REAL(KIND = 8), pointer, dimension(:,:) :: alirctmgat
      REAL(KIND = 8), pointer, dimension(:) :: cfluxcggat
      REAL(KIND = 8), pointer, dimension(:) :: cfluxcsgat
      REAL(KIND = 8), pointer, dimension(:) :: dstcemls3gat
      REAL(KIND = 8), pointer, dimension(:,:) :: anveggat
      REAL(KIND = 8), pointer, dimension(:,:) :: rmlveggat

      REAL(KIND = 8), pointer, dimension(:) :: twarmm    ! temperature of the warmest month (c)
      REAL(KIND = 8), pointer, dimension(:) :: tcoldm    ! temperature of the coldest month (c)
      REAL(KIND = 8), pointer, dimension(:) :: gdd5      ! growing degree days above 5 c
      REAL(KIND = 8), pointer, dimension(:) :: aridity   ! aridity index, ratio of potential evaporation to precipitation
      REAL(KIND = 8), pointer, dimension(:) :: srplsmon  ! number of months in a year with surplus water i.e.
                                                  !  precipitation more than potential evaporation
      REAL(KIND = 8), pointer, dimension(:) :: defctmon  ! number of months in a year with water deficit i.e.
                                                  ! precipitation less than potential evaporation
      REAL(KIND = 8), pointer, dimension(:) :: anndefct  ! annual water deficit (mm)
      REAL(KIND = 8), pointer, dimension(:) :: annsrpls  ! annual water surplus (mm)
      REAL(KIND = 8), pointer, dimension(:) :: annpcp    ! annual precipitation (mm)
      REAL(KIND = 8), pointer, dimension(:) :: dry_season_length  ! length of dry season (months)

      ! Mosaic level:

      REAL(KIND = 8), pointer, dimension(:,:) :: PREACC_M
      REAL(KIND = 8), pointer, dimension(:,:) :: GTACC_M
      REAL(KIND = 8), pointer, dimension(:,:) :: QEVPACC_M
      REAL(KIND = 8), pointer, dimension(:,:) :: HFSACC_M
      REAL(KIND = 8), pointer, dimension(:,:) :: HMFNACC_M
      REAL(KIND = 8), pointer, dimension(:,:) :: ROFACC_M
      REAL(KIND = 8), pointer, dimension(:,:) :: SNOACC_M
      REAL(KIND = 8), pointer, dimension(:,:) :: OVRACC_M
      REAL(KIND = 8), pointer, dimension(:,:) :: WTBLACC_M
      REAL(KIND = 8), pointer, dimension(:,:,:) :: TBARACC_M
      REAL(KIND = 8), pointer, dimension(:,:,:) :: THLQACC_M
      REAL(KIND = 8), pointer, dimension(:,:,:) :: THICACC_M
      REAL(KIND = 8), pointer, dimension(:,:,:) :: THALACC_M
      REAL(KIND = 8), pointer, dimension(:,:) :: ALVSACC_M
      REAL(KIND = 8), pointer, dimension(:,:) :: ALIRACC_M
      REAL(KIND = 8), pointer, dimension(:,:) :: RHOSACC_M
      REAL(KIND = 8), pointer, dimension(:,:) :: TSNOACC_M
      REAL(KIND = 8), pointer, dimension(:,:) :: WSNOACC_M
      REAL(KIND = 8), pointer, dimension(:,:) :: SNOARE_M
      REAL(KIND = 8), pointer, dimension(:,:) :: TCANACC_M
      REAL(KIND = 8), pointer, dimension(:,:) :: RCANACC_M
      REAL(KIND = 8), pointer, dimension(:,:) :: SCANACC_M
      REAL(KIND = 8), pointer, dimension(:,:) :: GROACC_M
      REAL(KIND = 8), pointer, dimension(:,:) :: FSINACC_M
      REAL(KIND = 8), pointer, dimension(:,:) :: FLINACC_M
      REAL(KIND = 8), pointer, dimension(:,:) :: TAACC_M
      REAL(KIND = 8), pointer, dimension(:,:) :: UVACC_M
      REAL(KIND = 8), pointer, dimension(:,:) :: PRESACC_M
      REAL(KIND = 8), pointer, dimension(:,:) :: QAACC_M
      REAL(KIND = 8), pointer, dimension(:,:) :: EVAPACC_M
      REAL(KIND = 8), pointer, dimension(:,:) :: FLUTACC_M

!      Outputs

       REAL(KIND = 8), pointer, dimension(:,:) :: tcanoaccrow_out
       REAL(KIND = 8), pointer, dimension(:) :: tcanoaccgat_out
       REAL(KIND = 8), pointer, dimension(:,:) :: qevpacc_m_save

!     -----------------------
!      Mosaic-level variables (denoted by an ending of "_m")

      REAL(KIND = 8) faregat(ilg)
c
      REAL(KIND = 8), pointer, dimension(:,:) :: leaflitr_m
      REAL(KIND = 8), pointer, dimension(:,:) :: tltrleaf_m
      REAL(KIND = 8), pointer, dimension(:,:) :: tltrstem_m
      REAL(KIND = 8), pointer, dimension(:,:) :: tltrroot_m
      REAL(KIND = 8), pointer, dimension(:,:) :: ailcg_m
      REAL(KIND = 8), pointer, dimension(:,:) :: ailcb_m
      REAL(KIND = 8), pointer, dimension(:,:,:) :: rmatctem_m
      REAL(KIND = 8), pointer, dimension(:,:) :: veghght_m
      REAL(KIND = 8), pointer, dimension(:,:) :: rootdpth_m
      REAL(KIND = 8), pointer, dimension(:,:) :: roottemp_m
      REAL(KIND = 8), pointer, dimension(:,:) :: slai_m
      REAL(KIND = 8), pointer, dimension(:,:) :: afrroot_m
      REAL(KIND = 8), pointer, dimension(:,:) :: afrleaf_m
      REAL(KIND = 8), pointer, dimension(:,:) :: afrstem_m
      REAL(KIND = 8), pointer, dimension(:,:) :: laimaxg_m
      REAL(KIND = 8), pointer, dimension(:,:) :: stemmass_m
      REAL(KIND = 8), pointer, dimension(:,:) :: rootmass_m
      REAL(KIND = 8), pointer, dimension(:,:) :: litrmass_m
      REAL(KIND = 8), pointer, dimension(:,:) :: gleafmas_m
      REAL(KIND = 8), pointer, dimension(:,:) :: bleafmas_m
      REAL(KIND = 8), pointer, dimension(:,:) :: soilcmas_m

      REAL(KIND = 8), pointer, dimension(:) :: fsnowacc_m
      REAL(KIND = 8), pointer, dimension(:) :: tcansacc_m
      REAL(KIND = 8), pointer, dimension(:) :: tcanoaccgat_m
      REAL(KIND = 8), pointer, dimension(:) :: taaccgat_m
      REAL(KIND = 8), pointer, dimension(:) :: uvaccgat_m
      REAL(KIND = 8), pointer, dimension(:) :: vvaccgat_m
      REAL(KIND = 8), pointer, dimension(:,:) :: tbaraccgat_m
      REAL(KIND = 8), pointer, dimension(:,:) :: tbarcacc_m
      REAL(KIND = 8), pointer, dimension(:,:) :: tbarcsacc_m
      REAL(KIND = 8), pointer, dimension(:,:) :: tbargacc_m
      REAL(KIND = 8), pointer, dimension(:,:) :: tbargsacc_m
      REAL(KIND = 8), pointer, dimension(:,:) :: thliqcacc_m
      REAL(KIND = 8), pointer, dimension(:,:) :: thliqgacc_m
      REAL(KIND = 8), pointer, dimension(:,:) :: thicecacc_m
      REAL(KIND = 8), pointer, dimension(:,:) :: ancsvgac_m
      REAL(KIND = 8), pointer, dimension(:,:) :: ancgvgac_m
      REAL(KIND = 8), pointer, dimension(:,:) :: rmlcsvga_m
      REAL(KIND = 8), pointer, dimension(:,:) :: rmlcgvga_m
      integer, pointer, dimension(:,:) :: ifcancmx_m

!     -----------------------
!     Grid-averaged variables (denoted with an ending of "_g")

      REAL(KIND = 8), pointer, dimension(:) :: WSNOROT_g
      REAL(KIND = 8), pointer, dimension(:) :: ROFSROT_g
      REAL(KIND = 8), pointer, dimension(:) :: SNOROT_g
      REAL(KIND = 8), pointer, dimension(:) :: RHOSROT_g
      REAL(KIND = 8), pointer, dimension(:) :: ROFROT_g
      REAL(KIND = 8), pointer, dimension(:) :: ZPNDROT_g
      REAL(KIND = 8), pointer, dimension(:) :: RCANROT_g
      REAL(KIND = 8), pointer, dimension(:) :: SCANROT_g
      REAL(KIND = 8), pointer, dimension(:) :: TROFROT_g
      REAL(KIND = 8), pointer, dimension(:) :: TROOROT_g
      REAL(KIND = 8), pointer, dimension(:) :: TROBROT_g
      REAL(KIND = 8), pointer, dimension(:) :: ROFOROT_g
      REAL(KIND = 8), pointer, dimension(:) :: ROFBROT_g
      REAL(KIND = 8), pointer, dimension(:) :: TROSROT_g
      REAL(KIND = 8), pointer, dimension(:) :: FSGVROT_g
      REAL(KIND = 8), pointer, dimension(:) :: FSGSROT_g
      REAL(KIND = 8), pointer, dimension(:) :: FLGVROT_g
      REAL(KIND = 8), pointer, dimension(:) :: FLGSROT_g
      REAL(KIND = 8), pointer, dimension(:) :: HFSCROT_g
      REAL(KIND = 8), pointer, dimension(:) :: HFSSROT_g
      REAL(KIND = 8), pointer, dimension(:) :: HEVCROT_g
      REAL(KIND = 8), pointer, dimension(:) :: HEVSROT_g
      REAL(KIND = 8), pointer, dimension(:) :: HMFCROT_g
      REAL(KIND = 8), pointer, dimension(:) :: HMFNROT_g
      REAL(KIND = 8), pointer, dimension(:) :: HTCSROT_g
      REAL(KIND = 8), pointer, dimension(:) :: HTCCROT_g
      REAL(KIND = 8), pointer, dimension(:) :: FSGGROT_g
      REAL(KIND = 8), pointer, dimension(:) :: FLGGROT_g
      REAL(KIND = 8), pointer, dimension(:) :: HFSGROT_g
      REAL(KIND = 8), pointer, dimension(:) :: HEVGROT_g
      REAL(KIND = 8), pointer, dimension(:) :: CDHROT_g
      REAL(KIND = 8), pointer, dimension(:) :: CDMROT_g
      REAL(KIND = 8), pointer, dimension(:) :: SFCUROT_g
      REAL(KIND = 8), pointer, dimension(:) :: SFCVROT_g
      REAL(KIND = 8), pointer, dimension(:) :: fc_g
      REAL(KIND = 8), pointer, dimension(:) :: fg_g
      REAL(KIND = 8), pointer, dimension(:) :: fcs_g
      REAL(KIND = 8), pointer, dimension(:) :: fgs_g
      REAL(KIND = 8), pointer, dimension(:) :: PCFCROT_g
      REAL(KIND = 8), pointer, dimension(:) :: PCLCROT_g
      REAL(KIND = 8), pointer, dimension(:) :: PCPGROT_g
      REAL(KIND = 8), pointer, dimension(:) :: QFCFROT_g
      REAL(KIND = 8), pointer, dimension(:) :: QFGROT_g
      REAL(KIND = 8), pointer, dimension(:,:) :: QFCROT_g
      REAL(KIND = 8), pointer, dimension(:) :: ROFCROT_g
      REAL(KIND = 8), pointer, dimension(:) :: ROFNROT_g
      REAL(KIND = 8), pointer, dimension(:) :: WTRSROT_g
      REAL(KIND = 8), pointer, dimension(:) :: WTRGROT_g
      REAL(KIND = 8), pointer, dimension(:) :: PCPNROT_g
      REAL(KIND = 8), pointer, dimension(:) :: QFCLROT_g
      REAL(KIND = 8), pointer, dimension(:) :: QFNROT_g
      REAL(KIND = 8), pointer, dimension(:) :: WTRCROT_g
      REAL(KIND = 8), pointer, dimension(:) :: gpp_g
      REAL(KIND = 8), pointer, dimension(:) :: npp_g
      REAL(KIND = 8), pointer, dimension(:) :: nbp_g
      REAL(KIND = 8), pointer, dimension(:) :: autores_g
      REAL(KIND = 8), pointer, dimension(:) :: socres_g
      REAL(KIND = 8), pointer, dimension(:) :: litres_g
      REAL(KIND = 8), pointer, dimension(:) :: dstcemls3_g
      REAL(KIND = 8), pointer, dimension(:) :: litrfall_g
      REAL(KIND = 8), pointer, dimension(:) :: rml_g
      REAL(KIND = 8), pointer, dimension(:) :: rms_g
      REAL(KIND = 8), pointer, dimension(:) :: rg_g
      REAL(KIND = 8), pointer, dimension(:) :: leaflitr_g
      REAL(KIND = 8), pointer, dimension(:) :: tltrstem_g
      REAL(KIND = 8), pointer, dimension(:) :: tltrroot_g
      REAL(KIND = 8), pointer, dimension(:) :: nep_g
      REAL(KIND = 8), pointer, dimension(:) :: hetrores_g
      REAL(KIND = 8), pointer, dimension(:) :: dstcemls_g
      REAL(KIND = 8), pointer, dimension(:) :: humiftrs_g
      REAL(KIND = 8), pointer, dimension(:) :: rmr_g
      REAL(KIND = 8), pointer, dimension(:) :: tltrleaf_g
      REAL(KIND = 8), pointer, dimension(:) :: gavgltms_g

      REAL(KIND = 8), pointer, dimension(:) :: vgbiomas_g
      REAL(KIND = 8), pointer, dimension(:) :: gavglai_g
      REAL(KIND = 8), pointer, dimension(:) :: gavgscms_g
      REAL(KIND = 8), pointer, dimension(:) :: gleafmas_g
      REAL(KIND = 8), pointer, dimension(:) :: bleafmas_g
      REAL(KIND = 8), pointer, dimension(:) :: stemmass_g
      REAL(KIND = 8), pointer, dimension(:) :: rootmass_g
      REAL(KIND = 8), pointer, dimension(:) :: litrmass_g
      REAL(KIND = 8), pointer, dimension(:) :: soilcmas_g
      REAL(KIND = 8), pointer, dimension(:) :: slai_g
      REAL(KIND = 8), pointer, dimension(:) :: ailcg_g
      REAL(KIND = 8), pointer, dimension(:) :: ailcb_g
      REAL(KIND = 8), pointer, dimension(:) :: veghght_g
      REAL(KIND = 8), pointer, dimension(:) :: rootdpth_g
      REAL(KIND = 8), pointer, dimension(:) :: roottemp_g
      REAL(KIND = 8), pointer, dimension(:) :: totcmass_g
      REAL(KIND = 8), pointer, dimension(:) :: tcanoacc_out_g
      REAL(KIND = 8), pointer, dimension(:) :: burnfrac_g
      REAL(KIND = 8), pointer, dimension(:) :: probfire_g
      REAL(KIND = 8), pointer, dimension(:) :: lucemcom_g
      REAL(KIND = 8), pointer, dimension(:) :: lucltrin_g
      REAL(KIND = 8), pointer, dimension(:) :: lucsocin_g
      REAL(KIND = 8), pointer, dimension(:) :: emit_co2_g
      REAL(KIND = 8), pointer, dimension(:) :: emit_co_g
      REAL(KIND = 8), pointer, dimension(:) :: emit_ch4_g
      REAL(KIND = 8), pointer, dimension(:) :: emit_nmhc_g
      REAL(KIND = 8), pointer, dimension(:) :: emit_h2_g
      REAL(KIND = 8), pointer, dimension(:) :: emit_nox_g
      REAL(KIND = 8), pointer, dimension(:) :: emit_n2o_g
      REAL(KIND = 8), pointer, dimension(:) :: emit_pm25_g
      REAL(KIND = 8), pointer, dimension(:) :: emit_tpm_g
      REAL(KIND = 8), pointer, dimension(:) :: emit_tc_g
      REAL(KIND = 8), pointer, dimension(:) :: emit_oc_g
      REAL(KIND = 8), pointer, dimension(:) :: emit_bc_g
      REAL(KIND = 8), pointer, dimension(:) :: bterm_g
      REAL(KIND = 8), pointer, dimension(:) :: lterm_g
      REAL(KIND = 8), pointer, dimension(:) :: mterm_g
      REAL(KIND = 8), pointer, dimension(:) :: ch4wet1_g
      REAL(KIND = 8), pointer, dimension(:) :: ch4wet2_g
      REAL(KIND = 8), pointer, dimension(:) :: wetfdyn_g
      REAL(KIND = 8), pointer, dimension(:) :: ch4dyn1_g
      REAL(KIND = 8), pointer, dimension(:) :: ch4dyn2_g
      REAL(KIND = 8), pointer, dimension(:,:) :: afrleaf_g
      REAL(KIND = 8), pointer, dimension(:,:) :: afrstem_g
      REAL(KIND = 8), pointer, dimension(:,:) :: afrroot_g
      REAL(KIND = 8), pointer, dimension(:,:) :: lfstatus_g
      REAL(KIND = 8), pointer, dimension(:,:) :: rmlvegrow_g
      REAL(KIND = 8), pointer, dimension(:,:) :: anvegrow_g
      REAL(KIND = 8), pointer, dimension(:,:) :: rmatctem_g
      REAL(KIND = 8), pointer, dimension(:,:) :: HMFGROT_g
      REAL(KIND = 8), pointer, dimension(:,:) :: HTCROT_g
      REAL(KIND = 8), pointer, dimension(:,:) :: TBARROT_g
      REAL(KIND = 8), pointer, dimension(:,:) :: THLQROT_g
      REAL(KIND = 8), pointer, dimension(:,:) :: THICROT_g
      REAL(KIND = 8), pointer, dimension(:,:) :: GFLXROT_g

!     -----------------------
!      Grid averaged monthly variables (denoted by name ending in "_mo_g")

        REAL(KIND = 8), pointer, dimension(:) :: laimaxg_mo_g
        REAL(KIND = 8), pointer, dimension(:) :: stemmass_mo_g
        REAL(KIND = 8), pointer, dimension(:) :: rootmass_mo_g
        REAL(KIND = 8), pointer, dimension(:) :: litrmass_mo_g
        REAL(KIND = 8), pointer, dimension(:) :: soilcmas_mo_g
        REAL(KIND = 8), pointer, dimension(:) :: npp_mo_g
        REAL(KIND = 8), pointer, dimension(:) :: gpp_mo_g
        REAL(KIND = 8), pointer, dimension(:) :: nep_mo_g
        REAL(KIND = 8), pointer, dimension(:) :: nbp_mo_g
        REAL(KIND = 8), pointer, dimension(:) :: hetrores_mo_g
        REAL(KIND = 8), pointer, dimension(:) :: autores_mo_g
        REAL(KIND = 8), pointer, dimension(:) :: litres_mo_g
        REAL(KIND = 8), pointer, dimension(:) :: soilcres_mo_g
        REAL(KIND = 8), pointer, dimension(:) :: vgbiomas_mo_g
        REAL(KIND = 8), pointer, dimension(:) :: totcmass_mo_g
        REAL(KIND = 8), pointer, dimension(:) :: emit_co2_mo_g
        REAL(KIND = 8), pointer, dimension(:) :: emit_co_mo_g
        REAL(KIND = 8), pointer, dimension(:) :: emit_ch4_mo_g
        REAL(KIND = 8), pointer, dimension(:) :: emit_nmhc_mo_g
        REAL(KIND = 8), pointer, dimension(:) :: emit_h2_mo_g
        REAL(KIND = 8), pointer, dimension(:) :: emit_nox_mo_g
        REAL(KIND = 8), pointer, dimension(:) :: emit_n2o_mo_g
        REAL(KIND = 8), pointer, dimension(:) :: emit_pm25_mo_g
        REAL(KIND = 8), pointer, dimension(:) :: emit_tpm_mo_g
        REAL(KIND = 8), pointer, dimension(:) :: emit_tc_mo_g
        REAL(KIND = 8), pointer, dimension(:) :: emit_oc_mo_g
        REAL(KIND = 8), pointer, dimension(:) :: emit_bc_mo_g
        REAL(KIND = 8), pointer, dimension(:) :: probfire_mo_g
        REAL(KIND = 8), pointer, dimension(:) :: luc_emc_mo_g
        REAL(KIND = 8), pointer, dimension(:) :: lucltrin_mo_g
        REAL(KIND = 8), pointer, dimension(:) :: lucsocin_mo_g
        REAL(KIND = 8), pointer, dimension(:) :: burnfrac_mo_g
        REAL(KIND = 8), pointer, dimension(:) :: bterm_mo_g
        REAL(KIND = 8), pointer, dimension(:) :: lterm_mo_g
        REAL(KIND = 8), pointer, dimension(:) :: mterm_mo_g
        REAL(KIND = 8), pointer, dimension(:) :: ch4wet1_mo_g
        REAL(KIND = 8), pointer, dimension(:) :: ch4wet2_mo_g
        REAL(KIND = 8), pointer, dimension(:) :: wetfdyn_mo_g
        REAL(KIND = 8), pointer, dimension(:) :: ch4dyn1_mo_g
        REAL(KIND = 8), pointer, dimension(:) :: ch4dyn2_mo_g

!      Mosaic monthly variables (denoted by name ending in "_mo_m")
c
      REAL(KIND = 8), pointer, dimension(:,:,:) :: laimaxg_mo_m
      REAL(KIND = 8), pointer, dimension(:,:,:) :: stemmass_mo_m
      REAL(KIND = 8), pointer, dimension(:,:,:) :: rootmass_mo_m
      REAL(KIND = 8), pointer, dimension(:,:,:) :: npp_mo_m
      REAL(KIND = 8), pointer, dimension(:,:,:) :: gpp_mo_m
      REAL(KIND = 8), pointer, dimension(:,:,:) :: vgbiomas_mo_m
      REAL(KIND = 8), pointer, dimension(:,:,:) :: autores_mo_m
      REAL(KIND = 8), pointer, dimension(:,:,:) :: totcmass_mo_m
      REAL(KIND = 8), pointer, dimension(:,:,:) :: litrmass_mo_m
      REAL(KIND = 8), pointer, dimension(:,:,:) :: soilcmas_mo_m
      REAL(KIND = 8), pointer, dimension(:,:,:) :: nep_mo_m
      REAL(KIND = 8), pointer, dimension(:,:,:) :: litres_mo_m
      REAL(KIND = 8), pointer, dimension(:,:,:) :: soilcres_mo_m
      REAL(KIND = 8), pointer, dimension(:,:,:) :: hetrores_mo_m
      REAL(KIND = 8), pointer, dimension(:,:,:) :: nbp_mo_m
      REAL(KIND = 8), pointer, dimension(:,:,:) :: emit_co2_mo_m
      REAL(KIND = 8), pointer, dimension(:,:,:) :: emit_co_mo_m
      REAL(KIND = 8), pointer, dimension(:,:,:) :: emit_ch4_mo_m
      REAL(KIND = 8), pointer, dimension(:,:,:) :: emit_nmhc_mo_m
      REAL(KIND = 8), pointer, dimension(:,:,:) :: emit_h2_mo_m
      REAL(KIND = 8), pointer, dimension(:,:,:) :: emit_nox_mo_m
      REAL(KIND = 8), pointer, dimension(:,:,:) :: emit_n2o_mo_m
      REAL(KIND = 8), pointer, dimension(:,:,:) :: emit_pm25_mo_m
      REAL(KIND = 8), pointer, dimension(:,:,:) :: emit_tpm_mo_m
      REAL(KIND = 8), pointer, dimension(:,:,:) :: emit_tc_mo_m
      REAL(KIND = 8), pointer, dimension(:,:,:) :: emit_oc_mo_m
      REAL(KIND = 8), pointer, dimension(:,:,:) :: emit_bc_mo_m
      REAL(KIND = 8), pointer, dimension(:,:,:) :: burnfrac_mo_m
      REAL(KIND = 8), pointer, dimension(:,:) :: probfire_mo_m
      REAL(KIND = 8), pointer, dimension(:,:) :: bterm_mo_m
      REAL(KIND = 8), pointer, dimension(:,:) :: luc_emc_mo_m
      REAL(KIND = 8), pointer, dimension(:,:) :: lterm_mo_m
      REAL(KIND = 8), pointer, dimension(:,:) :: lucsocin_mo_m
      REAL(KIND = 8), pointer, dimension(:,:) :: mterm_mo_m
      REAL(KIND = 8), pointer, dimension(:,:) :: lucltrin_mo_m
      REAL(KIND = 8), pointer, dimension(:,:) :: ch4wet1_mo_m
      REAL(KIND = 8), pointer, dimension(:,:) :: ch4wet2_mo_m
      REAL(KIND = 8), pointer, dimension(:,:) :: wetfdyn_mo_m
      REAL(KIND = 8), pointer, dimension(:,:) :: ch4dyn1_mo_m
      REAL(KIND = 8), pointer, dimension(:,:) :: ch4dyn2_mo_m

!     -----------------------
c      Annual output for CTEM grid-averaged variables:
c      (denoted by name ending in "_yr_g")

      REAL(KIND = 8), pointer, dimension(:) :: laimaxg_yr_g
      REAL(KIND = 8), pointer, dimension(:) :: stemmass_yr_g
      REAL(KIND = 8), pointer, dimension(:) :: rootmass_yr_g
      REAL(KIND = 8), pointer, dimension(:) :: litrmass_yr_g
      REAL(KIND = 8), pointer, dimension(:) :: soilcmas_yr_g
      REAL(KIND = 8), pointer, dimension(:) :: npp_yr_g
      REAL(KIND = 8), pointer, dimension(:) :: gpp_yr_g
      REAL(KIND = 8), pointer, dimension(:) :: nep_yr_g
      REAL(KIND = 8), pointer, dimension(:) :: nbp_yr_g
      REAL(KIND = 8), pointer, dimension(:) :: hetrores_yr_g
      REAL(KIND = 8), pointer, dimension(:) :: autores_yr_g
      REAL(KIND = 8), pointer, dimension(:) :: litres_yr_g
      REAL(KIND = 8), pointer, dimension(:) :: soilcres_yr_g
      REAL(KIND = 8), pointer, dimension(:) :: vgbiomas_yr_g
      REAL(KIND = 8), pointer, dimension(:) :: totcmass_yr_g
      REAL(KIND = 8), pointer, dimension(:) :: emit_co2_yr_g
      REAL(KIND = 8), pointer, dimension(:) :: emit_co_yr_g
      REAL(KIND = 8), pointer, dimension(:) :: emit_ch4_yr_g
      REAL(KIND = 8), pointer, dimension(:) :: emit_nmhc_yr_g
      REAL(KIND = 8), pointer, dimension(:) :: emit_h2_yr_g
      REAL(KIND = 8), pointer, dimension(:) :: emit_nox_yr_g
      REAL(KIND = 8), pointer, dimension(:) :: emit_n2o_yr_g
      REAL(KIND = 8), pointer, dimension(:) :: emit_pm25_yr_g
      REAL(KIND = 8), pointer, dimension(:) :: emit_tpm_yr_g
      REAL(KIND = 8), pointer, dimension(:) :: emit_tc_yr_g
      REAL(KIND = 8), pointer, dimension(:) :: emit_oc_yr_g
      REAL(KIND = 8), pointer, dimension(:) :: emit_bc_yr_g
      REAL(KIND = 8), pointer, dimension(:) :: probfire_yr_g
      REAL(KIND = 8), pointer, dimension(:) :: luc_emc_yr_g
      REAL(KIND = 8), pointer, dimension(:) :: lucltrin_yr_g
      REAL(KIND = 8), pointer, dimension(:) :: lucsocin_yr_g
      REAL(KIND = 8), pointer, dimension(:) :: burnfrac_yr_g
      REAL(KIND = 8), pointer, dimension(:) :: bterm_yr_g
      REAL(KIND = 8), pointer, dimension(:) :: lterm_yr_g
      REAL(KIND = 8), pointer, dimension(:) :: mterm_yr_g
      REAL(KIND = 8), pointer, dimension(:) :: ch4wet1_yr_g
      REAL(KIND = 8), pointer, dimension(:) :: ch4wet2_yr_g
      REAL(KIND = 8), pointer, dimension(:) :: wetfdyn_yr_g
      REAL(KIND = 8), pointer, dimension(:) :: ch4dyn1_yr_g
      REAL(KIND = 8), pointer, dimension(:) :: ch4dyn2_yr_g


! c      Annual output for CTEM mosaic variables:
! c      (denoted by name ending in "_yr_m")
!
      REAL(KIND = 8), pointer, dimension(:,:,:) :: laimaxg_yr_m
      REAL(KIND = 8), pointer, dimension(:,:,:) :: stemmass_yr_m
      REAL(KIND = 8), pointer, dimension(:,:,:) :: rootmass_yr_m
      REAL(KIND = 8), pointer, dimension(:,:,:) :: npp_yr_m
      REAL(KIND = 8), pointer, dimension(:,:,:) :: gpp_yr_m
      REAL(KIND = 8), pointer, dimension(:,:,:) :: vgbiomas_yr_m
      REAL(KIND = 8), pointer, dimension(:,:,:) :: autores_yr_m
      REAL(KIND = 8), pointer, dimension(:,:,:) :: totcmass_yr_m
      REAL(KIND = 8), pointer, dimension(:,:,:) :: litrmass_yr_m
      REAL(KIND = 8), pointer, dimension(:,:,:) :: soilcmas_yr_m
      REAL(KIND = 8), pointer, dimension(:,:,:) :: nep_yr_m
      REAL(KIND = 8), pointer, dimension(:,:,:) :: litres_yr_m
      REAL(KIND = 8), pointer, dimension(:,:,:) :: soilcres_yr_m
      REAL(KIND = 8), pointer, dimension(:,:,:) :: hetrores_yr_m
      REAL(KIND = 8), pointer, dimension(:,:,:) :: nbp_yr_m
      REAL(KIND = 8), pointer, dimension(:,:,:) :: emit_co2_yr_m
      REAL(KIND = 8), pointer, dimension(:,:,:) :: emit_co_yr_m
      REAL(KIND = 8), pointer, dimension(:,:,:) :: emit_ch4_yr_m
      REAL(KIND = 8), pointer, dimension(:,:,:) :: emit_nmhc_yr_m
      REAL(KIND = 8), pointer, dimension(:,:,:) :: emit_h2_yr_m
      REAL(KIND = 8), pointer, dimension(:,:,:) :: emit_nox_yr_m
      REAL(KIND = 8), pointer, dimension(:,:,:) :: emit_n2o_yr_m
      REAL(KIND = 8), pointer, dimension(:,:,:) :: emit_pm25_yr_m
      REAL(KIND = 8), pointer, dimension(:,:,:) :: emit_tpm_yr_m
      REAL(KIND = 8), pointer, dimension(:,:,:) :: emit_tc_yr_m
      REAL(KIND = 8), pointer, dimension(:,:,:) :: emit_oc_yr_m
      REAL(KIND = 8), pointer, dimension(:,:,:) :: emit_bc_yr_m
      REAL(KIND = 8), pointer, dimension(:,:,:) :: burnfrac_yr_m
      REAL(KIND = 8), pointer, dimension(:,:) :: probfire_yr_m
      REAL(KIND = 8), pointer, dimension(:,:) :: bterm_yr_m
      REAL(KIND = 8), pointer, dimension(:,:) :: luc_emc_yr_m
      REAL(KIND = 8), pointer, dimension(:,:) :: lterm_yr_m
      REAL(KIND = 8), pointer, dimension(:,:) :: lucsocin_yr_m
      REAL(KIND = 8), pointer, dimension(:,:) :: mterm_yr_m
      REAL(KIND = 8), pointer, dimension(:,:) :: lucltrin_yr_m
      REAL(KIND = 8), pointer, dimension(:,:) :: ch4wet1_yr_m
      REAL(KIND = 8), pointer, dimension(:,:) :: ch4wet2_yr_m
      REAL(KIND = 8), pointer, dimension(:,:) :: wetfdyn_yr_m
      REAL(KIND = 8), pointer, dimension(:,:) :: ch4dyn1_yr_m
      REAL(KIND = 8), pointer, dimension(:,:) :: ch4dyn2_yr_m

      logical, parameter :: obslght = .false.  ! if true the observed lightning will be used. False means you will use the
                                              ! lightning climatology from the CTM file.
c
c============= CTEM array declaration done =============================/
C
C=======================================================================
C     * PHYSICAL CONSTANTS.
C     * PARAMETERS IN THE FOLLOWING COMMON BLOCKS ARE NORMALLY DEFINED
C     * WITHIN THE GCM.

      COMMON /PARAMS/ X1,    X2,    X3,    X4,   G,GAS,   X5,
     1                X6,    CPRES, GASV,  X7
      COMMON /PARAM1/ CPI,   X8,    CELZRO,X9,    X10,    X11
      COMMON /PARAM3/ X12,   X13,   X14,   X15,   SIGMA,  X16
      COMMON  /TIMES/ DELTIM,K1,    K2,    K3,    K4,     K5,
     1                K6,    K7,    K8,    K9,    K10,    K11
C
C     * THE FOLLOWING COMMON BLOCKS ARE DEFINED SPECIFICALLY FOR USE
C     * IN CLASS, VIA BLOCK DATA AND THE SUBROUTINE "CLASSD".
C
      COMMON /CLASS1/ DELT,TFREZ
      COMMON /CLASS2/ RGAS,RGASV,GRAV,SBC,VKC,CT,VMIN
      COMMON /CLASS3/ TCW,TCICE,TCSAND,TCCLAY,TCOM,TCDRYS,
     1                RHOSOL,RHOOM
      COMMON /CLASS4/ HCPW,HCPICE,HCPSOL,HCPOM,HCPSND,HCPCLY,
     1                SPHW,SPHICE,SPHVEG,SPHAIR,RHOW,RHOICE,
     2                TCGLAC,CLHMLT,CLHVAP
      COMMON /CLASS5/ THPORG,THRORG,THMORG,BORG,PSISORG,GRKSORG
      COMMON /CLASS6/ PI,GROWYR,ZOLNG,ZOLNS,ZOLNI,ZORAT,ZORATG
      COMMON /CLASS7/ CANEXT,XLEAF
      COMMON /CLASS8/ ALVSI,ALIRI,ALVSO,ALIRO,ALBRCK
      COMMON /PHYCON/ DELTA,CGRAV,CKARM,CPD
      COMMON /CLASSD2/ AS,ASX,CI,BS,BETA,FACTN,HMIN,ANGMAX
C
      CALL CLASSD
C
      ZDMROW(1)=10.0
      ZDHROW(1)=2.0
      NTLD=NMOS
      CUMSNO = 0.0
C
C===================== CTEM ==============================================\

    ! Point pointers

      ALVSACC_MO        => class_out%ALVSACC_MO
      ALIRACC_MO        => class_out%ALIRACC_MO
      FLUTACC_MO        => class_out%FLUTACC_MO
      FSINACC_MO        => class_out%FSINACC_MO
      FLINACC_MO        => class_out%FLINACC_MO
      HFSACC_MO         => class_out%HFSACC_MO
      QEVPACC_MO        => class_out%QEVPACC_MO
      SNOACC_MO         => class_out%SNOACC_MO
      WSNOACC_MO        => class_out%WSNOACC_MO
      ROFACC_MO         => class_out%ROFACC_MO
      PREACC_MO         => class_out%PREACC_MO
      EVAPACC_MO        => class_out%EVAPACC_MO
      TAACC_MO          => class_out%TAACC_MO
      FSSTAR_MO         => class_out%FSSTAR_MO
      FLSTAR_MO         => class_out%FLSTAR_MO
      QH_MO             => class_out%QH_MO
      QE_MO             => class_out%QE_MO
      TBARACC_MO        => class_out%TBARACC_MO
      THLQACC_MO        => class_out%THLQACC_MO
      THICACC_MO        => class_out%THICACC_MO
      ALVSACC_YR        => class_out%ALVSACC_YR
      ALIRACC_YR        => class_out%ALIRACC_YR
      FLUTACC_YR        => class_out%FLUTACC_YR
      FSINACC_YR        => class_out%FSINACC_YR
      FLINACC_YR        => class_out%FLINACC_YR
      HFSACC_YR         => class_out%HFSACC_YR
      QEVPACC_YR        => class_out%QEVPACC_YR
      ROFACC_YR         => class_out%ROFACC_YR
      PREACC_YR         => class_out%PREACC_YR
      EVAPACC_YR        => class_out%EVAPACC_YR
      TAACC_YR          => class_out%TAACC_YR
      FSSTAR_YR         => class_out%FSSTAR_YR
      FLSTAR_YR         => class_out%FLSTAR_YR
      QH_YR             => class_out%QH_YR
      QE_YR             => class_out%QE_YR

      ctem_on           => c_switch%ctem_on
      parallelrun       => c_switch%parallelrun
      mosaic            => c_switch%mosaic
      cyclemet          => c_switch%cyclemet
      dofire            => c_switch%dofire
      run_model         => c_switch%run_model
      met_rewound       => c_switch%met_rewound
      reach_eof         => c_switch%reach_eof
      compete           => c_switch%compete
      start_bare        => c_switch%start_bare
      rsfile            => c_switch%rsfile
      lnduseon          => c_switch%lnduseon
      co2on             => c_switch%co2on
      popdon            => c_switch%popdon
      inibioclim        => c_switch%inibioclim
      start_from_rs     => c_switch%start_from_rs
      dowetlands        => c_switch%dowetlands
      obswetf           => c_switch%obswetf
      transient_run     => c_switch%transient_run

      tcanrs            => vrot%tcanrs
      tsnors            => vrot%tsnors
      tpndrs            => vrot%tpndrs
      csum              => vrot%csum
      tbaraccrow_m      => vrot%tbaraccrow_m
      tcanoaccrow_m     => vrot%tcanoaccrow_m
      uvaccrow_m        => vrot%uvaccrow_m
      vvaccrow_m        => vrot%vvaccrow_m

      ! ROW:
      ailcminrow        => vrot%ailcmin
      ailcmaxrow        => vrot%ailcmax
      dvdfcanrow        => vrot%dvdfcan
      gleafmasrow       => vrot%gleafmas
      bleafmasrow       => vrot%bleafmas
      stemmassrow       => vrot%stemmass
      rootmassrow       => vrot%rootmass
      pstemmassrow      => vrot%pstemmass
      pgleafmassrow     => vrot%pgleafmass
      fcancmxrow        => vrot%fcancmx
      gavglairow        => vrot%gavglai
      zolncrow          => vrot%zolnc
      ailcrow           => vrot%ailc
      ailcgrow          => vrot%ailcg
      ailcgsrow         => vrot%ailcgs
      fcancsrow         => vrot%fcancs
      fcancrow          => vrot%fcanc
      co2concrow        => vrot%co2conc
      co2i1cgrow        => vrot%co2i1cg
      co2i1csrow        => vrot%co2i1cs
      co2i2cgrow        => vrot%co2i2cg
      co2i2csrow        => vrot%co2i2cs
      ancsvegrow        => vrot%ancsveg
      ancgvegrow        => vrot%ancgveg
      rmlcsvegrow       => vrot%rmlcsveg
      rmlcgvegrow       => vrot%rmlcgveg
      slairow           => vrot%slai
      ailcbrow          => vrot%ailcb
      canresrow         => vrot%canres
      flhrlossrow       => vrot%flhrloss

      tcanoaccrow_out   => vrot%tcanoaccrow_out
      qevpacc_m_save    => vrot%qevpacc_m_save

      grwtheffrow       => vrot%grwtheff
      lystmmasrow       => vrot%lystmmas
      lyrotmasrow       => vrot%lyrotmas
      tymaxlairow       => vrot%tymaxlai
      vgbiomasrow       => vrot%vgbiomas
      gavgltmsrow       => vrot%gavgltms
      gavgscmsrow       => vrot%gavgscms
      stmhrlosrow       => vrot%stmhrlos
      rmatcrow          => vrot%rmatc
      rmatctemrow       => vrot%rmatctem
      litrmassrow       => vrot%litrmass
      soilcmasrow       => vrot%soilcmas
      vgbiomas_vegrow   => vrot%vgbiomas_veg

      emit_co2row       => vrot%emit_co2
      emit_corow        => vrot%emit_co
      emit_ch4row       => vrot%emit_ch4
      emit_nmhcrow      => vrot%emit_nmhc
      emit_h2row        => vrot%emit_h2
      emit_noxrow       => vrot%emit_nox
      emit_n2orow       => vrot%emit_n2o
      emit_pm25row      => vrot%emit_pm25
      emit_tpmrow       => vrot%emit_tpm
      emit_tcrow        => vrot%emit_tc
      emit_ocrow        => vrot%emit_oc
      emit_bcrow        => vrot%emit_bc
      burnfracrow       => vrot%burnfrac
      burnvegfrow       => vrot%burnvegf
      probfirerow       => vrot%probfire
      btermrow          => vrot%bterm
      ltermrow          => vrot%lterm
      mtermrow          => vrot%mterm

      extnprobgrd       => vrot%extnprob
      prbfrhucgrd       => vrot%prbfrhuc
      mlightnggrd       => vrot%mlightng

      bmasvegrow        => vrot%bmasveg
      cmasvegcrow       => vrot%cmasvegc
      veghghtrow        => vrot%veghght
      rootdpthrow       => vrot%rootdpth
      rmlrow            => vrot%rml
      rmsrow            => vrot%rms
      tltrleafrow       => vrot%tltrleaf
      tltrstemrow       => vrot%tltrstem
      tltrrootrow       => vrot%tltrroot
      leaflitrrow       => vrot%leaflitr
      roottemprow       => vrot%roottemp
      afrleafrow        => vrot%afrleaf
      afrstemrow        => vrot%afrstem
      afrrootrow        => vrot%afrroot
      wtstatusrow       => vrot%wtstatus
      ltstatusrow       => vrot%ltstatus
      rmrrow            => vrot%rmr

      wetfracgrd        => vrot%wetfrac
      ch4wet1row        => vrot%ch4wet1
      ch4wet2row        => vrot%ch4wet2
      wetfdynrow        => vrot%wetfdyn
      ch4dyn1row        => vrot%ch4dyn1
      ch4dyn2row        => vrot%ch4dyn2
      wetfrac_mon       => vrot%wetfrac_mon

      lucemcomrow       => vrot%lucemcom
      lucltrinrow       => vrot%lucltrin
      lucsocinrow       => vrot%lucsocin

      npprow            => vrot%npp
      neprow            => vrot%nep
      nbprow            => vrot%nbp
      gpprow            => vrot%gpp
      hetroresrow       => vrot%hetrores
      autoresrow        => vrot%autores
      soilcresprow      => vrot%soilcresp
      rmrow             => vrot%rm
      rgrow             => vrot%rg
      litresrow         => vrot%litres
      socresrow         => vrot%socres
      dstcemlsrow       => vrot%dstcemls
      litrfallrow       => vrot%litrfall
      humiftrsrow       => vrot%humiftrs

      gppvegrow         => vrot%gppveg
      nepvegrow         => vrot%nepveg
      nbpvegrow         => vrot%nbpveg
      nppvegrow         => vrot%nppveg
      hetroresvegrow    => vrot%hetroresveg
      autoresvegrow     => vrot%autoresveg
      litresvegrow      => vrot%litresveg
      soilcresvegrow    => vrot%soilcresveg
      rmlvegaccrow      => vrot%rmlvegacc
      rmsvegrow         => vrot%rmsveg
      rmrvegrow         => vrot%rmrveg
      rgvegrow          => vrot%rgveg

      rothrlosrow       => vrot%rothrlos
      pfcancmxrow       => vrot%pfcancmx
      nfcancmxrow       => vrot%nfcancmx
      alvsctmrow        => vrot%alvsctm
      paicrow           => vrot%paic
      slaicrow          => vrot%slaic
      alirctmrow        => vrot%alirctm
      cfluxcgrow        => vrot%cfluxcg
      cfluxcsrow        => vrot%cfluxcs
      dstcemls3row      => vrot%dstcemls3
      anvegrow          => vrot%anveg
      rmlvegrow         => vrot%rmlveg

      pftexistrow       => vrot%pftexist
      colddaysrow       => vrot%colddays
      icountrow         => vrot%icount
      lfstatusrow       => vrot%lfstatus
      pandaysrow        => vrot%pandays
      stdalngrd         => vrot%stdaln


      ! >>>>>>>>>>>>>>>>>>>>>>>>>>
      ! GAT:

      lightng           => vgat%lightng
      tcanoaccgat_out   => vgat%tcanoaccgat_out

      ailcmingat        => vgat%ailcmin
      ailcmaxgat        => vgat%ailcmax
      dvdfcangat        => vgat%dvdfcan
      gleafmasgat       => vgat%gleafmas
      bleafmasgat       => vgat%bleafmas
      stemmassgat       => vgat%stemmass
      rootmassgat       => vgat%rootmass
      pstemmassgat      => vgat%pstemmass
      pgleafmassgat     => vgat%pgleafmass
      fcancmxgat        => vgat%fcancmx
      gavglaigat        => vgat%gavglai
      zolncgat          => vgat%zolnc
      ailcgat           => vgat%ailc
      ailcggat          => vgat%ailcg
      ailcgsgat         => vgat%ailcgs
      fcancsgat         => vgat%fcancs
      fcancgat          => vgat%fcanc
      co2concgat        => vgat%co2conc
      co2i1cggat        => vgat%co2i1cg
      co2i1csgat        => vgat%co2i1cs
      co2i2cggat        => vgat%co2i2cg
      co2i2csgat        => vgat%co2i2cs
      ancsveggat        => vgat%ancsveg
      ancgveggat        => vgat%ancgveg
      rmlcsveggat       => vgat%rmlcsveg
      rmlcgveggat       => vgat%rmlcgveg
      slaigat           => vgat%slai
      ailcbgat          => vgat%ailcb
      canresgat         => vgat%canres
      flhrlossgat       => vgat%flhrloss

      grwtheffgat       => vgat%grwtheff
      lystmmasgat       => vgat%lystmmas
      lyrotmasgat       => vgat%lyrotmas
      tymaxlaigat       => vgat%tymaxlai
      vgbiomasgat       => vgat%vgbiomas
      gavgltmsgat       => vgat%gavgltms
      gavgscmsgat       => vgat%gavgscms
      stmhrlosgat       => vgat%stmhrlos
      rmatcgat          => vgat%rmatc
      rmatctemgat       => vgat%rmatctem
      litrmassgat       => vgat%litrmass
      soilcmasgat       => vgat%soilcmas
      vgbiomas_veggat   => vgat%vgbiomas_veg

      emit_co2gat       => vgat%emit_co2
      emit_cogat        => vgat%emit_co
      emit_ch4gat       => vgat%emit_ch4
      emit_nmhcgat      => vgat%emit_nmhc
      emit_h2gat        => vgat%emit_h2
      emit_noxgat       => vgat%emit_nox
      emit_n2ogat       => vgat%emit_n2o
      emit_pm25gat      => vgat%emit_pm25
      emit_tpmgat       => vgat%emit_tpm
      emit_tcgat        => vgat%emit_tc
      emit_ocgat        => vgat%emit_oc
      emit_bcgat        => vgat%emit_bc
      burnfracgat       => vgat%burnfrac
      burnvegfgat       => vgat%burnvegf
      probfiregat       => vgat%probfire
      btermgat          => vgat%bterm
      ltermgat          => vgat%lterm
      mtermgat          => vgat%mterm

      extnprobgat       => vgat%extnprob
      prbfrhucgat       => vgat%prbfrhuc
      mlightnggat       => vgat%mlightng

      bmasveggat        => vgat%bmasveg
      cmasvegcgat       => vgat%cmasvegc
      veghghtgat        => vgat%veghght
      rootdpthgat       => vgat%rootdpth
      rmlgat            => vgat%rml
      rmsgat            => vgat%rms
      tltrleafgat       => vgat%tltrleaf
      tltrstemgat       => vgat%tltrstem
      tltrrootgat       => vgat%tltrroot
      leaflitrgat       => vgat%leaflitr
      roottempgat       => vgat%roottemp
      afrleafgat        => vgat%afrleaf
      afrstemgat        => vgat%afrstem
      afrrootgat        => vgat%afrroot
      wtstatusgat       => vgat%wtstatus
      ltstatusgat       => vgat%ltstatus
      rmrgat            => vgat%rmr

      wetfrac_sgrd      => vgat%wetfrac_s
      ch4wet1gat        => vgat%ch4wet1
      ch4wet2gat        => vgat%ch4wet2
      wetfdyngat        => vgat%wetfdyn
      ch4dyn1gat        => vgat%ch4dyn1
      ch4dyn2gat        => vgat%ch4dyn2

      lucemcomgat       => vgat%lucemcom
      lucltringat       => vgat%lucltrin
      lucsocingat       => vgat%lucsocin

      nppgat            => vgat%npp
      nepgat            => vgat%nep
      nbpgat            => vgat%nbp
      gppgat            => vgat%gpp
      hetroresgat       => vgat%hetrores
      autoresgat        => vgat%autores
      soilcrespgat      => vgat%soilcresp
      rmgat             => vgat%rm
      rggat             => vgat%rg
      litresgat         => vgat%litres
      socresgat         => vgat%socres
      dstcemlsgat       => vgat%dstcemls
      litrfallgat       => vgat%litrfall
      humiftrsgat       => vgat%humiftrs

      gppveggat         => vgat%gppveg
      nepveggat         => vgat%nepveg
      nbpveggat         => vgat%nbpveg
      nppveggat         => vgat%nppveg
      hetroresveggat    => vgat%hetroresveg
      autoresveggat     => vgat%autoresveg
      litresveggat      => vgat%litresveg
      soilcresveggat    => vgat%soilcresveg
      rmlvegaccgat      => vgat%rmlvegacc
      rmsveggat         => vgat%rmsveg
      rmrveggat         => vgat%rmrveg
      rgveggat          => vgat%rgveg

      rothrlosgat       => vgat%rothrlos
      pfcancmxgat       => vgat%pfcancmx
      nfcancmxgat       => vgat%nfcancmx
      alvsctmgat        => vgat%alvsctm
      paicgat           => vgat%paic
      slaicgat          => vgat%slaic
      alirctmgat        => vgat%alirctm
      cfluxcggat        => vgat%cfluxcg
      cfluxcsgat        => vgat%cfluxcs
      dstcemls3gat      => vgat%dstcemls3
      anveggat          => vgat%anveg
      rmlveggat         => vgat%rmlveg

      twarmm            => vgat%twarmm
      tcoldm            => vgat%tcoldm
      gdd5              => vgat%gdd5
      aridity           => vgat%aridity
      srplsmon          => vgat%srplsmon
      defctmon          => vgat%defctmon
      anndefct          => vgat%anndefct
      annsrpls          => vgat%annsrpls
      annpcp            => vgat%annpcp
      dry_season_length => vgat%dry_season_length

      pftexistgat       => vgat%pftexist
      colddaysgat       => vgat%colddays
      icountgat         => vgat%icount
      lfstatusgat       => vgat%lfstatus
      pandaysgat        => vgat%pandays
      stdalngat         => vgat%stdaln

      ! Mosaic-level:

      PREACC_M          => vrot%PREACC_M
      GTACC_M           => vrot%GTACC_M
      QEVPACC_M         => vrot%QEVPACC_M
      HFSACC_M          => vrot%HFSACC_M
      HMFNACC_M         => vrot%HMFNACC_M
      ROFACC_M          => vrot%ROFACC_M
      SNOACC_M          => vrot%SNOACC_M
      OVRACC_M          => vrot%OVRACC_M
      WTBLACC_M         => vrot%WTBLACC_M
      TBARACC_M         => vrot%TBARACC_M
      THLQACC_M         => vrot%THLQACC_M
      THICACC_M         => vrot%THICACC_M
      THALACC_M         => vrot%THALACC_M
      ALVSACC_M         => vrot%ALVSACC_M
      ALIRACC_M         => vrot%ALIRACC_M
      RHOSACC_M         => vrot%RHOSACC_M
      TSNOACC_M         => vrot%TSNOACC_M
      WSNOACC_M         => vrot%WSNOACC_M
      SNOARE_M          => vrot%SNOARE_M
      TCANACC_M         => vrot%TCANACC_M
      RCANACC_M         => vrot%RCANACC_M
      SCANACC_M         => vrot%SCANACC_M
      GROACC_M          => vrot%GROACC_M
      FSINACC_M         => vrot%FSINACC_M
      FLINACC_M         => vrot%FLINACC_M
      TAACC_M           => vrot%TAACC_M
      UVACC_M           => vrot%UVACC_M
      PRESACC_M         => vrot%PRESACC_M
      QAACC_M           => vrot%QAACC_M
      EVAPACC_M         => vrot%EVAPACC_M
      FLUTACC_M         => vrot%FLUTACC_M

      ! grid-averaged

      WSNOROT_g         => ctem_grd%WSNOROT_g
      ROFSROT_g         => ctem_grd%ROFSROT_g
      SNOROT_g          => ctem_grd%SNOROT_g
      RHOSROT_g         => ctem_grd%RHOSROT_g
      ROFROT_g          => ctem_grd%ROFROT_g
      ZPNDROT_g         => ctem_grd%ZPNDROT_g
      RCANROT_g         => ctem_grd%RCANROT_g
      SCANROT_g         => ctem_grd%SCANROT_g
      TROFROT_g         => ctem_grd%TROFROT_g
      TROOROT_g         => ctem_grd%TROOROT_g
      TROBROT_g         => ctem_grd%TROBROT_g
      ROFOROT_g         => ctem_grd%ROFOROT_g
      ROFBROT_g         => ctem_grd%ROFBROT_g
      TROSROT_g         => ctem_grd%TROSROT_g
      FSGVROT_g         => ctem_grd%FSGVROT_g
      FSGSROT_g         => ctem_grd%FSGSROT_g
      FLGVROT_g         => ctem_grd%FLGVROT_g
      FLGSROT_g         => ctem_grd%FLGSROT_g
      HFSCROT_g         => ctem_grd%HFSCROT_g
      HFSSROT_g         => ctem_grd%HFSSROT_g
      HEVCROT_g         => ctem_grd%HEVCROT_g
      HEVSROT_g         => ctem_grd%HEVSROT_g
      HMFCROT_g         => ctem_grd%HMFCROT_g
      HMFNROT_g         => ctem_grd%HMFNROT_g
      HTCSROT_g         => ctem_grd%HTCSROT_g
      HTCCROT_g         => ctem_grd%HTCCROT_g
      FSGGROT_g         => ctem_grd%FSGGROT_g
      FLGGROT_g         => ctem_grd%FLGGROT_g
      HFSGROT_g         => ctem_grd%HFSGROT_g
      HEVGROT_g         => ctem_grd%HEVGROT_g
      CDHROT_g          => ctem_grd%CDHROT_g
      CDMROT_g          => ctem_grd%CDMROT_g
      SFCUROT_g         => ctem_grd%SFCUROT_g
      SFCVROT_g         => ctem_grd%SFCVROT_g
      fc_g              => ctem_grd%fc_g
      fg_g              => ctem_grd%fg_g
      fcs_g             => ctem_grd%fcs_g
      fgs_g             => ctem_grd%fgs_g
      PCFCROT_g         => ctem_grd%PCFCROT_g
      PCLCROT_g         => ctem_grd%PCLCROT_g
      PCPGROT_g         => ctem_grd%PCPGROT_g
      QFCFROT_g         => ctem_grd%QFCFROT_g
      QFGROT_g          => ctem_grd%QFGROT_g
      QFCROT_g          => ctem_grd%QFCROT_g
      ROFCROT_g         => ctem_grd%ROFCROT_g
      ROFNROT_g         => ctem_grd%ROFNROT_g
      WTRSROT_g         => ctem_grd%WTRSROT_g
      WTRGROT_g         => ctem_grd%WTRGROT_g
      PCPNROT_g         => ctem_grd%PCPNROT_g
      QFCLROT_g         => ctem_grd%QFCLROT_g
      QFNROT_g          => ctem_grd%QFNROT_g
      WTRCROT_g         => ctem_grd%WTRCROT_g
      gpp_g             => ctem_grd%gpp_g
      npp_g             => ctem_grd%npp_g
      nbp_g             => ctem_grd%nbp_g
      autores_g         => ctem_grd%autores_g
      socres_g          => ctem_grd%socres_g
      litres_g          => ctem_grd%litres_g
      dstcemls3_g       => ctem_grd%dstcemls3_g
      litrfall_g        => ctem_grd%litrfall_g
      rml_g             => ctem_grd%rml_g
      rms_g             => ctem_grd%rms_g
      rg_g              => ctem_grd%rg_g
      leaflitr_g        => ctem_grd%leaflitr_g
      tltrstem_g        => ctem_grd%tltrstem_g
      tltrroot_g        => ctem_grd%tltrroot_g
      nep_g             => ctem_grd%nep_g
      hetrores_g        => ctem_grd%hetrores_g
      dstcemls_g        => ctem_grd%dstcemls_g
      humiftrs_g        => ctem_grd%humiftrs_g
      rmr_g             => ctem_grd%rmr_g
      tltrleaf_g        => ctem_grd%tltrleaf_g
      gavgltms_g        => ctem_grd%gavgltms_g
      vgbiomas_g        => ctem_grd%vgbiomas_g
      gavglai_g         => ctem_grd%gavglai_g
      gavgscms_g        => ctem_grd%gavgscms_g
      gleafmas_g        => ctem_grd%gleafmas_g
      bleafmas_g        => ctem_grd%bleafmas_g
      stemmass_g        => ctem_grd%stemmass_g
      rootmass_g        => ctem_grd%rootmass_g
      litrmass_g        => ctem_grd%litrmass_g
      soilcmas_g        => ctem_grd%soilcmas_g
      slai_g            => ctem_grd%slai_g
      ailcg_g           => ctem_grd%ailcg_g
      ailcb_g           => ctem_grd%ailcb_g
      veghght_g         => ctem_grd%veghght_g
      rootdpth_g        => ctem_grd%rootdpth_g
      roottemp_g        => ctem_grd%roottemp_g
      totcmass_g        => ctem_grd%totcmass_g
      tcanoacc_out_g    => ctem_grd%tcanoacc_out_g
      burnfrac_g        => ctem_grd%burnfrac_g
      probfire_g        => ctem_grd%probfire_g
      lucemcom_g        => ctem_grd%lucemcom_g
      lucltrin_g        => ctem_grd%lucltrin_g
      lucsocin_g        => ctem_grd%lucsocin_g
      emit_co2_g        => ctem_grd%emit_co2_g
      emit_co_g         => ctem_grd%emit_co_g
      emit_ch4_g        => ctem_grd%emit_ch4_g
      emit_nmhc_g       => ctem_grd%emit_nmhc_g
      emit_h2_g         => ctem_grd%emit_h2_g
      emit_nox_g        => ctem_grd%emit_nox_g
      emit_n2o_g        => ctem_grd%emit_n2o_g
      emit_pm25_g       => ctem_grd%emit_pm25_g
      emit_tpm_g        => ctem_grd%emit_tpm_g
      emit_tc_g         => ctem_grd%emit_tc_g
      emit_oc_g         => ctem_grd%emit_oc_g
      emit_bc_g         => ctem_grd%emit_bc_g
      bterm_g           => ctem_grd%bterm_g
      lterm_g           => ctem_grd%lterm_g
      mterm_g           => ctem_grd%mterm_g
      ch4wet1_g         => ctem_grd%ch4wet1_g
      ch4wet2_g         => ctem_grd%ch4wet2_g
      wetfdyn_g         => ctem_grd%wetfdyn_g
      ch4dyn1_g         => ctem_grd%ch4dyn1_g
      ch4dyn2_g         => ctem_grd%ch4dyn2_g
      afrleaf_g         => ctem_grd%afrleaf_g
      afrstem_g         => ctem_grd%afrstem_g
      afrroot_g         => ctem_grd%afrroot_g
      lfstatus_g        => ctem_grd%lfstatus_g
      rmlvegrow_g       => ctem_grd%rmlvegrow_g
      anvegrow_g        => ctem_grd%anvegrow_g
      rmatctem_g        => ctem_grd%rmatctem_g
      HMFGROT_g         => ctem_grd%HMFGROT_g
      HTCROT_g          => ctem_grd%HTCROT_g
      TBARROT_g         => ctem_grd%TBARROT_g
      THLQROT_g         => ctem_grd%THLQROT_g
      THICROT_g         => ctem_grd%THICROT_g
      GFLXROT_g         => ctem_grd%GFLXROT_g

       fsstar_g         => ctem_grd%fsstar_g
       flstar_g         => ctem_grd%flstar_g
       qh_g             => ctem_grd%qh_g
       qe_g             => ctem_grd%qe_g
       snomlt_g         => ctem_grd%snomlt_g
       beg_g            => ctem_grd%beg_g
       gtout_g          => ctem_grd%gtout_g
       tpn_g            => ctem_grd%tpn_g
       altot_g          => ctem_grd%altot_g
       tcn_g            => ctem_grd%tcn_g
       tsn_g            => ctem_grd%tsn_g
       zsn_g            => ctem_grd%zsn_g

      ! mosaic level variables:

      leaflitr_m        => ctem_tile%leaflitr_m
      tltrleaf_m        => ctem_tile%tltrleaf_m
      tltrstem_m        => ctem_tile%tltrstem_m
      tltrroot_m        => ctem_tile%tltrroot_m
      ailcg_m           => ctem_tile%ailcg_m
      ailcb_m           => ctem_tile%ailcb_m
      rmatctem_m        => ctem_tile%rmatctem_m
      veghght_m         => ctem_tile%veghght_m
      rootdpth_m        => ctem_tile%rootdpth_m
      roottemp_m        => ctem_tile%roottemp_m
      slai_m            => ctem_tile%slai_m
      afrroot_m         => ctem_tile%afrroot_m
      afrleaf_m         => ctem_tile%afrleaf_m
      afrstem_m         => ctem_tile%afrstem_m
      laimaxg_m         => ctem_tile%laimaxg_m
      stemmass_m        => ctem_tile%stemmass_m
      rootmass_m        => ctem_tile%rootmass_m
      litrmass_m        => ctem_tile%litrmass_m
      gleafmas_m        => ctem_tile%gleafmas_m
      bleafmas_m        => ctem_tile%bleafmas_m
      soilcmas_m        => ctem_tile%soilcmas_m
      fsnowacc_m        => ctem_tile%fsnowacc_m
      tcansacc_m        => ctem_tile%tcansacc_m
      tcanoaccgat_m     => ctem_tile%tcanoaccgat_m
      taaccgat_m        => ctem_tile%taaccgat_m
      uvaccgat_m        => ctem_tile%uvaccgat_m
      vvaccgat_m        => ctem_tile%vvaccgat_m
      tbaraccgat_m      => ctem_tile%tbaraccgat_m
      tbarcacc_m        => ctem_tile%tbarcacc_m
      tbarcsacc_m       => ctem_tile%tbarcsacc_m
      tbargacc_m        => ctem_tile%tbargacc_m
      tbargsacc_m       => ctem_tile%tbargsacc_m
      thliqcacc_m       => ctem_tile%thliqcacc_m
      thliqgacc_m       => ctem_tile%thliqgacc_m
      thicecacc_m       => ctem_tile%thicecacc_m
      ancsvgac_m        => ctem_tile%ancsvgac_m
      ancgvgac_m        => ctem_tile%ancgvgac_m
      rmlcsvga_m        => ctem_tile%rmlcsvga_m
      rmlcgvga_m        => ctem_tile%rmlcgvga_m
      ifcancmx_m        => ctem_tile%ifcancmx_m


      ! grid level monthly outputs

        laimaxg_mo_g        =>ctem_grd_mo%laimaxg_mo_g
        stemmass_mo_g       =>ctem_grd_mo%stemmass_mo_g
        rootmass_mo_g       =>ctem_grd_mo%rootmass_mo_g
        litrmass_mo_g       =>ctem_grd_mo%litrmass_mo_g
        soilcmas_mo_g       =>ctem_grd_mo%soilcmas_mo_g
        npp_mo_g            =>ctem_grd_mo%npp_mo_g
        gpp_mo_g            =>ctem_grd_mo%gpp_mo_g
        nep_mo_g            =>ctem_grd_mo%nep_mo_g
        nbp_mo_g            =>ctem_grd_mo%nbp_mo_g
        hetrores_mo_g       =>ctem_grd_mo%hetrores_mo_g
        autores_mo_g        =>ctem_grd_mo%autores_mo_g
        litres_mo_g         =>ctem_grd_mo%litres_mo_g
        soilcres_mo_g       =>ctem_grd_mo%soilcres_mo_g
        vgbiomas_mo_g       =>ctem_grd_mo%vgbiomas_mo_g
        totcmass_mo_g       =>ctem_grd_mo%totcmass_mo_g
        emit_co2_mo_g       =>ctem_grd_mo%emit_co2_mo_g
        emit_co_mo_g        =>ctem_grd_mo%emit_co_mo_g
        emit_ch4_mo_g       =>ctem_grd_mo%emit_ch4_mo_g
        emit_nmhc_mo_g      =>ctem_grd_mo%emit_nmhc_mo_g
        emit_h2_mo_g        =>ctem_grd_mo%emit_h2_mo_g
        emit_nox_mo_g       =>ctem_grd_mo%emit_nox_mo_g
        emit_n2o_mo_g       =>ctem_grd_mo%emit_n2o_mo_g
        emit_pm25_mo_g      =>ctem_grd_mo%emit_pm25_mo_g
        emit_tpm_mo_g       =>ctem_grd_mo%emit_tpm_mo_g
        emit_tc_mo_g        =>ctem_grd_mo%emit_tc_mo_g
        emit_oc_mo_g        =>ctem_grd_mo%emit_oc_mo_g
        emit_bc_mo_g        =>ctem_grd_mo%emit_bc_mo_g
        probfire_mo_g       =>ctem_grd_mo%probfire_mo_g
        luc_emc_mo_g        =>ctem_grd_mo%luc_emc_mo_g
        lucltrin_mo_g       =>ctem_grd_mo%lucltrin_mo_g
        lucsocin_mo_g       =>ctem_grd_mo%lucsocin_mo_g
        burnfrac_mo_g       =>ctem_grd_mo%burnfrac_mo_g
        bterm_mo_g          =>ctem_grd_mo%bterm_mo_g
        lterm_mo_g          =>ctem_grd_mo%lterm_mo_g
        mterm_mo_g          =>ctem_grd_mo%mterm_mo_g
        ch4wet1_mo_g        =>ctem_grd_mo%ch4wet1_mo_g
        ch4wet2_mo_g        =>ctem_grd_mo%ch4wet2_mo_g
        wetfdyn_mo_g        =>ctem_grd_mo%wetfdyn_mo_g
        ch4dyn1_mo_g        =>ctem_grd_mo%ch4dyn1_mo_g
        ch4dyn2_mo_g        =>ctem_grd_mo%ch4dyn2_mo_g

      ! mosaic monthly outputs

      laimaxg_mo_m          =>ctem_tile_mo%laimaxg_mo_m
      stemmass_mo_m         =>ctem_tile_mo%stemmass_mo_m
      rootmass_mo_m         =>ctem_tile_mo%rootmass_mo_m
      npp_mo_m              =>ctem_tile_mo%npp_mo_m
      gpp_mo_m              =>ctem_tile_mo%gpp_mo_m
      vgbiomas_mo_m         =>ctem_tile_mo%vgbiomas_mo_m
      autores_mo_m          =>ctem_tile_mo%autores_mo_m
      totcmass_mo_m         =>ctem_tile_mo%totcmass_mo_m
      litrmass_mo_m         =>ctem_tile_mo%litrmass_mo_m
      soilcmas_mo_m         =>ctem_tile_mo%soilcmas_mo_m
      nep_mo_m              =>ctem_tile_mo%nep_mo_m
      litres_mo_m           =>ctem_tile_mo%litres_mo_m
      soilcres_mo_m         =>ctem_tile_mo%soilcres_mo_m
      hetrores_mo_m         =>ctem_tile_mo%hetrores_mo_m
      nbp_mo_m              =>ctem_tile_mo%nbp_mo_m
      emit_co2_mo_m         =>ctem_tile_mo%emit_co2_mo_m
      emit_co_mo_m          =>ctem_tile_mo%emit_co_mo_m
      emit_ch4_mo_m         =>ctem_tile_mo%emit_ch4_mo_m
      emit_nmhc_mo_m        =>ctem_tile_mo%emit_nmhc_mo_m
      emit_h2_mo_m          =>ctem_tile_mo%emit_h2_mo_m
      emit_nox_mo_m         =>ctem_tile_mo%emit_nox_mo_m
      emit_n2o_mo_m         =>ctem_tile_mo%emit_n2o_mo_m
      emit_pm25_mo_m        =>ctem_tile_mo%emit_pm25_mo_m
      emit_tpm_mo_m         =>ctem_tile_mo%emit_tpm_mo_m
      emit_tc_mo_m          =>ctem_tile_mo%emit_tc_mo_m
      emit_oc_mo_m          =>ctem_tile_mo%emit_oc_mo_m
      emit_bc_mo_m          =>ctem_tile_mo%emit_bc_mo_m
      burnfrac_mo_m         =>ctem_tile_mo%burnfrac_mo_m
      probfire_mo_m         =>ctem_tile_mo%probfire_mo_m
      bterm_mo_m            =>ctem_tile_mo%bterm_mo_m
      luc_emc_mo_m          =>ctem_tile_mo%luc_emc_mo_m
      lterm_mo_m            =>ctem_tile_mo%lterm_mo_m
      lucsocin_mo_m         =>ctem_tile_mo%lucsocin_mo_m
      mterm_mo_m            =>ctem_tile_mo%mterm_mo_m
      lucltrin_mo_m         =>ctem_tile_mo%lucltrin_mo_m
      ch4wet1_mo_m          =>ctem_tile_mo%ch4wet1_mo_m
      ch4wet2_mo_m          =>ctem_tile_mo%ch4wet2_mo_m
      wetfdyn_mo_m          =>ctem_tile_mo%wetfdyn_mo_m
      ch4dyn1_mo_m          =>ctem_tile_mo%ch4dyn1_mo_m
      ch4dyn2_mo_m          =>ctem_tile_mo%ch4dyn2_mo_m

      ! grid level annual outputs
      laimaxg_yr_g          =>ctem_grd_yr%laimaxg_yr_g
      stemmass_yr_g         =>ctem_grd_yr%stemmass_yr_g
      rootmass_yr_g         =>ctem_grd_yr%rootmass_yr_g
      litrmass_yr_g         =>ctem_grd_yr%litrmass_yr_g
      soilcmas_yr_g         =>ctem_grd_yr%soilcmas_yr_g
      npp_yr_g              =>ctem_grd_yr%npp_yr_g
      gpp_yr_g              =>ctem_grd_yr%gpp_yr_g
      nep_yr_g              =>ctem_grd_yr%nep_yr_g
      nbp_yr_g              =>ctem_grd_yr%nbp_yr_g
      hetrores_yr_g         =>ctem_grd_yr%hetrores_yr_g
      autores_yr_g          =>ctem_grd_yr%autores_yr_g
      litres_yr_g           =>ctem_grd_yr%litres_yr_g
      soilcres_yr_g         =>ctem_grd_yr%soilcres_yr_g
      vgbiomas_yr_g         =>ctem_grd_yr%vgbiomas_yr_g
      totcmass_yr_g         =>ctem_grd_yr%totcmass_yr_g
      emit_co2_yr_g         =>ctem_grd_yr%emit_co2_yr_g
      emit_co_yr_g          =>ctem_grd_yr%emit_co_yr_g
      emit_ch4_yr_g         =>ctem_grd_yr%emit_ch4_yr_g
      emit_nmhc_yr_g        =>ctem_grd_yr%emit_nmhc_yr_g
      emit_h2_yr_g          =>ctem_grd_yr%emit_h2_yr_g
      emit_nox_yr_g         =>ctem_grd_yr%emit_nox_yr_g
      emit_n2o_yr_g         =>ctem_grd_yr%emit_n2o_yr_g
      emit_pm25_yr_g        =>ctem_grd_yr%emit_pm25_yr_g
      emit_tpm_yr_g         =>ctem_grd_yr%emit_tpm_yr_g
      emit_tc_yr_g          =>ctem_grd_yr%emit_tc_yr_g
      emit_oc_yr_g          =>ctem_grd_yr%emit_oc_yr_g
      emit_bc_yr_g          =>ctem_grd_yr%emit_bc_yr_g
      probfire_yr_g         =>ctem_grd_yr%probfire_yr_g
      luc_emc_yr_g          =>ctem_grd_yr%luc_emc_yr_g
      lucltrin_yr_g         =>ctem_grd_yr%lucltrin_yr_g
      lucsocin_yr_g         =>ctem_grd_yr%lucsocin_yr_g
      burnfrac_yr_g         =>ctem_grd_yr%burnfrac_yr_g
      bterm_yr_g            =>ctem_grd_yr%bterm_yr_g
      lterm_yr_g            =>ctem_grd_yr%lterm_yr_g
      mterm_yr_g            =>ctem_grd_yr%mterm_yr_g
      ch4wet1_yr_g          =>ctem_grd_yr%ch4wet1_yr_g
      ch4wet2_yr_g          =>ctem_grd_yr%ch4wet2_yr_g
      wetfdyn_yr_g          =>ctem_grd_yr%wetfdyn_yr_g
      ch4dyn1_yr_g          =>ctem_grd_yr%ch4dyn1_yr_g
      ch4dyn2_yr_g          =>ctem_grd_yr%ch4dyn2_yr_g

      ! mosaic annual outputs

      laimaxg_yr_m          =>ctem_tile_yr%laimaxg_yr_m
      stemmass_yr_m         =>ctem_tile_yr%stemmass_yr_m
      rootmass_yr_m         =>ctem_tile_yr%rootmass_yr_m
      npp_yr_m              =>ctem_tile_yr%npp_yr_m
      gpp_yr_m              =>ctem_tile_yr%gpp_yr_m
      vgbiomas_yr_m         =>ctem_tile_yr%vgbiomas_yr_m
      autores_yr_m          =>ctem_tile_yr%autores_yr_m
      totcmass_yr_m         =>ctem_tile_yr%totcmass_yr_m
      litrmass_yr_m         =>ctem_tile_yr%litrmass_yr_m
      soilcmas_yr_m         =>ctem_tile_yr%soilcmas_yr_m
      nep_yr_m              =>ctem_tile_yr%nep_yr_m
      litres_yr_m           =>ctem_tile_yr%litres_yr_m
      soilcres_yr_m         =>ctem_tile_yr%soilcres_yr_m
      hetrores_yr_m         =>ctem_tile_yr%hetrores_yr_m
      nbp_yr_m              =>ctem_tile_yr%nbp_yr_m
      emit_co2_yr_m         =>ctem_tile_yr%emit_co2_yr_m
      emit_co_yr_m          =>ctem_tile_yr%emit_co_yr_m
      emit_ch4_yr_m         =>ctem_tile_yr%emit_ch4_yr_m
      emit_nmhc_yr_m        =>ctem_tile_yr%emit_nmhc_yr_m
      emit_h2_yr_m          =>ctem_tile_yr%emit_h2_yr_m
      emit_nox_yr_m         =>ctem_tile_yr%emit_nox_yr_m
      emit_n2o_yr_m         =>ctem_tile_yr%emit_n2o_yr_m
      emit_pm25_yr_m        =>ctem_tile_yr%emit_pm25_yr_m
      emit_tpm_yr_m         =>ctem_tile_yr%emit_tpm_yr_m
      emit_tc_yr_m          =>ctem_tile_yr%emit_tc_yr_m
      emit_oc_yr_m          =>ctem_tile_yr%emit_oc_yr_m
      emit_bc_yr_m          =>ctem_tile_yr%emit_bc_yr_m
      burnfrac_yr_m         =>ctem_tile_yr%burnfrac_yr_m
      probfire_yr_m         =>ctem_tile_yr%probfire_yr_m
      bterm_yr_m            =>ctem_tile_yr%bterm_yr_m
      luc_emc_yr_m          =>ctem_tile_yr%luc_emc_yr_m
      lterm_yr_m            =>ctem_tile_yr%lterm_yr_m
      lucsocin_yr_m         =>ctem_tile_yr%lucsocin_yr_m
      mterm_yr_m            =>ctem_tile_yr%mterm_yr_m
      lucltrin_yr_m         =>ctem_tile_yr%lucltrin_yr_m
      ch4wet1_yr_m          =>ctem_tile_yr%ch4wet1_yr_m
      ch4wet2_yr_m          =>ctem_tile_yr%ch4wet2_yr_m
      wetfdyn_yr_m          =>ctem_tile_yr%wetfdyn_yr_m
      ch4dyn1_yr_m          =>ctem_tile_yr%ch4dyn1_yr_m
      ch4dyn2_yr_m          =>ctem_tile_yr%ch4dyn2_yr_m


c     all model switches are read in from a namelist file
      call read_from_job_options(argbuff,mosaic,transient_run,
     1             trans_startyr,ctemloop,ctem_on,ncyear,lnduseon,
     2             spinfast,cyclemet,nummetcylyrs,metcylyrst,co2on,
     3             setco2conc,popdon,popcycleyr,parallelrun,dofire,
     4             dowetlands,obswetf,compete,inibioclim,start_bare,
     5             rsfile,start_from_rs,jmosty,idisp,izref,islfd,ipcp,
     6             itc,itcg,itg,iwf,ipai,ihgt,ialc,ials,ialg,isnoalb,
     7             igralb,jhhstd,jhhendd,jdstd,jdendd,jhhsty,jhhendy,
     8             jdsty,jdendy)

c     Initialize the CTEM parameters
      call initpftpars(compete)
c
c     set ictemmod, which is the class switch for coupling to ctem
c     either to 1 (ctem is coupled to class) or 0 (class runs alone)
c     this switch is set based on ctem_on that was set by read_from_job_options
c
      if (ctem_on) then
        ictemmod = 1
      else  !ctem_on is false
        ictemmod = 0
      end if
c
      lopcount = 1   ! initialize loop count to 1.
c
c     checking the time spent for running model
c
c      call idate(today)
c      call itime(now)
c      write(*,1000)   today(2), today(1), 2000+today(3), now
c 1000 format( 'start date: ', i2.2, '/', i2.2, '/', i4.4,
c     &      '; start time: ', i2.2, ':', i2.2, ':', i2.2 )
c
C     INITIALIZATION FOR COUPLING CLASS AND CTEM
C
       call initrowvars()
       call resetclassaccum(nltest,nmtest)

       IMONTH = 0

       do 11 i=1,nlat
        barf(i)                = 1.0
        do 11 m=1,nmos
         TCANOACCROW_M(I,M)       = 0.0
         UVACCROW_M(I,M)          = 0.0
         VVACCROW_M(I,M)          = 0.0
         TCANOACCROW_OUT(I,M)     = 0.0
11     continue
c
!     ==================================================================================================

c     do some initializations for the reading in of data from files. these
c     initializations primarily affect how the model does a spinup or transient
c     simulation and which years of the input data are being read.

      if (.not. cyclemet .and. transient_run) then !transient simulation, set to dummy values
        metcylyrst=-9999
        metcycendyr=9999
      else
c       find the final year of the cycling met
c       metcylyrst is defined in the joboptions file
        metcycendyr = metcylyrst + nummetcylyrs - 1
      endif

c     if cycling met (and not doing a transient run), find the popd and luc year to cycle with.
c     it is assumed that you always want to cycle the popd and luc
c     on the same year to be consistent. so if you are cycling the
c     met data, you can set a popd year (first case), or if cycling
c     the met data you can let the popcycleyr default to the met cycling
c     year by setting popcycleyr to -9999 (second case). if not cycling
c     the met data or you are doing a transient run that cycles the MET
c     at the start, cypopyr and cylucyr will default to a dummy value
c     (last case). (See example at bottom of read_from_job_options.f90
c     if confused)
c
      if (cyclemet .and. popcycleyr .ne. -9999 .and.
     &                                .not. transient_run) then
        cypopyr = popcycleyr
        cylucyr = popcycleyr
      else if (cyclemet .and. .not. transient_run) then
        cypopyr = metcylyrst
        cylucyr = metcylyrst
      else  ! give dummy value
        cypopyr = popcycleyr !-9999
        cylucyr = popcycleyr !-9999
      end if

c     ctem initialization done
c
c     open files for reading and writing.
c     these are for coupled model (class_ctem)
c     we added both grid and mosaic output files
c
c     * input files

c         If we wish to restart from the .CTM_RS and .INI_RS files, then
c         we move the original RS files into place and start from them.
          if (start_from_rs) then
             command='mv '//argbuff(1:strlen(argbuff))//'.INI_RS '
     &                    //argbuff(1:strlen(argbuff))//'.INI'
             call system(command)
             command='mv '//argbuff(1:strlen(argbuff))//'.CTM_RS '
     &                    //argbuff(1:strlen(argbuff))//'.CTM'
             call system(command)
          end if


        open(unit=10,file=argbuff(1:strlen(argbuff))//'.INI',
     &       status='old')

        open(unit=12,file=argbuff(1:strlen(argbuff))//'.MET',
     &      status='old')

        open(unit=18,file=argbuff(1:strlen(argbuff))//'.MET2',
     &      status='old')
c     luc file is opened in initialize_luc subroutine

      if (popdon) then
        open(unit=13,file=argbuff(1:strlen(argbuff))//'.POPD',
     &       status='old')
        read(13,*)  !Skip 3 lines of header
        read(13,*)
        read(13,*)
      endif
      if (co2on) then
        open(unit=14,file=argbuff(1:strlen(argbuff))//'.CO2',
     &         status='old')
      endif

c
      if (obswetf) then
        open(unit=16,file=argbuff(1:strlen(argbuff))//'.WET',
     &         status='old')
      endif 

      if (obslght) then ! this was brought in for FireMIP
        open(unit=17,file=argbuff(1:strlen(argbuff))//'.LGHT',
     &         status='old')
      endif
c
c     * CLASS output files
c
      if (.not. parallelrun) then ! stand alone mode, includes half-hourly and daily output
       OPEN(UNIT=61,FILE=ARGBUFF(1:STRLEN(ARGBUFF))//'.OF1_G')  ! GRID-LEVEL DAILY OUTPUT FROM CLASS
       OPEN(UNIT=62,FILE=ARGBUFF(1:STRLEN(ARGBUFF))//'.OF2_G')
       OPEN(UNIT=63,FILE=ARGBUFF(1:STRLEN(ARGBUFF))//'.OF3_G')

       OPEN(UNIT=611,FILE=ARGBUFF(1:STRLEN(ARGBUFF))//'.OF1_M') ! MOSAIC DAILY OUTPUT FROM CLASS
       OPEN(UNIT=621,FILE=ARGBUFF(1:STRLEN(ARGBUFF))//'.OF2_M')
       OPEN(UNIT=631,FILE=ARGBUFF(1:STRLEN(ARGBUFF))//'.OF3_M')

       OPEN(UNIT=64,FILE=ARGBUFF(1:STRLEN(ARGBUFF))//'.OF4_M')  ! MOSAIC HALF-HOURLY OUTPUT FROM CLASS
       OPEN(UNIT=65,FILE=ARGBUFF(1:STRLEN(ARGBUFF))//'.OF5_M')
       OPEN(UNIT=66,FILE=ARGBUFF(1:STRLEN(ARGBUFF))//'.OF6_M')
       OPEN(UNIT=67,FILE=ARGBUFF(1:STRLEN(ARGBUFF))//'.OF7_M')
       OPEN(UNIT=68,FILE=ARGBUFF(1:STRLEN(ARGBUFF))//'.OF8_M')
       OPEN(UNIT=69,FILE=ARGBUFF(1:STRLEN(ARGBUFF))//'.OF9_M')

       OPEN(UNIT=641,FILE=ARGBUFF(1:STRLEN(ARGBUFF))//'.OF4_G') ! GRID-LEVEL HALF-HOURLY OUTPUT FROM CLASS
       OPEN(UNIT=651,FILE=ARGBUFF(1:STRLEN(ARGBUFF))//'.OF5_G')
       OPEN(UNIT=661,FILE=ARGBUFF(1:STRLEN(ARGBUFF))//'.OF6_G')
       OPEN(UNIT=671,FILE=ARGBUFF(1:STRLEN(ARGBUFF))//'.OF7_G')
       OPEN(UNIT=681,FILE=ARGBUFF(1:STRLEN(ARGBUFF))//'.OF8_G')
       OPEN(UNIT=691,FILE=ARGBUFF(1:STRLEN(ARGBUFF))//'.OF9_G')
       end if

C     * READ AND PROCESS INITIALIZATION AND BACKGROUND INFORMATION.
C     * FIRST, MODEL RUN SPECIFICATIONS.

      READ (10,5010) TITLE1,TITLE2,TITLE3,TITLE4,TITLE5,TITLE6
      READ (10,5010) NAME1,NAME2,NAME3,NAME4,NAME5,NAME6
      READ (10,5010) PLACE1,PLACE2,PLACE3,PLACE4,PLACE5,PLACE6
!
!      the ctem output file suffix naming convention is as follows:
!                       ".CT##{time}_{mosaic/grid}"
!      where the ## is a numerical identifier, {time} is any of H, D, M,
!      or Y for half hourly, daily, monthly, or yearly, respectively.
!      after the underscore M or G is used to denote mosaic or grid
!      -averaged values, respectively. also possible is GM for competition
!      outputs since they are the same format in either composite or
!      mosaic modes.
!

       ! Set up the CTEM half-hourly, daily, monthly and yearly files (if all needed), also
       ! setup the CLASS monthly and annual output files:

       call create_outfiles(argbuff,title1, title2, title3, title4,
     1                     title5,title6,name1, name2, name3, name4,
     2                     name5, name6, place1,place2, place3,
     3                     place4, place5, place6)

      IF(CTEM_ON) THEN

        if(obswetf) then
         read(16,*) TITLEC1
        end if
       ENDIF
C
      IF (.NOT. PARALLELRUN) THEN ! STAND ALONE MODE, INCLUDES HALF-HOURLY AND DAILY OUTPUT
C
       WRITE(61,6001) TITLE1,TITLE2,TITLE3,TITLE4,TITLE5,TITLE6
       WRITE(61,6002) NAME1,NAME2,NAME3,NAME4,NAME5,NAME6
       WRITE(61,6011)
6011  FORMAT(2X,'DAY  YEAR  K*  L*  QH  QE  SM  QG  ',
     1          'TR  SWE  DS  WS  AL  ROF  CUMS')

       WRITE(62,6001) TITLE1,TITLE2,TITLE3,TITLE4,TITLE5,TITLE6
       WRITE(62,6002) NAME1,NAME2,NAME3,NAME4,NAME5,NAME6

       IF(IGND.GT.3) THEN
          WRITE(62,6012)
6012      FORMAT(2X,'DAY  YEAR  TG1  THL1  THI1  TG2  THL2  THI2  ',
     1              'TG3  THL3  THI3  TG4  THL4  THI4  TG5  THL5  ',
     2              'THI5')

       ELSE
          WRITE(62,6212)
6212      FORMAT(2X,'DAY  YEAR  TG1  THL1  THI1  TG2  THL2  THI2  ',
     1              'TG3  THL3  THI3  TCN  RCAN  SCAN  TSN  ZSN')

       ENDIF

       WRITE(63,6001) TITLE1,TITLE2,TITLE3,TITLE4,TITLE5,TITLE6
       WRITE(63,6002) NAME1,NAME2,NAME3,NAME4,NAME5,NAME6

       IF(IGND.GT.3) THEN
          WRITE(63,6013)
6013      FORMAT(2X,'DAY  YEAR  TG6  THL6  THI6  TG7  THL7  THI7  ',
     1              'TG8  THL8  THI8  TG9  THL9  THI9  TG10'  ,
     2              'THL10  THI10')

       ELSE
          WRITE(63,6313)
6313      FORMAT(2X,'DAY YEAR KIN LIN TA UV PRES QA PCP EVAP')

       ENDIF
C
       WRITE(64,6001) TITLE1,TITLE2,TITLE3,TITLE4,TITLE5,TITLE6
       WRITE(64,6002) NAME1,NAME2,NAME3,NAME4,NAME5,NAME6
       WRITE(64,6014)
6014  FORMAT(2X,'HOUR  MIN  DAY  YEAR  K*  L*  QH  QE  SM  QG  ',
     1          'TR  SWE  DS  WS  AL  ROF  TPN  ZPN  CDH  CDM  ',
     2          'SFCU  SFCV  UV')

       WRITE(65,6001) TITLE1,TITLE2,TITLE3,TITLE4,TITLE5,TITLE6
       WRITE(65,6002) NAME1,NAME2,NAME3,NAME4,NAME5,NAME6

       IF(IGND.GT.3) THEN
          WRITE(65,6015)
6015      FORMAT(2X,'HOUR  MIN  DAY  YEAR  TG1  THL1  THI1  TG2  ',
     1          'THL2  THI2  TG3  THL3  THI3  TG4  THL4  THI4  ',
     2          'TG5  THL5  THI5')

       ELSE
          WRITE(65,6515)
6515      FORMAT(2X,'HOUR  MIN  DAY  YEAR  TG1  THL1  THI1  TG2  ',
     1           'THL2  THI2  TG3  THL3  THI3  TCN  RCAN  SCAN  ',
     2           'TSN  ZSN  TCN-TA  TCANO  TAC  ACTLYR  FTABLE')

       ENDIF

       WRITE(66,6001) TITLE1,TITLE2,TITLE3,TITLE4,TITLE5,TITLE6
       WRITE(66,6002) NAME1,NAME2,NAME3,NAME4,NAME5,NAME6

       IF(IGND.GT.3) THEN
          WRITE(66,6016)
6016      FORMAT(2X,'HOUR  MIN  DAY  YEAR  TG6  THL6  THI6  TG7  ',
     1          'THL7  THI7  TG8  THL8  THI8  TG9  THL9  THI9  ',
     2          'TG10  THL10  THI10  G0  G1  G2  G3  G4  G5  G6  ',
     3          'G7  G8  G9')

       ELSE
          WRITE(66,6616)
          WRITE(66,6615)
6616  FORMAT(2X,'HOUR  MIN  DAY  SWIN  LWIN  PCP  TA  VA  PA  QA')
6615  FORMAT(2X,'IF IGND <= 3, THIS FILE IS EMPTY')
       ENDIF

       WRITE(67,6001) TITLE1,TITLE2,TITLE3,TITLE4,TITLE5,TITLE6
       WRITE(67,6002) NAME1,NAME2,NAME3,NAME4,NAME5,NAME6
       WRITE(67,6017)
!     6017  FORMAT(2X,'WCAN SCAN CWLCAP CWFCAP FC FG FCS FGS CDH ', !runclass formatted.
!     1          'TCANO TCANS ALBS')
6017  FORMAT(2X,'HOUR  MIN  DAY  YEAR  ',
     1  'TROF     TROO     TROS     TROB      ROF     ROFO   ',
     2  '  ROFS        ROFB         FCS        FGS        FC       FG')

       WRITE(68,6001) TITLE1,TITLE2,TITLE3,TITLE4,TITLE5,TITLE6
       WRITE(68,6002) NAME1,NAME2,NAME3,NAME4,NAME5,NAME6
       WRITE(68,6018)
6018  FORMAT(2X,'HOUR  MIN  DAY  YEAR  ',
     1          'FSGV FSGS FSGG FLGV FLGS FLGG HFSC HFSS HFSG ',
     2          'HEVC HEVS HEVG HMFC HMFS HMFG1 HMFG2 HMFG3 ',
     3          'HTCC HTCS HTC1 HTC2 HTC3')

       WRITE(69,6001) TITLE1,TITLE2,TITLE3,TITLE4,TITLE5,TITLE6
       WRITE(69,6002) NAME1,NAME2,NAME3,NAME4,NAME5,NAME6
       WRITE(69,6019)
6019  FORMAT(2X,'HOUR  MIN  DAY  YEAR  ',
     1   'PCFC PCLC PCPN PCPG QFCF QFCL QFN QFG QFC1 ',
     2          'QFC2 QFC3 ROFC ROFN ROFO ROF WTRC WTRS WTRG')
!       runclass also has: EVDF ','CTV CTS CT1 CT2 CT3')
C
       WRITE(611,6001) TITLE1,TITLE2,TITLE3,TITLE4,TITLE5,TITLE6
       WRITE(611,6002) NAME1,NAME2,NAME3,NAME4,NAME5,NAME6
       WRITE(611,6011)
       WRITE(621,6001) TITLE1,TITLE2,TITLE3,TITLE4,TITLE5,TITLE6
       WRITE(621,6002) NAME1,NAME2,NAME3,NAME4,NAME5,NAME6
C
       IF(IGND.GT.3) THEN
           WRITE(621,6012)
       ELSE
           WRITE(621,6212)
       ENDIF
C
       WRITE(631,6001) TITLE1,TITLE2,TITLE3,TITLE4,TITLE5,TITLE6
       WRITE(631,6002) NAME1,NAME2,NAME3,NAME4,NAME5,NAME6
C
       IF(IGND.GT.3) THEN
          WRITE(631,6013)
       ELSE
          WRITE(631,6313)
       ENDIF
C
       WRITE(641,6001) TITLE1,TITLE2,TITLE3,TITLE4,TITLE5,TITLE6
       WRITE(641,6002) NAME1,NAME2,NAME3,NAME4,NAME5,NAME6
       WRITE(641,6008)
       WRITE(651,6001) TITLE1,TITLE2,TITLE3,TITLE4,TITLE5,TITLE6
       WRITE(651,6002) NAME1,NAME2,NAME3,NAME4,NAME5,NAME6
C
       IF(IGND.GT.3) THEN
          WRITE(651,6015)
       ELSE
          WRITE(651,6515)
       ENDIF
C
       WRITE(661,6001) TITLE1,TITLE2,TITLE3,TITLE4,TITLE5,TITLE6
       WRITE(661,6002) NAME1,NAME2,NAME3,NAME4,NAME5,NAME6
C
       IF(IGND.GT.3) THEN
          WRITE(661,6016)
       ELSE
          WRITE(661,6616)
       ENDIF
C
       WRITE(671,6001) TITLE1,TITLE2,TITLE3,TITLE4,TITLE5,TITLE6
       WRITE(671,6002) NAME1,NAME2,NAME3,NAME4,NAME5,NAME6
       WRITE(671,6017)
       WRITE(681,6001) TITLE1,TITLE2,TITLE3,TITLE4,TITLE5,TITLE6
       WRITE(681,6002) NAME1,NAME2,NAME3,NAME4,NAME5,NAME6
       WRITE(681,6018)
       WRITE(691,6001) TITLE1,TITLE2,TITLE3,TITLE4,TITLE5,TITLE6
       WRITE(691,6002) NAME1,NAME2,NAME3,NAME4,NAME5,NAME6
       WRITE(691,6019)
C
6008  FORMAT(2X,'HOUR  MIN  DAY  YEAR  K*  L*  QH  QE  SM  QG  ',
     1          'TR  SWE  DS  WS  AL  ROF  TPN  ZPN  CDH  CDM  ',
     2          'SFCU  SFCV  UV')

C
       ENDIF !IF NOT PARALLELRUN

C     CTEM FILE TITLES DONE
C======================= CTEM ========================================== /
C
C=======================================================================

C     BEGIN READ IN OF THE .INI FILE

      READ(10,5020)DLATROW(1),DEGLON,ZRFMROW(1),ZRFHROW(1),ZBLDROW(1),
     1              GCROW(1),NLTEST,NMTEST
      JLAT=NINT(DLATROW(1))
      RADJROW(1)=DLATROW(1)*PI/180.
      DLONROW(1)=DEGLON
      Z0ORROW(1)=0.0
      GGEOROW(1)=0.0
C     GGEOROW(1)=-0.035

      DO 50 I=1,NLTEST
      DO 50 M=1,NMTEST
          READ(10,5040) (FCANROT(I,M,J),J=1,ICAN+1),(PAMXROT(I,M,J),
     1                  J=1,ICAN)
          READ(10,5040) (LNZ0ROT(I,M,J),J=1,ICAN+1),(PAMNROT(I,M,J),
     1                  J=1,ICAN)
          READ(10,5040) (ALVCROT(I,M,J),J=1,ICAN+1),(CMASROT(I,M,J),
     1                  J=1,ICAN)
          READ(10,5040) (ALICROT(I,M,J),J=1,ICAN+1),(ROOTROT(I,M,J),
     1                  J=1,ICAN)
          READ(10,5030) (RSMNROT(I,M,J),J=1,ICAN),
     1                  (QA50ROT(I,M,J),J=1,ICAN)
          READ(10,5030) (VPDAROT(I,M,J),J=1,ICAN),
     1                  (VPDBROT(I,M,J),J=1,ICAN)
          READ(10,5030) (PSGAROT(I,M,J),J=1,ICAN),
     1                  (PSGBROT(I,M,J),J=1,ICAN)
          READ(10,5040) DRNROT(I,M),SDEPROT(I,M),FAREROT(I,M)
          READ(10,5090) XSLPROT(I,M),GRKFROT(I,M),WFSFROT(I,M),
     1                  WFCIROT(I,M),MIDROT(I,M)
          READ(10,5080) (SANDROT(I,M,J),J=1,3)
          READ(10,5080) (CLAYROT(I,M,J),J=1,3)
          READ(10,5080) (ORGMROT(I,M,J),J=1,3)
          READ(10,5050) (TBARROT(I,M,J),J=1,3),TCANROT(I,M),
     1                  TSNOROT(I,M),TPNDROT(I,M)
          READ(10,5060) (THLQROT(I,M,J),J=1,3),(THICROT(I,M,J),
     1                  J=1,3),ZPNDROT(I,M)
          READ(10,5070) RCANROT(I,M),SCANROT(I,M),SNOROT(I,M),
     1                  ALBSROT(I,M),RHOSROT(I,M),GROROT(I,M)

          CMASROT(I,M,J)=CMASROT(I,M,J)*10.0
50    CONTINUE

C     ! In CLASS 3.6.2, we not include this soil info in the INI file.
      DO 25 J=1,IGND
          READ(10,*) DELZ(J),ZBOT(J) 
 25   CONTINUE
C
c     the output year ranges can be read in from the job options file or not.
c     if the values should be read in from the .ini file, and not
c     from the job options file, the job options file values are set to
c     -9999 thus triggering the read in of the .ini file values below
      if (jhhstd .eq. -9999) then
        read(10,5200) jhhstd,jhhendd,jdstd,jdendd
       read(10,5200) jhhsty,jhhendy,jdsty,jdendy
      end if
C======================= CTEM ========================================== /

      CLOSE(10)
C
C====================== CTEM =========================================== \
C
c     read from ctem initialization file (.CTM)

      if (ctem_on) then
      call read_from_ctm(nltest,nmtest,FCANROT,FAREROT,
     1                   RSMNROT,QA50ROT,VPDAROT,VPDBROT,PSGAROT,
     2                   PSGBROT,DRNROT,SDEPROT, XSLPROT,GRKFROT,
     3                   WFSFROT,WFCIROT,MIDROT,SANDROT, CLAYROT,
     4                   ORGMROT,TBARROT,THLQROT,THICROT,TCANROT,
     5                   TSNOROT,TPNDROT,ZPNDROT,RCANROT,SCANROT,
     6                   SNOROT, ALBSROT,RHOSROT,GROROT,argbuff)
      end if
c
C===================== CTEM =============================================== /
C
      DO 100 I=1,NLTEST
      DO 100 M=1,NMTEST

          TBARROT(I,M,1)=TBARROT(I,M,1)+TFREZ
          TBARROT(I,M,2)=TBARROT(I,M,2)+TFREZ
          TBARROT(I,M,3)=TBARROT(I,M,3)+TFREZ
          TSNOROT(I,M)=TSNOROT(I,M)+TFREZ
          TCANROT(I,M)=TCANROT(I,M)+TFREZ

          TPNDROT(I,M)=TPNDROT(I,M)+TFREZ
          TBASROT(I,M)=TBARROT(I,M,3)
          CMAIROT(I,M)=0.
          WSNOROT(I,M)=0.
          ZSNLROT(I,M)=0.10
          TSFSROT(I,M,1)=TFREZ
          TSFSROT(I,M,2)=TFREZ

          TSFSROT(I,M,3)=TBARROT(I,M,1)
          TSFSROT(I,M,4)=TBARROT(I,M,1)
          TACROT (I,M)=TCANROT(I,M)
          QACROT (I,M)=0.5E-2

          IF(IGND.GT.3)                                 THEN
              DO 65 J=4,IGND
                  TBARROT(I,M,J)=TBARROT(I,M,3)
                  IF(SDEPROT(I,M).LT.(ZBOT(J-1)+0.001) .AND.
     1                  SANDROT(I,M,3).GT.-2.5)     THEN
                      SANDROT(I,M,J)=-3.0
                      CLAYROT(I,M,J)=-3.0
                      ORGMROT(I,M,J)=-3.0
                      THLQROT(I,M,J)=0.0
                      THICROT(I,M,J)=0.0
                  ELSE
                      SANDROT(I,M,J)=SANDROT(I,M,3)
                      CLAYROT(I,M,J)=CLAYROT(I,M,3)
                      ORGMROT(I,M,J)=ORGMROT(I,M,3)
                      THLQROT(I,M,J)=THLQROT(I,M,3)
                      THICROT(I,M,J)=THICROT(I,M,3)
                  ENDIF
65            CONTINUE
          ENDIF

          DO 75 K=1,6
          DO 75 L=1,50
              ITCTROT(I,M,K,L)=0
75        CONTINUE
100   CONTINUE

      DO 150 I=1,NLTEST
          PREACC(I)=0.
          GTACC(I)=0.
          QEVPACC(I)=0.
          EVAPACC(I)=0.
          HFSACC(I)=0.
          HMFNACC(I)=0.
          ROFACC(I)=0.
          OVRACC(I)=0.
          WTBLACC(I)=0.
          ALVSACC(I)=0.
          ALIRACC(I)=0.
          RHOSACC(I)=0.
          SNOACC(I)=0.
          WSNOACC(I)=0.
          CANARE(I)=0.
          SNOARE(I)=0.
          TSNOACC(I)=0.
          TCANACC(I)=0.
          RCANACC(I)=0.
          SCANACC(I)=0.
          GROACC(I)=0.
          FSINACC(I)=0.
          FLINACC(I)=0.
          FLUTACC(I)=0.
          TAACC(I)=0.
          UVACC(I)=0.
          PRESACC(I)=0.
          QAACC(I)=0.
          DO 125 J=1,IGND
              TBARACC(I,J)=0.
              THLQACC(I,J)=0.
              THICACC(I,J)=0.
              THALACC(I,J)=0.
125       CONTINUE
150   CONTINUE
C
C===================== CTEM =============================================== \
c
c     initialize accumulated array for monthly & yearly output for class
c
      call resetclassmon(nltest)
      call resetclassyr(nltest)

C===================== CTEM =============================================== /

      DO 175 I=1,200
          TAHIST(I)=0.0
          TCHIST(I)=0.0
          TACHIST(I)=0.0
          TDHIST(I)=0.0
          TD2HIST(I)=0.0
          TD3HIST(I)=0.0
          TD4HIST(I)=0.0
          TSHIST(I)=0.0
          TSCRHIST(I)=0.0
175   CONTINUE
      ALAVG=0.0
      ALMAX=0.0
      ACTLYR=0.0
      FTAVG=0.0
      FTMAX=0.0
      FTABLE=0.0

      CALL CLASSB(THPROT,THRROT,THMROT,BIROT,PSISROT,GRKSROT,
     1            THRAROT,HCPSROT,TCSROT,THFCROT,THLWROT,PSIWROT,
     2            DLZWROT,ZBTWROT,ALGWROT,ALGDROT,
     +            ALGWVROT,ALGWNROT,ALGDVROT,ALGDNROT,
     3            SANDROT,CLAYROT,ORGMROT,SOCIROT,DELZ,ZBOT,
     4            SDEPROT,ISNDROT,IGDRROT,
     5            NLAT,NMOS,1,NLTEST,NMTEST,IGND,IGRALB)

5010  FORMAT(2X,6A4)
5020  FORMAT(5F10.2,F7.1,3I5)
5030  FORMAT(4F8.3,8X,4F8.3)
5040  FORMAT(9F8.3)
5050  FORMAT(6F10.2)
5060  FORMAT(7F10.3)
5070  FORMAT(2F10.4,F10.2,F10.3,F10.4,F10.3)
5080  FORMAT(3F10.1)
5090  FORMAT(4E8.1,I8)
5200  FORMAT(4I10)
5300  FORMAT(1X,I2,I3,I5,I6,2F9.2,E14.4,F9.2,E12.3,F8.2,F12.2,3F9.2,
     1       F9.4)
5301  FORMAT(I5,F10.4)
6001  FORMAT('CLASS TEST RUN:',6A4)
      WRITE(*, *)DLATROW(1),DEGLON,ZRFMROW(1),ZRFHROW(1),ZBLDROW(1),
     1              GCROW(1),NLTEST,NMTEST
6002  FORMAT('RESEARCHER:         ',6A4)
6003  FORMAT('INSTITUTION:        ',6A4)
C
C===================== CTEM =============================================== \
C
c     ctem initializations.
c
      if (ctem_on) then
c
c     calculate number of level 2 pfts using modelpft
c
      do 101 j = 1, ican
        isumc = 0
        k1c = (j-1)*l2max + 1
        k2c = k1c + (l2max - 1)
        do n = k1c, k2c
          if(modelpft(n).eq.1) isumc = isumc + 1
        enddo
        nol2pfts(j)=isumc  ! number of level 2 pfts
101   continue
c
      do 110 i=1,nltest
       do 110 m=1,nmtest
        do 111 j = 1, icc
          co2i1csrow(i,m,j)=0.0     !intercellular co2 concentrations
          co2i1cgrow(i,m,j)=0.0
          co2i2csrow(i,m,j)=0.0
          co2i2cgrow(i,m,j)=0.0
          slairow(i,m,j)=0.0        !if bio2str is not called we need to initialize this to zero
111     continue
110   continue
c
      do 123 i =1, ilg
         fsnowacc_m(i)=0.0         !daily accu. fraction of snow
         tcansacc_m(i)=0.0         !daily accu. canopy temp. over snow
         taaccgat_m(i)=0.0
c
         do 128 j = 1, icc
           ancsvgac_m(i,j)=0.0    !daily accu. net photosyn. for canopy over snow subarea
           ancgvgac_m(i,j)=0.0    !daily accu. net photosyn. for canopy over ground subarea
           rmlcsvga_m(i,j)=0.0    !daily accu. leaf respiration for canopy over snow subarea
           rmlcgvga_m(i,j)=0.0    !daily accu. leaf respiration for canopy over ground subarea
           todfrac(i,j)=0.0
128      continue
c
         do 112 j = 1,ignd       !soil temperature and moisture over different subareas
            tbarcacc_m (i,j)=0.0
            tbarcsacc_m(i,j)=0.0
            tbargacc_m (i,j)=0.0
            tbargsacc_m(i,j)=0.0
            thliqcacc_m(i,j)=0.0
            thliqgacc_m(i,j)=0.0
            thicecacc_m(i,j)=0.0
112      continue
123    continue
c
c     find fcancmx with class' fcanmxs and dvdfcans read from ctem's
c     initialization file. this is to divide needle leaf and broad leaf
c     into dcd and evg, and crops and grasses into c3 and c4.
c
      do 113 j = 1, ican
        do 114 i=1,nltest
        do 114 m=1,nmtest
c
          k1c = (j-1)*l2max + 1
          k2c = k1c + (l2max - 1)
c
          do n = k1c, k2c
            if(modelpft(n).eq.1)then
              icountrow(i,m) = icountrow(i,m) + 1
              csum(i,m,j) = csum(i,m,j) +
     &         dvdfcanrow(i,m,icountrow(i,m))

!              Added in seed here to prevent competition from getting
!              pfts with no seed fraction.  JM Feb 20 2014.
              if (compete .and. .not. mosaic) then
               fcancmxrow(i,m,icountrow(i,m))=max(seed,FCANROT(i,m,j)*
     &         dvdfcanrow(i,m,icountrow(i,m)))
               barf(i) = barf(i) - fcancmxrow(i,m,icountrow(i,m))
              else
               fcancmxrow(i,m,icountrow(i,m))=FCANROT(i,m,j)*
     &         dvdfcanrow(i,m,icountrow(i,m))
              end if
            endif
          enddo
c
          if( abs(csum(i,m,j)-1.0).gt.abszero ) then
           write(6,1130)i,m,j
1130       format('dvdfcans for (',i1,',',i1,',',i1,') must add to 1.0')
            call xit('runclass36ctem', -3)
          endif
c
114     continue
113   continue

!     Now make sure that you aren´t over 1.0 for a grid cell (i.e. with a negative
!     bare ground fraction due to the seed fractions being added in.) JM Mar 27 2014
      do i=1,nltest
       if (barf(i) .lt. 0.) then
        bigpftc=maxloc(fcancmxrow(i,:,:))
        ! reduce the most predominant PFT by barf and 1.0e-5,
        ! which ensures that our barefraction is non-zero to avoid
        ! problems later.
        fcancmxrow(i,bigpftc(1),bigpftc(2))=fcancmxrow
     &                (i,bigpftc(1),bigpftc(2))+barf(i) - 1.0e-5
       end if
      end do
c
c     ----------

c     preparation with the input datasets prior to launching run:

      iyear=-99999  ! initialization, forces entry to loop below
      obswetyr=-99999
      obslghtyr=-99999

c     find the first year of met data

       do while (iyear .lt. metcylyrst)
c
        do i=1,nltest  ! formatting was 5300
          read(12,*) ihour,imin,iday,iyear,FSSROW(I),FDLROW(i),
     1         PREROW(i),TAROW(i),QAROW(i),UVROW(i),PRESROW(i)
        enddo
       enddo

c      back up one space in the met file so it is ready for the next readin
       backspace(12)

       if(obswetf) then
         do while (obswetyr .lt. metcylyrst)
            do i=1,nltest
              read(16,*) obswetyr,(wetfrac_mon(i,j),j=1,12)
            end do
         end do
         backspace(16)
       else
           do i=1,nltest
             do j = 1,12
               wetfrac_mon(i,j) = 0.0
             enddo
           enddo

       end if

       if(obslght) then
        do while (obslghtyr .lt. metcylyrst)
            do i=1,nltest
              read(17,*) obslghtyr,(mlightnggrd(i,j),j=1,12)
            end do
         end do
         backspace(17)
       end if

c      If you are not cycling over the MET, you can still specify to end on a
c      year that is shorter than the total climate file length.
       if (.not. cyclemet) endyr = iyear + ncyear

c      find the popd data to cycle over, popd is only cycled over when the met is cycled.
       popyr=-99999  ! initialization, forces entry to loop below

       if (cyclemet .and. popdon) then
        do while (popyr .lt. cypopyr)
         do i = 1, nltest
          read(13,5301) popyr,popdin(i)
         enddo
        enddo
       endif
c
c     if land use change switch is on then read the fractional coverages
c     of ctem's 9 pfts for the first year.
c
      if (lnduseon .and. transient_run) then

         reach_eof=.false.  !flag for when read to end of luc input file

         call initialize_luc(iyear,argbuff,nmtest,nltest,
     1                     mosaic,nol2pfts,cyclemet,
     2                     cylucyr,lucyr,FCANROT,FAREROT,nfcancmxrow,
     3                     pfcancmxrow,fcancmxrow,reach_eof,start_bare,
     4                     compete)

         if (reach_eof) goto 999

      endif ! if (lnduseon)
c
c     with fcancmx calculated above and initialized values of all ctem pools,
c     find mosaic tile (grid) average vegetation biomass, litter mass, and soil c mass.
c     also initialize additional variables which are used by ctem.
c
      do 115 i = 1,nltest
        do 115 m = 1,nmtest
          vgbiomasrow(i,m)=0.0
          gavglairow(i,m)=0.0
          gavgltmsrow(i,m)=0.0
          gavgscmsrow(i,m)=0.0
          lucemcomrow(i,m)=0.0      !land use change combustion emission losses
          lucltrinrow(i,m)=0.0      !land use change inputs to litter pool
          lucsocinrow(i,m)=0.0      !land use change inputs to soil c pool
          colddaysrow(i,m,1)=0      !cold days counter for ndl dcd
          colddaysrow(i,m,2)=0      !cold days counter for crops

          do 116 j = 1, icc
            vgbiomasrow(i,m)=vgbiomasrow(i,m)+fcancmxrow(i,m,j)*
     &        (gleafmasrow(i,m,j)+stemmassrow(i,m,j)+
     &         rootmassrow(i,m,j)+bleafmasrow(i,m,j))
            gavgltmsrow(i,m)=gavgltmsrow(i,m)+fcancmxrow(i,m,j)*
     &                       litrmassrow(i,m,j)
            gavgscmsrow(i,m)=gavgscmsrow(i,m)+fcancmxrow(i,m,j)*
     &         soilcmasrow(i,m,j)
            grwtheffrow(i,m,j)=100.0   !set growth efficiency to some large number
c                                      !so that no growth related mortality occurs in
c                                      !first year
            flhrlossrow(i,m,j)=0.0     !fall/harvest loss
            stmhrlosrow(i,m,j)=0.0     !stem harvest loss for crops
            rothrlosrow(i,m,j)=0.0     !root death for crops
            lystmmasrow(i,m,j)=stemmassrow(i,m,j)
            lyrotmasrow(i,m,j)=rootmassrow(i,m,j)
            tymaxlairow(i,m,j)=0.0

116      continue

c
c *     initialize accumulated array for monthly and yearly output for ctem
c

         call resetmonthend_m(nltest,nmtest)
         call resetyearend_m(nltest,nmtest)
c
115   continue
c
      do 117 i = 1,nltest
        do 117 m = 1,nmtest
         gavgltmsrow(i,m)=gavgltmsrow(i,m)+ (1.0-FCANROT(i,m,1)-
     &       FCANROT(i,m,2)-
     &    FCANROT(i,m,3)-FCANROT(i,m,4))*litrmassrow(i,m,icc+1)
         gavgscmsrow(i,m)=gavgscmsrow(i,m)+ (1.0-FCANROT(i,m,1)-
     &   FCANROT(i,m,2)-
     &    FCANROT(i,m,3)-FCANROT(i,m,4))*soilcmasrow(i,m,icc+1)
c
117   continue
c

      CALL GATPREP(ILMOS,JLMOS,IWMOS,JWMOS,
     1             NML,NMW,GCROW,FAREROT,MIDROT,
     2             NLAT,NMOS,ILG,1,NLTEST,NMTEST)
c
      call ctemg1(gleafmasgat,bleafmasgat,stemmassgat,rootmassgat,
     1      fcancmxgat,zbtwgat,dlzwgat,sdepgat,ailcggat,ailcbgat,
     2      ailcgat,zolncgat,rmatcgat,rmatctemgat,slaigat,
     3      bmasveggat,cmasvegcgat,veghghtgat,
     4      rootdpthgat,alvsctmgat,alirctmgat,
     5      paicgat,    slaicgat,
     6      ilmos,jlmos,iwmos,jwmos,
     7      nml,
     8      gleafmasrow,bleafmasrow,stemmassrow,rootmassrow,
     9      fcancmxrow,ZBTWROT,DLZWROT,SDEPROT,ailcgrow,ailcbrow,
     a      ailcrow,zolncrow,rmatcrow,rmatctemrow,slairow,
     b      bmasvegrow,cmasvegcrow,veghghtrow,
     c      rootdpthrow,alvsctmrow,alirctmrow,
     d      paicrow,    slaicrow)
c
c
      call bio2str( gleafmasgat,bleafmasgat,stemmassgat,rootmassgat,
     1                           1,      nml,    fcancmxgat, zbtwgat,
     2                        dlzwgat, nol2pfts,   sdepgat,
     4                       ailcggat, ailcbgat,  ailcgat, zolncgat,
     5                       rmatcgat, rmatctemgat,slaigat,bmasveggat,
     6                 cmasvegcgat,veghghtgat, rootdpthgat,alvsctmgat,
     7                     alirctmgat, paicgat,  slaicgat )
c
      call ctems1(gleafmasrow,bleafmasrow,stemmassrow,rootmassrow,
     1      fcancmxrow,ZBTWROT,DLZWROT,SDEPROT,ailcgrow,ailcbrow,
     2      ailcrow,zolncrow,rmatcrow,rmatctemrow,slairow,
     3      bmasvegrow,cmasvegcrow,veghghtrow,
     4      rootdpthrow,alvsctmrow,alirctmrow,
     5      paicrow,    slaicrow,
     6      ilmos,jlmos,iwmos,jwmos,
     7      nml,
     8      gleafmasgat,bleafmasgat,stemmassgat,rootmassgat,
     9      fcancmxgat,zbtwgat,dlzwgat,sdepgat,ailcggat,ailcbgat,
     a      ailcgat,zolncgat,rmatcgat,rmatctemgat,slaigat,
     b      bmasveggat,cmasvegcgat,veghghtgat,
     c      rootdpthgat,alvsctmgat,alirctmgat,
     d      paicgat,    slaicgat)
c
      endif   ! if (ctem_on)
c
c     ctem initial preparation done

C===================== CTEM ============================================ /
C
C     **** LAUNCH RUN. ****

      N=0
      NCOUNT=1
      NDAY=86400/NINT(DELT)

C===================== CTEM ============================================ \

      run_model=.true.
      met_rewound=.false.

200   continue

c     start up the main model loop

      do while (run_model)

c     if the met file has been rewound (due to cycling over the met data)
c     then we need to find the proper year in the file before we continue
c     on with the run
      if (met_rewound) then
        do while (iyear .lt. metcylyrst)
         do i=1,nltest
c         this reads in one 30 min slice of met data, when it reaches
c         the end of file it will go to label 999.  !formatting was 5300
          read(12,*,end=999) ihour,imin,iday,iyear,FSSROW(I),FDLROW(i),
     1         PREROW(i),TAROW(i),QAROW(i),UVROW(i),PRESROW(i)
         enddo
        enddo

c       back up one space in the met file so it is ready for the next readin
c       but only if it was read in during the loop above.
        if (metcylyrst .ne. -9999) backspace(12)

      if (ctem_on) then
        if (obswetf) then
          do while (obswetyr .lt. metcylyrst)
              do i = 1,nltest
                read(16,*) obswetyr,(wetfrac_mon(i,j),j=1,12)
              enddo
          enddo
         if (metcylyrst .ne. -9999) backspace(16)
        else
           do i=1,nltest
             do j = 1,12
               wetfrac_mon(i,j) = 0.0
             enddo
           enddo
        endif !obswetf

       if(obslght) then
        do while (obslghtyr .lt. metcylyrst)
            do i=1,nltest
              read(17,*) obslghtyr,(mlightnggrd(i,j),j=1,12)
            end do
         end do
         if (metcylyrst .ne. -9999) backspace(17)
       end if

       endif ! ctem_on 

      met_rewound = .false.

      endif

C===================== CTEM ============================================ /
C
C========================================================================
C     * READ IN METEOROLOGICAL FORCING DATA FOR CURRENT TIME STEP;
C     * CALCULATE SOLAR ZENITH ANGLE AND COMPONENTS OF INCOMING SHORT-
C     * WAVE RADIATION FLUX; ESTIMATE FLUX PARTITIONS IF NECESSARY.
C
      N=N+1

      DO 250 I=1,NLTEST
C         THIS READS IN ONE 30 MIN SLICE OF MET DATA, WHEN IT REACHES
C         THE END OF FILE IT WILL GO TO 999. !formatting was 5300
          READ(12,*,END=999) IHOUR,IMIN,IDAY,IYEAR,FSSROW(I),FDLROW(I),
     1         PREROW(I),TAROW(I),QAROW(I),UVROW(I),PRESROW(I)

C===================== CTEM ============================================ \
          if (iday.eq.1.and.ihour.eq.0.and.imin.eq.0) then
            if (ctem_on) then
              if (obswetf) then
                  read(16,*,end=1001) obswetyr,(wetfrac_mon(i,j),j=1,12)
              else
                   do j = 1,12
                     wetfrac_mon(i,j) = 0.0
                   enddo
              endif !obswetf

              if(obslght) then
                read(17,*,end=212) obslghtyr,(mlightnggrd(i,j),j=1,12)
212       continue !if end of file, just keep using the last year of lighting data.
              end if
            endif ! ctem_on 
 
          endif 

c         assign the met climate year to climiyear
          climiyear = iyear

!         If in a transient_run that has to cycle over MET then change
!         the iyear here:
          if (transient_run .and. cyclemet) then
            iyear = iyear - (metcylyrst - trans_startyr)
          end if
c
          if(lopcount .gt. 1) then
            if (cyclemet) then
              iyear=iyear + nummetcylyrs*(lopcount-1)
            else
              iyear=iyear + ncyear*(lopcount-1)
            end if
          endif   ! lopcount .gt. 1

c
c         write(*,*)'year=',iyear,'day=',iday,' hour=',ihour,' min=',imin
c
C===================== CTEM ============================================ /

          FSVHROW(I)=0.5*FSSROW(I)
          FSIHROW(I)=0.5*FSSROW(I)
          TAROW(I)=TAROW(I)+TFREZ
          ULROW(I)=UVROW(I)
          VLROW(I)=0.0
          VMODROW(I)=UVROW(I)
250   CONTINUE

      DO 251 I=1,NLTEST
C         THIS READS IN ONE 30 MIN SLICE OF MET DATA, WHEN IT REACHES
C         THE END OF FILE IT WILL GO TO 999. !formatting was 5300
          READ(18,*,END=999) IHOUR2,IMIN2,IDAY2,IYEAR2,FSSROW2(I),FDLROW2(I),
     1         PREROW2(I),TAROW2(I),QAROW2(I),UVROW2(I),PRESROW2(I)
C
251   CONTINUE

      DAY=REAL(IDAY)+(REAL(IHOUR)+REAL(IMIN)/60.)/24.
      DECL=SIN(2.*PI*(284.+DAY)/365.)*23.45*PI/180.
      HOUR=(REAL(IHOUR)+REAL(IMIN)/60.)*PI/12.-PI
      COSZ=SIN(RADJROW(1))*SIN(DECL)+COS(RADJROW(1))*COS(DECL)*COS(HOUR)

      DO 300 I=1,NLTEST
          CSZROW(I)=SIGN(MAX(ABS(COSZ),1.0E-3),COSZ)
          IF(PREROW(I).GT.0.) THEN
              XDIFFUS(I)=1.0
          ELSE
              XDIFFUS(I)=MAX(0.0,MIN(1.0-0.9*COSZ,1.0))
          ENDIF
          FCLOROW(I)=XDIFFUS(I)
300   CONTINUE
C
C===================== CTEM ============================================ \
C
      if (iday.eq.1.and.ihour.eq.0.and.imin.eq.0) then
c
c      if popdon=true
c      calculate fire extinguishing probability and
c      probability of fire due to human causes
c      from population density input data. In disturb.f90 this will
c      overwrite extnprobgrd(i) and prbfrhucgrd(i) that are
c      read in from the .ctm file. Set
c      cypopyr = -9999 when we don't want to cycle over the popd data
c      so this allows us to grab a new value each year.

       if(popdon .and. transient_run) then
         do while (popyr .lt. iyear) 
          do i=1,nltest
           read(13,5301,end=999) popyr,popdin(i)
          enddo
         enddo
       endif
c
c      if co2on is true
c      read co2concin from input datafile and
c      overwrite co2concrow, otherwise set to constant value
c
       if(co2on) then

        do while (co2yr .lt. iyear)
          do i=1,nltest
           read(14,*,end=999) co2yr,co2concin
           do m=1,nmtest
            co2concrow(i,m)=co2concin
           enddo !nmtest
          enddo !nltest
        enddo !co2yr < iyear

       else !constant co2

         do i=1,nltest
          do m=1,nmtest
           co2concrow(i,m)=setco2conc
          enddo
         enddo

       endif !co2on

c      if lnduseon is true, read in the luc data now

       if (ctem_on .and. lnduseon .and. transient_run) then

         call readin_luc(iyear,nmtest,nltest,mosaic,lucyr,
     &                   nfcancmxrow,pfcancmxrow,reach_eof,compete)
         if (reach_eof) goto 999

       else ! lnduseon = false or met is cycling in a spin up run

c          land use is not on or the met data is being cycled, so the
c          pfcancmx value is also the nfcancmx value.
c
           nfcancmxrow=pfcancmxrow

       endif ! lnduseon/cyclemet
c
      endif   ! at the first day of each year i.e.
c             ! if (iday.eq.1.and.ihour.eq.0.and.imin.eq.0)
c

C===================== CTEM ============================================ /
C
      CALL CLASSI(VPDROW,TADPROW,PADRROW,RHOAROW,RHSIROW,
     1            RPCPROW,TRPCROW,SPCPROW,TSPCROW,TAROW,QAROW,
     2            PREROW,RPREROW,SPREROW,PRESROW,
     3            IPCP,NLAT,1,NLTEST)

C
      CUMSNO=CUMSNO+SPCPROW(1)*RHSIROW(1)*DELT
C
      CALL GATPREP(ILMOS,JLMOS,IWMOS,JWMOS,
     1             NML,NMW,GCROW,FAREROT,MIDROT,
     2             NLAT,NMOS,ILG,1,NLTEST,NMTEST)
C
      CALL CLASSG (TBARGAT,THLQGAT,THICGAT,TPNDGAT,ZPNDGAT,
     1             TBASGAT,ALBSGAT,TSNOGAT,RHOSGAT,SNOGAT,
     2             TCANGAT,RCANGAT,SCANGAT,GROGAT, CMAIGAT,
     3             FCANGAT,LNZ0GAT,ALVCGAT,ALICGAT,PAMXGAT,
     4             PAMNGAT,CMASGAT,ROOTGAT,RSMNGAT,QA50GAT,
     5             VPDAGAT,VPDBGAT,PSGAGAT,PSGBGAT,PAIDGAT,
     6             HGTDGAT,ACVDGAT,ACIDGAT,TSFSGAT,WSNOGAT,
     7             THPGAT, THRGAT, THMGAT, BIGAT,  PSISGAT,
     8             GRKSGAT,THRAGAT,HCPSGAT,TCSGAT, IGDRGAT,
     9             THFCGAT,THLWGAT,PSIWGAT,DLZWGAT,ZBTWGAT,
     A             VMODGAT,ZSNLGAT,ZPLGGAT,ZPLSGAT,TACGAT,
     B             QACGAT,DRNGAT, XSLPGAT,GRKFGAT,WFSFGAT,
     C             WFCIGAT,ALGWVGAT,ALGWNGAT,ALGDVGAT,ALGDNGAT,
     +             ALGWGAT,ALGDGAT,ASVDGAT,ASIDGAT,AGVDGAT,
     D             AGIDGAT,ISNDGAT,RADJGAT,ZBLDGAT,Z0ORGAT,
     E             ZRFMGAT,ZRFHGAT,ZDMGAT, ZDHGAT, FSVHGAT,
     F             FSIHGAT,FSDBGAT,FSFBGAT,FSSBGAT,CSZGAT,
     +             FSGGAT, FLGGAT, FDLGAT, ULGAT,  VLGAT,
     G             TAGAT,  QAGAT,  PRESGAT,PREGAT, PADRGAT,
     H             VPDGAT, TADPGAT,RHOAGAT,RPCPGAT,TRPCGAT,
     I             SPCPGAT,TSPCGAT,RHSIGAT,FCLOGAT,DLONGAT,
     J             GGEOGAT,GUSTGAT,REFGAT, BCSNGAT,DEPBGAT,
     K             ILMOS,JLMOS,
     L             NML,NLAT,NTLD,NMOS,ILG,IGND,ICAN,ICAN+1,NBS,
     M             TBARROT,THLQROT,THICROT,TPNDROT,ZPNDROT,
     N             TBASROT,ALBSROT,TSNOROT,RHOSROT,SNOROT,
     O             TCANROT,RCANROT,SCANROT,GROROT, CMAIROT,
     P             FCANROT,LNZ0ROT,ALVCROT,ALICROT,PAMXROT,
     Q             PAMNROT,CMASROT,ROOTROT,RSMNROT,QA50ROT,
     R             VPDAROT,VPDBROT,PSGAROT,PSGBROT,PAIDROT,
     S             HGTDROT,ACVDROT,ACIDROT,TSFSROT,WSNOROT,
     T             THPROT, THRROT, THMROT, BIROT,  PSISROT,
     U             GRKSROT,THRAROT,HCPSROT,TCSROT, IGDRROT,
     V             THFCROT,THLWROT,PSIWROT,DLZWROT,ZBTWROT,
     W             VMODROW,ZSNLROT,ZPLGROT,ZPLSROT,TACROT,
     X             QACROT,DRNROT, XSLPROT,GRKFROT,WFSFROT,
     Y             WFCIROT,ALGWVROT,ALGWNROT,ALGDVROT,ALGDNROT,
     +             ALGWROT,ALGDROT,ASVDROT,ASIDROT,AGVDROT,
     Z             AGIDROT,ISNDROT,RADJROW,ZBLDROW,Z0ORROW,
     +             ZRFMROW,ZRFHROW,ZDMROW, ZDHROW, FSVHROW,
     +             FSIHROW,FSDBROL,FSFBROL,FSSBROL,CSZROW,
     +             FSGROL, FLGROL, FDLROW, ULROW,  VLROW,
     +             TAROW,  QAROW,  PRESROW,PREROW, PADRROW,
     +             VPDROW, TADPROW,RHOAROW,RPCPROW,TRPCROW,
     +             SPCPROW,TSPCROW,RHSIROW,FCLOROW,DLONROW,
     +             GGEOROW,GUSTROL,REFROT, BCSNROT,DEPBROW )
C
C    * INITIALIZATION OF DIAGNOSTIC VARIABLES SPLIT OUT OF CLASSG
C    * FOR CONSISTENCY WITH GCM APPLICATIONS.
C
      DO 330 K=1,ILG
          CDHGAT (K)=0.0
          CDMGAT (K)=0.0
          HFSGAT (K)=0.0
          TFXGAT (K)=0.0
          QEVPGAT(K)=0.0
          QFSGAT (K)=0.0
          QFXGAT (K)=0.0
          PETGAT (K)=0.0
          GAGAT  (K)=0.0
          EFGAT  (K)=0.0
          GTGAT  (K)=0.0
          QGGAT  (K)=0.0
          ALVSGAT(K)=0.0
          ALIRGAT(K)=0.0
          SFCTGAT(K)=0.0
          SFCUGAT(K)=0.0
          SFCVGAT(K)=0.0
          SFCQGAT(K)=0.0
          FSNOGAT(K)=0.0
          FSGVGAT(K)=0.0
          FSGSGAT(K)=0.0
          FSGGGAT(K)=0.0
          FLGVGAT(K)=0.0
          FLGSGAT(K)=0.0
          FLGGGAT(K)=0.0
          HFSCGAT(K)=0.0
          HFSSGAT(K)=0.0
          HFSGGAT(K)=0.0
          HEVCGAT(K)=0.0
          HEVSGAT(K)=0.0
          HEVGGAT(K)=0.0
          HMFCGAT(K)=0.0
          HMFNGAT(K)=0.0
          HTCCGAT(K)=0.0
          HTCSGAT(K)=0.0
          PCFCGAT(K)=0.0
          PCLCGAT(K)=0.0
          PCPNGAT(K)=0.0
          PCPGGAT(K)=0.0
          QFGGAT (K)=0.0
          QFNGAT (K)=0.0
          QFCFGAT(K)=0.0
          QFCLGAT(K)=0.0
          ROFGAT (K)=0.0
          ROFOGAT(K)=0.0
          ROFSGAT(K)=0.0
          ROFBGAT(K)=0.0
          TROFGAT(K)=0.0
          TROOGAT(K)=0.0
          TROSGAT(K)=0.0
          TROBGAT(K)=0.0
          ROFCGAT(K)=0.0
          ROFNGAT(K)=0.0
          ROVGGAT(K)=0.0
          WTRCGAT(K)=0.0
          WTRSGAT(K)=0.0
          WTRGGAT(K)=0.0
          DRGAT  (K)=0.0
330   CONTINUE
C
      DO 334 L=1,IGND
      DO 332 K=1,ILG
          HMFGGAT(K,L)=0.0
          HTCGAT (K,L)=0.0
          QFCGAT (K,L)=0.0
          GFLXGAT(K,L)=0.0
332   CONTINUE
334   CONTINUE
C
      DO 340 M=1,50
          DO 338 L=1,6
              DO 336 K=1,NML
                  ITCTGAT(K,L,M)=0
336           CONTINUE
338       CONTINUE
340   CONTINUE
C
C========================================================================
C
      CALL CLASSZ (0,      CTVSTP, CTSSTP, CT1STP, CT2STP, CT3STP,
     1             WTVSTP, WTSSTP, WTGSTP,
     2             FSGVGAT,FLGVGAT,HFSCGAT,HEVCGAT,HMFCGAT,HTCCGAT,
     3             FSGSGAT,FLGSGAT,HFSSGAT,HEVSGAT,HMFNGAT,HTCSGAT,
     4             FSGGGAT,FLGGGAT,HFSGGAT,HEVGGAT,HMFGGAT,HTCGAT,
     5             PCFCGAT,PCLCGAT,QFCFGAT,QFCLGAT,ROFCGAT,WTRCGAT,
     6             PCPNGAT,QFNGAT, ROFNGAT,WTRSGAT,PCPGGAT,QFGGAT,
     7             QFCGAT, ROFGAT, WTRGGAT,CMAIGAT,RCANGAT,SCANGAT,
     8             TCANGAT,SNOGAT, WSNOGAT,TSNOGAT,THLQGAT,THICGAT,
     9             HCPSGAT,THPGAT, DLZWGAT,TBARGAT,ZPNDGAT,TPNDGAT,
     A             DELZ,   FCS,    FGS,    FC,     FG,
     B             1,      NML,    ILG,    IGND,   N    )
C
C========================================================================
C
C===================== CTEM ============================================ \
C
      call ctemg2(fcancmxgat,rmatcgat,zolncgat,paicgat,
     1      ailcgat,     ailcggat,    cmasvegcgat,  slaicgat,
     2      ailcgsgat,   fcancsgat,   fcancgat,     rmatctemgat,
     3      co2concgat,  co2i1cggat,  co2i1csgat,   co2i2cggat,
     4      co2i2csgat,  xdiffusgat,  slaigat,      cfluxcggat,
     5      cfluxcsgat,  ancsveggat,  ancgveggat,   rmlcsveggat,
     6      rmlcgveggat, canresgat,   sdepgat,
     7      sandgat,     claygat,     orgmgat,
     8      anveggat,    rmlveggat,   tcanoaccgat_m,tbaraccgat_m,
     9      uvaccgat_m,  vvaccgat_m,  mlightnggat,  prbfrhucgat,
     a      extnprobgat, stdalngat,   pfcancmxgat,  nfcancmxgat,
     b      stemmassgat, rootmassgat, litrmassgat,  gleafmasgat,
     c      bleafmasgat, soilcmasgat, ailcbgat,     flhrlossgat,
     d      pandaysgat,  lfstatusgat, grwtheffgat,  lystmmasgat,
     e      lyrotmasgat, tymaxlaigat, vgbiomasgat,  gavgltmsgat,
     f      stmhrlosgat, bmasveggat,  colddaysgat,  rothrlosgat,
     g      alvsctmgat,  alirctmgat,  gavglaigat,   nppgat,
     h      nepgat,      hetroresgat, autoresgat,   soilcrespgat,
     i      rmgat,       rggat,       nbpgat,       litresgat,
     j      socresgat,   gppgat,      dstcemlsgat,  litrfallgat,
     k      humiftrsgat, veghghtgat,  rootdpthgat,  rmlgat,
     l      rmsgat,      rmrgat,      tltrleafgat,  tltrstemgat,
     m      tltrrootgat, leaflitrgat, roottempgat,  afrleafgat,
     n      afrstemgat,  afrrootgat,  wtstatusgat,  ltstatusgat,
     o      burnfracgat, probfiregat, lucemcomgat,  lucltringat,
     p      lucsocingat, nppveggat,   dstcemls3gat,
     q      faregat,     gavgscmsgat, rmlvegaccgat, pftexistgat,
     &      rmsveggat,   rmrveggat,   rgveggat,    vgbiomas_veggat,
     &      gppveggat,   nepveggat,   ailcmingat,   ailcmaxgat,
     &      emit_co2gat,  emit_cogat, emit_ch4gat,  emit_nmhcgat,
     &      emit_h2gat,   emit_noxgat,emit_n2ogat,  emit_pm25gat,
     &      emit_tpmgat,  emit_tcgat, emit_ocgat,   emit_bcgat,
     &      btermgat,     ltermgat,   mtermgat,
     &      nbpveggat,    hetroresveggat, autoresveggat,litresveggat,
     &      soilcresveggat, burnvegfgat, pstemmassgat, pgleafmassgat,
     &      CH4WET1GAT, CH4WET2GAT,
     &      WETFDYNGAT, CH4DYN1GAT,  CH4DYN2GAT,
c
     r      ilmos,       jlmos,       iwmos,        jwmos,
     s      nml,      fcancmxrow,  rmatcrow,    zolncrow,  paicrow,
     v      ailcrow,     ailcgrow,    cmasvegcrow,  slaicrow,
     w      ailcgsrow,   fcancsrow,   fcancrow,     rmatctemrow,
     x      co2concrow,  co2i1cgrow,  co2i1csrow,   co2i2cgrow,
     y      co2i2csrow,  xdiffus,     slairow,      cfluxcgrow,
     z      cfluxcsrow,  ancsvegrow,  ancgvegrow,   rmlcsvegrow,
     1      rmlcgvegrow, canresrow,   SDEPROT,
     2      SANDROT,     CLAYROT,     ORGMROT,
     3      anvegrow,    rmlvegrow,   tcanoaccrow_m,tbaraccrow_m,
     4      uvaccrow_m,  vvaccrow_m,  mlightnggrd,  prbfrhucgrd,
     5      extnprobgrd, stdalngrd,   pfcancmxrow,  nfcancmxrow,
     6      stemmassrow, rootmassrow, litrmassrow,  gleafmasrow,
     7      bleafmasrow, soilcmasrow, ailcbrow,     flhrlossrow,
     8      pandaysrow,  lfstatusrow, grwtheffrow,  lystmmasrow,
     9      lyrotmasrow, tymaxlairow, vgbiomasrow,  gavgltmsrow,
     a      stmhrlosrow, bmasvegrow,  colddaysrow,  rothrlosrow,
     b      alvsctmrow,  alirctmrow,  gavglairow,   npprow,
     c      neprow,      hetroresrow, autoresrow,   soilcresprow,
     d      rmrow,       rgrow,       nbprow,       litresrow,
     e      socresrow,   gpprow,      dstcemlsrow,  litrfallrow,
     f      humiftrsrow, veghghtrow,  rootdpthrow,  rmlrow,
     g      rmsrow,      rmrrow,      tltrleafrow,  tltrstemrow,
     h      tltrrootrow, leaflitrrow, roottemprow,  afrleafrow,
     i      afrstemrow,  afrrootrow,  wtstatusrow,  ltstatusrow,
     j      burnfracrow, probfirerow, lucemcomrow,  lucltrinrow,
     k      lucsocinrow, nppvegrow,   dstcemls3row,
     l      FAREROT,     gavgscmsrow, rmlvegaccrow, pftexistrow,
     &      rmsvegrow,   rmrvegrow,   rgvegrow,    vgbiomas_vegrow,
     &      gppvegrow,   nepvegrow,   ailcminrow,   ailcmaxrow,
     &      emit_co2row,  emit_corow, emit_ch4row,  emit_nmhcrow,
     &      emit_h2row,   emit_noxrow,emit_n2orow,  emit_pm25row,
     &      emit_tpmrow,  emit_tcrow, emit_ocrow,   emit_bcrow,
     &      btermrow,     ltermrow,   mtermrow,
     &      nbpvegrow,    hetroresvegrow, autoresvegrow,litresvegrow,
     &      soilcresvegrow, burnvegfrow, pstemmassrow, pgleafmassrow,
     &      CH4WET1ROW, CH4WET2ROW,
     &      WETFDYNROW, CH4DYN1ROW, CH4DYN2ROW)
c
C===================== CTEM ============================================ /
C
C-----------------------------------------------------------------------
C     * ALBEDO AND TRANSMISSIVITY CALCULATIONS; GENERAL VEGETATION
C     * CHARACTERISTICS.
C     * ADAPTED TO COUPLING OF CLASS3.6 AND CTEM by including zolnc,
!     * cmasvegc, alvsctm, alirctm in the arguments.
C
      CALL CLASSA    (FC,     FG,     FCS,    FGS,    ALVSCN, ALIRCN,
     1                ALVSG,  ALIRG,  ALVSCS, ALIRCS, ALVSSN, ALIRSN,
     2                ALVSGC, ALIRGC, ALVSSC, ALIRSC, TRVSCN, TRIRCN,
     3                TRVSCS, TRIRCS, FSVF,   FSVFS,
     4                RAICAN, RAICNS, SNOCAN, SNOCNS, FRAINC, FSNOWC,
     5                FRAICS, FSNOCS, DISP,   DISPS,  ZOMLNC, ZOMLCS,
     6                ZOELNC, ZOELCS, ZOMLNG, ZOMLNS, ZOELNG, ZOELNS,
     7                CHCAP,  CHCAPS, CMASSC, CMASCS, CWLCAP, CWFCAP,
     8                CWLCPS, CWFCPS, RC,     RCS,    RBCOEF, FROOT,
     9                FROOTS, ZPLIMC, ZPLIMG, ZPLMCS, ZPLMGS, ZSNOW,
     A                WSNOGAT,ALVSGAT,ALIRGAT,HTCCGAT,HTCSGAT,HTCGAT,
     +                ALTG,   ALSNO,  TRSNOWC,TRSNOWG,
     B                WTRCGAT,WTRSGAT,WTRGGAT,CMAIGAT,FSNOGAT,
     C                FCANGAT,LNZ0GAT,ALVCGAT,ALICGAT,PAMXGAT,PAMNGAT,
     D                CMASGAT,ROOTGAT,RSMNGAT,QA50GAT,VPDAGAT,VPDBGAT,
     E                PSGAGAT,PSGBGAT,PAIDGAT,HGTDGAT,ACVDGAT,ACIDGAT,
     F                ASVDGAT,ASIDGAT,AGVDGAT,AGIDGAT,ALGWGAT,ALGDGAT,
     +                ALGWVGAT,ALGWNGAT,ALGDVGAT,ALGDNGAT,
     G                THLQGAT,THICGAT,TBARGAT,RCANGAT,SCANGAT,TCANGAT,
     H                GROGAT, SNOGAT, TSNOGAT,RHOSGAT,ALBSGAT,ZBLDGAT,
     I                Z0ORGAT,ZSNLGAT,ZPLGGAT,ZPLSGAT,
     J                FCLOGAT,TAGAT,  VPDGAT, RHOAGAT,CSZGAT,
     +                FSDBGAT,FSFBGAT,REFGAT, BCSNGAT,
     K                FSVHGAT,RADJGAT,DLONGAT,RHSIGAT,DELZ,   DLZWGAT,
     L                ZBTWGAT,THPGAT, THMGAT, PSISGAT,BIGAT,  PSIWGAT,
     M                HCPSGAT,ISNDGAT,
     P                FCANCMXGAT,ICC,ICTEMMOD,RMATCGAT,ZOLNCGAT,
     Q                CMASVEGCGAT,AILCGAT,PAICGAT,L2MAX, NOL2PFTS,
     R                SLAICGAT,AILCGGAT,AILCGSGAT,FCANCGAT,FCANCSGAT,
     R                IDAY,   ILG,    1,      NML,  NBS,
     N                JLAT,N, ICAN,   ICAN+1, IGND,   IDISP,  IZREF,
     O                IWF,    IPAI,   IHGT,   IALC,   IALS,   IALG,
     P                ISNOALB,IGRALB, alvsctmgat,alirctmgat )
C
C-----------------------------------------------------------------------
C          * SURFACE TEMPERATURE AND FLUX CALCULATIONS.
C          * ADAPTED TO COUPLING OF CLASS3.6 AND CTEM
!          * by including in the arguments lfstatus
C
      CALL CLASST     (TBARC,  TBARG,  TBARCS, TBARGS, THLIQC, THLIQG,
     1  THICEC, THICEG, HCPC,   HCPG,   TCTOPC, TCBOTC, TCTOPG, TCBOTG,
     2  GZEROC, GZEROG, GZROCS, GZROGS, G12C,   G12G,   G12CS,  G12GS,
     3  G23C,   G23G,   G23CS,  G23GS,  QFREZC, QFREZG, QMELTC, QMELTG,
     4  EVAPC,  EVAPCG, EVAPG,  EVAPCS, EVPCSG, EVAPGS, TCANO,  TCANS,
     5  RAICAN, SNOCAN, RAICNS, SNOCNS, CHCAP,  CHCAPS, TPONDC, TPONDG,
     6  TPNDCS, TPNDGS, TSNOCS, TSNOGS, WSNOCS, WSNOGS, RHOSCS, RHOSGS,
     7  ITCTGAT,CDHGAT, CDMGAT, HFSGAT, TFXGAT, QEVPGAT,QFSGAT,
     8  PETGAT, GAGAT,  EFGAT,  GTGAT,  QGGAT,
     +  SFCTGAT,SFCUGAT,SFCVGAT,SFCQGAT,SFRHGAT,
     +  GTBS,   SFCUBS, SFCVBS, USTARBS,
     9  FSGVGAT,FSGSGAT,FSGGGAT,FLGVGAT,FLGSGAT,FLGGGAT,
     A  HFSCGAT,HFSSGAT,HFSGGAT,HEVCGAT,HEVSGAT,HEVGGAT,HMFCGAT,HMFNGAT,
     B  HTCCGAT,HTCSGAT,HTCGAT, QFCFGAT,QFCLGAT,DRGAT,  WTABGAT,ILMOGAT,
     C  UEGAT,  HBLGAT, TACGAT, QACGAT, ZRFMGAT,ZRFHGAT,ZDMGAT, ZDHGAT,
     D  VPDGAT, TADPGAT,RHOAGAT,FSVHGAT,FSIHGAT,FDLGAT, ULGAT,  VLGAT,
     E  TAGAT,  QAGAT,  PADRGAT,FC,     FG,     FCS,    FGS,    RBCOEF,
     F  FSVF,   FSVFS,  PRESGAT,VMODGAT,ALVSCN, ALIRCN, ALVSG,  ALIRG,
     G  ALVSCS, ALIRCS, ALVSSN, ALIRSN, ALVSGC, ALIRGC, ALVSSC, ALIRSC,
     H  TRVSCN, TRIRCN, TRVSCS, TRIRCS, RC,     RCS,    WTRGGAT,QLWOGAT,
     I  FRAINC, FSNOWC, FRAICS, FSNOCS, CMASSC, CMASCS, DISP,   DISPS,
     J  ZOMLNC, ZOELNC, ZOMLNG, ZOELNG, ZOMLCS, ZOELCS, ZOMLNS, ZOELNS,
     K  TBARGAT,THLQGAT,THICGAT,TPNDGAT,ZPNDGAT,TBASGAT,TCANGAT,TSNOGAT,
     L  ZSNOW,  RHOSGAT,WSNOGAT,THPGAT, THRGAT, THMGAT, THFCGAT,THLWGAT,
     +  TRSNOWC,TRSNOWG,ALSNO,  FSSBGAT, FROOT, FROOTS,
     M  RADJGAT,PREGAT, HCPSGAT,TCSGAT, TSFSGAT,DELZ,   DLZWGAT,ZBTWGAT,
     N  FTEMP,  FVAP,   RIB,    ISNDGAT,
     O  AILCGGAT,  AILCGSGAT, FCANCGAT,FCANCSGAT,CO2CONCGAT,CO2I1CGGAT,
     P  CO2I1CSGAT,CO2I2CGGAT,CO2I2CSGAT,CSZGAT,XDIFFUSGAT,SLAIGAT,ICC,
     Q  ICTEMMOD,RMATCTEMGAT,FCANCMXGAT,L2MAX,  NOL2PFTS,CFLUXCGGAT,
     R  CFLUXCSGAT,ANCSVEGGAT,ANCGVEGGAT,RMLCSVEGGAT,RMLCGVEGGAT,
     S  TCSNOW,GSNOW,ITC,ITCG,ITG,    ILG,    1,NML,  JLAT,N, ICAN,
     T  IGND,   IZREF,  ISLFD,  NLANDCS,NLANDGS,NLANDC, NLANDG, NLANDI,
     U  NBS,    ISNOALB,lfstatusgat, RN, H_MEP, LE_MEP, G_MEP)  


C
C-----------------------------------------------------------------------
C          * WATER BUDGET CALCULATIONS.
C

          CALL CLASSW  (THLQGAT,THICGAT,TBARGAT,TCANGAT,RCANGAT,SCANGAT,
     1                  ROFGAT, TROFGAT,SNOGAT, TSNOGAT,RHOSGAT,ALBSGAT,
     2                  WSNOGAT,ZPNDGAT,TPNDGAT,GROGAT, TBASGAT,GFLXGAT,
     3                  PCFCGAT,PCLCGAT,PCPNGAT,PCPGGAT,QFCFGAT,QFCLGAT,
     4                  QFNGAT, QFGGAT, QFCGAT, HMFCGAT,HMFGGAT,HMFNGAT,
     5                  HTCCGAT,HTCSGAT,HTCGAT, ROFCGAT,ROFNGAT,ROVGGAT,
     6                  WTRSGAT,WTRGGAT,ROFOGAT,ROFSGAT,ROFBGAT,
     7                  TROOGAT,TROSGAT,TROBGAT,QFSGAT, QFXGAT, RHOAGAT,
     8                  TBARC,  TBARG,  TBARCS, TBARGS, THLIQC, THLIQG,
     9                  THICEC, THICEG, HCPC,   HCPG,   RPCPGAT,TRPCGAT,
     A                  SPCPGAT,TSPCGAT,PREGAT, TAGAT,  RHSIGAT,GGEOGAT,
     B                  FC,     FG,     FCS,    FGS,    TPONDC, TPONDG,
     C                  TPNDCS, TPNDGS, EVAPC,  EVAPCG, EVAPG,  EVAPCS,
     D                  EVPCSG, EVAPGS, QFREZC, QFREZG, QMELTC, QMELTG,
     E                  RAICAN, SNOCAN, RAICNS, SNOCNS, FSVF,    FSVFS,
     F                  CWLCAP, CWFCAP, CWLCPS, CWFCPS, TCANO,
     G                  TCANS,  CHCAP,  CHCAPS, CMASSC, CMASCS, ZSNOW,
     H                  GZEROC, GZEROG, GZROCS, GZROGS, G12C,   G12G,
     I                  G12CS,  G12GS,  G23C,   G23G,   G23CS,  G23GS,
     J                  TSNOCS, TSNOGS, WSNOCS, WSNOGS, RHOSCS, RHOSGS,
     K                  ZPLIMC, ZPLIMG, ZPLMCS, ZPLMGS, TSFSGAT,
     L                  TCTOPC, TCBOTC, TCTOPG, TCBOTG, FROOT,   FROOTS,
     M                  THPGAT, THRGAT, THMGAT, BIGAT,  PSISGAT,GRKSGAT,
     N                  THRAGAT,THFCGAT,DRNGAT, HCPSGAT,DELZ,
     O                  DLZWGAT,ZBTWGAT,XSLPGAT,GRKFGAT,WFSFGAT,WFCIGAT,
     P                  ISNDGAT,IGDRGAT,
     Q                  IWF,    ILG,    1,      NML,    N,
     R                  JLAT,   ICAN,   IGND,   IGND+1, IGND+2,
     S                  NLANDCS,NLANDGS,NLANDC, NLANDG, NLANDI )
C

C========================================================================
C
      CALL CLASSZ (1,      CTVSTP, CTSSTP, CT1STP, CT2STP, CT3STP,
     1             WTVSTP, WTSSTP, WTGSTP,
     2             FSGVGAT,FLGVGAT,HFSCGAT,HEVCGAT,HMFCGAT,HTCCGAT,
     3             FSGSGAT,FLGSGAT,HFSSGAT,HEVSGAT,HMFNGAT,HTCSGAT,
     4             FSGGGAT,FLGGGAT,HFSGGAT,HEVGGAT,HMFGGAT,HTCGAT,
     5             PCFCGAT,PCLCGAT,QFCFGAT,QFCLGAT,ROFCGAT,WTRCGAT,
     6             PCPNGAT,QFNGAT, ROFNGAT,WTRSGAT,PCPGGAT,QFGGAT,
     7             QFCGAT, ROFGAT, WTRGGAT,CMAIGAT,RCANGAT,SCANGAT,
     8             TCANGAT,SNOGAT, WSNOGAT,TSNOGAT,THLQGAT,THICGAT,
     9             HCPSGAT,THPGAT, DLZWGAT,TBARGAT,ZPNDGAT,TPNDGAT,
     A             DELZ,   FCS,    FGS,    FC,     FG,
     B             1,      NML,    ILG,    IGND,   N    )
C
C=======================================================================

C===================== CTEM ============================================ \
c     * accumulate output data for running ctem.
c
      do 660 i=1,nml
         uvaccgat_m(i)=uvaccgat_m(i)+ulgat(i)
660   continue
c
c     accumulate variables not already accumulated but which are required by
c     ctem.
c
      if (ctem_on) then
        do 700 i = 1, nml
c
          alswacc_gat(i)=alswacc_gat(i)+alvsgat(i)*fsvhgat(i)
          allwacc_gat(i)=allwacc_gat(i)+alirgat(i)*fsihgat(i)
          fsinacc_gat(i)=fsinacc_gat(i)+FSSROW(I)
          flinacc_gat(i)=flinacc_gat(i)+fdlgat(i)
          flutacc_gat(i)=flutacc_gat(i)+sbc*gtgat(i)**4
          pregacc_gat(i)=pregacc_gat(i)+pregat(i)*delt
c
          fsnowacc_m(i)=fsnowacc_m(i)+fsnogat(i)
          tcanoaccgat_m(i)=tcanoaccgat_m(i)+tcano(i)
          tcansacc_m(i)=tcansacc_m(i)+tcans(i)
          taaccgat_m(i)=taaccgat_m(i)+tagat(i)
          vvaccgat_m(i)=vvaccgat_m(i)+ vlgat(i)
c
          do 710 j=1,ignd
             tbaraccgat_m(i,j)=tbaraccgat_m(i,j)+tbargat(i,j)
             tbarcacc_m(i,j)=tbarcacc_m(i,j)+tbarc(i,j)
             tbarcsacc_m(i,j)=tbarcsacc_m(i,j)+tbarcs(i,j)
             tbargacc_m(i,j)=tbargacc_m(i,j)+tbarg(i,j)
             tbargsacc_m(i,j)=tbargsacc_m(i,j)+tbargs(i,j)
             thliqcacc_m(i,j)=thliqcacc_m(i,j)+thliqc(i,j)
             thliqgacc_m(i,j)=thliqgacc_m(i,j)+thliqg(i,j)
             thicecacc_m(i,j)=thicecacc_m(i,j)+thicec(i,j)
710       continue
c
          do 713 j = 1, icc
            ancsvgac_m(i,j)=ancsvgac_m(i,j)+ancsveggat(i,j)
            ancgvgac_m(i,j)=ancgvgac_m(i,j)+ancgveggat(i,j)
            rmlcsvga_m(i,j)=rmlcsvga_m(i,j)+rmlcsveggat(i,j)
            rmlcgvga_m(i,j)=rmlcgvga_m(i,j)+rmlcgveggat(i,j)
713       continue
c
700     continue
      endif !if (ctem_on)
c
      if(ncount.eq.nday) then
c
      do 855 i=1,nml
          uvaccgat_m(i)=uvaccgat_m(i)/real(nday)
          vvaccgat_m(i)=vvaccgat_m(i)/real(nday)
c
c         daily averages of accumulated variables for ctem
c
          if (ctem_on) then
c
c           net radiation and precipitation estimates for ctem's bioclim
c
            if(fsinacc_gat(i).gt.0.0) then
              alswacc_gat(i)=alswacc_gat(i)/(fsinacc_gat(i)*0.5)
              allwacc_gat(i)=allwacc_gat(i)/(fsinacc_gat(i)*0.5)
            else
              alswacc_gat(i)=0.0
              allwacc_gat(i)=0.0
            endif
c
            fsinacc_gat(i)=fsinacc_gat(i)/real(nday)
            flinacc_gat(i)=flinacc_gat(i)/real(nday)
            flutacc_gat(i)=flutacc_gat(i)/real(nday)
c
            altot_gat=(alswacc_gat(i)+allwacc_gat(i))/2.0
            fsstar_gat=fsinacc_gat(i)*(1.-altot_gat)
            flstar_gat=flinacc_gat(i)-flutacc_gat(i)
            netrad_gat(i)=fsstar_gat+flstar_gat
            preacc_gat(i)=pregacc_gat(i)
c
            fsnowacc_m(i)=fsnowacc_m(i)/real(nday)
            tcanoaccgat_m(i)=tcanoaccgat_m(i)/real(nday)
            tcansacc_m(i)=tcansacc_m(i)/real(nday)
            taaccgat_m(i)=taaccgat_m(i)/real(nday)
c
            do 831 j=1,ignd
              tbaraccgat_m(i,j)=tbaraccgat_m(i,j)/real(nday)
              tbarcacc_m(i,j) = tbaraccgat_m(i,j)
              tbarcsacc_m(i,j) = tbaraccgat_m(i,j)
              tbargacc_m(i,j) = tbaraccgat_m(i,j)
              tbargsacc_m(i,j) = tbaraccgat_m(i,j)
c
              thliqcacc_m(i,j)=thliqcacc_m(i,j)/real(nday)
              thliqgacc_m(i,j)=thliqgacc_m(i,j)/real(nday)
              thicecacc_m(i,j)=thicecacc_m(i,j)/real(nday)
831         continue
c
            do 832 j = 1, icc
              ancsvgac_m(i,j)=ancsvgac_m(i,j)/real(nday)
              ancgvgac_m(i,j)=ancgvgac_m(i,j)/real(nday)
              rmlcsvga_m(i,j)=rmlcsvga_m(i,j)/real(nday)
              rmlcgvga_m(i,j)=rmlcgvga_m(i,j)/real(nday)
832         continue
c
c           pass on mean monthly lightning for the current month to ctem
c           lightng(i)=mlightng(i,month)
c
c           in a very simple way try to interpolate monthly lightning to
c           daily lightning
c
            if(iday.ge.15.and.iday.le.45)then ! mid jan - mid feb
              month1=1
              month2=2
              xday=iday-15
            else if(iday.ge.46.and.iday.le.74)then ! mid feb - mid mar
              month1=2
              month2=3
              xday=iday-46
            else if(iday.ge.75.and.iday.le.105)then ! mid mar - mid apr
              month1=3
              month2=4
              xday=iday-75
            else if(iday.ge.106.and.iday.le.135)then ! mid apr - mid may
              month1=4
              month2=5
              xday=iday-106
            else if(iday.ge.136.and.iday.le.165)then ! mid may - mid june
              month1=5
              month2=6
              xday=iday-136
            else if(iday.ge.166.and.iday.le.196)then ! mid june - mid july
              month1=6
              month2=7
              xday=iday-166
            else if(iday.ge.197.and.iday.le.227)then ! mid july - mid aug
              month1=7
              month2=8
              xday=iday-197
            else if(iday.ge.228.and.iday.le.258)then ! mid aug - mid sep
              month1=8
              month2=9
              xday=iday-228
            else if(iday.ge.259.and.iday.le.288)then ! mid sep - mid oct
              month1=9
              month2=10
              xday=iday-259
            else if(iday.ge.289.and.iday.le.319)then ! mid oct - mid nov
              month1=10
              month2=11
              xday=iday-289
            else if(iday.ge.320.and.iday.le.349)then ! mid nov - mid dec
              month1=11
              month2=12
              xday=iday-320
            else if(iday.ge.350.or.iday.lt.14)then ! mid dec - mid jan
              month1=12
              month2=1
              xday=iday-350
              if(xday.lt.0)xday=iday
            endif
c
            lightng(i)=mlightnggat(i,month1)+(real(xday)/30.0)*
     &                 (mlightnggat(i,month2)-mlightnggat(i,month1))
c
            if (obswetf) then
              wetfracgrd(i)=wetfrac_mon(i,month1)+(real(xday)/30.0)*
     &                 (wetfrac_mon(i,month2)-wetfrac_mon(i,month1))
            endif !obswetf

          endif ! if(ctem_on)
c
855   continue
c
c     call canadian terrestrial ecosystem model which operates at a
c     daily time step, and uses daily accumulated values of variables
c     simulated by class.
c
      if (ctem_on) then
c
        call ctem ( fcancmxgat, fsnowacc_m,    sandgat,    claygat,
     2                      1,        nml,        iday,    radjgat,
     4          tcanoaccgat_m,  tcansacc_m, tbarcacc_m,tbarcsacc_m,
     5             tbargacc_m, tbargsacc_m, taaccgat_m,    dlzwgat,
     6             ancsvgac_m,  ancgvgac_m, rmlcsvga_m, rmlcgvga_m,
     7                zbtwgat, thliqcacc_m,thliqgacc_m,     deltat,
     8             uvaccgat_m,  vvaccgat_m,    lightng,prbfrhucgat,
     9            extnprobgat,   stdalngat,tbaraccgat_m,  popdon,
     a               nol2pfts, pfcancmxgat, nfcancmxgat,  lnduseon,
     b            thicecacc_m,     sdepgat,    spinfast,   todfrac,
     &                compete,  netrad_gat,  preacc_gat,
     &                 popdin,  dofire, dowetlands,obswetf, isndgat,
     &                faregat,      mosaic, WETFRACGRD, wetfrac_sgrd,
c    -------------- inputs used by ctem are above this line ---------
     c            stemmassgat, rootmassgat, litrmassgat, gleafmasgat,
     d            bleafmasgat, soilcmasgat,    ailcggat,    ailcgat,
     e               zolncgat,  rmatctemgat,   rmatcgat,  ailcbgat,
     f            flhrlossgat,  pandaysgat, lfstatusgat, grwtheffgat,
     g            lystmmasgat, lyrotmasgat, tymaxlaigat, vgbiomasgat,
     h            gavgltmsgat, gavgscmsgat, stmhrlosgat,     slaigat,
     i             bmasveggat, cmasvegcgat,  colddaysgat, rothrlosgat,
     j                fcangat,  alvsctmgat,   alirctmgat,  gavglaigat,
     &                  tcurm,    srpcuryr,     dftcuryr,  inibioclim,
     &                 tmonth,    anpcpcur,      anpecur,     gdd5cur,
     &               surmncur,    defmncur,     srplscur,    defctcur,
     &            geremortgat, intrmortgat,    lambdagat, lyglfmasgat,
     &            pftexistgat,      twarmm,       tcoldm,        gdd5,
     1                aridity,    srplsmon,     defctmon,    anndefct,
     2               annsrpls,      annpcp,  dry_season_length,
     &              burnvegfgat, pstemmassgat, pgleafmassgat,
c    -------------- inputs updated by ctem are above this line ------
     k                 nppgat,      nepgat, hetroresgat, autoresgat,
     l            soilcrespgat,       rmgat,       rggat,      nbpgat,
     m              litresgat,    socresgat,     gppgat, dstcemlsgat,
     n            litrfallgat,  humiftrsgat, veghghtgat, rootdpthgat,
     o                 rmlgat,      rmsgat,     rmrgat,  tltrleafgat,
     p            tltrstemgat, tltrrootgat, leaflitrgat, roottempgat,
     q             afrleafgat,  afrstemgat,  afrrootgat, wtstatusgat,
     r            ltstatusgat, burnfracgat, probfiregat, lucemcomgat,
     s            lucltringat, lucsocingat,   nppveggat, grclarea,
     t            dstcemls3gat,    paicgat,    slaicgat,
     u            emit_co2gat,  emit_cogat,  emit_ch4gat, emit_nmhcgat,
     v             emit_h2gat, emit_noxgat,  emit_n2ogat, emit_pm25gat,
     w            emit_tpmgat,  emit_tcgat,   emit_ocgat,   emit_bcgat,
     &               btermgat,    ltermgat,     mtermgat,
     &            ccgat,             mmgat,
     &          rmlvegaccgat,    rmsveggat,  rmrveggat,  rgveggat,
     &       vgbiomas_veggat, gppveggat,  nepveggat, nbpveggat,
     &        hetroresveggat, autoresveggat, litresveggat,
     &           soilcresveggat, nml, ilmos, jlmos, CH4WET1GAT,
     &          CH4WET2GAT, WETFDYNGAT, CH4DYN1GAT, CH4DYN2GAT)
c    ---------------- outputs are listed above this line ------------
c
      endif  !if(ctem_on)
c
c     reset mosaic accumulator arrays.
c
      do 655 i=1,nml
         uvaccgat_m(i)=0.0
655   continue
c
      if (ctem_on) then
        do 705 i = 1, nml
c
c         competitition related variables added by y. peng \\
          fsinacc_gat(i)=0.
          flinacc_gat(i)=0.
          flutacc_gat(i)=0.
          alswacc_gat(i)=0.
          allwacc_gat(i)=0.
          pregacc_gat(i)=0.
c         competitition related variables added by y. peng //
c
          fsnowacc_m(i)=0.0
          tcanoaccgat_out(i)=tcanoaccgat_m(i)
          tcanoaccgat_m(i)=0.0
c
          tcansacc_m(i)=0.0
          taaccgat_m(i)=0.0
          vvaccgat_m(i)=0.0
c
          do 715 j=1,ignd
             tbaraccgat_m(i,j)=0.0
             tbarcacc_m(i,j)=0.0
             tbarcsacc_m(i,j)=0.0
             tbargacc_m(i,j)=0.0
             tbargsacc_m(i,j)=0.0
             thliqcacc_m(i,j)=0.0
             thliqgacc_m(i,j)=0.0
             thicecacc_m(i,j)=0.0
715       continue
c
          do 716 j = 1, icc
            ancsvgac_m(i,j)=0.0
            ancgvgac_m(i,j)=0.0
            rmlcsvga_m(i,j)=0.0
            rmlcgvga_m(i,j)=0.0
716       continue
c
705     continue
      endif  ! if(ctem_on)
      endif  ! if(ncount.eq.nday)
C===================== CTEM ============================================ /
C
      CALL CLASSS (TBARROT,THLQROT,THICROT,TSFSROT,TPNDROT,
     1             ZPNDROT,TBASROT,ALBSROT,TSNOROT,RHOSROT,
     2             SNOROT, GTROT, TCANROT,RCANROT,SCANROT,
     3             GROROT, CMAIROT,TACROT, QACROT, WSNOROT,
     +             REFROT, BCSNROT,EMISROT,SALBROT,CSALROT,
     4             ILMOS,JLMOS,NML,NLAT,NTLD,NMOS,
     5             ILG,IGND,ICAN,ICAN+1,NBS,
     6             TBARGAT,THLQGAT,THICGAT,TSFSGAT,TPNDGAT,
     7             ZPNDGAT,TBASGAT,ALBSGAT,TSNOGAT,RHOSGAT,
     8             SNOGAT, GTGAT, TCANGAT,RCANGAT,SCANGAT,
     9             GROGAT, CMAIGAT,TACGAT, QACGAT, WSNOGAT,
     +             REFGAT, BCSNGAT,EMISGAT,SALBGAT,CSALGAT)

C
C    * SCATTER OPERATION ON DIAGNOSTIC VARIABLES SPLIT OUT OF
C    * CLASSS FOR CONSISTENCY WITH GCM APPLICATIONS.
C
      DO 380 K=1,NML
          CDHROT (ILMOS(K),JLMOS(K))=CDHGAT (K)
          CDMROT (ILMOS(K),JLMOS(K))=CDMGAT (K)
          HFSROT (ILMOS(K),JLMOS(K))=HFSGAT (K)
          TFXROT (ILMOS(K),JLMOS(K))=TFXGAT (K)
          QEVPROT(ILMOS(K),JLMOS(K))=QEVPGAT(K)
          QFSROT (ILMOS(K),JLMOS(K))=QFSGAT (K)
          QFXROT (ILMOS(K),JLMOS(K))=QFXGAT (K)
          PETROT (ILMOS(K),JLMOS(K))=PETGAT (K)
          GAROT  (ILMOS(K),JLMOS(K))=GAGAT  (K)
          EFROT  (ILMOS(K),JLMOS(K))=EFGAT  (K)
          QGROT  (ILMOS(K),JLMOS(K))=QGGAT  (K)
          ALVSROT(ILMOS(K),JLMOS(K))=ALVSGAT(K)
          ALIRROT(ILMOS(K),JLMOS(K))=ALIRGAT(K)
          SFCTROT(ILMOS(K),JLMOS(K))=SFCTGAT(K)
          SFCUROT(ILMOS(K),JLMOS(K))=SFCUGAT(K)
          SFCVROT(ILMOS(K),JLMOS(K))=SFCVGAT(K)
          SFCQROT(ILMOS(K),JLMOS(K))=SFCQGAT(K)
          FSNOROT(ILMOS(K),JLMOS(K))=FSNOGAT(K)
          FSGVROT(ILMOS(K),JLMOS(K))=FSGVGAT(K)
          FSGSROT(ILMOS(K),JLMOS(K))=FSGSGAT(K)
          FSGGROT(ILMOS(K),JLMOS(K))=FSGGGAT(K)
          FLGVROT(ILMOS(K),JLMOS(K))=FLGVGAT(K)
          FLGSROT(ILMOS(K),JLMOS(K))=FLGSGAT(K)
          FLGGROT(ILMOS(K),JLMOS(K))=FLGGGAT(K)
          HFSCROT(ILMOS(K),JLMOS(K))=HFSCGAT(K)
          HFSSROT(ILMOS(K),JLMOS(K))=HFSSGAT(K)
          HFSGROT(ILMOS(K),JLMOS(K))=HFSGGAT(K)
          HEVCROT(ILMOS(K),JLMOS(K))=HEVCGAT(K)
          HEVSROT(ILMOS(K),JLMOS(K))=HEVSGAT(K)
          HEVGROT(ILMOS(K),JLMOS(K))=HEVGGAT(K)
          HMFCROT(ILMOS(K),JLMOS(K))=HMFCGAT(K)
          HMFNROT(ILMOS(K),JLMOS(K))=HMFNGAT(K)
          HTCCROT(ILMOS(K),JLMOS(K))=HTCCGAT(K)
          HTCSROT(ILMOS(K),JLMOS(K))=HTCSGAT(K)
          PCFCROT(ILMOS(K),JLMOS(K))=PCFCGAT(K)
          PCLCROT(ILMOS(K),JLMOS(K))=PCLCGAT(K)
          PCPNROT(ILMOS(K),JLMOS(K))=PCPNGAT(K)
          PCPGROT(ILMOS(K),JLMOS(K))=PCPGGAT(K)
          QFGROT (ILMOS(K),JLMOS(K))=QFGGAT (K)
          QFNROT (ILMOS(K),JLMOS(K))=QFNGAT (K)
          QFCLROT(ILMOS(K),JLMOS(K))=QFCLGAT(K)
          QFCFROT(ILMOS(K),JLMOS(K))=QFCFGAT(K)
          ROFROT (ILMOS(K),JLMOS(K))=ROFGAT (K)
          ROFOROT(ILMOS(K),JLMOS(K))=ROFOGAT(K)
          ROFSROT(ILMOS(K),JLMOS(K))=ROFSGAT(K)
          ROFBROT(ILMOS(K),JLMOS(K))=ROFBGAT(K)
          TROFROT(ILMOS(K),JLMOS(K))=TROFGAT(K)
          TROOROT(ILMOS(K),JLMOS(K))=TROOGAT(K)
          TROSROT(ILMOS(K),JLMOS(K))=TROSGAT(K)
          TROBROT(ILMOS(K),JLMOS(K))=TROBGAT(K)
          ROFCROT(ILMOS(K),JLMOS(K))=ROFCGAT(K)
          ROFNROT(ILMOS(K),JLMOS(K))=ROFNGAT(K)
          ROVGROT(ILMOS(K),JLMOS(K))=ROVGGAT(K)
          WTRCROT(ILMOS(K),JLMOS(K))=WTRCGAT(K)
          WTRSROT(ILMOS(K),JLMOS(K))=WTRSGAT(K)
          WTRGROT(ILMOS(K),JLMOS(K))=WTRGGAT(K)
          DRROT  (ILMOS(K),JLMOS(K))=DRGAT  (K)
          WTABROT(ILMOS(K),JLMOS(K))=WTABGAT(K)
          ILMOROT(ILMOS(K),JLMOS(K))=ILMOGAT(K)
          UEROT  (ILMOS(K),JLMOS(K))=UEGAT(K)
          HBLROT (ILMOS(K),JLMOS(K))=HBLGAT(K)
380   CONTINUE
C
      DO 390 L=1,IGND
      DO 390 K=1,NML
          HMFGROT(ILMOS(K),JLMOS(K),L)=HMFGGAT(K,L)
          HTCROT (ILMOS(K),JLMOS(K),L)=HTCGAT (K,L)
          QFCROT (ILMOS(K),JLMOS(K),L)=QFCGAT (K,L)
          GFLXROT(ILMOS(K),JLMOS(K),L)=GFLXGAT(K,L)
390   CONTINUE
C
      DO 430 M=1,50
          DO 420 L=1,6
              DO 410 K=1,NML
                  ITCTROT(ILMOS(K),JLMOS(K),L,M)=ITCTGAT(K,L,M)
410           CONTINUE
420       CONTINUE
430   CONTINUE
C

C
C===================== CTEM ============================================ \
C
      call ctems2(fcancmxrow,rmatcrow,zolncrow,paicrow,
     1      ailcrow,     ailcgrow,    cmasvegcrow,  slaicrow,
     2      ailcgsrow,   fcancsrow,   fcancrow,     rmatctemrow,
     3      co2concrow,  co2i1cgrow,  co2i1csrow,   co2i2cgrow,
     4      co2i2csrow,  xdiffus,     slairow,      cfluxcgrow,
     5      cfluxcsrow,  ancsvegrow,  ancgvegrow,   rmlcsvegrow,
     6      rmlcgvegrow, canresrow,   SDEPROT,
     7      SANDROT,     CLAYROT,     ORGMROT,
     8      anvegrow,    rmlvegrow,   tcanoaccrow_m,tbaraccrow_m,
     9      uvaccrow_m,  vvaccrow_m,  mlightnggrd,  prbfrhucgrd,
     a      extnprobgrd, stdalngrd,   pfcancmxrow,  nfcancmxrow,
     b      stemmassrow, rootmassrow, litrmassrow,  gleafmasrow,
     c      bleafmasrow, soilcmasrow, ailcbrow,     flhrlossrow,
     d      pandaysrow,  lfstatusrow, grwtheffrow,  lystmmasrow,
     e      lyrotmasrow, tymaxlairow, vgbiomasrow,  gavgltmsrow,
     f      stmhrlosrow, bmasvegrow,  colddaysrow,  rothrlosrow,
     g      alvsctmrow,  alirctmrow,  gavglairow,   npprow,
     h      neprow,      hetroresrow, autoresrow,   soilcresprow,
     i      rmrow,       rgrow,       nbprow,       litresrow,
     j      socresrow,   gpprow,      dstcemlsrow,  litrfallrow,
     k      humiftrsrow, veghghtrow,  rootdpthrow,  rmlrow,
     l      rmsrow,      rmrrow,      tltrleafrow,  tltrstemrow,
     m      tltrrootrow, leaflitrrow, roottemprow,  afrleafrow,
     n      afrstemrow,  afrrootrow,  wtstatusrow,  ltstatusrow,
     o      burnfracrow, probfirerow, lucemcomrow,  lucltrinrow,
     p      lucsocinrow, nppvegrow,   dstcemls3row,
     q      FAREROT,     gavgscmsrow, tcanoaccrow_out,
     &      rmlvegaccrow, rmsvegrow,  rmrvegrow,    rgvegrow,
     &      vgbiomas_vegrow,gppvegrow,nepvegrow,ailcminrow,ailcmaxrow,
     &      FCANROT,      pftexistrow,
     &      emit_co2row,  emit_corow, emit_ch4row,  emit_nmhcrow,
     &      emit_h2row,   emit_noxrow,emit_n2orow,  emit_pm25row,
     &      emit_tpmrow,  emit_tcrow, emit_ocrow,   emit_bcrow,
     &      btermrow,     ltermrow,   mtermrow,
     &      nbpvegrow,   hetroresvegrow, autoresvegrow,litresvegrow,
     &      soilcresvegrow, burnvegfrow, pstemmassrow, pgleafmassrow,
     &      CH4WET1ROW, CH4WET2ROW,
     &      WETFDYNROW, CH4DYN1ROW, CH4DYN2ROW,
c    ----
     r      ilmos,       jlmos,       iwmos,        jwmos,
     s      nml,     fcancmxgat,  rmatcgat,    zolncgat,     paicgat,
     v      ailcgat,     ailcggat,    cmasvegcgat,  slaicgat,
     w      ailcgsgat,   fcancsgat,   fcancgat,     rmatctemgat,
     x      co2concgat,  co2i1cggat,  co2i1csgat,   co2i2cggat,
     y      co2i2csgat,  xdiffusgat,  slaigat,      cfluxcggat,
     z      cfluxcsgat,  ancsveggat,  ancgveggat,   rmlcsveggat,
     1      rmlcgveggat, canresgat,   sdepgat,
     2      sandgat,     claygat,     orgmgat,
     3      anveggat,    rmlveggat,   tcanoaccgat_m,tbaraccgat_m,
     4      uvaccgat_m,  vvaccgat_m,  mlightnggat,  prbfrhucgat,
     5      extnprobgat, stdalngat,   pfcancmxgat,  nfcancmxgat,
     6      stemmassgat, rootmassgat, litrmassgat,  gleafmasgat,
     7      bleafmasgat, soilcmasgat, ailcbgat,     flhrlossgat,
     8      pandaysgat,  lfstatusgat, grwtheffgat,  lystmmasgat,
     9      lyrotmasgat, tymaxlaigat, vgbiomasgat,  gavgltmsgat,
     a      stmhrlosgat, bmasveggat,  colddaysgat,  rothrlosgat,
     b      alvsctmgat,  alirctmgat,  gavglaigat,   nppgat,
     c      nepgat,      hetroresgat, autoresgat,   soilcrespgat,
     d      rmgat,       rggat,       nbpgat,       litresgat,
     e      socresgat,   gppgat,      dstcemlsgat,  litrfallgat,
     f      humiftrsgat, veghghtgat,  rootdpthgat,  rmlgat,
     g      rmsgat,      rmrgat,      tltrleafgat,  tltrstemgat,
     h      tltrrootgat, leaflitrgat, roottempgat,  afrleafgat,
     i      afrstemgat,  afrrootgat,  wtstatusgat,  ltstatusgat,
     j      burnfracgat, probfiregat, lucemcomgat,  lucltringat,
     k      lucsocingat, nppveggat,   dstcemls3gat,
     l      faregat,     gavgscmsgat, tcanoaccgat_out,
     &      rmlvegaccgat, rmsveggat,  rmrveggat,    rgveggat,
     &      vgbiomas_veggat,gppveggat,nepveggat,ailcmingat,ailcmaxgat,
     &      fcangat,      pftexistgat,
     &      emit_co2gat,  emit_cogat, emit_ch4gat,  emit_nmhcgat,
     &      emit_h2gat,   emit_noxgat,emit_n2ogat,  emit_pm25gat,
     &      emit_tpmgat,  emit_tcgat, emit_ocgat,   emit_bcgat,
     &      btermgat,     ltermgat,   mtermgat,
     &      nbpveggat, hetroresveggat, autoresveggat,litresveggat,
     &      soilcresveggat, burnvegfgat, pstemmassgat, pgleafmassgat,
     &      CH4WET1GAT, CH4WET2GAT,
     &      WETFDYNGAT, CH4DYN1GAT, CH4DYN2GAT)
c
C===================== CTEM ============================================ /
C
C=======================================================================
C     * WRITE FIELDS FROM CURRENT TIME STEP TO OUTPUT FILES.

6100  FORMAT(1X,I4,I5,9F8.2,2F8.3,F12.4,F8.2,2(A6,I2))
6200  FORMAT(1X,I4,I5,3(F8.2,2F6.3),F8.2,2F8.4,F8.2,F8.3,2(A6,I2))
6201  FORMAT(1X,I4,I5,5(F7.2,2F6.3),2(A6,I2))
6300  FORMAT(1X,I4,I5,3F9.2,F8.2,F10.2,E12.3,2F12.3,A6,I2)
6400  FORMAT(I2,2X,I2,2X,I3,2X,I4,11(2X,F9.2),E12.3,7(2X,F7.3))
!6400  FORMAT(1X,I2,I3,I5,I6,9F8.2,2F7.3,E11.3,F8.2,F12.4,5F9.5,2(A6,I2))
6500  FORMAT(1X,I2,I3,I5,I6,3(F7.2,2F6.3),F8.2,2F8.4,F8.2,4F8.3,
     &       2F7.3,2(A6,I2))
6600  FORMAT(1X,I2,I3,I5,2F10.2,E12.3,F10.2,F8.2,F10.2,E12.3,2(A6,I2))
6501  FORMAT(1X,I2,I3,I5,I6,5(F7.2,2F6.3),2(A6,I2))
6601  FORMAT(1X,I2,I3,I5,I6,7(F7.2,2F6.3),10F9.4,2(A6,I2))
6700  FORMAT(1X,I2,I3,I5,I6,2X,12E11.4,2(A6,I2))
6800  FORMAT(1X,I2,I3,I5,I6,2X,22(F10.4,2X),2(A6,I2))
6900  FORMAT(1X,I2,I3,I5,I6,2X,18(E12.4,2X),2(A6,I2))
C
C===================== CTEM ============================================ \
c
c  fc,fg,fcs and fgs are one_dimensional in class subroutines
c  the transformations here to grid_cell mean fc_g,fg_g,fcs_g and fgs_g
c  are only applicable when nltest=1 (e.g., one grid cell)
c
      do i=1,nltest
        fc_g(i)=0.0
        fg_g(i)=0.0
        fcs_g(i)=0.0
        fgs_g(i)=0.0
        do m=1,nmtest
          fc_g(i)=fc_g(i)+fc(m)
          fg_g(i)=fg_g(i)+fg(m)
          fcs_g(i)=fcs_g(i)+fcs(m)
          fgs_g(i)=fgs_g(i)+fgs(m)
        enddo
      enddo
c
C===================== CTEM =====================================/
C

      ACTLYR=0.0
      FTABLE=0.0
      DO 440 J=1,IGND
          IF(ABS(TBARGAT(1,J)-TFREZ).LT.0.0001) THEN
              IF(ISNDGAT(1,J).GT.-3) THEN
                  ACTLYR=ACTLYR+(THLQGAT(1,J)/(THLQGAT(1,J)+
     1                THICGAT(1,J)))*DLZWGAT(1,J)
              ELSEIF(ISNDGAT(1,J).EQ.-3) THEN
                  ACTLYR=ACTLYR+DELZ(J)
              ENDIF
          ELSEIF(TBARGAT(1,J).GT.TFREZ) THEN
              ACTLYR=ACTLYR+DELZ(J)
          ENDIF
          IF(ABS(TBARGAT(1,J)-TFREZ).LT.0.0001) THEN
              IF(ISNDGAT(1,J).GT.-3) THEN
                  FTABLE=FTABLE+(THICGAT(1,J)/(THLQGAT(1,J)+
     1                THICGAT(1,J)-THMGAT(1,J)))*DLZWGAT(1,J)
              ELSE
                  FTABLE=FTABLE+DELZ(J)
              ENDIF
          ELSEIF(TBARGAT(1,J).LT.TFREZ) THEN
              FTABLE=FTABLE+DELZ(J)
          ENDIF
440   CONTINUE
C
      IF(IDAY.GE.182 .AND. IDAY.LE.243)  THEN
          ALAVG=ALAVG+ACTLYR
          NAL=NAL+1
          IF(ACTLYR.GT.ALMAX) ALMAX=ACTLYR
      ENDIF
C
      IF(IDAY.GE.1 .AND. IDAY.LE.59)   THEN
          FTAVG=FTAVG+FTABLE
          NFT=NFT+1
          IF(FTABLE.GT.FTMAX) FTMAX=FTABLE
      ENDIF

      if (.not. parallelrun) then ! stand alone mode, include half-hourly
c                                 ! output for class & ctem
C
      DO 450 I=1,NLTEST

c       initialization of various grid-averaged variables
        call resetgridavg(nltest)

       DO 425 M=1,NMTEST
          IF(FSSROW(I).GT.0.0) THEN
C              ALTOT=(ALVSROT(I,M)+ALIRROT(I,M))/2.0
              ALTOT=(FSSROW(I)-FSGGGAT(1))/FSSROW(I)  !FLAG I adopt the runclass approach of using 1 for index here. JM Jul 2015.
          ELSE
              ALTOT=0.0
          ENDIF
          FSSTAR=FSSROW(I)*(1.0-ALTOT)
          FLSTAR=FDLROW(I)-SBC*GTROT(I,M)**4
          QH=HFSROT(I,M)
          QE=QEVPROT(I,M)
C          BEG=FSSTAR+FLSTAR-QH-QE !(commented out in runclass.fieldsite)
          BEG=GFLXGAT(1,1)  !FLAG!
C          USTARBS=UVROW(1)*SQRT(CDMROT(I,M)) !FLAG (commented out in runclass.fieldsite)
          SNOMLT=HMFNROT(I,M)
          IF(RHOSROT(I,M).GT.0.0) THEN
              ZSN=SNOROT(I,M)/RHOSROT(I,M)
          ELSE
              ZSN=0.0
          ENDIF
          IF(TCANROT(I,M).GT.0.01) THEN
              TCN=TCANROT(I,M)-TFREZ
          ELSE
              TCN=0.0
          ENDIF
          TSURF=FCS(I)*TSFSGAT(I,1)+FGS(I)*TSFSGAT(I,2)+
     1           FC(I)*TSFSGAT(I,3)+FG(I)*TSFSGAT(I,4)
C          IF(FSSROW(I).GT.0.0 .AND. (FCS(I)+FC(I)).GT.0.0) THEN
C          IF(FSSROW(I).GT.0.0) THEN
              NFS=NFS+1
              ITA=NINT(TAROW(I)-TFREZ)
              ITCAN=NINT(TCN)
              ITAC=NINT(TACGAT(I)-TFREZ)
              ITSCR=NINT(SFCTGAT(I)-TFREZ)
              ITS=NINT(TSURF-TFREZ)
C              ITD=ITS-ITA
              ITD=ITCAN-ITA
              ITD2=ITCAN-ITSCR
              ITD3=ITCAN-ITAC
              ITD4=ITAC-ITA
C              IF(ITA.GT.0.0) THEN
                  TAHIST(ITA+100)=TAHIST(ITA+100)+1.0
                  TCHIST(ITCAN+100)=TCHIST(ITCAN+100)+1.0
                  TSHIST(ITS+100)=TSHIST(ITS+100)+1.0
                  TACHIST(ITAC+100)=TACHIST(ITAC+100)+1.0
                  TDHIST(ITD+100)=TDHIST(ITD+100)+1.0
                  TD2HIST(ITD2+100)=TD2HIST(ITD2+100)+1.0
                  TD3HIST(ITD3+100)=TD3HIST(ITD3+100)+1.0
                  TD4HIST(ITD4+100)=TD4HIST(ITD4+100)+1.0
                  TSCRHIST(ITSCR+100)=TSCRHIST(ITSCR+100)+1.0
C              ENDIF
C          ENDIF     
          IF(FC(I).GT.0.1 .AND. RC(I).GT.1.0E5) NDRY=NDRY+1
!           IF((ITCAN-ITA).GE.10) THEN
!               WRITE(6,6070) IHOUR,IMIN,IDAY,IYEAR,FSSTAR,FLSTAR,QH,QE,
!      1                      BEG,TAROW(I)-TFREZ,TCN,TCN-(TAROW(I)-TFREZ),
!      2                      PAICAN(I),FSVF(I),UVROW(I),RC(I)
! 6070          FORMAT(2X,2I2,I4,I5,9F6.1,F6.3,F6.1,F8.1)
!           ENDIF
C
          IF(TSNOROT(I,M).GT.0.01) THEN
              TSN=TSNOROT(I,M)-TFREZ
          ELSE
              TSN=0.0
          ENDIF
          IF(TPNDROT(I,M).GT.0.01) THEN
              TPN=TPNDROT(I,M)-TFREZ
          ELSE
              TPN=0.0
          ENDIF
          GTOUT=GTROT(I,M)-TFREZ
          EVAPSUM=QFCFROT(I,M)+QFCLROT(I,M)+QFNROT(I,M)+QFGROT(I,M)+
     1                   QFCROT(I,M,1)+QFCROT(I,M,2)+QFCROT(I,M,3)
C
C===================== CTEM =====================================\
c         start writing output
c
          if ((iyear .ge. jhhsty) .and. (iyear .le. jhhendy)) then
           if ((iday .ge. jhhstd) .and. (iday .le. jhhendd)) then
C===================== CTEM =====================================/
          WRITE(64,6400) IHOUR,IMIN,IDAY,IYEAR,FSSTAR,FLSTAR,QH,QE,
     1                   SNOMLT,BEG,GTOUT,SNOROT(I,M),RHOSROT(I,M),
     2                   WSNOROT(I,M),ALTOT,ROFROT(I,M),
     3                   TPN,ZPNDROT(I,M),CDHROT(I,M),CDMROT(I,M), 
     4                   SFCUROT(I,M),SFCVROT(I,M),UVROW(I)
               open(unit=22,file='MEP_FLUXES.MET',status='unknown')
        WRITE(22,'(1X,I2,I3,I5,I6,4(3X,F15.5))') IHOUR,IMIN,IDAY,IYEAR,
     1                               RN(I),H_MEP(I),LE_MEP(I),G_MEP(I)
                                         
          IF(IGND.GT.3) THEN
C===================== CTEM =====================================\

              write(65,6500) ihour,imin,iday,iyear,(TBARROT(i,m,j)-
     1                   tfrez,THLQROT(i,m,j),THICROT(i,m,j),j=1,3),
     2                  tcn,RCANROT(i,m),SCANROT(i,m),tsn,zsn,
     3                   TCN-(TAROW(I)-TFREZ),TCANO(I)-TFREZ,
     4                   TACGAT(I)-TFREZ,ACTLYR,FTABLE,' TILE ',m
              write(66,6601) ihour,imin,iday,iyear,(TBARROT(i,m,j)-
     1                   tfrez,THLQROT(i,m,j),THICROT(i,m,j),j=4,10),
     2                   (GFLXROT(i,m,j),j=1,10),
     3                   ' TILE ',m
          else
              write(65,6500) ihour,imin,iday,iyear,(TBARROT(i,m,j)-
     1                   tfrez,THLQROT(i,m,j),THICROT(i,m,j),j=1,3),
     2                  tcn,RCANROT(i,m),SCANROT(i,m),tsn,zsn,
     3                   TCN-(TAROW(I)-TFREZ),TCANO(I)-TFREZ,
     4                   TACGAT(I)-TFREZ,ACTLYR,FTABLE,' TILE ',m

C===================== CTEM =====================================/
          ENDIF
C
          WRITE(67,6700) IHOUR,IMIN,IDAY,IYEAR,
     1                   TROFROT(I,M),TROOROT(I,M),TROSROT(I,M),
     2                   TROBROT(I,M),ROFROT(I,M),ROFOROT(I,M),
     3                   ROFSROT(I,M),ROFBROT(I,M),
     4                   FCS(M),FGS(M),FC(M),FG(M),' TILE ',M
          WRITE(68,6800) IHOUR,IMIN,IDAY,IYEAR,
     1                   FSGVROT(I,M),FSGSROT(I,M),FSGGROT(I,M),
     2                   FLGVROT(I,M),FLGSROT(I,M),FLGGROT(I,M),
     3                   HFSCROT(I,M),HFSSROT(I,M),HFSGROT(I,M),
     4                   HEVCROT(I,M),HEVSROT(I,M),HEVGROT(I,M),
     5                   HMFCROT(I,M),HMFNROT(I,M),
     6                   (HMFGROT(I,M,J),J=1,3),
     7                   HTCCROT(I,M),HTCSROT(I,M),
     8                   (HTCROT(I,M,J),J=1,3),' TILE ',M
          WRITE(69,6900) IHOUR,IMIN,IDAY,IYEAR,
     1                   PCFCROT(I,M),PCLCROT(I,M),PCPNROT(I,M),
     2                   PCPGROT(I,M),QFCFROT(I,M),QFCLROT(I,M),
     3                   QFNROT(I,M),QFGROT(I,M),(QFCROT(I,M,J),J=1,3),
     4                   ROFCROT(I,M),ROFNROT(I,M),ROFOROT(I,M),
     5                   ROFROT(I,M),WTRCROT(I,M),WTRSROT(I,M),
     6                   WTRGROT(I,M),' TILE ',M
C===================== CTEM =====================================\
C
         endif
        endif ! half hourly output loop.
c
c         Write half-hourly CTEM results to file *.CT01H
c
c         Net photosynthetic rates and leaf maintenance respiration for
c         each pft. however, if ctem_on then physyn subroutine
c         is using storage lai while actual lai is zero. if actual lai is
c         zero then we make anveg and rmlveg zero as well because these
c         are imaginary just like storage lai. note that anveg and rmlveg
c         are not passed to ctem. rather ancsveg, ancgveg, rmlcsveg, and
c         rmlcgveg are passed.
c
          if (ctem_on) then

            do 760 j = 1,icc
             if(ailcgrow(i,m,j).le.0.0) then
                anvegrow(i,m,j)=0.0
                rmlvegrow(i,m,j)=0.0
              else
                anvegrow(i,m,j)=ancsvegrow(i,m,j)*FSNOROT(i,m) +
     &                          ancgvegrow(i,m,j)*(1. - FSNOROT(i,m))
                rmlvegrow(i,m,j)=rmlcsvegrow(i,m,j)*FSNOROT(i,m) +
     &                         rmlcgvegrow(i,m,j)*(1. - FSNOROT(i,m))
              endif
760         continue
c
          if ((iyear .ge. jhhsty) .and. (iyear .le. jhhendy)) then
           if ((iday .ge. jhhstd) .and. (iday .le. jhhendd)) then

              write(71,7200)ihour,imin,iday,iyear,(anvegrow(i,m,j),
     1                    j=1,icc),(rmlvegrow(i,m,j),j=1,icc),
     2                    ' TILE ',m
            endif
           end if

           do j = 1,icc
              anvegrow_g(i,j)=anvegrow_g(i,j)+anvegrow(i,m,j)
     1                                        *FAREROT(i,m)
              rmlvegrow_g(i,j)=rmlvegrow_g(i,j)+rmlvegrow(i,m,j)
     1                                         *FAREROT(i,m)
            enddo

          endif   ! ctem_on

7200      format(1x,i2,1x,i2,i5,i5,9f11.3,9f11.3,2(a6,i2))
c
          fsstar_g    =fsstar_g + fsstar*FAREROT(i,m)
          flstar_g    =flstar_g + flstar*FAREROT(i,m)
          qh_g        =qh_g     + qh*FAREROT(i,m)
          qe_g        =qe_g     + qe*FAREROT(i,m)
          snomlt_g    =snomlt_g + snomlt*FAREROT(i,m)
          beg_g       =beg_g    + beg*FAREROT(i,m)
          gtout_g     =gtout_g  + gtout*FAREROT(i,m)
          tcn_g=tcn_g + tcn*FAREROT(i,m)
          tsn_g=tsn_g + tsn*FAREROT(i,m)
          zsn_g=zsn_g + zsn*FAREROT(i,m)
          altot_g     =altot_g   + altot*FAREROT(i,m)
          tpn_g       =tpn_g       + tpn*FAREROT(i,m)
c
          do j=1,ignd
            TBARROT_g(i,j)=TBARROT_g(i,j) + TBARROT(i,m,j)*FAREROT(i,m)
            THLQROT_g(i,j)=THLQROT_g(i,j) + THLQROT(i,m,j)*FAREROT(i,m)
            THICROT_g(i,j)=THICROT_g(i,j) + THICROT(i,m,j)*FAREROT(i,m)
            GFLXROT_g(i,j)=GFLXROT_g(i,j) + GFLXROT(i,m,j)*FAREROT(i,m)
            HMFGROT_g(i,j)=HMFGROT_g(i,j) + HMFGROT(i,m,j)*FAREROT(i,m)
            HTCROT_g(i,j)=HTCROT_g(i,j) + HTCROT(i,m,j)*FAREROT(i,m)
            QFCROT_g(i,j)=QFCROT_g(i,j) + QFCROT(i,m,j)*FAREROT(i,m)
          enddo
c
          ZPNDROT_g(i)=ZPNDROT_g(i) + ZPNDROT(i,m)*FAREROT(i,m)
          RHOSROT_g(i)=RHOSROT_g(i) + RHOSROT(i,m)*FAREROT(i,m)
          WSNOROT_g(i)=WSNOROT_g(i) + WSNOROT(i,m)*FAREROT(i,m)
          RCANROT_g(i)=RCANROT_g(i) + RCANROT(i,m)*FAREROT(i,m)
          SCANROT_g(i)=SCANROT_g(i) + SCANROT(i,m)*FAREROT(i,m)
          TROFROT_g(i)=TROFROT_g(i) + TROFROT(i,m)*FAREROT(i,m)
          TROOROT_g(i)=TROOROT_g(i) + TROOROT(i,m)*FAREROT(i,m)
          TROSROT_g(i)=TROSROT_g(i) + TROSROT(i,m)*FAREROT(i,m)
          TROBROT_g(i)=TROBROT_g(i) + TROBROT(i,m)*FAREROT(i,m)
          ROFOROT_g(i)=ROFOROT_g(i) + ROFOROT(i,m)*FAREROT(i,m)
          ROFSROT_g(i)=ROFSROT_g(i) + ROFSROT(i,m)*FAREROT(i,m)
          ROFBROT_g(i)=ROFBROT_g(i) + ROFBROT(i,m)*FAREROT(i,m)
          FSGVROT_g(i)=FSGVROT_g(i) + FSGVROT(i,m)*FAREROT(i,m)
          FSGSROT_g(i)=FSGSROT_g(i) + FSGSROT(i,m)*FAREROT(i,m)
          FSGGROT_g(i)=FSGGROT_g(i) + FSGGROT(i,m)*FAREROT(i,m)
          FLGVROT_g(i)=FLGVROT_g(i) + FLGVROT(i,m)*FAREROT(i,m)
          FLGSROT_g(i)=FLGSROT_g(i) + FLGSROT(i,m)*FAREROT(i,m)
          FLGGROT_g(i)=FLGGROT_g(i) + FLGGROT(i,m)*FAREROT(i,m)
          HFSCROT_g(i)=HFSCROT_g(i) + HFSCROT(i,m)*FAREROT(i,m)
          HFSSROT_g(i)=HFSSROT_g(i) + HFSSROT(i,m)*FAREROT(i,m)
          HFSGROT_g(i)=HFSGROT_g(i) + HFSGROT(i,m)*FAREROT(i,m)
          HEVCROT_g(i)=HEVCROT_g(i) + HEVCROT(i,m)*FAREROT(i,m)
          HEVSROT_g(i)=HEVSROT_g(i) + HEVSROT(i,m)*FAREROT(i,m)
          HEVGROT_g(i)=HEVGROT_g(i) + HEVGROT(i,m)*FAREROT(i,m)
          HMFCROT_g(i)=HMFCROT_g(i) + HMFCROT(i,m)*FAREROT(i,m)
          HMFNROT_g(i)=HMFNROT_g(i) + HMFNROT(i,m)*FAREROT(i,m)
          HTCCROT_g(i)=HTCCROT_g(i) + HTCCROT(i,m)*FAREROT(i,m)
          HTCSROT_g(i)=HTCSROT_g(i) + HTCSROT(i,m)*FAREROT(i,m)
          PCFCROT_g(i)=PCFCROT_g(i) + PCFCROT(i,m)*FAREROT(i,m)
          PCLCROT_g(i)=PCLCROT_g(i) + PCLCROT(i,m)*FAREROT(i,m)
          PCPNROT_g(i)=PCPNROT_g(i) + PCPNROT(i,m)*FAREROT(i,m)
          PCPGROT_g(i)=PCPGROT_g(i) + PCPGROT(i,m)*FAREROT(i,m)
          QFCFROT_g(i)=QFCFROT_g(i) + QFCFROT(i,m)*FAREROT(i,m)
          QFCLROT_g(i)=QFCLROT_g(i) + QFCLROT(i,m)*FAREROT(i,m)
          ROFCROT_g(i)=ROFCROT_g(i) + ROFCROT(i,m)*FAREROT(i,m)
          ROFNROT_g(i)=ROFNROT_g(i) + ROFNROT(i,m)*FAREROT(i,m)
          WTRCROT_g(i)=WTRCROT_g(i) + WTRCROT(i,m)*FAREROT(i,m)
          WTRSROT_g(i)=WTRSROT_g(i) + WTRSROT(i,m)*FAREROT(i,m)
          WTRGROT_g(i)=WTRGROT_g(i) + WTRGROT(i,m)*FAREROT(i,m)
          QFNROT_g(i) =QFNROT_g(i) + QFNROT(i,m)*FAREROT(i,m)
          QFGROT_g(i) =QFGROT_g(i) + QFGROT(i,m)*FAREROT(i,m)
          ROFROT_g(i) =ROFROT_g(i) + ROFROT(i,m)*FAREROT(i,m)
          SNOROT_g(i) =SNOROT_g(i) + SNOROT(i,m)*FAREROT(i,m)
          CDHROT_g(i) =CDHROT_g(i) + CDHROT(i,m)*FAREROT(i,m)
          CDMROT_g(i) =CDMROT_g(i) + CDMROT(i,m)*FAREROT(i,m)
          SFCUROT_g(i) =SFCUROT_g(i) + SFCUROT(i,m)*FAREROT(i,m)
          SFCVROT_g(i) =SFCVROT_g(i) + SFCVROT(i,m)*FAREROT(i,m)
C
C======================== CTEM =====================================/
425    CONTINUE
C===================== CTEM =====================================\
C      WRITE CTEM OUTPUT FILES
C
      if ((iyear .ge. jhhsty) .and. (iyear .le. jhhendy)) then
       if ((iday .ge. jhhstd) .and. (iday .le. jhhendd)) then

       IF (CTEM_ON) THEN
           WRITE(711,7200)IHOUR,IMIN,IDAY,IYEAR,(ANVEGROW_G(I,J),
     1                 J=1,ICC),(RMLVEGROW_G(I,J),J=1,ICC)
       ENDIF !CTEM_ON

       WRITE(641,6400) IHOUR,IMIN,IDAY,IYEAR,FSSTAR_G,FLSTAR_G,QH_G,
     1      QE_G,SNOMLT_G,BEG_G,GTOUT_G,SNOROT_G(I),RHOSROT_G(I),
     2                   WSNOROT_G(I),ALTOT_G,ROFROT_G(I),
     3                   TPN_G,ZPNDROT_G(I),CDHROT_G(I),CDMROT_G(I),
     4                   SFCUROT_G(I),SFCVROT_G(I),UVROW(I)
         WRITE(651,6500) IHOUR,IMIN,IDAY,IYEAR,(TBARROT_G(I,J)-
     1                   TFREZ,THLQROT_G(I,J),THICROT_G(I,J),J=1,3),
     2                   TCN_G,RCANROT_G(I),SCANROT_G(I),TSN_G,ZSN_G,
     3                   TCN_G-(TAROW(I)-TFREZ),TCANO(I)-TFREZ,
     4                   TACGAT(I)-TFREZ,ACTLYR,FTABLE
C
         IF(IGND.GT.3) THEN
          WRITE(661,6601) IHOUR,IMIN,IDAY,IYEAR,(TBARROT_G(I,J)-
     1                   TFREZ,THLQROT_G(I,J),THICROT_G(I,J),J=4,10),
     2                   (GFLXROT_G(I,J),J=1,10)
         ELSE
          WRITE(661,6600) IHOUR,IMIN,IDAY,FSSROW(I),FDLROW(I),PREROW(I),
     1                   TAROW(I)-TFREZ,UVROW(I),PRESROW(I),QAROW(I)
         ENDIF
C
         WRITE(671,6700) IHOUR,IMIN,IDAY,IYEAR,
     &                   TROFROT_G(I),TROOROT_G(I),TROSROT_G(I),
     1                   TROBROT_G(I),ROFROT_G(I),ROFOROT_G(I),
     2                   ROFSROT_G(I),ROFBROT_G(I),
     3                   FCS_G(I),FGS_G(I),FC_G(I),FG_G(I)
         WRITE(681,6800) IHOUR,IMIN,IDAY,IYEAR,
     &                   FSGVROT_G(I),FSGSROT_G(I),FSGGROT_G(I),
     1                   FLGVROT_G(I),FLGSROT_G(I),FLGGROT_G(I),
     2                   HFSCROT_G(I),HFSSROT_G(I),HFSGROT_G(I),
     3                   HEVCROT_G(I),HEVSROT_G(I),HEVGROT_G(I),
     4                   HMFCROT_G(I),HMFNROT_G(I),
     5                   (HMFGROT_G(I,J),J=1,3),
     6                   HTCCROT_G(I),HTCSROT_G(I),
     7                   (HTCROT_G(I,J),J=1,3)
         WRITE(691,6900) IHOUR,IMIN,IDAY,IYEAR,
     &                   PCFCROT_G(I),PCLCROT_G(I),PCPNROT_G(I),
     1                   PCPGROT_G(I),QFCFROT_G(I),QFCLROT_G(I),
     2                   QFNROT_G(I),QFGROT_G(I),(QFCROT_G(I,J),J=1,3),
     3                   ROFCROT_G(I),ROFNROT_G(I),ROFOROT_G(I),
     4                   ROFROT_G(I),WTRCROT_G(I),WTRSROT_G(I),
     5                   WTRGROT_G(I)
C
        endif
       ENDIF ! if write half-hourly
C===================== CTEM =====================================/
450   CONTINUE
C
C===================== CTEM =====================================\

      endif ! not parallelrun
C===================== CTEM =====================================/
C
C=======================================================================
C     * CALCULATE GRID CELL AVERAGE DIAGNOSTIC FIELDS.
C
C===================== CTEM =====================================\

      if(.not.parallelrun) then ! stand alone mode, includes
c                               ! diagnostic fields
C===================== CTEM =====================================/
C
      DO 525 I=1,NLTEST
          CDHROW(I)=0.
          CDMROW(I)=0.
          HFSROW(I)=0.
          TFXROW(I)=0.
          QEVPROW(I)=0.
          QFSROW(I)=0.
          QFXROW(I)=0.
          PETROW(I)=0.
          GAROW(I)=0.
          EFROW(I)=0.
          GTROW(I)=0.
          QGROW(I)=0.
          ALVSROW(I)=0.
          ALIRROW(I)=0.
          SFCTROW(I)=0.
          SFCUROW(I)=0.
          SFCVROW(I)=0.
          SFCQROW(I)=0.
          SFRHROW(I)=0.
          FSNOROW(I)=0.
          FSGVROW(I)=0.
          FSGSROW(I)=0.
          FSGGROW(I)=0.
          FLGVROW(I)=0.
          FLGSROW(I)=0.
          FLGGROW(I)=0.
          HFSCROW(I)=0.
          HFSSROW(I)=0.
          HFSGROW(I)=0.
          HEVCROW(I)=0.
          HEVSROW(I)=0.
          HEVGROW(I)=0.
          HMFCROW(I)=0.
          HMFNROW(I)=0.
          HTCCROW(I)=0.
          HTCSROW(I)=0.
          PCFCROW(I)=0.
          PCLCROW(I)=0.
          PCPNROW(I)=0.
          PCPGROW(I)=0.
          QFGROW(I)=0.
          QFNROW(I)=0.
          QFCLROW(I)=0.
          QFCFROW(I)=0.
          ROFROW(I)=0.
          ROFOROW(I)=0.
          ROFSROW(I)=0.
          ROFBROW(I)=0.
          ROFCROW(I)=0.
          ROFNROW(I)=0.
          ROVGROW(I)=0.
          WTRCROW(I)=0.
          WTRSROW(I)=0.
          WTRGROW(I)=0.
          DRROW(I)=0.
          WTABROW(I)=0.
          ILMOROW(I)=0.
          UEROW(I)=0.
          HBLROW(I)=0.
          DO 500 J=1,IGND
              HMFGROW(I,J)=0.
              HTCROW(I,J)=0.
              QFCROW(I,J)=0.
              GFLXROW(I,J)=0.
500       CONTINUE
525   CONTINUE
C
      DO 600 I=1,NLTEST
      DO 575 M=1,NMTEST
          CDHROW(I)=CDHROW(I)+CDHROT(I,M)*FAREROT(I,M)
          CDMROW(I)=CDMROW(I)+CDMROT(I,M)*FAREROT(I,M)
          HFSROW(I)=HFSROW(I)+HFSROT(I,M)*FAREROT(I,M)
          TFXROW(I)=TFXROW(I)+TFXROT(I,M)*FAREROT(I,M)
          QEVPROW(I)=QEVPROW(I)+QEVPROT(I,M)*FAREROT(I,M)
          QFSROW(I)=QFSROW(I)+QFSROT(I,M)*FAREROT(I,M)
          QFXROW(I)=QFXROW(I)+QFXROT(I,M)*FAREROT(I,M)
          PETROW(I)=PETROW(I)+PETROT(I,M)*FAREROT(I,M)
          GAROW(I)=GAROW(I)+GAROT(I,M)*FAREROT(I,M)
          EFROW(I)=EFROW(I)+EFROT(I,M)*FAREROT(I,M)
          GTROW(I)=GTROW(I)+GTROT(I,M)*FAREROT(I,M)
          QGROW(I)=QGROW(I)+QGROT(I,M)*FAREROT(I,M)
          ALVSROW(I)=ALVSROW(I)+ALVSROT(I,M)*FAREROT(I,M)
          ALIRROW(I)=ALIRROW(I)+ALIRROT(I,M)*FAREROT(I,M)
          SFCTROW(I)=SFCTROW(I)+SFCTROT(I,M)*FAREROT(I,M)
          SFCUROW(I)=SFCUROW(I)+SFCUROT(I,M)*FAREROT(I,M)
          SFCVROW(I)=SFCVROW(I)+SFCVROT(I,M)*FAREROT(I,M)
          SFCQROW(I)=SFCQROW(I)+SFCQROT(I,M)*FAREROT(I,M)
          SFRHROW(I)=SFRHROW(I)+SFRHROT(I,M)*FAREROT(I,M)
          FSNOROW(I)=FSNOROW(I)+FSNOROT(I,M)*FAREROT(I,M)
          FSGVROW(I)=FSGVROW(I)+FSGVROT(I,M)*FAREROT(I,M)
          FSGSROW(I)=FSGSROW(I)+FSGSROT(I,M)*FAREROT(I,M)
          FSGGROW(I)=FSGGROW(I)+FSGGROT(I,M)*FAREROT(I,M)
          FLGVROW(I)=FLGVROW(I)+FLGVROT(I,M)*FAREROT(I,M)
          FLGSROW(I)=FLGSROW(I)+FLGSROT(I,M)*FAREROT(I,M)
          FLGGROW(I)=FLGGROW(I)+FLGGROT(I,M)*FAREROT(I,M)
          HFSCROW(I)=HFSCROW(I)+HFSCROT(I,M)*FAREROT(I,M)
          HFSSROW(I)=HFSSROW(I)+HFSSROT(I,M)*FAREROT(I,M)
          HFSGROW(I)=HFSGROW(I)+HFSGROT(I,M)*FAREROT(I,M)
          HEVCROW(I)=HEVCROW(I)+HEVCROT(I,M)*FAREROT(I,M)
          HEVSROW(I)=HEVSROW(I)+HEVSROT(I,M)*FAREROT(I,M)
          HEVGROW(I)=HEVGROW(I)+HEVGROT(I,M)*FAREROT(I,M)
          HMFCROW(I)=HMFCROW(I)+HMFCROT(I,M)*FAREROT(I,M)
          HMFNROW(I)=HMFNROW(I)+HMFNROT(I,M)*FAREROT(I,M)
          HTCCROW(I)=HTCCROW(I)+HTCCROT(I,M)*FAREROT(I,M)
          HTCSROW(I)=HTCSROW(I)+HTCSROT(I,M)*FAREROT(I,M)
          PCFCROW(I)=PCFCROW(I)+PCFCROT(I,M)*FAREROT(I,M)
          PCLCROW(I)=PCLCROW(I)+PCLCROT(I,M)*FAREROT(I,M)
          PCPNROW(I)=PCPNROW(I)+PCPNROT(I,M)*FAREROT(I,M)
          PCPGROW(I)=PCPGROW(I)+PCPGROT(I,M)*FAREROT(I,M)
          QFGROW(I)=QFGROW(I)+QFGROT(I,M)*FAREROT(I,M)
          QFNROW(I)=QFNROW(I)+QFNROT(I,M)*FAREROT(I,M)
          QFCLROW(I)=QFCLROW(I)+QFCLROT(I,M)*FAREROT(I,M)
          QFCFROW(I)=QFCFROW(I)+QFCFROT(I,M)*FAREROT(I,M)
          ROFROW(I)=ROFROW(I)+ROFROT(I,M)*FAREROT(I,M)
          ROFOROW(I)=ROFOROW(I)+ROFOROT(I,M)*FAREROT(I,M)
          ROFSROW(I)=ROFSROW(I)+ROFSROT(I,M)*FAREROT(I,M)
          ROFBROW(I)=ROFBROW(I)+ROFBROT(I,M)*FAREROT(I,M)
          ROFCROW(I)=ROFCROW(I)+ROFCROT(I,M)*FAREROT(I,M)
          ROFNROW(I)=ROFNROW(I)+ROFNROT(I,M)*FAREROT(I,M)
          ROVGROW(I)=ROVGROW(I)+ROVGROT(I,M)*FAREROT(I,M)
          WTRCROW(I)=WTRCROW(I)+WTRCROT(I,M)*FAREROT(I,M)
          WTRSROW(I)=WTRSROW(I)+WTRSROT(I,M)*FAREROT(I,M)
          WTRGROW(I)=WTRGROW(I)+WTRGROT(I,M)*FAREROT(I,M)
          DRROW(I)=DRROW(I)+DRROT(I,M)*FAREROT(I,M)
          WTABROW(I)=WTABROW(I)+WTABROT(I,M)*FAREROT(I,M)
          ILMOROW(I)=ILMOROW(I)+ILMOROT(I,M)*FAREROT(I,M)
          UEROW(I)=UEROW(I)+UEROT(I,M)*FAREROT(I,M)
          HBLROW(I)=HBLROW(I)+HBLROT(I,M)*FAREROT(I,M)
          DO 550 J=1,IGND
              HMFGROW(I,J)=HMFGROW(I,J)+HMFGROT(I,M,J)*FAREROT(I,M)
              HTCROW(I,J)=HTCROW(I,J)+HTCROT(I,M,J)*FAREROT(I,M)
              QFCROW(I,J)=QFCROW(I,J)+QFCROT(I,M,J)*FAREROT(I,M)
              GFLXROW(I,J)=GFLXROW(I,J)+GFLXROT(I,M,J)*FAREROT(I,M)
550       CONTINUE
575   CONTINUE
600   CONTINUE
C
C===================== CTEM =====================================\

      endif ! not parallelrun, for diagnostic fields
c
      if(.not.parallelrun) then ! stand alone mode, includes daily output for class
C===================== CTEM =====================================/
C
C     * ACCUMULATE OUTPUT DATA FOR DIURNALLY AVERAGED FIELDS. BOTH GRID
C       MEAN AND MOSAIC MEAN
C
      DO 675 I=1,NLTEST
      DO 650 M=1,NMTEST
          PREACC(I)=PREACC(I)+PREROW(I)*FAREROT(I,M)*DELT
          GTACC(I)=GTACC(I)+GTROT(I,M)*FAREROT(I,M)
          QEVPACC(I)=QEVPACC(I)+QEVPROT(I,M)*FAREROT(I,M)
          EVAPACC(I)=EVAPACC(I)+QFSROT(I,M)*FAREROT(I,M)*DELT
          HFSACC(I)=HFSACC(I)+HFSROT(I,M)*FAREROT(I,M)
          HMFNACC(I)=HMFNACC(I)+HMFNROT(I,M)*FAREROT(I,M)
          ROFACC(I)=ROFACC(I)+ROFROT(I,M)*FAREROT(I,M)*DELT
          OVRACC(I)=OVRACC(I)+ROFOROT(I,M)*FAREROT(I,M)*DELT
          WTBLACC(I)=WTBLACC(I)+WTABROT(I,M)*FAREROT(I,M)
          DO 625 J=1,IGND
              TBARACC(I,J)=TBARACC(I,J)+TBARROT(I,M,J)*FAREROT(I,M)
              THLQACC(I,J)=THLQACC(I,J)+THLQROT(I,M,J)*FAREROT(I,M)
              THICACC(I,J)=THICACC(I,J)+THICROT(I,M,J)*FAREROT(I,M)
              THALACC(I,J)=THALACC(I,J)+(THLQROT(I,M,J)+THICROT(I,M,J))
     1                    *FAREROT(I,M)
625       CONTINUE
          ALVSACC(I)=ALVSACC(I)+ALVSROT(I,M)*FAREROT(I,M)*FSVHROW(I)
          ALIRACC(I)=ALIRACC(I)+ALIRROT(I,M)*FAREROT(I,M)*FSIHROW(I)
          IF(SNOROT(I,M).GT.0.0) THEN
              RHOSACC(I)=RHOSACC(I)+RHOSROT(I,M)*FAREROT(I,M)
              TSNOACC(I)=TSNOACC(I)+TSNOROT(I,M)*FAREROT(I,M)
              WSNOACC(I)=WSNOACC(I)+WSNOROT(I,M)*FAREROT(I,M)
              SNOARE(I)=SNOARE(I)+FAREROT(I,M)
          ENDIF
          IF(TCANROT(I,M).GT.0.5) THEN
              TCANACC(I)=TCANACC(I)+TCANROT(I,M)*FAREROT(I,M)
              CANARE(I)=CANARE(I)+FAREROT(I,M)
          ENDIF
          SNOACC(I)=SNOACC(I)+SNOROT(I,M)*FAREROT(I,M)
          RCANACC(I)=RCANACC(I)+RCANROT(I,M)*FAREROT(I,M)
          SCANACC(I)=SCANACC(I)+SCANROT(I,M)*FAREROT(I,M)
          GROACC(I)=GROACC(I)+GROROT(I,M)*FAREROT(I,M)
          FSINACC(I)=FSINACC(I)+FSSROW(I)*FAREROT(I,M)
          FLINACC(I)=FLINACC(I)+FDLROW(I)*FAREROT(I,M)
          FLUTACC(I)=FLUTACC(I)+SBC*GTROT(I,M)**4*FAREROT(I,M)
          TAACC(I)=TAACC(I)+TAROW(I)*FAREROT(I,M)
          UVACC(I)=UVACC(I)+UVROW(I)*FAREROT(I,M)
          PRESACC(I)=PRESACC(I)+PRESROW(I)*FAREROT(I,M)
          QAACC(I)=QAACC(I)+QAROW(I)*FAREROT(I,M)
650   CONTINUE
675   CONTINUE
C
C     * CALCULATE AND PRINT DAILY AVERAGES.
C
      IF(NCOUNT.EQ.NDAY) THEN

      DO 800 I=1,NLTEST
          PREACC(I)=PREACC(I)
          GTACC(I)=GTACC(I)/REAL(NDAY)
          QEVPACC(I)=QEVPACC(I)/REAL(NDAY)
          EVAPACC(I)=EVAPACC(I)
          HFSACC(I)=HFSACC(I)/REAL(NDAY)
          HMFNACC(I)=HMFNACC(I)/REAL(NDAY)
          ROFACC(I)=ROFACC(I)
          OVRACC(I)=OVRACC(I)
          WTBLACC(I)=WTBLACC(I)/REAL(NDAY)
          DO 725 J=1,IGND
              TBARACC(I,J)=TBARACC(I,J)/REAL(NDAY)
              THLQACC(I,J)=THLQACC(I,J)/REAL(NDAY)
              THICACC(I,J)=THICACC(I,J)/REAL(NDAY)
              THALACC(I,J)=THALACC(I,J)/REAL(NDAY)
725       CONTINUE
          IF(FSINACC(I).GT.0.0) THEN
              ALVSACC(I)=ALVSACC(I)/(FSINACC(I)*0.5)
              ALIRACC(I)=ALIRACC(I)/(FSINACC(I)*0.5)
          ELSE
              ALVSACC(I)=0.0
              ALIRACC(I)=0.0
          ENDIF
          IF(SNOARE(I).GT.0.0) THEN
              RHOSACC(I)=RHOSACC(I)/SNOARE(I)
              TSNOACC(I)=TSNOACC(I)/SNOARE(I)
              WSNOACC(I)=WSNOACC(I)/SNOARE(I)
          ENDIF
          IF(CANARE(I).GT.0.0) THEN
              TCANACC(I)=TCANACC(I)/CANARE(I)
          ENDIF
          SNOACC(I)=SNOACC(I)/REAL(NDAY)
          RCANACC(I)=RCANACC(I)/REAL(NDAY)
          SCANACC(I)=SCANACC(I)/REAL(NDAY)
          GROACC(I)=GROACC(I)/REAL(NDAY)
          FSINACC(I)=FSINACC(I)/REAL(NDAY)
          FLINACC(I)=FLINACC(I)/REAL(NDAY)
          FLUTACC(I)=FLUTACC(I)/REAL(NDAY)
          TAACC(I)=TAACC(I)/REAL(NDAY)
          UVACC(I)=UVACC(I)/REAL(NDAY)
          PRESACC(I)=PRESACC(I)/REAL(NDAY)
          QAACC(I)=QAACC(I)/REAL(NDAY)

              ALTOT=(ALVSACC(I)+ALIRACC(I))/2.0
              FSSTAR=FSINACC(I)*(1.-ALTOT)
              FLSTAR=FLINACC(I)-FLUTACC(I)
              QH=HFSACC(I)
              QE=QEVPACC(I)
              BEG=FSSTAR+FLSTAR-QH-QE
              SNOMLT=HMFNACC(I)
              IF(RHOSACC(I).GT.0.0) THEN
                  ZSN=SNOACC(I)/RHOSACC(I)
              ELSE
                  ZSN=0.0
              ENDIF
              IF(TCANACC(I).GT.0.01) THEN
                  TCN=TCANACC(I)-TFREZ
              ELSE
                  TCN=0.0
              ENDIF
              IF(TSNOACC(I).GT.0.01) THEN
                  TSN=TSNOACC(I)-TFREZ
              ELSE
                  TSN=0.0
              ENDIF
              GTOUT=GTACC(I)-TFREZ
C
             if ((iyear .ge. jdsty) .and. (iyear .le. jdendy)) then
              if ((iday .ge. jdstd) .and. (iday .le. jdendd)) then

              WRITE(61,6100) IDAY,IYEAR,FSSTAR,FLSTAR,QH,QE,SNOMLT,
     1                       BEG,GTOUT,SNOACC(I),RHOSACC(I),
     2                       WSNOACC(I),ALTOT,ROFACC(I),CUMSNO
              IF(IGND.GT.3) THEN
                  WRITE(62,6201) IDAY,IYEAR,(TBARACC(I,J)-TFREZ,
     1                       THLQACC(I,J),THICACC(I,J),J=1,5)
                  WRITE(63,6201) IDAY,IYEAR,(TBARACC(I,J)-TFREZ,
     1                       THLQACC(I,J),THICACC(I,J),J=6,10)
              ELSE
                  WRITE(62,6200) IDAY,IYEAR,(TBARACC(I,J)-TFREZ,
     1                       THLQACC(I,J),THICACC(I,J),J=1,3),
     2                       TCN,RCANACC(I),SCANACC(I),TSN,ZSN
                  WRITE(63,6300) IDAY,IYEAR,FSINACC(I),FLINACC(I),
     1                       TAACC(I)-TFREZ,UVACC(I),PRESACC(I),
     2                       QAACC(I),PREACC(I),EVAPACC(I)
              ENDIF
             endif
            ENDIF
C
C     * RESET ACCUMULATOR ARRAYS.
C
          PREACC(I)=0.
          GTACC(I)=0.
          QEVPACC(I)=0.
          HFSACC(I)=0.
          HMFNACC(I)=0.
          ROFACC(I)=0.
          SNOACC(I)=0.
          CANARE(I)=0.
          SNOARE(I)=0.
          OVRACC(I)=0.
          WTBLACC(I)=0.
          DO 750 J=1,IGND
              TBARACC(I,J)=0.
              THLQACC(I,J)=0.
              THICACC(I,J)=0.
              THALACC(I,J)=0.
750       CONTINUE
          ALVSACC(I)=0.
          ALIRACC(I)=0.
          RHOSACC(I)=0.
          TSNOACC(I)=0.
          WSNOACC(I)=0.
          TCANACC(I)=0.
          RCANACC(I)=0.
          SCANACC(I)=0.
          GROACC(I)=0.
          FSINACC(I)=0.
          FLINACC(I)=0.
          TAACC(I)=0.
          UVACC(I)=0.
          PRESACC(I)=0.
          QAACC(I)=0.
          EVAPACC(I)=0.
          FLUTACC(I)=0.
800   CONTINUE

      ENDIF ! IF(NCOUNT.EQ.NDAY)

C===================== CTEM =====================================\
C
C     CALCULATE AND PRINT MOSAIC DAILY AVERAGES.
C
!       start -> FLAG JM
      DO 676 I=1,NLTEST
      DO 658 M=1,NMTEST
          PREACC_M(I,M)=PREACC_M(I,M)+PREROW(I)*DELT
          GTACC_M(I,M)=GTACC_M(I,M)+GTROT(I,M)
          QEVPACC_M(I,M)=QEVPACC_M(I,M)+QEVPROT(I,M)
          EVAPACC_M(I,M)=EVAPACC_M(I,M)+QFSROT(I,M)*DELT
          HFSACC_M(I,M)=HFSACC_M(I,M)+HFSROT(I,M)
          HMFNACC_M(I,M)=HMFNACC_M(I,M)+HMFNROT(I,M)
          ROFACC_M(I,M)=ROFACC_M(I,M)+ROFROT(I,M)*DELT
          OVRACC_M(I,M)=OVRACC_M(I,M)+ROFOROT(I,M)*DELT
          WTBLACC_M(I,M)=WTBLACC_M(I,M)+WTABROT(I,M)
          DO 626 J=1,IGND
              TBARACC_M(I,M,J)=TBARACC_M(I,M,J)+TBARROT(I,M,J)
              THLQACC_M(I,M,J)=THLQACC_M(I,M,J)+THLQROT(I,M,J)
              THICACC_M(I,M,J)=THICACC_M(I,M,J)+THICROT(I,M,J)
              THALACC_M(I,M,J)=THALACC_M(I,M,J)+(THLQROT(I,M,J)+
     1           THICROT(I,M,J))
626       CONTINUE
          ALVSACC_M(I,M)=ALVSACC_M(I,M)+ALVSROT(I,M)*FSVHROW(I)
          ALIRACC_M(I,M)=ALIRACC_M(I,M)+ALIRROT(I,M)*FSIHROW(I)
          IF(SNOROT(I,M).GT.0.0) THEN
              RHOSACC_M(I,M)=RHOSACC_M(I,M)+RHOSROT(I,M)
              TSNOACC_M(I,M)=TSNOACC_M(I,M)+TSNOROT(I,M)
              WSNOACC_M(I,M)=WSNOACC_M(I,M)+WSNOROT(I,M)
              SNOARE_M(I,M) = SNOARE_M(I,M) + 1.0 !FLAG test.
          ENDIF
          IF(TCANROT(I,M).GT.0.5) THEN
              TCANACC_M(I,M)=TCANACC_M(I,M)+TCANROT(I,M)
C              CANARE(I)=CANARE(I)+FAREROT(I,M)
          ENDIF
          SNOACC_M(I,M)=SNOACC_M(I,M)+SNOROT(I,M)
          RCANACC_M(I,M)=RCANACC_M(I,M)+RCANROT(I,M)
          SCANACC_M(I,M)=SCANACC_M(I,M)+SCANROT(I,M)
          GROACC_M(I,M)=GROACC_M(I,M)+GROROT(I,M)
          FSINACC_M(I,M)=FSINACC_M(I,M)+FSSROW(I)
          FLINACC_M(I,M)=FLINACC_M(I,M)+FDLROW(I)
          FLUTACC_M(I,M)=FLUTACC_M(I,M)+SBC*GTROT(I,M)**4
          TAACC_M(I,M)=TAACC_M(I,M)+TAROW(I)
          UVACC_M(I,M)=UVACC_M(I,M)+UVROW(I)
          PRESACC_M(I,M)=PRESACC_M(I,M)+PRESROW(I)
          QAACC_M(I,M)=QAACC_M(I,M)+QAROW(I)
658   CONTINUE
676   CONTINUE
C
C     CALCULATE AND PRINT DAILY AVERAGES.
C
      IF(NCOUNT.EQ.NDAY) THEN

      DO 808 I=1,NLTEST
        DO 809 M=1,NMTEST
          PREACC_M(I,M)=PREACC_M(I,M)     !became [kg m-2 day-1] instead of [kg m-2 s-1]
          GTACC_M(I,M)=GTACC_M(I,M)/REAL(NDAY)
          QEVPACC_M(I,M)=QEVPACC_M(I,M)/REAL(NDAY)
          EVAPACC_M(I,M)=EVAPACC_M(I,M)   !became [kg m-2 day-1] instead of [kg m-2 s-1]
          HFSACC_M(I,M)=HFSACC_M(I,M)/REAL(NDAY)
          HMFNACC_M(I,M)=HMFNACC_M(I,M)/REAL(NDAY)
          ROFACC_M(I,M)=ROFACC_M(I,M)   !became [kg m-2 day-1] instead of [kg m-2 s-1
          OVRACC_M(I,M)=OVRACC_M(I,M)   !became [kg m-2 day-1] instead of [kg m-2 s-1]
          WTBLACC_M(I,M)=WTBLACC_M(I,M)/REAL(NDAY)
          DO 726 J=1,IGND
            TBARACC_M(I,M,J)=TBARACC_M(I,M,J)/REAL(NDAY)
            THLQACC_M(I,M,J)=THLQACC_M(I,M,J)/REAL(NDAY)
            THICACC_M(I,M,J)=THICACC_M(I,M,J)/REAL(NDAY)
            THALACC_M(I,M,J)=THALACC_M(I,M,J)/REAL(NDAY)
726       CONTINUE
C
          IF(FSINACC_M(I,M).GT.0.0) THEN
            ALVSACC_M(I,M)=ALVSACC_M(I,M)/(FSINACC_M(I,M)*0.5)
            ALIRACC_M(I,M)=ALIRACC_M(I,M)/(FSINACC_M(I,M)*0.5)
          ELSE
            ALVSACC_M(I,M)=0.0
            ALIRACC_M(I,M)=0.0
          ENDIF
C
          SNOACC_M(I,M)=SNOACC_M(I,M)/REAL(NDAY)
          if (SNOARE_M(I,M) .GT. 0.) THEN
             RHOSACC_M(I,M)=RHOSACC_M(I,M)/SNOARE_M(I,M)
             TSNOACC_M(I,M)=TSNOACC_M(I,M)/SNOARE_M(I,M)
             WSNOACC_M(I,M)=WSNOACC_M(I,M)/SNOARE_M(I,M)
          END IF
          TCANACC_M(I,M)=TCANACC_M(I,M)/REAL(NDAY)
          RCANACC_M(I,M)=RCANACC_M(I,M)/REAL(NDAY)
          SCANACC_M(I,M)=SCANACC_M(I,M)/REAL(NDAY)
          GROACC_M(I,M)=GROACC_M(I,M)/REAL(NDAY)
          FSINACC_M(I,M)=FSINACC_M(I,M)/REAL(NDAY)
          FLINACC_M(I,M)=FLINACC_M(I,M)/REAL(NDAY)
          FLUTACC_M(I,M)=FLUTACC_M(I,M)/REAL(NDAY)
          TAACC_M(I,M)=TAACC_M(I,M)/REAL(NDAY)
          UVACC_M(I,M)=UVACC_M(I,M)/REAL(NDAY)
          PRESACC_M(I,M)=PRESACC_M(I,M)/REAL(NDAY)
          QAACC_M(I,M)=QAACC_M(I,M)/REAL(NDAY)
          ALTOT=(ALVSACC_M(I,M)+ALIRACC_M(I,M))/2.0
          FSSTAR=FSINACC_M(I,M)*(1.-ALTOT)
          FLSTAR=FLINACC_M(I,M)-FLUTACC_M(I,M)
          QH=HFSACC_M(I,M)
          QE=QEVPACC_M(I,M)
          QEVPACC_M_SAVE(I,M)=QEVPACC_M(I,M)   !FLAG!! What is the point of this? JM Apr 1 2015
          BEG=FSSTAR+FLSTAR-QH-QE
          SNOMLT=HMFNACC_M(I,M)
C
          IF(RHOSACC_M(I,M).GT.0.0) THEN
              ZSN=SNOACC_M(I,M)/RHOSACC_M(I,M)
          ELSE
              ZSN=0.0
          ENDIF
C
          IF(TCANACC_M(I,M).GT.0.01) THEN
              TCN=TCANACC_M(I,M)-TFREZ
          ELSE
              TCN=0.0
          ENDIF
C
          IF(TSNOACC_M(I,M).GT.0.01) THEN
              TSN=TSNOACC_M(I,M)-TFREZ
          ELSE
              TSN=0.0
          ENDIF
C
          GTOUT=GTACC_M(I,M)-TFREZ
C 
          if ((iyear .ge. jdsty) .and. (iyear .le. jdendy)) then
           if ((iday .ge. jdstd) .and. (iday .le. jdendd)) then



C
C         WRITE TO OUTPUT FILES
C
          WRITE(611,6100) IDAY,IYEAR,FSSTAR,FLSTAR,QH,QE,SNOMLT,
     1                    BEG,GTOUT,SNOACC_M(I,M),RHOSACC_M(I,M),
     2                    WSNOACC_M(I,M),ALTOT,ROFACC_M(I,M),
     3                    CUMSNO,' TILE ',M
            IF(IGND.GT.3) THEN
               WRITE(621,6201) IDAY,IYEAR,(TBARACC_M(I,M,J)-TFREZ,
     1                  THLQACC_M(I,M,J),THICACC_M(I,M,J),J=1,5),
     2                  ' TILE ',M
               WRITE(631,6201) IDAY,IYEAR,(TBARACC_M(I,M,J)-TFREZ,
     1                  THLQACC_M(I,M,J),THICACC_M(I,M,J),J=6,10),
     2                  ' TILE ',M
            ELSE
               WRITE(621,6200) IDAY,IYEAR,(TBARACC_M(I,M,J)-TFREZ,
     1                  THLQACC_M(I,M,J),THICACC_M(I,M,J),J=1,3),
     2                  TCN,RCANACC_M(I,M),SCANACC_M(I,M),TSN,ZSN,
     3                  ' TILE ',M
               WRITE(631,6300) IDAY,IYEAR,FSINACC_M(I,M),FLINACC_M(I,M),
     1                  TAACC_M(I,M)-TFREZ,UVACC_M(I,M),PRESACC_M(I,M),
     2                  QAACC_M(I,M),PREACC_M(I,M),EVAPACC_M(I,M),
     3                  ' TILE ',M
            ENDIF
C
           endif
          ENDIF ! IF write daily
C
C          INITIALIZTION FOR MOSAIC TILE AND GRID VARIABLES
C
            call resetclassaccum(nltest,nmtest)
C
809   CONTINUE
808   CONTINUE

      ENDIF ! IF(NCOUNT.EQ.NDAY)
C
      ENDIF !  IF(.NOT.PARALLELRUN)
C
C=======================================================================
C
!       CTEM--------------\
!     Only bother with monthly calculations if we desire those outputs to be written out.
      if (iyear .ge. jmosty) then

        call class_monthly_aw(IDAY,IYEAR,NCOUNT,NDAY,SBC,DELT,
     1                       nltest,nmtest,ALVSROT,FAREROT,FSVHROW,
     2                       ALIRROT,FSIHROW,GTROT,FSSROW,FDLROW,
     3                       HFSROT,ROFROT,PREROW,QFSROT,QEVPROT,
     4                       SNOROT,TAROW,WSNOROT,TBARROT,THLQROT,
     5                       THICROT,TFREZ)
!       CTEM--------------/

       DO NT=1,NMON
        IF(IDAY.EQ.monthend(NT+1).AND.NCOUNT.EQ.NDAY)THEN
         IMONTH=NT
        ENDIF ! IF(IDAY.EQ.monthend(NT+1).AND.NCOUNT.EQ.NDAY)
       ENDDO ! NMON

!       CTEM--------------\

      end if !skip the monthly calculations/writing unless iyear>=jmosty
!       CTEM--------------/
C
C     ACCUMULATE OUTPUT DATA FOR YEARLY AVERAGED FIELDS FOR CLASS GRID-MEAN.
C     FOR BOTH PARALLEL MODE AND STAND ALONE MODE
C
      FSSTAR_YR   =0.0
      FLSTAR_YR   =0.0
      QH_YR       =0.0
      QE_YR       =0.0
      ALTOT_YR    =0.0
C
      DO 827 I=1,NLTEST
       DO 828 M=1,NMTEST
          ALVSACC_YR(I)=ALVSACC_YR(I)+ALVSROT(I,M)*FAREROT(I,M)
     1                  *FSVHROW(I)
          ALIRACC_YR(I)=ALIRACC_YR(I)+ALIRROT(I,M)*FAREROT(I,M)
     1                  *FSIHROW(I)
          FLUTACC_YR(I)=FLUTACC_YR(I)+SBC*GTROT(I,M)**4*FAREROT(I,M)
          FSINACC_YR(I)=FSINACC_YR(I)+FSSROW(I)*FAREROT(I,M)
          FLINACC_YR(I)=FLINACC_YR(I)+FDLROW(I)*FAREROT(I,M)
          HFSACC_YR(I) =HFSACC_YR(I)+HFSROT(I,M)*FAREROT(I,M)
          QEVPACC_YR(I)=QEVPACC_YR(I)+QEVPROT(I,M)*FAREROT(I,M)
          TAACC_YR(I)=TAACC_YR(I)+TAROW(I)*FAREROT(I,M)
          ROFACC_YR(I) =ROFACC_YR(I)+ROFROT(I,M)*FAREROT(I,M)*DELT
          PREACC_YR(I) =PREACC_YR(I)+PREROW(I)*FAREROT(I,M)*DELT
          EVAPACC_YR(I)=EVAPACC_YR(I)+QFSROT(I,M)*FAREROT(I,M)*DELT
828    CONTINUE
827   CONTINUE
C
      IF (IDAY.EQ.365.AND.NCOUNT.EQ.NDAY) THEN
C
       DO 829 I=1,NLTEST
         IF(FSINACC_YR(I).GT.0.0) THEN
          ALVSACC_YR(I)=ALVSACC_YR(I)/(FSINACC_YR(I)*0.5)
          ALIRACC_YR(I)=ALIRACC_YR(I)/(FSINACC_YR(I)*0.5)
         ELSE
          ALVSACC_YR(I)=0.0
          ALIRACC_YR(I)=0.0
         ENDIF
         FLUTACC_YR(I)=FLUTACC_YR(I)/(REAL(NDAY)*365.)
         FSINACC_YR(I)=FSINACC_YR(I)/(REAL(NDAY)*365.)
         FLINACC_YR(I)=FLINACC_YR(I)/(REAL(NDAY)*365.)
         HFSACC_YR(I) =HFSACC_YR(I)/(REAL(NDAY)*365.)
         QEVPACC_YR(I)=QEVPACC_YR(I)/(REAL(NDAY)*365.)
         ROFACC_YR(I) =ROFACC_YR(I)
         PREACC_YR(I) =PREACC_YR(I)
         EVAPACC_YR(I)=EVAPACC_YR(I)
         TAACC_YR(I)=TAACC_YR(I)/(REAL(NDAY)*365.)
C
         ALTOT_YR=(ALVSACC_YR(I)+ALIRACC_YR(I))/2.0
         FSSTAR_YR=FSINACC_YR(I)*(1.-ALTOT_YR)
         FLSTAR_YR=FLINACC_YR(I)-FLUTACC_YR(I)
         QH_YR=HFSACC_YR(I)
         QE_YR=QEVPACC_YR(I)
C
         WRITE(*,*) 'IYEAR=',IYEAR,' CLIMATE YEAR=',CLIMIYEAR

         WRITE(83,8103)IYEAR,FSSTAR_YR,FLSTAR_YR,QH_YR,
     1                  QE_YR,ROFACC_YR(I),PREACC_YR(I),
     2                  EVAPACC_YR(I)
C
C ADD INITIALIZTION FOR YEARLY ACCUMULATED ARRAYS
C
      call resetclassyr(nltest)
C
829    CONTINUE ! I
C
      ENDIF ! IDAY.EQ.365 .AND. NDAY
C
8103  FORMAT(1X,I5,4(F8.2,1X),F12.4,1X,2(F12.3,1X),2(A5,I1))
C
c     CTEM output and write out
c
      if(.not.parallelrun) then ! stand alone mode, includes daily and yearly mosaic-mean output for ctem
c
c     calculate daily outputs from ctem
c
       if (ctem_on) then
         if(ncount.eq.nday) then
          call ctem_daily_aw(nltest,nmtest,iday,FAREROT,
     1                      iyear,jdstd,jdsty,jdendd,jdendy,grclarea)
         endif ! if(ncount.eq.nday)
       endif ! if(ctem_on)
! c
       endif ! if(not.parallelrun)
c
c=======================================================================
c     Calculate monthly & yearly output for ctem
c

c     First initialize some output variables
c     initialization is done just before use.

      if (ctem_on) then
       if(ncount.eq.nday) then
c
        !do 861 i=1,nltest

c
          do nt=1,nmon
           if (iday.eq.mmday(nt)) then

            call resetmidmonth(nltest)

           endif
          enddo
c

          if(iday.eq.monthend(imonth+1))then

            call resetmonthend_g(nltest)

          endif
c
          if (iday .eq. 365) then

            call resetyearend_g(nltest)

          endif

861     continue
c

!       CTEM--------------\
!     Only bother with monthly calculations if we desire those outputs to be written out.
      if (iyear .ge. jmosty) then
!       CTEM--------------/

        call ctem_monthly_aw(nltest,nmtest,iday,FAREROT,iyear,nday)

        end if !to write out the monthly outputs or not
c
c       accumulate yearly outputs
            call ctem_annual_aw(nltest,nmtest,iday,FAREROT,iyear)
c
      endif ! if(ncount.eq.nday)
      endif ! if(ctem_on)
C
C     OPEN AND WRITE TO THE RESTART FILES
C
      IF (RSFILE) THEN
       IF (IDAY.EQ.365.AND.NCOUNT.EQ.NDAY) THEN
C
C       WRITE .INI_RS FOR CLASS RESTART DATA
C
        OPEN(UNIT=100,FILE=ARGBUFF(1:STRLEN(ARGBUFF))//'.INI_RS')
C
        WRITE(100,5010) TITLE1,TITLE2,TITLE3,TITLE4,TITLE5,TITLE6
        WRITE(100,5010) NAME1,NAME2,NAME3,NAME4,NAME5,NAME6
        WRITE(100,5010) PLACE1,PLACE2,PLACE3,PLACE4,PLACE5,PLACE6
C
        WRITE(100,5020)DLATROW(1),DEGLON,ZRFMROW(1),ZRFHROW(1),
     1                 ZBLDROW(1),GCROW(1),NLTEST,NMTEST
        DO I=1,NLTEST
          DO M=1,NMTEST

C         IF START_BARE (SO EITHER COMPETE OR LNDUSEON), THEN WE NEED TO CREATE
C         THE FCANROT FOR THE RS FILE.
          IF (START_BARE .AND. MOSAIC) THEN
           IF (M .LE. 2) THEN                     !NDL
            FCANROT(I,M,1)=1.0
           ELSEIF (M .GE. 3 .AND. M .LE. 5) THEN  !BDL
            FCANROT(I,M,2)=1.0
           ELSEIF (M .EQ. 6 .OR. M .EQ. 7) THEN  !CROP
            FCANROT(I,M,3)=1.0
           ELSEIF (M .EQ. 8 .OR. M .EQ. 9) THEN  !GRASSES
            FCANROT(I,M,4)=1.0
           ELSE                                  !BARE
            FCANROT(I,M,5)=1.0
           ENDIF
          ENDIF !START_BARE/MOSAIC

            WRITE(100,5040) (FCANROT(I,M,J),J=1,ICAN+1),(PAMXROT(I,M,J),
     1                      J=1,ICAN)
            WRITE(100,5040) (LNZ0ROT(I,M,J),J=1,ICAN+1),(PAMNROT(I,M,J),
     1                      J=1,ICAN)
            WRITE(100,5040) (ALVCROT(I,M,J),J=1,ICAN+1),(CMASROT(I,M,J),
     1                      J=1,ICAN)
            WRITE(100,5040) (ALICROT(I,M,J),J=1,ICAN+1),(ROOTROT(I,M,J),
     1                      J=1,ICAN)
            WRITE(100,5030) (RSMNROT(I,M,J),J=1,ICAN),
     1                      (QA50ROT(I,M,J),J=1,ICAN)
            WRITE(100,5030) (VPDAROT(I,M,J),J=1,ICAN),
     1                      (VPDBROT(I,M,J),J=1,ICAN)
            WRITE(100,5030) (PSGAROT(I,M,J),J=1,ICAN),
     1                      (PSGBROT(I,M,J),J=1,ICAN)
            WRITE(100,5040) DRNROT(I,M),SDEPROT(I,M),FAREROT(I,M)
            WRITE(100,5090) XSLPROT(I,M),GRKFROT(I,M),WFSFROT(I,M),
     1                      WFCIROT(I,M),MIDROT(I,M)
            WRITE(100,5080) (SANDROT(I,M,J),J=1,3)
            WRITE(100,5080) (CLAYROT(I,M,J),J=1,3)
            WRITE(100,5080) (ORGMROT(I,M,J),J=1,3)
C           Temperatures are in degree C
            IF (TCANROT(I,M).NE.0.0) TCANRS(I,M)=TCANROT(I,M)-273.16
            IF (TSNOROT(I,M).NE.0.0) TSNORS(I,M)=TSNOROT(I,M)-273.16
            IF (TPNDROT(I,M).NE.0.0) TPNDRS(I,M)=TPNDROT(I,M)-273.16
            WRITE(100,5050) (TBARROT(I,M,J)-273.16,J=1,3),TCANRS(I,M),
     2                      TSNORS(I,M),TPNDRS(I,M)
            WRITE(100,5060) (THLQROT(I,M,J),J=1,3),(THICROT(I,M,J),
     1                      J=1,3),ZPNDROT(I,M)
C
            WRITE(100,5070) RCANROT(I,M),SCANROT(I,M),SNOROT(I,M),
     1                      ALBSROT(I,M),RHOSROT(I,M),GROROT(I,M)
C           WRITE(100,5070) 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
          ENDDO
        ENDDO
C
        DO J=1,IGND
          WRITE(100,5002) DELZ(J),ZBOT(J)
        ENDDO
5002  FORMAT(2F8.3)
C
        WRITE(100,5200) JHHSTD,JHHENDD,JDSTD,JDENDD
        WRITE(100,5200) JHHSTY,JHHENDY,JDSTY,JDENDY
        CLOSE(100)
C
c       write .CTM_RS for ctem restart data
c
         if (ctem_on) then
            call write_ctm_rs(nltest,nmtest,FCANROT,argbuff)
        endif ! ctem_on
c
       endif ! if iday=365
      endif ! if generate restart files
c
c      check if the model is done running.
       if (iday.eq.365.and.ncount.eq.nday) then

          if (cyclemet .and. climiyear .ge. metcycendyr) then

            lopcount = lopcount+1

             if(lopcount.le.ctemloop .and. .not. transient_run)then

              rewind(12)   ! rewind met file

               if(obswetf) then
                rewind(16) !rewind obswetf file
                read(16,*) ! read in the header
               endif
              if (obslght) then ! FLAG
                 obslghtyr=-9999
                 rewind(17)
              endif

              met_rewound = .true.
              iyear=-9999
              obswetyr=-9999     !Rudra

               if(popdon) then
                 rewind(13) !rewind popd file
                 read(13,*) ! skip header (3 lines)
                 read(13,*) ! skip header (3 lines)
                 read(13,*) ! skip header (3 lines)
               endif
               if(co2on) then
                 rewind(14) !rewind co2 file
               endif

             else if (lopcount.le.ctemloop .and. transient_run)then
             ! rewind only the MET file (since we are looping over the MET  while
             ! the other inputs continue on.
               rewind(12)   ! rewind met file

               if (obslght) then ! FLAG
                 obslghtyr=-999
                 rewind(17)
                 do while (obslghtyr .lt. metcylyrst)
                   do i=1,nltest
                    read(17,*) obslghtyr,(mlightnggrd(i,j),j=1,12)
                   end do
                 end do
                 backspace(17)
               endif
             else

              if (transient_run .and. cyclemet) then
              ! Now switch from cycling over the MET to running through the file
              rewind(12)   ! rewind met file
              if (obslght) then !FLAG
                 obslghtyr=-999
                 rewind(17)
                 do while (obslghtyr .lt. metcylyrst)
                   do i=1,nltest
                    read(17,*) obslghtyr,(mlightnggrd(i,j),j=1,12)
                   end do
                 end do
               backspace(17)
              endif
              cyclemet = .false.
              lopcount = 1
              endyr = iyear + ncyear  !set the new end year

              else
               run_model = .false.
              endif

             endif

          else if (iyear .eq. endyr .and. .not. cyclemet) then

             run_model = .false.

          endif !if cyclemet and iyear > metcycendyr
       endif !last day of year check

C===================== CTEM =====================================/
C
        NCOUNT=NCOUNT+1
        IF(NCOUNT.GT.NDAY) THEN
            NCOUNT=1
        ENDIF

      ENDDO !MAIN MODEL LOOP

C     MODEL RUN HAS COMPLETED SO NOW CLOSE OUTPUT FILES AND EXIT
C==================================================================
c
c      checking the time spent for running model
c
c      call idate(today)
c      call itime(now)
c      write(*,1001) today(2), today(1), 2000+today(3), now
c 1001 format( 'end date: ', i2.2, '/', i2.2, '/', i4.4,
c     &      '; end time: ', i2.2, ':', i2.2, ':', i2.2 )
c
      IF (.NOT. PARALLELRUN) THEN
C       FIRST ANY CLASS OUTPUT FILES
        CLOSE(61)
        CLOSE(62)
        CLOSE(63)
        CLOSE(64)
        CLOSE(65)
        CLOSE(66)
        CLOSE(67)
        CLOSE(68)
        CLOSE(69)
        CLOSE(611)
        CLOSE(621)
        CLOSE(631)
        CLOSE(641)
        CLOSE(651)
        CLOSE(661)
        CLOSE(671)
        CLOSE(681)
        CLOSE(691)
        end if ! moved this up from below so it calls the close subroutine. JRM.

c       then ctem ones

        call close_outfiles()
c
c     close the input files too
      close(12)
      close(13)
      close(14)
      if (obswetf) then
        close(16)  !*.WET
      end if
      if (obslght) then
         close(17)
      end if
      call exit
C
c         the 999 label below is hit when an input file reaches its end.
999       continue

            lopcount = lopcount+1

             if(lopcount.le.ctemloop)then

              rewind(12)   ! rewind met file
c /-----------Rudra-----------------/
                if(obswetf) then
                  rewind(16) !rewind obswetf file
                  read(16,*) ! read in the header
                endif
c \------------Rudra---------------\
              met_rewound = .true.
              iyear=-9999
              obswetyr=-9999   !Rudra

               if(popdon) then
                 rewind(13) !rewind popd file
                 read(13,*) ! skip header (3 lines)
                 read(13,*) ! skip header
                 read(13,*) ! skip header
               endif
               if(co2on) then
                 rewind(14) !rewind co2 file
               endif
              if (obslght) then
                 rewind(17)
              endif

             else

              run_model = .false.

             endif

c     return to the time stepping loop
      if (run_model) then
         goto 200
      else

c     close the output files
C
      IF (.NOT. PARALLELRUN) THEN
C       FIRST ANY CLASS OUTPUT FILES
        CLOSE(61)
        CLOSE(62)
        CLOSE(63)
        CLOSE(64)
        CLOSE(65)
        CLOSE(66)
        CLOSE(67)
        CLOSE(68)
        CLOSE(69)
        CLOSE(611)
        CLOSE(621)
        CLOSE(631)
        CLOSE(641)
        CLOSE(651)
        CLOSE(661)
        CLOSE(671)
        CLOSE(681)
        CLOSE(691)
        end if ! moved this up from below so it calls the close subroutine. JRM.
c       then ctem ones
        call close_outfiles()
C
C     CLOSE THE INPUT FILES TOO
      CLOSE(12)
      CLOSE(13)
      CLOSE(14)

         CALL EXIT
      END IF

1001  continue

       write(*,*)'Error while reading WETF file'
       run_model=.false.



C ============================= CTEM =========================/

      END

      INTEGER FUNCTION STRLEN(ST)
      INTEGER       I
      CHARACTER     ST*(*)
      I = LEN(ST)
      DO WHILE (ST(I:I) .EQ. ' ')
        I = I - 1
      ENDDO
      STRLEN = I
      RETURN
      END
