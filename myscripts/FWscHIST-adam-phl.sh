#!/bin/tcsh
setenv short "T"
#setenv proj "P93300642"
setenv proj "P03010039"
unset src
setenv src "bug-N2O"
unset res
#setenv res "ne30_ne30_mg17"
#setenv res "ne30pg3_ne30pg3_mg17"
setenv res "f09_f09_mg17"
unset comp
setenv comp "HIST_CAM60%WCSC_CLM50%BGC-CROP_CICE%PRES_DOCN%DOM_MOSART_SGLC_SWAV"
unset wall

unset pes
if ($short == "T") then
  setenv pes "225"
  setenv wall "00:12:00"
else
  if ($res == "ne30_ne30_mg17" || $res == "ne30pg3_ne30pg3_mg17") then
    setenv pes "1800"
  else
    setenv pes "900"
  endif
  setenv wall "02:30:00"
endif
unset drv
setenv drv "nuopc"

# for nuopc
module load python/3.7.9
ncar_pylib

if ($short == "T") then
  echo "Short run"
endif


#setenv caze ${src}_${comp}_${res}_L58dev_${drv}_${pes}pes_`date '+%y%m%d'`_test
unset caze
if ($short == "T") then
  setenv caze N20_TT_${res}_${pes} #${src}_FWscHIST_${res}_L58dev_${drv}_${pes}pes_`date '+%y%m%d'`_wexit
else
  setenv caze N20_TT_${res}_${pes}_long #${src}_FWscHIST_${res}_L58dev_${drv}_${pes}pes_`date '+%y%m%d'`_wexit
endif

/glade/scratch/pel/$src/cime/scripts/create_newcase --case /glade/scratch/$USER/$caze --compset $comp --res $res --driver $drv --walltime $wall --mach cheyenne --pecount $pes --project $proj --compiler intel --queue regular --run-unsupported
#/glade/u/home/$USER/src/$src/cime/scripts/create_newcase --case /glade/scratch/$USER/$caze --compset $comp --res $res --driver $drv --walltime $wall --mach cheyenne --pecount $pes --project $proj --compiler intel --queue regular --run-unsupported

cd /glade/scratch/$USER/$caze 
#./xmlchange REST_N=12
#./xmlchange REST_OPTION=nmonths
./xmlchange NTHRDS=1
if ($short == "T") then
  ./xmlchange STOP_OPTION=ndays
  ./xmlchange STOP_N=1
else
  ./xmlchange STOP_OPTION=ndays
  ./xmlchange STOP_N=60
endif


#xxx ./xmlchange STOP_OPTION=nmonths
#xxx ./xmlchange STOP_N=2
./xmlchange RESUBMIT=0
./xmlchange DOUT_S=FALSE
./xmlchange RUN_STARTDATE=2011-01-01

#./xmlchange CAM_CONFIG_OPTS='-phys cam_dev -nlev 58' --append
./xmlchange CAM_CONFIG_OPTS='-phys cam_dev -nlev 58 -chem waccm_sc_mam4 -microphys mg2'
./xmlchange --append CAM_CONFIG_OPTS="-cppdefs -DN2O_diag -nadv_tt=5" 

if ($res == "ne30_ne30_mg17" || $res == "ne30pg3_ne30pg3_mg17") then
  echo "ncdata = '/glade/p/cesm/amwg_dev/juliob/FWsc_ne30pg3_58L_GRID_48_taperstart10km_lowtop_BL10_v3_beta1p75_Top_43km.nc'">>user_nl_cam
##echo "ncdata = '/glade/p/acom/acom-climate/jzhan166/acclip/ncdata/FCnudged_f09.mam.mar27.2000_2021.001.cam.i.2011-01-01-00000.ncne0np4_58L_cdf5.nc'">>user_nl_cam
#echo "ncdata         = '/glade/scratch/jzhan166/archive/f.e21.FCnudged_58L.ne30_ne30_mg17.cam6_3_058.dtsens.2011_fvIC_wetexit/atm/hist/f.e21.FCnudged_58L.ne30_ne30_mg17.cam6_3_058.dtsens.2011_fvIC_wetexit.cam.i.2011-05-01-00000.nc'">>user_nl_cam
else
  echo "ncdata = '/glade/p/cesm/amwg_dev/juliob/cam_ic_files/fv_dycore/cami_test_FV_1_0.9x1.25_58L_GRID_48_taperstart10km_lowtop_BL10_v3p1_beta1p75_Top_43km.nc'">>user_nl_cam
endif

#echo "ncdata = '/glade/p/cesm/amwg_dev/juliob/FWsc_ne30pg3_58L_GRID_48_taperstart10km_lowtop_BL10_v3_beta1p75_Top_43km.nc'">>user_nl_cam
#echo "ncdata = '/glade/work/aherring/tmp/f.e22.FCnudged.ne30_ne30_mg17.release-cesm2.2.0_spinup.2010_2020.001.map_TO_L58.cam.i.2011-06-01-00000.nc'">>user_nl_cam
echo "inithist = 'MONTHLY'" >>user_nl_cam
echo "empty_htapes      = .true.                                        ">>user_nl_cam
if ($short == "T") then
  echo "mfilt  = 1,1,1,1,1,73,365                                                                    ">>user_nl_cam
else
  echo "mfilt  = 1,1,1,1,1,1                                                        ">>user_nl_cam
endif

echo "fexcl1 = ' '                                                      ">>user_nl_cam

echo "fincl1 = 'N2O',               " >>user_nl_cam                
echo "         'N2O_AC1','N2O_AC2','N2O_AC3','N2O_AC4','N2O_AC5', ">>user_nl_cam     
echo "         'N2O_AC6','N2O_AC7','N2O_AC8','N2O_AC9','N2O_AC10',">>user_nl_cam  
echo "	       'N2O_AC11b','N2O_AC11c','N2O_AC11d','N2O_AC11e',   ">>user_nl_cam
echo "	       'N2O_AC11f','N2O_AC11g','N2O_AC11h',               ">>user_nl_cam
echo "         'N2O_AC11','N2O_AC12', 'N2O_AC13',                 ">>user_nl_cam                                                              
echo "	       'N2O_BC1','N2O_BC2'                                ">>user_nl_cam



if ($res == "ne30_ne30_mg17" || $res == "ne30pg3_ne30pg3_mg17") then
  echo " ,'N2O_DP','N2O_dyn1','N2O_dyn2'                         ">>user_nl_cam
  echo " ,'N2O_dyn_remap1','N2O_dyn_remap2'                      ">>user_nl_cam

  echo "se_statediag_numtrac = 300                                         ">>user_nl_cam
  if ($short == "T") then
    echo "se_statefreq      = 1                                            ">>user_nl_cam
  else
    echo "se_statefreq      = 144                                            ">>user_nl_cam
  endif
endif


echo "fincl4 = 'TT1_BC1',       'TT2_BC1',       'TT3_BC1',       'TT4_BC1',       'TT5_BC1',       'Q_BC1',">>user_nl_cam
echo "         'TT1_AC4',       'TT2_AC4',       'TT3_AC4',       'TT4_AC4',       'TT5_AC4',       'Q_AC4',">>user_nl_cam
echo "         'TT1_AC5',       'TT2_AC5',       'TT3_AC5',       'TT4_AC5',       'TT5_AC5',       'Q_AC5',">>user_nl_cam
echo "         'TT1_AC11d',     'TT2_AC11d',     'TT3_AC11d',     'TT4_AC11d',     'TT5_AC11d',     'Q_AC11d',">>user_nl_cam
echo "         'TT1_AC11e',     'TT2_AC11e',     'TT3_AC11e',     'TT4_AC11e',     'TT5_AC11e',     'Q_AC11e',">>user_nl_cam
echo "         'TT1_AC12',      'TT2_AC12',      'TT3_AC12',      'TT4_AC12',      'TT5_AC12',      'Q_AC12' ">>user_nl_cam
if ($res == "ne30_ne30_mg17" || $res == "ne30pg3_ne30pg3_mg17") then
echo "        ,'TT1_DP',        'TT2_DP',        'TT3_DP',        'TT4_DP',        'TT5_DP',        'Q_DP',">>user_nl_cam
echo "         'TT1_dyn1',      'TT2_dyn1',      'TT3_dyn1',      'TT4_dyn1',      'TT5_dyn1',      'Q_dyn1',">>user_nl_cam
echo "         'TT1_dyn2',      'TT2_dyn2',      'TT3_dyn2',      'TT4_dyn2',      'TT5_dyn2',      'Q_dyn2',">>user_nl_cam
echo "         'TT1_dyn_remap1','TT1_dyn_remap1','TT3_dyn_remap1','TT4_dyn_remap1','TT5_dyn_remap1','Q_dyn_remap1',">>user_nl_cam
echo "         'TT1_dyn_remap2','TT1_dyn_remap2','TT3_dyn_remap2','TT4_dyn_remap2','TT5_dyn_remap2','Q_dyn_remap2'">>user_nl_cam
endif


if ($res == "ne30_ne30_mg17" || $res == "ne30pg3_ne30pg3_mg17") then
  echo "fincl3 = 'N2O',               " >>user_nl_cam                
  echo "         'N2O_AC1','N2O_AC2','N2O_AC3','N2O_AC4','N2O_AC5', ">>user_nl_cam     
  echo "         'N2O_AC6','N2O_AC7','N2O_AC8','N2O_AC9','N2O_AC10',">>user_nl_cam  
  echo "	       'N2O_AC11b','N2O_AC11c','N2O_AC11d','N2O_AC11e',   ">>user_nl_cam
  echo "	       'N2O_AC11f','N2O_AC11g','N2O_AC11h',               ">>user_nl_cam
  echo "         'N2O_AC11','N2O_AC12', 'N2O_AC13',                 ">>user_nl_cam                                                              
  echo "	       'N2O_BC1','N2O_BC2'                                ">>user_nl_cam
  echo " ,'N2O_DP','N2O_dyn1','N2O_dyn2'                         ">>user_nl_cam
endif

echo "fincl2 =   'WV_phBF','WL_phBF','WI_phBF','SE_phBF','KE_phBF','TT_phBF',  ">> user_nl_cam
echo "           'WV_phBP','WL_phBP','WI_phBP','SE_phBP','KE_phBP','TT_phBP',  ">> user_nl_cam
echo "           'WV_phAP','WL_phAP','WI_phAP','SE_phAP','KE_phAP','TT_phAP',  ">> user_nl_cam
echo "           'WV_phAM','WL_phAM','WI_phAM','SE_phAM','KE_phAM','TT_phAM',  ">> user_nl_cam
echo "           'WV_phB1','WL_phB1','WI_phB1','SE_phB1','KE_phB1','TT_phB1',  ">> user_nl_cam
echo "           'WV_phA1','WL_phA1','WI_phA1','SE_phA1','KE_phA1','TT_phA1'   ">> user_nl_cam

#  echo "         'WV_dED','WL_dED','WI_dED','SE_dED','KE_dED','TT_dED',  ">> user_nl_cam
#  echo "         'WV_dAF','WL_dAF','WI_dAF','SE_dAF','KE_dAF','TT_dAF',  ">> user_nl_cam
#  echo "         'WV_dBD','WL_dBD','WI_dBD','SE_dBD','KE_dBD','TT_dBD',  ">> user_nl_cam
#  echo "         'WV_dAD','WL_dAD','WI_dAD','SE_dAD','KE_dAD','TT_dAD',  ">> user_nl_cam
#  echo "         'WV_dAR','WL_dAR','WI_dAR','SE_dAR','KE_dAR','TT_dAR',  ">> user_nl_cam
#  echo "         'WV_dBF','WL_dBF','WI_dBF','SE_dBF','KE_dBF','TT_dBF', ">> user_nl_cam
#  echo "         'WV_dBH','WL_dBH','WI_dBH','SE_dBH','KE_dBH','TT_dBH',  ">> user_nl_cam
#  echo "         'WV_dCH','WL_dCH','WI_dCH','SE_dCH','KE_dCH','TT_dCH',  ">> user_nl_cam
#  echo "         'WV_dAH','WL_dAH','WI_dAH','SE_dAH','KE_dAH','TT_dAH',  ">> user_nl_cam
#  echo "         'WV_dBS','WL_dBS','WI_dBS','SE_dBS','KE_dBS','TT_dBS',  ">> user_nl_cam
#  echo "         'WV_dAS','WL_dAS','WI_dAS','SE_dAS','KE_dAS','TT_dAS'   ">> user_nl_cam

if ($res == "ne30_ne30_mg17" || $res == "ne30pg3_ne30pg3_mg17") then
  echo "interpolate_output = .false.,.false.,.true.,.false.                                                         ">>user_nl_cam
  echo "interpolate_nlat = 192,192,192,192                                                                       ">>user_nl_cam
  echo "interpolate_nlon = 288,288,288,288                                                                       ">>user_nl_cam
endif

if ($short == "T") then
  echo "avgflag_pertape(1) = 'A'                                   ">>user_nl_cam
  echo "avgflag_pertape(2) = 'A'                                   ">>user_nl_cam
  echo "avgflag_pertape(3) = 'A'                                   ">>user_nl_cam
  echo "avgflag_pertape(4) = 'A'                                   ">>user_nl_cam
  echo "nhtfrq             =1,1,1,1                                ">>user_nl_cam
else
  echo "avgflag_pertape(1) = 'A'                                   ">>user_nl_cam
  echo "avgflag_pertape(2) = 'A'                                   ">>user_nl_cam
  echo "avgflag_pertape(3) = 'A'                                   ">>user_nl_cam
  echo "avgflag_pertape(4) = 'A'                                   ">>user_nl_cam
  echo "nhtfrq             =-24,-24,-24,-24                        ">>user_nl_cam
endif

#phl cp /glade/u/home/aherring/src/cam6_3_058.dtsens/usr_src/n2o/wetexit/* /glade/scratch/$USER/$caze/SourceMods/src.cam/

./case.setup

#cp /glade/scratch/hannay/archive/f.e21.FWscHIST.ne30_L48_BL10_cam6_3_035.tphysac_reorder_zm2.001.hf2/rest/1990-01-01-00000/* /glade/scratch/$USER/$caze/run/

qcmd -- ./case.build
#qcmd -- ./case.build --skip-provenance-check
./case.submit
