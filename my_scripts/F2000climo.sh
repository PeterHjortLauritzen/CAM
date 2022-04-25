#!/bin/tcsh
module load  python/3.7.12
if ( "$#argv" != 3) then
  echo "Wrong number of arguments specified:"
  echo "  -arg 1 res"
  echo "supported resolutions/dycores are"
  echo "se-cslam: ne30pg3_ne30pg3_mg17"
  echo "se      : ne30_ne30_mg17"
  echo "fv3     : C96_C96_mg17"
  echo "fv      : f09_f09_mg17"
  echo "mpas    : mpasa120_mpasa120"  
  echo " "
  echo "  -arg 2 duration"
  echo "debug   : 3 time-step run"
  echo "year    : 1 year run"
  echo " "
  echo "  -arg 3 src"
  echo "full path of source code"
endif
set n = 1
unset res
set res = "$argv[$n]" 

set n = 2
unset duration
set duration = "$argv[$n]" 

set n = 3
unset src
set src = "$argv[$n]" 

set cset="F2000climo"
set NTHRDS="1"
#
# set pecount (low number of procs for short run; larger for long run)
#
if ($duration == "debug") then
  echo "setting up for short debug run"
  set pecount = 225
  if ($res == "mpasa120_mpasa120") then
    set pecount = "192x1" 
  endif
  set stopoption="nsteps"
  set steps="3"
  set walltime="00:15:00"
else
  echo "setting up for 1 year run"
  if ($res == "C96_C96_mg17") then
    set pecount="384"
  else
    set pecount="900"
  endif
  if ($res == "mpasa120_mpasa120") then
    set pecount = "192x1" 
  endif
  set stopoption="nmonths"
  set steps="12"
  set walltime="03:30:00"
endif

source machine_settings.sh startup
set PBS_ACCOUNT="P93300642"

set pw=`pwd`

set caze=${cset}_${res}_ebudgets

$src/cime/scripts/create_newcase --case $scratch/$USER/$caze --compset $cset --res $res  --q $queue --walltime $walltime --pecount $pecount  --project $PBS_ACCOUNT --compiler $compiler --machine $machine --run#-unsupported

cd $scratch/$USER/$caze
#
# xmlchange commands
#
./xmlchange STOP_OPTION=$stopoption,STOP_N=$steps
./xmlchange DOUT_S=FALSE
./xmlchange NTHRDS=$NTHRDS
./xmlquery CASEROOT
./xmlchange --append CAM_CONFIG_OPTS="-cppdefs -DTRACER_CHECK -nadv_tt=5"
./xmlchange ROF_NCPL=48
#./xmlchange REST_N=12
./xmlchange REST_OPTION=nmonths
./case.setup
#
# user_nl_cam settings
#
if ($res == "ne30_ne30_mg17" || $res == "ne30pg3_ne30pg3_mg17") then
  if ($duration == "debug") then
    echo "se_statefreq       = 1"                                   >> user_nl_cam
  else
    echo "se_statefreq       = 256"                                 >> user_nl_cam
  endif
  echo "interpolate_output = .true.,.false.,.false.,.false.,.true." >> user_nl_cam      
  echo "interpolate_nlat = 192,192,192"                             >> user_nl_cam
  echo "interpolate_nlon = 288,288,288"                             >> user_nl_cam

# Energy consistent configuration
#
  echo " se_lcp_moist           = .false." >> user_nl_cam
#  echo " water_species_in_air   = 'Q'"     >> user_nl_cam
  echo " se_ftype =1"                      >> user_nl_cam
endif


endif

echo "fincl1 = 'PMID','PS','T','U','V','OMEGA500','OMEGA850'  " >> user_nl_cam
if ($duration == "debug") then
  echo "avgflag_pertape(1) = 'I'"                           >> user_nl_cam
  echo "avgflag_pertape(2) = 'A'"                           >> user_nl_cam
  echo "avgflag_pertape(3) = 'I'"                           >> user_nl_cam
  echo "nhtfrq             = 1,   1, 1, -720,-720  	" >> user_nl_cam
else
  echo "avgflag_pertape(1) = 'A'"                           >> user_nl_cam
  echo "avgflag_pertape(2) = 'A'"                           >> user_nl_cam
  echo "avgflag_pertape(3) = 'A'"                           >> user_nl_cam
  echo "nhtfrq             = 0,   0, 0, -720,-720  	" >> user_nl_cam
  echo "mfilt		 = 1,   1, 1,  30  		" >> user_nl_cam
endif
echo "inithist           =  'YEARLY'"                     >> user_nl_cam
echo "ndens              = 2,1,2,2                              ">> user_nl_cam

echo "fincl2 =   'WV_phBF','WL_phBF','WI_phBF','SE_phBF','KE_phBF','TT_phBF',  ">> user_nl_cam
echo "           'WV_phBP','WL_phBP','WI_phBP','SE_phBP','KE_phBP','TT_phBP',  ">> user_nl_cam
echo "           'WV_phAP','WL_phAP','WI_phAP','SE_phAP','KE_phAP','TT_phAP',  ">> user_nl_cam
echo "           'WV_phAM','WL_phAM','WI_phAM','SE_phAM','KE_phAM','TT_phAM',  ">> user_nl_cam
if ($res == "mpasa120_mpasa120") then
  echo "fsurdat = '/glade/p/cesmdata/inputdata/lnd/clm2/surfdata_map/release-clm5.0.34/surfdata_mpasa120_hist_78pfts_CMIP6_simyr2000_c201215.nc'" >> user_nl_clm

  echo "         'WV_dyBF','WL_dyBF','WI_dyBF','SE_dyBF','KE_dyBF','TT_dyBF',  ">> user_nl_cam
  echo "         'WV_dyBP','WL_dyBP','WI_dyBP','SE_dyBP','KE_dyBP','TT_dyBP',  ">> user_nl_cam
  echo "         'WV_dyAP','WL_dyAP','WI_dyAP','SE_dyAP','KE_dyAP','TT_dyAP',  ">> user_nl_cam
  echo "         'WV_dyAM','WL_dyAM','WI_dyAM','SE_dyAM','KE_dyAM','TT_dyAM',  ">> user_nl_cam

  echo "         'WV_dBF','WL_dBF','WI_dBF','SE_dBF','KE_dBF','TT_dBF',  ">> user_nl_cam
  echo "         'WV_dAP','WL_dAP','WI_dAP','SE_dAP','KE_dAP','TT_dAP',  ">> user_nl_cam
  echo "         'WV_dAM','WL_dAM','WI_dAM','SE_dAM','KE_dAM','TT_dAM'   ">> user_nl_cam
endif

if ($res == "ne30_ne30_mg17" || $res == "ne30pg3_ne30pg3_mg17") then
  echo "         'WV_dED','WL_dED','WI_dED','SE_dED','KE_dED','TT_dED',  ">> user_nl_cam
  echo "         'WV_dAF','WL_dAF','WI_dAF','SE_dAF','KE_dAF','TT_dAF',  ">> user_nl_cam
  echo "         'WV_dBD','WL_dBD','WI_dBD','SE_dBD','KE_dBD','TT_dBD',  ">> user_nl_cam
  echo "         'WV_dAD','WL_dAD','WI_dAD','SE_dAD','KE_dAD','TT_dAD',  ">> user_nl_cam
  echo "         'WV_dAR','WL_dAR','WI_dAR','SE_dAR','KE_dAR','TT_dAR',  ">> user_nl_cam
  echo "         'WV_dBF','WL_dBF','WI_dBF','SE_dBF','KE_dBF','TT_dBF', ">> user_nl_cam
  echo "         'WV_dBH','WL_dBH','WI_dBH','SE_dBH','KE_dBH','TT_dBH',  ">> user_nl_cam
  echo "         'WV_dCH','WL_dCH','WI_dCH','SE_dCH','KE_dCH','TT_dCH',  ">> user_nl_cam
  echo "         'WV_dAH','WL_dAH','WI_dAH','SE_dAH','KE_dAH','TT_dAH',  ">> user_nl_cam
  echo "         'WV_dBS','WL_dBS','WI_dBS','SE_dBS','KE_dBS','TT_dBS',  ">> user_nl_cam
  echo "         'WV_dAS','WL_dAS','WI_dAS','SE_dAS','KE_dAS','TT_dAS',  ">> user_nl_cam
  echo "         'WV_p2d','WL_p2d','WI_p2d','SE_p2d','KE_p2d','TT_p2d'   ">> user_nl_cam
endif


cd SourceMods/src.cam/
ln -s $pw/src.cam .
cd ../../
source $pw/machine_settings.sh cleanup
if(`hostname` == 'hobart.cgd.ucar.edu') then
  ./case.build
endif
if(`hostname` == 'izumi.unified.ucar.edu') then
 ./case.build
endif
if(`hostname` == 'cheyenne1' || `hostname` == 'cheyenne2' || `hostname` == 'cheyenne3' || `hostname` == 'cheyenne4' || `hostname` == 'cheyenne5' || `hostname` == 'cheyenne6') then
  qcmd -- ./case.build
endif
./case.submit
