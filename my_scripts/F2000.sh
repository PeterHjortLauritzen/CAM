#!/bin/tcsh
setenv short "T"
unset proj
setenv proj "P93300642"
#setenv proj "P03010039"
unset src
setenv src "cam-mpas-pdc-paper"
unset res
#setenv res "ne30_ne30_mg17"
#setenv res "ne30pg3_ne30pg3_mg17"
setenv res "mpasa120_mpasa120"
#setenv res "f09_f09_mg17"
unset comp
setenv comp "F2000climo"
unset wall

unset pes
if ($short == "T") then
  echo "Short run"
  if ($res == "mpasa120_mpasa120") then
    unset pes
    setenv pes "192"
  else
    setenv pes "225"
  endif
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

unset caze
if ($short == "T") then
  setenv caze ${src}_${res}_${pes}_short
else
  setenv caze ${src}_${res}_${pes}_long
endif
/glade/scratch/pel/$src/cime/scripts/create_newcase --case /glade/scratch/$USER/$caze --compset $comp --res $res --driver $drv --walltime $wall --mach cheyenne --pecount $pes --project $proj --compiler intel --queue regular --run-unsupported

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

./xmlchange DOUT_S=FALSE
echo "empty_htapes       = .true." >> user_nl_cam

echo "fincl1 = 'PS','zint_mpas','zmid_mpas','zm_phBF','zm_phAP','zm_phAM','zi_phBF','zi_phAP','zi_phAM','PHIS'," >> user_nl_cam
echo "         'zm_diff_phAP_phBF','zm_diff_phAM_phAP','dPSdt'" >> user_nl_cam
#
# Physics energy diagnostics
#
echo "fincl2 =   'WV_phBF','WL_phBF','WI_phBF','SE_phBF','KE_phBF','TT_phBF',  ">> user_nl_cam
echo "           'WV_phBP','WL_phBP','WI_phBP','SE_phBP','KE_phBP','TT_phBP',  ">> user_nl_cam
echo "           'WV_phAP','WL_phAP','WI_phAP','SE_phAP','KE_phAP','TT_phAP',  ">> user_nl_cam
echo "           'WV_phAM','WL_phAM','WI_phAM','SE_phAM','KE_phAM','TT_phAM'   ">> user_nl_cam
if ($res == "mpasa120_mpasa120") then
#
# dycore energy in physics
#
echo "          ,'WV_dyBF','WL_dyBF','WI_dyBF','SE_dyBF','KE_dyBF','TT_dyBF',  ">> user_nl_cam
echo "           'WV_dyBP','WL_dyBP','WI_dyBP','SE_dyBP','KE_dyBP','TT_dyBP',  ">> user_nl_cam
echo "           'WV_dyAP','WL_dyAP','WI_dyAP','SE_dyAP','KE_dyAP','TT_dyAP',  ">> user_nl_cam
echo "           'WV_dyAM','WL_dyAM','WI_dyAM','SE_dyAM','KE_dyAM','TT_dyAM'   ">> user_nl_cam
#
# MPAS energy diags in dp_coupling
#
echo "          ,'WV_dBF','WL_dBF','WI_dBF','SE_dBF','KE_dBF',  ">> user_nl_cam
echo "           'WV_dAP','WL_dAP','WI_dAP','SE_dAP','KE_dAP',  ">> user_nl_cam
echo "           'WV_dAM','WL_dAM','WI_dAM','SE_dAM','KE_dAM'   ">> user_nl_cam
endif

if ($short == "T") then
  echo "avgflag_pertape(1) = 'I'                                         ">>user_nl_cam
  echo "avgflag_pertape(2) = 'I'                                         ">>user_nl_cam
  echo "avgflag_pertape(3) = 'I'                                         ">>user_nl_cam
  echo "avgflag_pertape(4) = 'I'                                         ">>user_nl_cam
  echo "nhtfrq             =-12,-12,-12,-12                              ">>user_nl_cam
else
  echo "avgflag_pertape(1) = 'A'                                         ">>user_nl_cam
  echo "avgflag_pertape(2) = 'A'                                         ">>user_nl_cam
  echo "avgflag_pertape(3) = 'A'                                         ">>user_nl_cam
  echo "avgflag_pertape(4) = 'A'                                         ">>user_nl_cam
  echo "nhtfrq             = 0,0,0,0,0                                   ">>user_nl_cam
endif

if ($res == "mpasa120_mpasa120") then
#  echo "flanduse_timeseries = '/glade/p/cesmdata/inputdata/lnd/clm2/surfdata_map/landuse.timeseries_mpasa120_hist_78pfts_CMIP6_simyr1850-2015_c211108.nc'" >> user_nl_clm
  echo "fsurdat = '/glade/p/cesmdata/inputdata/lnd/clm2/surfdata_map/surfdata_mpasa120_hist_78pfts_CMIP6_simyr2000_c211108.nc'" >> user_nl_clm
endif

./case.setup


#qcmd -- ./case.build
#./case.submit
