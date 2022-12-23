#!/bin/tcsh
set short="T"
#set res = ne30_ne30_mg17
set res=ne16_ne16_mg17
#
# source code (assumed to be in /glade/u/home/$USER/src)
#
set src="cam-mars"
if ($short == "T") then
  set stopoption="ndays"
  set steps="1"
  set wall="00:05:00"
  set caze=mars_${res}_L49_short
else
  set stopoption="nmonths"
  set steps="12"
  set wall="00:40:00"
  set caze=mars_${res}_L49
endif
set cset="FHS94"
set pecount="225"
set queue="regular"
set PBS_ACCOUNT="P93300642"
echo $PBS_ACCOUNT
set scratch="/glade/scratch"


/glade/u/home/$USER/src/$src/cime/scripts/create_newcase --case $scratch/$USER/$caze --compset $cset --res $res  --q $queue --walltime $wall --pecount $pecount  --project $PBS_ACCOUNT --run-unsupported

cd $scratch/$USER/$caze
./xmlchange STOP_OPTION=$stopoption,STOP_N=$steps
./xmlchange DOUT_S=FALSE
./xmlchange DEBUG=FALSE
./xmlchange TIMER_LEVEL=10
#./xmlchange CAM_CONFIG_OPTS="-phys held_suarez -analytic_ic -cppdefs -Dplanet_mars -nlev 49"
#./xmlchange CAM_CONFIG_OPTS="-phys held_suarez -analytic_ic -cppdefs -Dplanet_mars -nlev 32"
./xmlchange CAM_CONFIG_OPTS="-phys held_suarez -cppdefs -Dplanet_mars -nlev 49"
./xmlquery CAM_CONFIG_OPTS
./case.setup
if ($short == "T") then
echo "se_statefreq       = 1"                       >> user_nl_cam
echo "nhtfrq = -24,-24"                             >> user_nl_cam
else
echo "se_statefreq       = 144"                     >> user_nl_cam
echo "nhtfrq = 0,0"                                 >> user_nl_cam
endif
echo "interpolate_output = .true.,.true.,.true."    >> user_nl_cam
echo "interpolate_nlat   = 192,192,192"             >> user_nl_cam
echo "interpolate_nlon   = 288,288,288"             >> user_nl_cam
echo "empty_htapes       = .true." >> user_nl_cam
#held_suarez_1994,moist_baroclinic_wave_dcmip2016,dry_baroclinic_wave_dcmip2016,dry_baroclinic_wave_jw2006,us_standard_atmosphere
#echo "analytic_ic_type='dry_baroclinic_wave_dcmip2016'"  >> user_nl_cam
echo "ncdata = 'ncdata = '/glade/u/home/pel/src/cam-mars/scripts/mars_cam_vcoords_L49_c221221.cdf5.nc'" >> user_nl_cam
#
# spun-up Held-Suarez initial condition
#
echo "ncdata = 'ncdata = 'FHS94_ne16_ne16_mg17.cam.i.0001-03-01-00000.nc'" >> user_nl_cam
echo "analytic_ic_type='us_standard_atmosphere'" >> user_nl_cam
echo "mfilt = 144" >> user_nl_cam
echo "avgflag_pertape(1) = 'I'" >> user_nl_cam
echo "fincl1 = 'PS','U','V','T','OMEGA'" >> user_nl_cam
echo "fincl2 = 'PS','U','V','T','OMEGA'" >> user_nl_cam
echo "inithist          =  'MONTHLY'" >> user_nl_cam
#
# time-steps have not been optimized for Mars at this point
#
echo "se_nsplit = 12" >> user_nl_cam
echo "se_rsplit = 1"  >> user_nl_cam
echo "se_hypervis_subcycle           = 3" >> user_nl_cam
echo "se_hypervis_subcycle_q         = 1" >> user_nl_cam
#echo "se_nsplit = 6" >> user_nl_cam

qcmd -- ./case.build
./case.submit
