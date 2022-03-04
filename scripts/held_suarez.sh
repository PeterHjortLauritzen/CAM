#!/bin/tcsh
set res = ne30_ne30_mg17
#
# source code (assumed to be in /glade/u/home/$USER/src)
#
set src="cam-mars"
set stopoption="ndays"
set steps="1"
set cset="FHS94"
set pecount="225"
set queue="regular"
set PBS_ACCOUNT="P93300642"
echo $PBS_ACCOUNT
set scratch="/glade/scratch"

set caze=${cset}_${res}
/glade/u/home/$USER/src/$src/cime/scripts/create_newcase --case $scratch/$USER/$caze --compset $cset --res $res  --q $queue --walltime 00:10:00 --pecount $pecount  --project $PBS_ACCOUNT --run-unsupported

cd $scratch/$USER/$caze
./xmlchange STOP_OPTION=$stopoption,STOP_N=$steps
./xmlchange DOUT_S=FALSE
./xmlchange DEBUG=FALSE
./xmlchange TIMER_LEVEL=10
./xmlchange CAM_CONFIG_OPTS="-phys held_suarez -analytic_ic -cppdefs -Dplanet_mars -nlev 32"
./xmlquery CAM_CONFIG_OPTS
./case.setup
echo "se_statefreq       = 1"                     >> user_nl_cam

echo "interpolate_output = .true.,.true.,.true."    >> user_nl_cam
echo "interpolate_nlat   = 192,192,192"             >> user_nl_cam
echo "interpolate_nlon   = 288,288,288"             >> user_nl_cam
echo "empty_htapes       = .true." >> user_nl_cam
#held_suarez_1994,moist_baroclinic_wave_dcmip2016,dry_baroclinic_wave_dcmip2016,dry_baroclinic_wave_jw2006,us_standard_atmosphere
echo "analytic_ic_type='dry_baroclinic_wave_dcmip2016'"  >> user_nl_cam
echo "nhtfrq = -24" >> user_nl_cam
echo "mfilt = 144" >> user_nl_cam
echo "avgflag_pertape(1) = 'I'" >> user_nl_cam
echo "fincl1 = 'PS','U','V','T','OMEGA'" >> user_nl_cam
#
# time-steps have not been optimized for Mars at this point
#
echo "se_nsplit = 6" >> user_nl_cam

qcmd -- ./case.build
./case.submit
