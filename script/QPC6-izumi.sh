#!/bin/tcsh
setenv PBS_ACCOUNT P93300642
#
# source code (assumed to be in /glade/u/home/$USER/src)
#
set src="cam_enth_simple"
set NTHRDS="1"
#
# run with CSLAM or without
#
set res="ne16pg3_ne16pg3_mg17" #
set stopoption="nsteps"
set steps="10"
set cset="QPC6"
set homedir="/home"
set scratch="/scratch/cluster"
set queue="monster"
set pecount="192"
set compiler="nag"

set caze=new_${src}_${cset}_${res}

/project/amp/pel/src/$src/cime/scripts/create_newcase --case $scratch/$USER/$caze --compset $cset --res $res  --q $queue --walltime 00:45:00 --pecount $pecount   --compiler $compiler --run-unsupported

cd $scratch/$USER/$caze
./xmlchange CAM_CONFIG_OPTS='-phys cam7 -aquaplanet'
./xmlchange STOP_OPTION=$stopoption,STOP_N=$steps
./xmlchange DOUT_S=FALSE
./xmlchange NTHRDS=$NTHRDS

#./xmlchange --append CAM_CONFIG_OPTS="-nlev 30" #-nadv_tt=5 -cppdefs -Dwdc_debug" #there are already 6 tracers in FKESSLER


./case.setup

#./case.build
#./case.submit

