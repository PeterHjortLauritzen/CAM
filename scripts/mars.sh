#!/bin/tcsh
setenv PBS_ACCOUNT P93300642
#
# source code (assumed to be in /glade/u/home/$USER/src)
#
set src="cam-mars"
set NTHRDS="1"
#
# run with CSLAM or without
#
#set res="ne30pg3_ne30pg3_mg17" #cslam
set res="ne30_ne30_mg17"        #no cslam

set stopoption="nsteps"
set steps="5"
#
# DO NOT MODIFY BELOW THIS LINE
#
set cset="FADIAB"

echo "setting up for Cheyenne"
set inic="/glade/p/cgd/amp/pel/inic"
set homedir="/glade/u/home"
set scratch="/glade/scratch"
set queue="regular"
#
# 900, 1800, 2700, 5400 (pecount should divide 6*30*30 evenly)
#
set pecount="450"
set walltime="00:05:00"  
set compiler="intel"

set caze=test_${src}_${cset}_CAM_${res}_${pecount}_NTHRDS${NTHRDS}_${steps}${stopoption}
$homedir/$USER/src/$src/cime/scripts/create_newcase --case $scratch/$USER/$caze --compset $cset --res $res  --q $queue --walltime $walltime --pecount $pecount  --project $PBS_ACCOUNT --compiler $compiler --run-unsupported

cd $scratch/$USER/$caze
./xmlchange STOP_OPTION=$stopoption,STOP_N=$steps
./xmlchange DOUT_S=FALSE
./xmlchange --append CAM_CONFIG_OPTS="-cppdefs -Dplanet_mars"
#./xmlchange DEBUG=TRUE
./xmlchange NTHRDS=$NTHRDS

./case.setup

#echo "ncdata = '$inic/mars.i.spunup-from-Mars-data-base-i.nc'"     >>  user_nl_cam
echo "ncdata = '$inic/mars.i.3month-spinup.nc'"     >>  user_nl_cam
echo "bnd_topo = '$inic/mars_ne30np4_nc3000_Co092_Fi001_PF_nullRR_Nsw065_20170928.nc'" >>  user_nl_cam
echo "fincl1             = 'PS:I','U:I','V:I','T:I','PHIS'" >>  user_nl_cam
echo "se_statefreq       = 1"        >> user_nl_cam
echo "avgflag_pertape(1) = 'I'" >> user_nl_cam
echo "nhtfrq             = -24,-24 " >> user_nl_cam

#echo "se_nsplit          = 10"   >> user_nl_cam
#echo "inithist           = 'DAILY'"   >> user_nl_cam
echo "se_hypervis_subcycle = 1"   >> user_nl_cam
echo "interpolate_output   = .true.,.true.,.true.,.true.,.true.,.true.,.true."   >> user_nl_cam


#
# set dynamics time-steps
#
echo "se_nsplit = 6" >> user_nl_cam
echo "se_hypervis_subcycle=8">> user_nl_cam
echo "se_nu     =  0.30E+15">> user_nl_cam
echo "se_nu_div  =  0.15E+16">> user_nl_cam
echo "se_nu_p   =  0.15E+16">> user_nl_cam
echo "se_nu_top =  0.25E+06">> user_nl_cam
echo "se_hypervis_on_plevs           = .false.">> user_nl_cam


#
# set physical constants
#
echo "gravit = 3.72"      >> user_nl_cam
echo "sday   = 88642.0"   >> user_nl_cam
echo "mwdry  = 43.34"     >> user_nl_cam #0.04334 - original value
echo "cpair  = 735.0"     >> user_nl_cam
echo "rearth = 3.38992e6" >> user_nl_cam

if(`hostname` == 'hobart.cgd.ucar.edu') then
  ./case.build
else
qcmd -- ./case.build
endif
./case.submit
