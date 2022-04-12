#!/bin/tcsh
module load  python/3.7.12
if ( "$#argv" != 1) then
  echo "Wrong number of arguments specified:"
  echo "  -arg 1 res"
  echo "supported resolutions/dycores are"
  echo "se-cslam: ne30pg3_ne30pg3_mg17"
  echo "se      : ne30_ne30_mg17"
  echo "fv3     : C96_C96_mg17"
  echo "fv      : f09_f09_mg17"
  echo "mpas    : mpasa120_mpasa120"
endif
set n = 1
unset res
set res = "$argv[$n]" 
setenv PBS_ACCOUNT "P93300642"
#
# source code (assumed to be in /glade/u/home/$USER/src)
#
set src="topo-mods"
#
# number of test tracers
#
set NTHRDS="1"
set pw=`pwd`
set stopoption="nsteps"
set steps="2"
set cset="F2000climo"
if ($res == "C96_C96_mg17") then
  set pecount="1536"
else
  set pecount="450"
endif

if ($res == "f09_f09_mg17") then
  set pecount="450"
endif
if ($res == "mpasa120_mpasa120") then
  set src="mpas" # remove when MPAS on cam_development
#  set pecount = "36x1"
  set pecount = "576x1"
endif
set pw=`pwd`
source machine_settings.sh startup
echo $PBS_ACCOUNT
unset caze
set caze=${cset}_${src}_${res}_debug
set pw=`pwd`
/glade/scratch/$USER/$src/cime/scripts/create_newcase --case $scratch/$USER/$caze --compset $cset --res $res  --q $queue --walltime 00:10:00 --pecount $pecount  --project $PBS_ACCOUNT --compiler $compiler --run-unsupported

cd $scratch/$USER/$caze
./xmlchange STOP_OPTION=$stopoption,STOP_N=$steps
./xmlchange DOUT_S=FALSE
./xmlchange DEBUG=FALSE
./xmlchange NTHRDS=$NTHRDS
./xmlchange TIMER_LEVEL=10
./xmlchange ROF_NCPL=48
./case.setup

cd SourceMods/src.cam/
ln -s $pw/src.cam/*.F90 .
cd ../../
echo "nhtfrq=1,0" >> user_nl_cam
if ( $res == "ne30pg3_ne30pg3_mg17" || $res == "ne30_ne30_mg17" ) then
  echo "se_statefreq       = 1"                     >> user_nl_cam
  echo "interpolate_output = .true.,.true.,.true."    >> user_nl_cam
  echo "interpolate_nlat   = 192,192,192"             >> user_nl_cam
  echo "interpolate_nlon   = 288,288,288"             >> user_nl_cam
  echo "fincl1='ABS_dPSdt','OMEGA500','OMEGA850','Z500','Z300','Z200','Z100','Z3' "    >> user_nl_cam
else
  echo "fincl1='PSL','OMEGA500','OMEGA850','Z500','Z300','Z200','Z100','Z3' "    >> user_nl_cam
endif

#echo "se_nu    =1E15" >> user_nl_cam
#echo "se_nu_div=1E15" >> user_nl_cam

#echo "empty_htapes       = .true." >> user_nl_cam
#echo "avgflag_pertape(1) = 'I'" >> user_nl_cam
#echo "avgflag_pertape(2) = 'I'" >> user_nl_cam
#echo "avgflag_pertape(5) = 'I'" >> user_nl_cam
#echo "nhtfrq             = -24,-24,-24,-24,-24" >> user_nl_cam

#cd SourceMods/src.cam/
#ln -s $pw/src.cam/*.F90 .
#cd ../../
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

