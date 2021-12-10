
#!/bin/tcsh
if ( "$#argv" != 3) then
  echo "Wrong number of arguments specified:"
  echo "  -arg 1 res"
  echo "supported resolutions/dycores are"
  echo "se-cslam: ne30pg3_ne30pg3_mg17"
  echo "se      : ne30_ne30_mg17"
  echo "fv3     : C96_C96_mg17"
  echo "fv      : f09_f09_mg17"
  echo "mpas    : mpasa120_mpasa120"
  echo ""
  echo "  -arg 2 src"
  echo "full path to source code"
  echo ""
  echo "  -arg 3 src short name"
  echo " "
endif
set n = 1
unset res
set res = "$argv[$n]" 
set n = 2
unset src
set src = "$argv[$n]" 
set n = 3
unset source_short_name
set source_short_name = "$argv[$n]" 
#
# source code (assumed to be in /glade/u/home/$USER/src)
#
#set src="mpas_fixes"
#endif
#
# number of test tracers
#
set NTHRDS="1"
set pw=`pwd`
#set stopoption="nsteps"
#set steps="3"
#set walltime = "00:03:00"
set stopoption="ndays"
set steps="10"
set walltime = "00:10:00"
set cset="FKESSLER"

if ($res == "C96_C96_mg17") then
  set pecount="384"
else
  set pecount="450"
endif
if ($res == "mpasa120_mpasa120") then
  set pecount = "360x1" #"540x1" #"36x1"
endif
set pw=`pwd`
source machine_settings.sh startup
set PBS_ACCOUNT="P93300642"
echo $PBS_ACCOUNT

set caze=${cset}_${res}_${source_short_name}
$src/cime/scripts/create_newcase --case $scratch/$USER/$caze --compset $cset --res $res  --q $queue --walltime $walltime --pecount $pecount  --project $PBS_ACCOUNT --compiler $compiler --run-unsupported

cd $scratch/$USER/$caze
./xmlchange STOP_OPTION=$stopoption,STOP_N=$steps
./xmlchange DOUT_S=FALSE
./xmlchange NTHRDS=$NTHRDS
./xmlchange TIMER_LEVEL=10
./xmlchange --append CAM_CONFIG_OPTS='-cppdefs -Dphl'
./case.setup

if ( $res == "ne30pg3_ne30pg3_mg17" || $res == "ne30_ne30_mg17" ) then
  echo "se_statefreq       = 144"                     >> user_nl_cam
  echo "interpolate_output = .true.,.false.,.true."   >> user_nl_cam
  echo "interpolate_nlat   = 192, 192, 192"           >> user_nl_cam
  echo "interpolate_nlon   = 288,288,288"             >> user_nl_cam
  echo "interpolate_type = 1,1,1"                     >> user_nl_cam
endif
#if ($res == "mpasa120_mpasa120") then
#  if ($cset == "FKESSLER") then
#    echo "ncdata='/project/amp/pel/inic/x1.40962.init.umjs.dry.32levels.nc'" >> user_nl_cam
#  endif
#endif

echo "avgflag_pertape(1) = 'I'" >> user_nl_cam
echo "avgflag_pertape(2) = 'I'" >> user_nl_cam
echo "avgflag_pertape(3) = 'I'" >> user_nl_cam
if ( $stopoption == "nsteps") then
    echo "nhtfrq             = 1,1,0,0,0" >> user_nl_cam
else
    echo "nhtfrq             = -24,-24,0,0,0" >> user_nl_cam
endif
echo "use_gw_front = .true." >> user_nl_cam
echo "fincl1 = 'frontgf_phys','frontga_phys'" >> user_nl_cam

source $pw/machine_settings.sh cleanup
qcmd -- ./case.build
./case.submit

