#!/bin/tcsh
#
# cd /glade/work/hannay/cesm_tags/cam6_3_119/cime/script
#
# ./create_clone --clone /glade/p/cesmdata/cseg/runs/cesm2_0/f.cam6_3_119.FLTHIST_ne30.r328.opt.001/ --case /glade/derecho/scratch/pel/cloneA --cime-output-root /glade/derecho/scratch/pel/
#
setenv PBS_ACCOUNT P93300642
#
# source code (assumed to be in /glade/u/home/$USER/src)
#
set short = "T"
#set src="cam_development"
set src="cam_enthalpy"
#set src="cam6_3_132"
#set src="/glade/work/hannay/cesm_tags/cam6_3_119/"
set cset="FLTHIST"
#set cset="FMTHIST"
set dycore = "cslam"
#set dycore = "mpas"
set proj="P93300642"
if ($dycore == "cslam") then
      	set res="ne30pg3_ne30pg3_mg17"
	#set res="ne30_ne30_mg17"
  if ($short == "T") then
    set pecount="384" 
    set stopoption="ndays"
    set steps="1"
    set wall="00:15:00"
  else
#    set pecount="900" #takes 1h20m for 3 months approximately
    set pecount="1280" 
    set stopoption="nmonths"
    set steps="3"
    set wall="04:00:00" #takes about 24 minutes
  endif
else
  set res="mpasa480_mpasa480"
  if ($short == "T") then
    set pecount="360"
    set stopoption="nsteps"
    set steps="3"
    set wall="00:10:00"
  else
#    set pecount="900" #takes 1h20m for 3 months approximately                                                                                                                         
    set pecount="1280"
    set stopoption="nmonths"
    set steps="3"
    set wall="04:00:00" #takes about 24 minutes                                                                                                                                        
  endif

endif  

#
# DO NOT MODIFY BELOW THIS LINE
#
set homedir="/glade/u/home"
set scratch="/glade/derecho/scratch"
#set scratch="/glade/scratch"
set queue="regular" #  set queue="short

set caze=${src}_${pecount}_${cset}_${res}
#/glade/work/pel/cesm_tags/cam6_3_132/cime/scripts/create_newcase --case $scratch/$USER/$caze --compset $cset --res $res  --q $queue --walltime $wall --pecount $pecount    --project $proj --walltime $wall  --run-unsupported

$home/src/$src/cime/scripts/create_newcase --case $scratch/$USER/$caze --compset $cset --res $res  --q $queue --walltime $wall    --project $proj --walltime $wall --pecount $pecount --run-unsupported

cd $scratch/$USER/$caze
./xmlchange STOP_OPTION=$stopoption,STOP_N=$steps
./xmlchange DOUT_S=FALSE
./xmlchange TIMER_LEVEL=10
#./xmlchange DEBUG=true
if ($short == "T") then
  ./xmlchange ROF_NCPL=48 #fix nuopc bug when only running a couple of time-steps
endif

./xmlchange RUN_STARTDATE=1995-01-01
#./xmlchange RUN_TYPE=hybrid
#./xmlchange RUN_REFCASE=f.cam6_3_107.FLTHIST_v0a.ne30.clm5_1.001
#./xmlchange RUN_REFDATE=1994-01-01
#./xmlchange GET_REFCASE=TRUE
#./xmlchange RUN_REFDIR=cesm2_init
./case.setup
#
# my personal io settings
#
if ($short == "T") then
  echo "avgflag_pertape(1) = 'I'"                                            >> user_nl_cam
  echo "nhtfrq = -24,0,0"                                                      >> user_nl_cam
else
  echo "avgflag_pertape(1) = 'I'"                                            >> user_nl_cam
  echo "avgflag_pertape(2) = 'A'"                                            >> user_nl_cam
  echo "nhtfrq = -96,0,0"                                                      >> user_nl_cam
endif
echo "empty_htapes = .true." >> user_nl_cam
#echo "history_amwg = .false."                                              >> user_nl_cam
echo "fincl1 = 'PS','PSDRY','PSL','OMEGA500','OMEGA850','PRECT',"     >> user_nl_cam
echo "         'TFIX','PTTEND','DTCORE','KVH','KVH_CLUBB','ABS_dPSdt',"     >> user_nl_cam
echo "fincl2 = 'PS','PSDRY','PSL','OMEGA500','OMEGA850','PRECT',"     >> user_nl_cam
#echo "         'TFIX','PTTEND','DTCORE','KVH','KVH_CLUBB','ABS_dPSdt','dT_dadadj',"     >> user_nl_cam
#echo "         'dU_dadadj','dV_dadadj'"     >> user_nl_cam

if ($dycore == "se" || $dycore == "cslam") then
    echo "interpolate_output = .true.,.true.,.true."                       >> user_nl_cam
    echo "se_statefreq       = 48"                                        >> user_nl_cam
    echo "se_statediag_numtrac = 99"                                       >> user_nl_cam
endif
if ($cset == "FLTHIST") then
#  echo "se_nsplit = 1" >> user_nl_cam
#  echo "se_rsplit = 6" >> user_nl_cam
#  echo "se_sponge_del4_nu_div_fac  = 1" >> user_nl_cam
#  echo "se_hypervis_subcycle = 1" >> user_nl_cam
#  echo "se_nu_div = 1.0E15" >> user_nl_cam
#  echo "se_nu     = 1.0E15" >> user_nl_cam
  
#  echo "se_hyperivis_subcycle_sponge = 2" >> user_nl_cam
#  echo "se_nu_top = 5E5" >> user_nl_cam
endif
if ($cset == "FMTHIST") then

#  echo "se_nsplit = 2" >> user_nl_cam
#  echo "se_rsplit = 5" >> user_nl_cam
#  echo "se_qsplit = 1" >> user_nl_cam
#  echo "se_hypervis_subcycle = 1" >> user_nl_cam
#  echo "se_nu_div = 1E15" >> user_nl_cam
#  echo "se_nu = 1E15" >> user_nl_cam
#  echo "se_sponge_del4_nu_div_fac  = 1" >> user_nl_cam
#  echo "se_sponge_del4_nu_fac  = 1" >> user_nl_cam
#  echo "se_sponge_del4_lev = 1" >> user_nl_cam
endif


#cp ~/src/for-hb-nohack-freeatm_hb-cam6_3_110/*.F90 SourceMods/src.cam/
#cp  ~/src/for-science-opt/*.F90 SourceMods/src.cam/
#qcmd -A $proj -- ./case.build 
#./case.submit
