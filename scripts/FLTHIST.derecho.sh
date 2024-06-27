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
set src="cam_enth_simple"   #cam_enthalpy"
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

/glade/campaign/cgd/amp/pel/src/$src/cime/scripts/create_newcase --case $scratch/$USER/$caze --compset $cset --res $res   --walltime $wall    --project $proj --walltime $wall --pecount $pecount --run-unsupported

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
echo "interpolate_output  =  .true.,  .false., .true., .false., .false., .true.,  .true." >> user_nl_cam
echo "interpolate_nlat    =     192,     192,    192,     192,     192,     192,   192" >> user_nl_cam
echo "interpolate_nlon    =     288,     288,    288,     288,     288,     288,   288" >> user_nl_cam

echo "print_energy_errors = .true." >> user_nl_cam
echo "se_statefreq = 1" >> user_nl_cam
echo "avgflag_pertape(1) = 'I'" >> user_nl_cam
echo "fincl1 = 'U','V','T','Q','PRECT','PRECC','PRECL','OMEGA500'" >> user_nl_cam
echo "se_ftype=1" >> user_nl_cam
echo "thermo_budget_history = .true." >> user_nl_cam
echo "thermo_budget_histfile_num = 2" >> user_nl_cam
echo "avgflag_pertape(2) = 'N'" >> user_nl_cam
echo "ndens(2) = 1" >> user_nl_cam
echo "nhtfrq =-24,1,-24,1,-24" >> user_nl_cam

echo "!avgflag_pertape(3) = 'I'" >> user_nl_cam
echo "!fincl3 = 'FRAIN_coupler','FSNOW_coupler','FEVAP_coupler','FRAIN_AC','FSNOW_AC','FEVAP_AC','FRAIN_BC','FSNOW_BC','FEVAP_BC','FRAIN_tot','FSNOW_tot','FEVAP_tot','dtot_wv','dtot_ice','dtot_liq','dtot_wv_coupler','dtot_ice_coupler','dtot_liq_coupler','rliqbc','HRAIN_AC','HSNOW_AC','HEVAP_AC','HRAIN_BC','HSNOW_BC','HEVAP_BC','HRAIN','HSNOW','HEVAP','SHFLX','FRHS_FLX','FRHS_FLXA','FRHS_FLXB','FRHS_FLXC','te_tnd','heating','radiation','te_sen','te_lat','T','TS','imbalance','enth_flx','PRECT','PRECC','PRECL','zmQ','prect_diagnosed','evap_te_pver','evap_se_pver','evap_ke_pver','evap_po_pver','evap_te_pver_cpdry','evap_se_pver_cpdry','evap_ke_pver_cpdry','evap_po_pver_cpdry','evap_lat_pver','evap_cpwv_pver','zm_te','zm_se','zm_ke','zm_po','zm_enth','zm_lat','cp_dycore_init','cp_dycore','zm_te_dme','heating_zm','heating_zm_cpice','zm_te_constantP','te_tnd_cnst_lat','cpice_dme_surf','zm_se_dme','mwv_lsT_minus_ls00_T','mliq_lfT_minus_lf00_T','Fliq_lfT_minus_lf00_TS','dmwv_dt','dmliq_dt','zm_se_constantP'" >> user_nl_cam

echo "!avgflag_pertape(4) = 'I'" >> user_nl_cam
echo "!fincl4 = 'FRAIN_coupler','FSNOW_coupler','FEVAP_coupler','FRAIN_AC','FSNOW_AC','FEVAP_AC','FRAIN_BC','FSNOW_BC','FEVAP_BC','FRAIN_tot','FSNOW_tot','FEVAP_tot','dtot_wv','dtot_ice','dtot_liq','dtot_wv_coupler','dtot_ice_coupler','dtot_liq_coupler','rliqbc','HRAIN_AC','HSNOW_AC','HEVAP_AC','HRAIN_BC','HSNOW_BC','HEVAP_BC','HRAIN','HSNOW','HEVAP','SHFLX','FRHS_FLX','FRHS_FLXA','FRHS_FLXB','FRHS_FLXC','te_tnd','heating','radiation','te_sen','te_lat','T','TS','imbalance','enth_flx','PRECT','PRECC','PRECL','zmQ','prect_diagnosed','evap_te_pver','evap_se_pver','evap_ke_pver','evap_po_pver','evap_te_pver_cpdry','evap_se_pver_cpdry','evap_ke_pver_cpdry','evap_po_pver_cpdry','evap_lat_pver','evap_cpwv_pver','zm_te','zm_se','zm_ke','zm_po','zm_enth','zm_lat','cp_dycore_init','cp_dycore','zm_te_dme','heating_zm','heating_zm_cpice','zm_te_constantP','te_tnd_cnst_lat','cpice_dme_surf','zm_se_dme','mwv_lsT_minus_ls00_T','mliq_lfT_minus_lf00_T','Fliq_lfT_minus_lf00_TS','dmwv_dt','dmliq_dt','zm_se_constantP'" >> user_nl_cam


qcmd -A $proj -- ./case.build 
./case.submit
