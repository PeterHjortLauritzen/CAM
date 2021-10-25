#!/bin/tcsh
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

#
# source code (assumed to be in /glade/u/home/$USER/src)
#
set src="cam_energy_analysis"
#
# number of test tracers
#
set NTHRDS="1"
set pw=`pwd`

set cset="F2000climo"

#set walltime = "01:00:00"
#set stopoption="nmonths"
#set steps="13"
#set pecount="2700"

#set pecount="900"

set pecount="1800"
set walltime = "00:45:00"
set stopoption="nmonths"
set steps="2"

#set walltime = "00:15:00"
#set stopoption="nsteps"
#set steps="2"
#set pecount="360"

set pw=`pwd`
source machine_settings.sh startup
set PBS_ACCOUNT="P93300642"
set queue="premium"
echo $PBS_ACCOUNT

set caze=debug_cam_energy_pe1800
#set caze=cam_energy_clubb_scale_shf
$homedir/$USER/src/$src/cime/scripts/create_newcase --case $scratch/$USER/$caze --compset $cset --res $res  --q $queue --walltime $walltime --pecount $pecount  --project $PBS_ACCOUNT --compiler $compiler --run-unsupported
#set caze=energy_${cset}_adam_${res}
#/glade/u/home/aherring/src/cam6_2_017/cime/scripts/create_newcase --case $scratch/$USER/$caze --compset $cset --res $res  --q $queue --walltime $walltime --pecount $pecount  --project $PBS_ACCOUNT --compiler $compiler --run-unsupported


cd $scratch/$USER/$caze
./xmlchange STOP_OPTION=$stopoption,STOP_N=$steps
./xmlchange DOUT_S=FALSE
#./xmlchange DEBUG=TRUE
./xmlchange NTHRDS=$NTHRDS
./xmlchange TIMER_LEVEL=10
#
# CLUBB mods
#
#./xmlchange --append CAM_CONFIG_OPTS="-cppdefs -Dscale_shf"
#./xmlchange --append CAM_CONFIG_OPTS="-cppdefs -Dfriction_heating_clubb" #NOT WORKING!

./case.setup
  echo "se_statefreq       = 144"                     >> user_nl_cam
  echo "interpolate_output = .true.,.false.,.true."    >> user_nl_cam
  echo "interpolate_nlat   = 192,192,192"             >> user_nl_cam
  echo "interpolate_nlon   = 288,288,288"             >> user_nl_cam

if ($stopoption == "nmonths") then
  echo "avgflag_pertape(1) = 'A'" >> user_nl_cam
  echo "avgflag_pertape(2) = 'A'" >> user_nl_cam
  echo "nhtfrq             = 0,0,0,0,0" >> user_nl_cam
else
  echo "avgflag_pertape(1) = 'A'" >> user_nl_cam
  echo "avgflag_pertape(2) = 'A'" >> user_nl_cam
  echo "avgflag_pertape(3) = 'A'" >> user_nl_cam
  echo "nhtfrq             = 0,1,0,0,0" >> user_nl_cam
endif
#echo "fincl4 = 'TEINP','TEOUT','TEFIX','EFIX','DTCORE'" >> user_nl_cam
#echo "fincl1 = 'PS', 'PRECL', 'Q', 'CLDLIQ', 'RAINQM', 'T', 'U', 'V', 'iCLy', 'iCL', 'iCL2', 'OMEGA'" >> user_nl_cam
echo "ndens              = 2,1,2,2                                            ">> user_nl_cam

echo "fincl2 ='WV_pBF','WL_pBF','WI_pBF','SE_pBF','KE_pBF',  ">> user_nl_cam    
echo "        'WV_pBP','WL_pBP','WI_pBP','SE_pBP','KE_pBP',  ">> user_nl_cam  
echo "        'WV_pAP','WL_pAP','WI_pAP','SE_pAP','KE_pAP',  ">> user_nl_cam  
echo "        'WV_pAM','WL_pAM','WI_pAM','SE_pAM','KE_pAM',  ">> user_nl_cam  
echo "        'WV_dED','WL_dED','WI_dED','SE_dED','KE_dED',  ">> user_nl_cam  
echo "        'WV_dAF','WL_dAF','WI_dAF','SE_dAF','KE_dAF',  ">> user_nl_cam
echo "        'WV_dBD','WL_dBD','WI_dBD','SE_dBD','KE_dBD',  ">> user_nl_cam
echo "        'WV_dAD','WL_dAD','WI_dAD','SE_dAD','KE_dAD',  ">> user_nl_cam  
echo "        'WV_dAR','WL_dAR','WI_dAR','SE_dAR','KE_dAR',  ">> user_nl_cam  
echo "        'WV_dBF','WL_dBF','WI_dBF','SE_dBF','KE_dBF',  ">> user_nl_cam  
echo "        'WV_dBH','WL_dBH','WI_dBH','SE_dBH','KE_dBH',  ">> user_nl_cam  
echo "        'WV_dCH','WL_dCH','WI_dCH','SE_dCH','KE_dCH',  ">> user_nl_cam  
echo "        'WV_dAH','WL_dAH','WI_dAH','SE_dAH','KE_dAH',  ">> user_nl_cam  
echo "        'WV_dBS','WL_dBS','WI_dBS','SE_dBS','KE_dBS',  ">> user_nl_cam  
echo "        'WV_dAS','WL_dAS','WI_dAS','SE_dAS','KE_dAS',  ">> user_nl_cam  
echo "        'WV_p2d','WL_p2d','WI_p2d','SE_p2d','KE_p2d',  ">> user_nl_cam  
echo "        'WV_PDC','WL_PDC','WI_PDC',                    ">> user_nl_cam  
echo "        'WV_AP1','WL_AP1','WI_AP1','SE_AP1','KE_AP1',  ">> user_nl_cam  
echo "        'WV_AP2','WL_AP2','WI_AP2','SE_AP2','KE_AP2',  ">> user_nl_cam  
echo "        'WV_AP3','WL_AP3','WI_AP3','SE_AP3','KE_AP3',  ">> user_nl_cam  
echo "        'WV_AP4','WL_AP4','WI_AP4','SE_AP4','KE_AP4',  ">> user_nl_cam  
echo "        'WV_AP5','WL_AP5','WI_AP5','SE_AP5','KE_AP5',  ">> user_nl_cam  
echo "        'WV_AP6','WL_AP6','WI_AP6','SE_AP6','KE_AP6',  ">> user_nl_cam  
echo "        'WV_AP7','WL_AP7','WI_AP7','SE_AP7','KE_AP7',  ">> user_nl_cam  
echo "        'WV_AP8','WL_AP8','WI_AP8','SE_AP8','KE_AP8',  ">> user_nl_cam  
echo "        'WV_AP9','WL_AP9','WI_AP9','SE_AP9','KE_AP9',  ">> user_nl_cam  
echo "        'WV_AP10','WL_AP10','WI_AP10','SE_AP10','KE_AP10',  ">> user_nl_cam  
echo "        'WV_AP11','WL_AP11','WI_AP11','SE_AP11','KE_AP11',  ">> user_nl_cam  
echo "        'WV_AP12','WL_AP12','WI_AP12','SE_AP12','KE_AP12',  ">> user_nl_cam  
echo "        'WV_BP1','WL_BP1','WI_BP1','SE_BP1','KE_BP1',  ">> user_nl_cam  
echo "        'WV_BP2','WL_BP2','WI_BP2','SE_BP2','KE_BP2',  ">> user_nl_cam  
echo "        'WV_BP3','WL_BP3','WI_BP3','SE_BP3','KE_BP3',  ">> user_nl_cam  
echo "        'WV_BP4','WL_BP4','WI_BP4','SE_BP4','KE_BP4',  ">> user_nl_cam  
echo "        'WV_BP5','WL_BP5','WI_BP5','SE_BP5','KE_BP5',  ">> user_nl_cam  
echo "        'WV_BP6','WL_BP6','WI_BP6','SE_BP6','KE_BP6',  ">> user_nl_cam  
echo "        'WV_BP7','WL_BP7','WI_BP7','SE_BP7','KE_BP7',  ">> user_nl_cam  
echo "        'WV_BP8','WL_BP8','WI_BP8','SE_BP8','KE_BP8',  ">> user_nl_cam  
echo "        'WV_BP9','WL_BP9','WI_BP9','SE_BP9','KE_BP9',  ">> user_nl_cam 
echo "        'WV_BP10','WL_BP10','WI_BP10','SE_BP10','KE_BP10',  ">> user_nl_cam 
echo "        'WV_BP11','WL_BP11','WI_BP11','SE_BP11','KE_BP11',  ">> user_nl_cam 
echo "        'WV_BP12','WL_BP12','WI_BP12','SE_BP12','KE_BP12',  ">> user_nl_cam 
echo "        'TS','FTURB','FLAT','FNH2O','FNET_TBOT' ,'T','TBOT','FNET_TS','SST','PRECT',">> user_nl_cam 
echo "        'FLATE','FLATP','FTAU','FKE','FTAU','FPHIS'">> user_nl_cam 
echo "        'FNWV','FNLIQ','FNICE','FLAT_T'">> user_nl_cam 
echo "        'ELEAK_CLUBB','KLEAK_CLUBB','SLEAK_CLUBB','TFIX_CLUBB'" >> user_nl_cam 
 


echo "inithist           = 'YEARLY'"   >> user_nl_cam




#
# Energy consistent configuration
#

#  echo " se_lcp_moist           = .false." >> user_nl_cam
#  echo " water_species_in_air   = 'Q'"     >> user_nl_cam
#  echo " se_ftype =1"                      >> user_nl_cam

#cd SourceMods/src.cam/
#ln -s $pw/src.cam/*.F90 .
#cd ../../

source $pw/machine_settings.sh cleanup
#if(`hostname` == 'cheyenne2.ib0.cheyenne.ucar.edu' || `hostname` == 'cheyenne1.ib0.cheyenne.ucar.edu' ||`hostname` == 'cheyenne3.ib0.cheyenne.ucar.edu' ) then
qcmd -- ./case.build
#endif
./case.submit

