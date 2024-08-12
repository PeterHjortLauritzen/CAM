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
#set src="cam_development"

set src="cam_enth_simple"   #cam_enthalpy"
#set caze=cam_enth_simple
#set caze=cam_enth_simple_ocnfrc_zmconv_ke_1E-5
#set caze=cam_enth_simple_ocnfrc_efix
set caze=cam_enth_simple_ocnfrc_temp_ave
#set caze=cam_enth_simple_ocnfrc_fix2
#set src="cam6_4_015"
#set caze=cam6_4_015

set homedir="/glade/u/home"
set scratch="/glade/derecho/scratch"
#set scratch="/glade/scratch"
set queue="regular" #  set queue="short


/glade/u/home/pel/src/$src/cime/scripts/create_newcase --compset FLTHIST --res ne30pg3_ne30pg3_mg17 --case $scratch/$USER/$caze --run-unsupported --pecount 2160 --project  P93300042
cd $scratch/$USER/$caze
./xmlchange --append CAM_CONFIG_OPTS="-rad rrtmgp"
./xmlchange RUN_STARTDATE=1995-01-01
./xmlchange STOP_N=2
./xmlchange STOP_OPTION=nyears
./xmlchange RESUBMIT=5

./xmlchange TIMER_LEVEL=10
./xmlchange RUN_STARTDATE=1995-01-01
./case.setup

#
# specialized experiment
#
#echo "clubb_l_do_expldiff_rtm_thlm=.false." >> user_nl_cam

#echo "microp_aero_wsubi_scale = 0.5" >> user_nl_cam #set caze=cam_enth_simple_ocnfrc_microp_aero_wsubi_scale_0.5

#set caze=cam_enth_simple_ocnfrc_zmconv_num_cin_5
#echo "zmconv_num_cin = 5" >> user_nl_cam

#set caze=cam_enth_simple_ocnfrc_zmconv_tau7200
#echo "zmconv_tau = 7200" >> user_nl_cam

#set caze=cam_enth_simple_ocnfrc_micro_mg_dcs250
#echo "micro_mg_dcs = 250e-6" >> user_nl_cam

#set caze=cam_enth_simple_ocnfrc_micro_mg_vtrmi_factor_1.5
#echo "micro_mg_vtrmi_factor = 1.5" >> user_nl_cam

#set caze=cam_enth_simple_ocnfrc_cldfrc_dp1_0.05
#echo "cldfrc_dp1 =0.05" >> user_nl_cam

#set caze=cam_enth_simple_ocnfrc_zmconv_ke_1E-5
#echo "zmconv_ke = 1.E-5" >> user_nl_cam

#
#
#
echo "compute_enthalpy_flux = .true."                >> user_nl_cam
echo "mfilt    =       0,       5,     20,      40,      12,       120,      1,   1,12"                >> user_nl_cam
echo "nhtfrq              =       0,     -24,    -24,      -3,       0,       -2,      0,  -8760,0,0"   >> user_nl_cam
echo "ndens               =       2,       2,      2,       2,       2,       1,      2,   1,1"       >> user_nl_cam
echo "interpolate_output  =  .true.,  .true., .true., .false., .false., .true.,  .true., .false.,.false."   >> user_nl_cam
echo "interpolate_nlat    =     192,     192,    192,     192,     192,     192,   192,   192"      >> user_nl_cam
echo "interpolate_nlon    =     288,     288,    288,     288,     288,     288,   288 ,   288"     >> user_nl_cam

echo "empty_htapes = .true."              >> user_nl_cam

echo "fincl1 = 'ACTNI', 'ACTNL', 'ACTREI', 'ACTREL', 'AODDUST', 'AODVIS', 'AODVISdn','BURDENBC', 'BURDENDUST', 'BURDENPOM', 'BURDENSEASALT', "                   >> user_nl_cam
echo "'BURDENSO4', 'BURDENSOA', 'CAPE', 'CCN3', 'CDNUMC', 'CH4', 'CLDHGH', 'CLDICE', 'CLDLIQ', 'CLDLOW', 'CLDMED', 'CLDTOT', 'CLOUD', 'CMFMC_DP', "              >> user_nl_cam
echo "'CT_H2O', 'DCQ', 'DQCORE', 'DTCOND', 'DTCORE', 'DTV', 'EVAPPREC', 'EVAPSNOW', 'FCTI', 'FCTL', 'FICE', 'FLDS', 'FLNS', 'FLNSC', 'FLNT', 'FLNTC', 'FLUT', "  >> user_nl_cam
echo "'FREQZM', 'FSDS', 'FSDSC', 'FSNS', 'FSNSC', 'FSNT', 'FSNTC', 'FSNTOA', 'ICEFRAC', 'LANDFRAC', 'LHFLX', 'LWCF', 'MPDICE', 'MPDLIQ', 'MPDQ', 'MPDT', "       >> user_nl_cam
echo "'OCNFRAC', 'OMEGA', 'OMEGA500', 'PBLH', 'PHIS', 'PINT', 'PMID', 'PRECC', 'PRECL', 'PRECSC', 'PRECSL', 'PRECT', 'PS', 'PSL', 'PTEQ', 'PTTEND', 'Q', "       >> user_nl_cam
echo "'QFLX', 'QRL', 'QRS', 'QTGW', 'RCMTEND_CLUBB', 'RELHUM', 'RVMTEND_CLUBB', 'SHFLX', 'SOLIN', 'SST', 'STEND_CLUBB', 'SWCF', "                                >> user_nl_cam
echo "'T', 'TAUX', 'TAUY', 'TFIX', 'TGCLDIWP', 'TGCLDLWP', 'TMQ', 'TREFHT', 'TS', 'TTGW', 'U', 'U10', 'UBOT', 'UTGWORO', 'UTGW_TOTAL', "                         >> user_nl_cam
echo "'V', 'VBOT', 'VTGWORO', 'VTGW_TOTAL', 'WPRTP_CLUBB', 'WPTHLP_CLUBB', 'Z3', 'ZMDQ', 'ZMDT', 'N2O', 'CO2','CFC11','CFC12', "                                 >> user_nl_cam
echo "'AODVISdn','CCN3', 'CDNUMC', 'H2O', 'NUMICE', 'NUMLIQ','OMEGA500'"                                                                                         >> user_nl_cam
echo " 'AQSO4_H2O2','AQSO4_O3', 'bc_a1', 'bc_a4', 'dst_a1', 'dst_a2', 'dst_a3', 'ncl_a1',"                                                                       >> user_nl_cam
echo "'ncl_a1', 'ncl_a2', 'ncl_a3', 'pom_a1', 'pom_a4', 'so4_a1', 'so4_a2', 'so4_a3', 'soa_a2' ,"                                                                >> user_nl_cam
echo "'soa_a1', 'num_a1', 'num_a2', 'num_a3', 'num_a4',"                                                                                                         >> user_nl_cam
echo "'bc_a1SFWET', 'bc_a4SFWET', 'dst_a1SFWET', 'dst_a2SFWET', 'dst_a3SFWET', 'ncl_a1SFWET',"                                                                   >> user_nl_cam
echo "'ncl_a2SFWET', 'ncl_a3SFWET', 'pom_a1SFWET', 'pom_a4SFWET', 'so4_a1SFWET', 'so4_a2SFWET', 'so4_a3SFWET', 'soa_a1SFWET',"                                   >> user_nl_cam
echo "'soa_a2SFWET', 'bc_c1SFWET', 'bc_c4SFWET', 'dst_c1SFWET', 'dst_c2SFWET', 'dst_c3SFWET', 'ncl_c1SFWET', 'ncl_c2SFWET',"                                     >> user_nl_cam
echo "'ncl_c3SFWET', 'pom_c1SFWET', 'pom_c4SFWET', 'so4_c1SFWET', 'so4_c2SFWET', 'so4_c3SFWET', 'soa_c1SFWET', 'soa_c2SFWET',"                                   >> user_nl_cam
echo "'bc_a1DDF', 'bc_a4DDF', 'dst_a1DDF', 'dst_a2DDF', 'dst_a3DDF', 'ncl_a1DDF', 'ncl_a2DDF', 'ncl_a3DDF',"                                                     >> user_nl_cam
echo "'pom_a1DDF', 'pom_a4DDF', 'so4_a1DDF', 'so4_a2DDF', 'so4_a3DDF', 'soa_a1DDF', 'soa_a2DDF',"                                                                >> user_nl_cam
echo "'so4_a1_CLXF', 'so4_a2_CLXF', 'SFbc_a4', 'SFpom_a4', 'SFso4_a1', 'SFso4_a2',"                                                                              >> user_nl_cam
echo "'so4_a1_sfgaex1', 'so4_a2_sfgaex1', 'so4_a3_sfgaex1', 'soa_a1_sfgaex1', 'soa_a2_sfgaex1',"                                                                 >> user_nl_cam
echo "'SFdst_a1','SFdst_a2', 'SFdst_a3', 'SFncl_a1', 'SFncl_a2', 'SFncl_a3',"                                                                                    >> user_nl_cam
echo "'num_a2_sfnnuc1', 'SFSO2', 'OCN_FLUX_DMS', 'SAD_SULFC', 'SAD_TROP', 'SAD_AERO',"                                                                           >> user_nl_cam
echo "'TTEND_HFIX','WP2_CLUBB','enth_fix_fct_bc','enth_fix_fct_ac'"      >> user_nl_cam
echo "fincl3 = 'PRECT', 'PRECC', 'FLUT', 'U850', 'U200', 'V850', 'V200', 'OMEGA500', 'TS', 'SST', 'PSL'"                                                         >> user_nl_cam

echo "fincl4 =  'PRECC','PRECL'"                                                                                                                                 >> user_nl_cam

echo "fincl5 = 'Uzm','Vzm','Wzm','THzm', 'VTHzm','WTHzm','UVzm','UWzm'"                                                                                          >> user_nl_cam
echo "phys_grid_ctem_nfreq=-6"                                                                                                                                   >> user_nl_cam
echo "phys_grid_ctem_zm_nbas=120"                                                                                                                                >> user_nl_cam
echo "phys_grid_ctem_za_nlat=90"                                                                                                                                 >> user_nl_cam

echo "clubb_c8 = 4.35 "                                                                                                                                          >> user_nl_cam

#echo "fincl9 = 'EFIX'"                  >> user_nl_cam
echo "fincl9 =  'enth_prec_ac_hice:A','enth_prec_ac_hliq:A','enth_prec_bc_hice:A','enth_prec_bc_hliq:A','enth_prec_ac_fice:A','enth_prec_ac_fliq:A','enth_prec_bc_fice:A',"  >> user_nl_cam
echo " 'enth_prec_bc_fliq:A','enth_evap_hevap:A','cpice_srf:A','te_tnd:A','te_lat:A','ls_srf:A','lf_srf:A','dEdt_dme:A','dEdt_physics:A',         " >> user_nl_cam
echo "'dEdt_cpdycore:A','residual:A','dEdt_enth_fix:A','enth_fix_fct_bc_tot:A','enth_fix_fct_ac_tot:A','enthalpy_heating_fix_bc:A','enthalpy_heating_fix_ac:A',   " >> user_nl_cam
echo "'dEdt_efix_physics:A','EFIX:A'"    >> user_nl_cam
echo "avgflag_pertape(9) = 'A'" >> user_nl_cam
#qcmd -A $proj -- ./case.build 
#./case.submit
