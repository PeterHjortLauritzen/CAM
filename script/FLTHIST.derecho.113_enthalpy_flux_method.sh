#!/bin/tcsh
#
# cd /glade/work/hannay/cesm_tags/cam6_3_119/cime/script
#
# ./create_clone --clone /glade/p/cesmdata/cseg/runs/cesm2_0/f.cam6_3_119.FLTHIST_ne30.r328.opt.001/ --case /glade/derecho/scratch/pel/cloneA --cime-output-root /glade/derecho/scratch/pel/
#
setenv PBS_ACCOUNT P93300642
set caze="113_enthalpy_flux_method1"

set homedir="/glade/u/home"
set scratch="/glade/derecho/scratch"
#set scratch="/glade/scratch"
set queue="regular" #  set queue="short


 /glade/work/gmarques/cesm.sandboxes/cesm3_0_enthalpy_branch/cime/scripts/create_newcase --compset FLTHIST --res ne30pg3_ne30pg3_mg17 --case $scratch/$USER/$caze --run-unsupported --pecount 2160 --project  P93300042
cd $scratch/$USER/$caze
./xmlchange --append CAM_CONFIG_OPTS="-rad rrtmgp"
./xmlchange RUN_STARTDATE=1995-01-01
./xmlchange STOP_N=2
./xmlchange STOP_OPTION=nyears
./xmlchange RESUBMIT=5

./xmlchange TIMER_LEVEL=10
./xmlchange RUN_STARTDATE=1995-01-01
./case.setup

cp /glade/u/home/pel/src/cam_enth_simple/src/physics/cam7/physpkg.F90 SourceMods/src.cam/
cp /glade/u/home/pel/src/cam_enth_simple/src/physics/cam/check_energy.F90 SourceMods/src.cam/

echo "mfilt               =       0,       5,     20,      40,      12,    120,       1,      1,      1,    1" >> user_nl_cam
echo "nhtfrq              =       0,     -24,    -24,      -3,       0,     -2,       0,      0,      0,    0" >> user_nl_cam
echo "ndens               =       2,       2,      2,       2,       2,      1,       2,      1,      1,    2" >> user_nl_cam
echo "interpolate_output  =  .true.,  .true., .true., .false., .false., .true.,  .true.,.false.,.false.,.true." >> user_nl_cam
echo "interpolate_nlat    =     192,     192,    192,     192,     192,     192,   192 ,   192 ,   192 ,   192" >> user_nl_cam
echo "interpolate_nlon    =     288,     288,    288,     288,     288,     288,   288 ,   288 ,   288 ,   288" >> user_nl_cam
echo "" >> user_nl_cam
echo "" >> user_nl_cam
echo "fincl1 = 'ACTNI', 'ACTNL', 'ACTREI', 'ACTREL', 'AODDUST', 'AODVIS', 'AODVISdn','BURDENBC', 'BURDENDUST', 'BURDENPOM', 'BURDENSEASALT', " >> user_nl_cam
echo "'BURDENSO4', 'BURDENSOA', 'CAPE', 'CCN3', 'CDNUMC', 'CH4', 'CLDHGH', 'CLDICE', 'CLDLIQ', 'CLDLOW', 'CLDMED', 'CLDTOT', 'CLOUD', 'CMFMC_DP', " >> user_nl_cam
echo "'CT_H2O', 'DCQ', 'DQCORE', 'DTCOND', 'DTCORE', 'DTV', 'EVAPPREC', 'EVAPSNOW', 'FCTI', 'FCTL', 'FICE', 'FLDS', 'FLNS', 'FLNSC', 'FLNT', 'FLNTC', 'FLUT', " >> user_nl_cam
echo "'FREQZM', 'FSDS', 'FSDSC', 'FSNS', 'FSNSC', 'FSNT', 'FSNTC', 'FSNTOA', 'ICEFRAC', 'LANDFRAC', 'LHFLX', 'LWCF', 'MPDICE', 'MPDLIQ', 'MPDQ', 'MPDT', " >> user_nl_cam
echo "'OCNFRAC', 'OMEGA', 'OMEGA500', 'PBLH', 'PHIS', 'PINT', 'PMID', 'PRECC', 'PRECL', 'PRECSC', 'PRECSL', 'PRECT', 'PS', 'PSL', 'PTEQ', 'PTTEND', 'Q', " >> user_nl_cam
echo "'QFLX', 'QRL', 'QRS', 'QTGW', 'RCMTEND_CLUBB', 'RELHUM', 'RVMTEND_CLUBB', 'SHFLX', 'SOLIN', 'SST', 'STEND_CLUBB', 'SWCF', " >> user_nl_cam
echo "'T', 'TAUX', 'TAUY', 'TFIX', 'TGCLDIWP', 'TGCLDLWP', 'TMQ', 'TREFHT', 'TS', 'TTGW', 'U', 'U10', 'UBOT', 'UTGWORO', 'UTGW_TOTAL', " >> user_nl_cam
echo "'V', 'VBOT', 'VTGWORO', 'VTGW_TOTAL', 'WPRTP_CLUBB', 'WPTHLP_CLUBB', 'Z3', 'ZMDQ', 'ZMDT', 'N2O', 'CO2','CFC11','CFC12', " >> user_nl_cam
echo "'AODVISdn','CCN3', 'CDNUMC', 'H2O', 'NUMICE', 'NUMLIQ','OMEGA500'" >> user_nl_cam
echo " 'AQSO4_H2O2','AQSO4_O3', 'bc_a1', 'bc_a4', 'dst_a1', 'dst_a2', 'dst_a3', 'ncl_a1'," >> user_nl_cam
echo "'ncl_a1', 'ncl_a2', 'ncl_a3', 'pom_a1', 'pom_a4', 'so4_a1', 'so4_a2', 'so4_a3', 'soa_a2' ," >> user_nl_cam
echo "'soa_a1', 'num_a1', 'num_a2', 'num_a3', 'num_a4'," >> user_nl_cam
echo "'bc_a1SFWET', 'bc_a4SFWET', 'dst_a1SFWET', 'dst_a2SFWET', 'dst_a3SFWET', 'ncl_a1SFWET'," >> user_nl_cam
echo "'ncl_a2SFWET', 'ncl_a3SFWET', 'pom_a1SFWET', 'pom_a4SFWET', 'so4_a1SFWET', 'so4_a2SFWET', 'so4_a3SFWET', 'soa_a1SFWET'," >> user_nl_cam
echo "'soa_a2SFWET', 'bc_c1SFWET', 'bc_c4SFWET', 'dst_c1SFWET', 'dst_c2SFWET', 'dst_c3SFWET', 'ncl_c1SFWET', 'ncl_c2SFWET'," >> user_nl_cam
echo "'ncl_c3SFWET', 'pom_c1SFWET', 'pom_c4SFWET', 'so4_c1SFWET', 'so4_c2SFWET', 'so4_c3SFWET', 'soa_c1SFWET', 'soa_c2SFWET'," >> user_nl_cam
echo "'bc_a1DDF', 'bc_a4DDF', 'dst_a1DDF', 'dst_a2DDF', 'dst_a3DDF', 'ncl_a1DDF', 'ncl_a2DDF', 'ncl_a3DDF'," >> user_nl_cam
echo "'pom_a1DDF', 'pom_a4DDF', 'so4_a1DDF', 'so4_a2DDF', 'so4_a3DDF', 'soa_a1DDF', 'soa_a2DDF'," >> user_nl_cam
echo "'so4_a1_CLXF', 'so4_a2_CLXF', 'SFbc_a4', 'SFpom_a4', 'SFso4_a1', 'SFso4_a2'," >> user_nl_cam
echo "'so4_a1_sfgaex1', 'so4_a2_sfgaex1', 'so4_a3_sfgaex1', 'soa_a1_sfgaex1', 'soa_a2_sfgaex1'," >> user_nl_cam
echo "'SFdst_a1','SFdst_a2', 'SFdst_a3', 'SFncl_a1', 'SFncl_a2', 'SFncl_a3'," >> user_nl_cam
echo "'num_a2_sfnnuc1', 'SFSO2', 'OCN_FLUX_DMS', 'SAD_SULFC', 'SAD_TROP', 'SAD_AERO'" >> user_nl_cam
echo "" >> user_nl_cam
echo "fincl3 = 'PRECT', 'PRECC', 'FLUT', 'U850', 'U200', 'V850', 'V200', 'OMEGA500', 'TS', 'SST', 'PSL'" >> user_nl_cam
echo "" >> user_nl_cam
echo "fincl4 =  'PRECC','PRECL'" >> user_nl_cam
echo "" >> user_nl_cam
echo "fincl5 = 'Uzm','Vzm','Wzm','THzm', 'VTHzm','WTHzm','UVzm','UWzm'" >> user_nl_cam
echo "phys_grid_ctem_nfreq=-6" >> user_nl_cam
echo "phys_grid_ctem_zm_nbas=120" >> user_nl_cam
echo "phys_grid_ctem_za_nlat=90" >> user_nl_cam
echo "" >> user_nl_cam
echo " rad_climate= 'A:Q:H2O', 'N:O2:O2'," >> user_nl_cam
echo "         'A:CO2:CO2', 'N:ozone:O3'," >> user_nl_cam
echo "         'A:N2O:N2O', 'A:CH4:CH4'," >> user_nl_cam
echo "         'N:CFC11STAR:CFC11', 'A:CFC12:CFC12'," >> user_nl_cam
echo "         'M:mam4_mode1:/glade/campaign/cesm/cesmdata/inputdata/atm/cam/physprops/mam4_mode1_rrtmg_aeronetdust_sig1.6_dgnh.48_c140304.nc', 'M:mam4_mode2:/glade/campaign/cesm/cesmdata/inputdata/atm/cam/physprops/mam4_mode2_rrtmg_aitkendust_c141106.nc'," >> user_nl_cam
echo "         'M:mam4_mode3:/glade/campaign/cesm/cesmdata/inputdata/atm/cam/physprops/mam4_mode3_rrtmg_aeronetdust_c141106.nc', 'M:mam4_mode4:/glade/campaign/cesm/cesmdata/inputdata/atm/cam/physprops/mam4_mode4_rrtmg_c130628.nc'," >> user_nl_cam
echo "         'N:VOLC_MMR1:/glade/campaign/cesm/cesmdata/inputdata/atm/cam/physprops/volc_camRRTMG_byradius_sigma1.6_mode1_c210211.nc', 'N:VOLC_MMR2:/glade/campaign/cesm/cesmdata/inputdata/atm/cam/physprops/volc_camRRTMG_byradius_sigma1.6_mode2_c210211.nc'," >> user_nl_cam
echo "         'N:VOLC_MMR3:/glade/campaign/cesm/cesmdata/inputdata/atm/cam/physprops/volc_camRRTMG_byradius_sigma1.2_mode3_c210211.nc'" >> user_nl_cam
echo "" >> user_nl_cam
echo " seasalt_emis_scale             = 0.75D0" >> user_nl_cam
echo "" >> user_nl_cam
echo "clubb_c8=  4.95 " >> user_nl_cam
echo "" >> user_nl_cam
echo "" >> user_nl_cam
echo "enthalpy_flux_method = 1" >> user_nl_cam
echo "" >> user_nl_cam
echo "thermo_budget_histfile_num             = 8" >> user_nl_cam
echo "thermo_budget_history          = .true." >> user_nl_cam
echo "avgflag_pertape(8) = 'N'" >> user_nl_cam
echo "" >> user_nl_cam
echo "avgflag_pertape(9) = 'A'" >> user_nl_cam
echo "fincl9 =  'enth_prec_ac_hice','enth_prec_ac_hliq','enth_prec_bc_hice','enth_prec_bc_hliq','enth_prec_ac_fice','enth_prec_ac_fliq','enth_prec_bc_fice'," >> user_nl_cam
echo " 'enth_prec_bc_fliq','enth_evap_hevap','cpice_srf','te_tnd','te_lat','ls_srf','lf_srf','dEdt_dme','dEdt_physics',         " >> user_nl_cam
echo "'dEdt_cpdycore','residual','dEdt_enth_fix','enth_fix_fct_bc_tot','enth_fix_fct_ac_tot','enthalpy_heating_fix_bc','enthalpy_heating_fix_ac',   " >> user_nl_cam
echo "'dEdt_efix_physics','EFIX','enth_flux_to_not_ocn','enth_flux_to_ocn','hsnow_liq_ref','hrain_liq_ref','hevap_liq_ref','OCNFRAC','ICEFRAC','LANDFRAC'" >> user_nl_cam
echo "avgflag_pertape(10) = 'A'" >> user_nl_cam
echo "fincl10 =  'enth_prec_ac_hice','enth_prec_ac_hliq','enth_prec_bc_hice','enth_prec_bc_hliq','enth_prec_ac_fice','enth_prec_ac_fliq','enth_prec_bc_fice'," >> user_nl_cam
echo " 'enth_prec_bc_fliq','enth_evap_hevap','cpice_srf','te_tnd','te_lat','ls_srf','lf_srf','dEdt_dme','dEdt_physics',         " >> user_nl_cam
echo "'dEdt_cpdycore','residual','dEdt_enth_fix','enth_fix_fct_bc_tot','enth_fix_fct_ac_tot','enthalpy_heating_fix_bc','enthalpy_heating_fix_ac',   " >> user_nl_cam
echo "'dEdt_efix_physics','EFIX','enth_flux_to_not_ocn','enth_flux_to_ocn','hsnow_liq_ref','hrain_liq_ref','hevap_liq_ref','OCNFRAC','ICEFRAC','LANDFRAC'" >> user_nl_cam


#qcmd -A $proj -- ./case.build 
#./case.submit
