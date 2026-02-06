#!/bin/csh -fv
set echo
#**********************************************************************
# Run SCAM with a single IOP
#    This script will build and run one IOP
#    If a user wishes to run more than one IOP, use create_scam6_iop_multi
#
# Usage:
#   ./create_scam6_iop <IOP>    # where IOP name is from list below
#      - or -
#   ./create_scam6_iop          # IOP is specified in the script below
#**********************************************************************

#### BEWARE OF ###
# - Have to blow away case directory to rebuild.

#------------------
# User sets options in this section
#------------------
### Case Name

set CASETITLE=cesm3_planets.020526.L32.4node.smtopo.113

### Full path of cesm source code and case (output) directories (see examples)

###set CESMDIR=/project/amp/jet/collections/CESM3-planets-gwdiff.010826
set CESMDIR=/project/amp/jet/collections/testgw/CESM3-planets-gwdiff.020526
set CAMDIR=${CESMDIR}/components/cam
set CASEDIR=/project/amp/$USER/cases

### Set location of user source mods (if any)
setenv this_dir  `pwd`
setenv usrsrc  ${this_dir}/mods/$CASETITLE

### Standard Run Settings
set RES=ne16_ne16
set COMPSET=FMARSEXORT
set COMPILER=intel
set DOUT_S=FALSE
set PES='--pecount 192'
#set PES=''
#------------------
# User should not need to set any options in this section
#------------------

cd  $CASEDIR


##set MODSDIR = $CESMDIR/components/cam/cime_config/usermods_dirs
#Create full casename
set CASENAME=${COMPSET}.${RES}.${CASETITLE}

#------------------
# create case
#------------------

#$CESMDIR/cime/scripts/create_newcase --compset $COMPSET  --res $RES --compiler $COMPILER --case $CASEDIR/$CASENAME  --user-mods-dir $MODSDIR --run-unsupported
$CESMDIR/cime/scripts/create_newcase --compset $COMPSET  --res $RES --compiler $COMPILER --case $CASEDIR/$CASENAME --run-unsupported $PES

cd  $CASEDIR/$CASENAME

### Set build and run directories to be under case directory.

set RUNDIR=${CASEDIR}/${CASENAME}/run
./xmlchange RUNDIR=$RUNDIR
./xmlchange EXEROOT=${CASEDIR}/${CASENAME}/bld

#------------------
# XMLCHANGE OPTIONS HERE
#------------------

# DEBUG
#./xmlchange DEBUG='TRUE'
./xmlchange STOP_OPTION=ndays
./xmlchange STOP_N=10 # Total timesteps for run f(dt)
./xmlchange DOUT_S=${DOUT_S}              ### Dont run archive script leave files in run directory
./xmlchange GMAKE_J=48
# nadv is set to the number of advected constituents we added for simple mars - right now co2, n2 and q
#./xmlchange --append CAM_CONFIG_OPTS="-chem none -nadv=2"
./xmlchange --append CAM_CONFIG_OPTS='-nlev 32'
./xmlchange --append CLM_BLDNML_OPTS=-ignore_warnings
./xmlchange CLM_FORCE_COLDSTART=on
#------------------
# Setup Case
#------------------

#./case.setup
./case.setup

#------------------
#  source mods: copy them into case directory
#------------------

/bin/cp  /project/amp/jet/cases/FMARSEXORT.ne16_ne16.cesm3_planets.020526.L32.4node.smtopo.112/SourceMods/src.cam/* SourceMods/src.cam/
#sed -i "s/0.7_r8\*phis_tmp/0.7_r8\*phis_tmp/" SourceMods/src.cam/dyn_comp.F90
#------------------
# Add all user specific cam namelist changes here
#
# Users should add all user specific namelist changes below in the form of
#    namelist_var = new_namelist_value
# Namelist settings which appear in usermods_dir and here will use the values
#    specified below
# Other namelist settings from usermods_dirs will be unchanged
# Output can also be specified here (e.g. fincl1)
#------------------
#ncdata = '/glade/p/cgd/amp/jet/mars_cam_vcoords_L32_c220317.nc'
#analytic_ic_type='dry_baroclinic_wave_dcmip2016'
#se_nsplit = 6
#

cat >> user_nl_cam << EOF
ncdata = '/project/amp/jet/cases/FMARSEXORT.ne16_ne16.cesm3_planets.020526.L32.4node.smtopo.112/run/FMARSEXORT.ne16_ne16.cesm3_planets.020526.L32.4node.smtopo.112.cam.i.0001-01-10-00000.nc'
bnd_topo = '/fs/cgd/csm/inputdata/atm/cam/topo/se/ne16np4_mars_nc3000_Laplace0300_20260129.nc'
use_gw_oro        = .true.
effgw_oro         = 1.0D0
use_topo_file		= .true.
se_nsplit         = 8
use_gw_rdg_beta          = .true.
effgw_rdg_beta           = 1.0      ! Start here, increase if jet grows
rdg_beta_cd_llb          = 1.0      ! Low level blocking
se_sponge_del4_lev            = 20
se_sponge_del4_nu_div_fac = 7.5
se_hypervis_subcycle_sponge = 4
se_statefreq       = 1
se_nu_top                = 0.50E+07 ! Increased from 0.1E7 to handle 270 m/s jet
do_exo_rt_spectral		= .false.
inithist = 'DAILY'
nhtfrq = -24
interpolate_output = .false.,.false.,.false.
interpolate_nlat   = 96,96,96
interpolate_nlon   = 144,144,144
mfilt = 10
write_nstep0 = .true.
avgflag_pertape(1) = 'I'
fincl1='TS','PS','USTAR','KVH','KVM','KVT','CGS','DTVKE','DTV','DUV','DVV',
        'qt_pre_PBL','qt_pre_PBL',   'sl_pre_PBL',   'slv_pre_PBL',
        'u_pre_PBL',    'v_pre_PBL',    'qv_pre_PBL',   'ql_pre_PBL',
        'qi_pre_PBL',   't_pre_PBL',    'rh_pre_PBL',   'qt_aft_PBL',
        'sl_aft_PBL',   'slv_aft_PBL',  'u_aft_PBL',    'v_aft_PBL',
        'qv_aft_PBL',   'ql_aft_PBL',   'qi_aft_PBL',   't_aft_PBL',
        'rh_aft_PBL',   'slflx_PBL',    'qtflx_PBL',    'uflx_PBL',
        'vflx_PBL',     'slflx_cg_PBL', 'qtflx_cg_PBL', 'uflx_cg_PBL',
        'vflx_cg_PBL',  'qtten_PBL',    'slten_PBL',    'uten_PBL',
        'vten_PBL',     'qvten_PBL',    'qlten_PBL',    'qiten_PBL',
        'tten_PBL',     'rhten_PBL'
EOF

cat >> user_nl_clm << EOFclm
! --- 1. Surface Data & Initialization ---
! Point to the file generated by your Python script
fsurdat = "/fs/cgd/csm/inputdata/lnd/clm2/surfdata_esmf/ctsm5.4.0/surfdata_ne16np4_mars_desert.c250204.nc"

! Disable interpolation since we are cold-starting
use_init_interp = .false.

! --- 2. Disable Earth-Centric Physics ---
! Turn off dynamic vegetation and photosynthesis
use_luna        = .false.  ! Leaf Utilization of Nitrogen
use_hydrstress  = .false.  ! Plant hydraulic stress
use_bedrock     = .true.   ! Keep true (Mars has regolith depth)

! Turn off Urban models
urban_hac       = 'OFF'    ! No AC/Heating
urban_traffic   = .false.

! Turn off Hydrology features irrelevant for desert
irrigate        = .false.
use_excess_ice  = .true.   ! Keep TRUE. Allows ground ice to form/persist.

! --- 4. Tuning ---
! Reduce roughness length (z0) to bare soil defaults
! (Optional, helps stability in high winds)
z0param_method = 'ZengWang2007'
dust_emis_method= 'Zender_2003'
calc_human_stress_indices='NONE'

! --- Turn off Earth-specific Hydrology/Glaciology ---
use_excess_ice          = .false.   ! Disables reading Earth 'exice_init' file [cite: 1083]
use_excess_ice_streams  = .false.   ! Stop the stream reading [cite: 1099]

! --- Turn off Earth Human Systems ---
!create_crop_landunit    = .false.   ! No agriculture on Mars [cite: 1076]
use_crop                = .false.   ! Ensure crop model is off [cite: 1083]
urban_hac               = 'OFF'     ! Disable urban air conditioning/heating [cite: 1093]
urban_traffic           = .false.   ! No traffic flux [cite: 1093]
collapse_urban          = .true.    ! Collapse urban patches to avoid reading urban data

! --- Disable Earth Data Streams ---
! Setting these to empty strings prevents the model from trying to open Earth files
cropcals_rx = .false.
cropcals_rx_adapt = .false.
! --- Disable Plant/Canopy Physics ---
! Even if PFTs=0, these flags consume compute or memory
use_clm5_fpi            = .false.   ! Turn off foliage interception [cite: 1089]
use_biomass_heat_storage = .false.  ! No biomass to store heat [cite: 1092]
modifyphoto_and_lmr_forcrop = .false. !  [cite: 1098]

! --- Disable Earth Aerosols on Snow ---
! SNICAR expects Earth dust/black carbon.
!snicar_use_aerosol      = .false.   ! Disable aerosol deposition on snow/ice [cite: 1082]
do_sno_oc               = .false.   ! No organic carbon in snow [cite: 1077]
snicar_snobc_intmix     = .false.   ! No black carbon mixing [cite: 1081]
EOFclm

cat >> user_nl_cpl << EOFcpl
orb_mode   = 'fixed_parameters'
orb_obliq  = 25.19
orb_eccen  = 0.0934
orb_mvelp  = 251.0
EOFcpl

#------------------
# Build
#------------------

./case.build

### make timing dir kludge  [REMOVE WHEN FIXED]
mkdir -p $RUNDIR/timing/checkpoints

#------------------
# Choose type of job submission (batch or interactive)
#------------------

### Submit to Queue (If you have one)
#./case.submit
