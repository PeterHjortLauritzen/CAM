#!/usr/bin/env python3
      
import os, sys
import subprocess
from netCDF4 import Dataset
from itertools import islice
# This script must be run in an environment with netCDF4 python libraries 
# active and the ability to build cesm cases.
# On Cheyenne, the best way to do this is from an interactive session.
# > qsub -X -I -l select=1:ncpus=36:mpiprocs=36 -l walltime=02:30:00 -q regular -A P93300642
# > conda activate my_env
# > ./build_ppe_CAM_cases.py

# Edit below to set up your cases
#cesmroot = '/glade/u/home/trude/src/cam6_3_026_ppe/'
cesmroot = '/glade/work/katec/mg3work/peterLau_CAM_fork'
basecasename = "PPE7_test_4"
baseroot = os.path.join("/glade/work/katec/mg3work","ppe_cases",basecasename)
res = "ne30pg3_ne30pg3_mg17"
compset = "FLTHIST"

paramfile = "./CESM_parameter_229_231.nc"


# ++TE: do not change ensemble_startval
ensemble_startval = "001" # The startval strings should be the same length, or else. 
basecase_startval = "000"
project = "P93300642"
wall_time = "12:00:00"
ntasks = 2160

specify_user_num_sims = False # If True, use "user_num_sims" number of sims. This will
                             # use the first X values of each parameter, and will 
                             # probably crash if you specify more simulations than
                             # parameter values in your paramfile. 
                             # If False, automatically create the same number of cases
                             # as param values in the paramfile. This is the safest 
                             # option.
user_num_sims = 2

specify_sstice_file = True   # If true, the base case will be updated with the following
                             # SSTICE file settings, if false, these will not be used.
sstice_filename = "$DIN_LOC_ROOT/atm/cam/sst/sst_HadOIBl_bc_1x1_2000climo_c180511.nc"
sstice_yr_align = 1
sstice_yr_start = 0
sstice_yr_end   = 0

# Currently only number of years set (no months or days)
# Will be the same in every case

stop_yrs = 3
rest_yrs = 1
run_type = "hybrid"
start_date = "1999-01-01"

# REFCASE settings below are optional
# Will be the same in every case
ref_case_get = True
ref_case_name = "f.cam6_3_107.FLTHIST_v0a.ne30.clm5_1.001"
ref_case_path = "cesm2_init"
ref_case_date = "1994-01-01"

## CAM_CONFIG_OPTS
cam_conopts = "-cosp"

## CAM SourceMods
camSourceModsPath = "/glade/work/pel/cesm_tags/PPE/SourceMods/src.cam/*.F90"
clmSourceModsPath = "/glade/campaign/cesm/cesmdata/cseg/runs/cesm2_0/f.cam6_3_119.FLTHIST_ne30.r328_gamma0.33_soae.001/SourceMods/src.clm/*"

## User_nl_cam text that will be cloned in every ensemble case
## Text in user_nl_cam will appear exactly as in this string specifier
user_nl_cam_string = """
dust_emis_fact = 1.3
mfilt    =       0,       5,     20,      40,      12,       120,      1,   1
nhtfrq              =       0,     -24,    -24,      -3,       0,       -2,      0,  -8760
ndens               =       2,       2,      2,       2,       2,       1,      2,   1
interpolate_output  =  .true.,  .true., .true., .false., .false., .true.,  .true.
interpolate_nlat    =     192,     192,    192,     192,     192,     192,   192
interpolate_nlon    =     288,     288,    288,     288,     288,     288,   288

empty_htapes = .true.

fincl1 = 'ACTNI', 'ACTNL', 'ACTREI', 'ACTREL', 'AODDUST', 'AODVIS', 'AODVISdn','BURDENBC', 'BURDENDUST', 'BURDENPOM', 'BURDENSEASALT',
'BURDENSO4', 'BURDENSOA', 'CAPE', 'CCN3', 'CDNUMC', 'CH4', 'CLDHGH', 'CLDICE', 'CLDLIQ', 'CLDLOW', 'CLDMED', 'CLDTOT', 'CLOUD', 'CMFMC_DP',
'CT_H2O', 'DCQ', 'DQCORE', 'DTCOND', 'DTCORE', 'DTV', 'EVAPPREC', 'EVAPSNOW', 'FCTI', 'FCTL', 'FICE', 'FLDS', 'FLNS', 'FLNSC', 'FLNT', 'FLNTC', 'FLUT',
'FREQZM', 'FSDS', 'FSDSC', 'FSNS', 'FSNSC', 'FSNT', 'FSNTC', 'FSNTOA', 'ICEFRAC', 'LANDFRAC', 'LHFLX', 'LWCF', 'MPDICE', 'MPDLIQ', 'MPDQ', 'MPDT',
'OCNFRAC', 'OMEGA', 'OMEGA500', 'PBLH', 'PHIS', 'PINT', 'PMID', 'PRECC', 'PRECL', 'PRECSC', 'PRECSL', 'PRECT', 'PS', 'PSL', 'PTEQ', 'PTTEND', 'Q',
'QFLX', 'QRL', 'QRS', 'QTGW', 'RCMTEND_CLUBB', 'RELHUM', 'RVMTEND_CLUBB', 'SHFLX', 'SOLIN', 'SST', 'STEND_CLUBB', 'SWCF',
'T', 'TAUX', 'TAUY', 'TFIX', 'TGCLDIWP', 'TGCLDLWP', 'TMQ', 'TREFHT', 'TS', 'TTGW', 'U', 'U10', 'UBOT', 'UTGWORO', 'UTGW_TOTAL',
'V', 'VBOT', 'VTGWORO', 'VTGW_TOTAL', 'WPRTP_CLUBB', 'WPTHLP_CLUBB', 'Z3', 'ZMDQ', 'ZMDT', 'N2O', 'CO2','CFC11','CFC12',
'CLD_MISR','FISCCP1_COSP','CLD_CAL','CLD_MISR','CLDTOT_CAL','CLDHGH_CAL', 'CLDMED_CAL','CLDLOW_CAL','CLMODIS', 'AODVISdn', 'AODDUSTdn',
'CCN3', 'CDNUMC', 'H2O', 'NUMICE', 'NUMLIQ'

fincl3 = 'PRECT', 'PRECC', 'FLUT', 'U850', 'U200', 'V850', 'V200', 'OMEGA', 'PSL'

fincl4 =  'PRECC','PRECL'

fincl5 = 'Uzm','Vzm','Wzm','THzm', 'VTHzm','WTHzm','UVzm','UWzm'
phys_grid_ctem_nfreq=-6
phys_grid_ctem_zm_nbas=120
phys_grid_ctem_za_nlat=90

fincl7= 'AQSO4_H2O2','AQSO4_O3', 'bc_a1', 'bc_a4', 'dst_a1', 'dst_a2', 'dst_a3', 'ncl_a1',
'ncl_a1', 'ncl_a2', 'ncl_a3', 'pom_a1', 'pom_a4', 'so4_a1', 'so4_a2', 'so4_a3',
'soa_a1', 'num_a1', 'num_a2', 'num_a3', 'num_a4',
'bc_a1SFWET', 'bc_a4SFWET', 'dst_a1SFWET', 'dst_a2SFWET', 'dst_a3SFWET', 'ncl_a1SFWET',
'ncl_a2SFWET', 'ncl_a3SFWET', 'pom_a1SFWET', 'pom_a4SFWET', 'so4_a1SFWET', 'so4_a2SFWET', 'so4_a3SFWET', 'soa_a1SFWET',
'soa_a2SFWET', 'bc_c1SFWET', 'bc_c4SFWET', 'dst_c1SFWET', 'dst_c2SFWET', 'dst_c3SFWET', 'ncl_c1SFWET', 'ncl_c2SFWET',
'ncl_c3SFWET', 'pom_c1SFWET', 'pom_c4SFWET', 'so4_c1SFWET', 'so4_c2SFWET', 'so4_c3SFWET', 'soa_c1SFWET', 'soa_c2SFWET',
'bc_a1DDF', 'bc_a4DDF', 'dst_a1DDF', 'dst_a2DDF', 'dst_a3DDF', 'ncl_a1DDF', 'ncl_a2DDF', 'ncl_a3DDF',
'pom_a1DDF', 'pom_a4DDF', 'so4_a1DDF', 'so4_a2DDF', 'so4_a3DDF', 'soa_a1DDF', 'soa_a2DDF',
'so4_a1_CLXF', 'so4_a2_CLXF', 'SFbc_a4', 'SFpom_a4', 'SFso4_a1', 'SFso4_a2',
'so4_a1_sfgaex1', 'so4_a2_sfgaex1', 'so4_a3_sfgaex1', 'soa_a1_sfgaex1', 'soa_a2_sfgaex1',
'SFdst_a1','SFdst_a2', 'SFdst_a3', 'SFncl_a1', 'SFncl_a2', 'SFncl_a3',
'num_a2_sfnnuc1', 'SFSO2', 'OCN_FLUX_DMS', 'SAD_SULFC', 'SAD_TROP', 'SAD_AERO'

clubb_C1                           =   0.800000000000000
clubb_C1b                          =   0.800000000000000
clubb_C1c                          =   0.750000000000000                                                      
clubb_C4                           =    1.00000000000000                                           
clubb_C_uu_shr                     =   0.100000000000000                                               
clubb_C_uu_buoy                    =   0.000000000000000                                           
clubb_C6rt                         =    1.00000000000000                                                  
clubb_C6rtb                        =    1.00000000000000                                                 
clubb_C6rtc                        =   0.500000000000000                                          
clubb_C6thl                        =    1.00000000000000                                           
clubb_C6thlb                       =    1.00000000000000                                           
clubb_C6thlc                       =   0.500000000000000                                  
clubb_C7b                          =   0.700000000000000                                       
clubb_C8                           =   0.800000000000000                                        
clubb_C11                          =   0.500000000000000                                        
clubb_C11b                         =   0.800000000000000                                    
clubb_C11c                         =   0.850000000000000                                     
clubb_C14                          =   0.800000000000000                                   
clubb_C_wp3_pr_turb                =    1.00000000000000                                   
clubb_wpxp_L_thresh                =    100.000000000000                                  
clubb_c_K1                         =    1.00000000000000                                     
clubb_nu1                          =    10.0000000000000                                       
clubb_c_K2                         =   0.100000000000000                               
clubb_nu2                          =    1.00000000000000                                    
clubb_c_K8                         =    10.0000000000000                                  
clubb_gamma_coef                   =   0.300000000000000                                      
clubb_gamma_coefb                  =   0.300000000000000                                      
clubb_gamma_coefc                  =    1.20000000000000                                
clubb_mu                           =   5.000000000000000E-004                              
clubb_beta                         =    2.00000000000000
clubb_lambda0_stability_coef       =   3.000000000000000E-002                                    
clubb_c_K10                        =   0.300000000000000                                         
clubb_c_K10h                       =   0.350000000000000                                       
clubb_C_invrs_tau_bkgnd            =    1.10000000000000                                      
clubb_C_invrs_tau_sfc              =   5.000000000000000E-002                                        
clubb_C_invrs_tau_shear            =   0.220000000000000                               
clubb_C_invrs_tau_N2               =   0.300000000000000                                      
clubb_C_invrs_tau_N2_xp2           =   0.000000000000000E+000                            
clubb_C_invrs_tau_N2_wpxp          =    2.00000000000000                                       
clubb_C_invrs_tau_N2_clear_wp3     =    2.00000000000000                                          
clubb_C_invrs_tau_wpxp_Ri          =    3.00000000000000                                           
clubb_C_invrs_tau_wpxp_N2_thresh   =   3.500000000000000E-004                                     
clubb_altitude_threshold           =    150.000000000000                                            
clubb_a3_coef_min                  =    1.60000000000000                                          
clubb_l_diag_Lscale_from_tau=.true.                                                      
clubb_l_stability_correct_tau_zm=.false.                                       
clubb_l_min_wp2_from_corr_wx=.false.                                           
clubb_l_min_xp2_from_corr_wx=.false.                               
clubb_l_vert_avg_closure=.false.                                     
clubb_l_use_cloud_cover=.false.                                      
clubb_l_damp_wp2_using_em=.true.                                  
clubb_l_use_C7_Richardson=.true.                                   
clubb_l_rcm_supersat_adj=.true.                                       
clubb_l_damp_wp3_Skw_squared=.true.                           
clubb_l_enable_relaxed_clipping=.true. 

"""
## No /EOF needed for python

## User_nl_clm text that will be cloned in every ensemble case
## Text in user_nl_cam will appear exactly as in this string specifier
user_nl_clm_string = """
use_init_interp = .true.
check_finidat_year_consistency = .false.
fsurdat = '/glade/work/slevis/git/mksurfdata_toolchain/tools/mksurfdata_esmf/surfdata_ne30np4.pg3_SSP5-8.5_78pfts_CMIP6_1850-2100_c230227.nc'
flanduse_timeseries = '/glade/work/slevis/git/mksurfdata_toolchain/tools/mksurfdata_esmf/landuse.timeseries_ne30np4.pg3_SSP5-8.5_78_CMIP6_1850-2100_c230227.nc'
"""
## No /EOF needed for python

## Should not be a need to edit below this line ##
if cesmroot is None:
    raise SystemExit("ERROR: CESM_ROOT must be defined in environment")
#_LIBDIR = os.path.join(cesmroot,"cime","scripts","Tools")
#sys.path.append(_LIBDIR)
#_LIBDIR = os.path.join(cesmroot,"cime","scripts","lib")
#sys.path.append(_LIBDIR)


_LIBDIR = os.path.join(cesmroot,"cime","CIME","Tools")
sys.path.append(_LIBDIR)
_LIBDIR = os.path.join(cesmroot,"cime")
sys.path.append(_LIBDIR)




import datetime, glob, shutil
import CIME.build as build
from standard_script_setup import *
from CIME.case             import Case
from CIME.utils            import safe_copy
from argparse              import RawTextHelpFormatter
from CIME.locked_files          import lock_file, unlock_file

def per_run_case_updates(case, ensemble_str, nint, paramdict, ensemble_idx):
    print(">>>>> BUILDING CLONE CASE...")
    caseroot = case.get_value("CASEROOT")
    basecasename = os.path.basename(caseroot)[:-nint]

    unlock_file("env_case.xml",caseroot=caseroot)
    casename = basecasename+ensemble_str
    case.set_value("CASE",casename)
    rundir = case.get_value("RUNDIR")
    rundir = rundir[:-nint]+ensemble_str
    case.set_value("RUNDIR",rundir)
    case.flush()
    lock_file("env_case.xml",caseroot=caseroot)
    print("...Casename is {}".format(casename))
    print("...Caseroot is {}".format(caseroot))
    print("...Rundir is {}".format(rundir))

    # Add user_nl updates for each run                                                        

    paramLines = []
#    ens_idx = int(ensemble_str)-int(ensemble_startval)
    ens_idx = int(ensemble_idx)-int(ensemble_startval)
    print('ensemble_index')
    print(ens_idx)
    for var in paramdict.keys():
        paramLines.append("{} = {}\n".format(var,paramdict[var][ens_idx]))

    usernlfile = os.path.join(caseroot,"user_nl_cam")
    print("...Writing to user_nl file: "+usernlfile)
    file1 = open(usernlfile, "a")
    file1.writelines(paramLines)
    file1.close()

    print(">> Clone {} case_setup".format(ensemble_str))
    case.case_setup()
    print(">> Clone {} create_namelists".format(ensemble_str))
    case.create_namelists()
    print(">> Clone {} submit".format(ensemble_str))
    case.submit()


def build_base_case(baseroot, basecasename, res, compset, overwrite):
    print(">>>>> BUILDING BASE CASE...")
    caseroot = os.path.join(baseroot,basecasename+'.'+basecase_startval)
    if overwrite and os.path.isdir(caseroot):
        shutil.rmtree(caseroot)
    with Case(caseroot, read_only=False) as case:
        if not os.path.isdir(caseroot):
            case.create(os.path.basename(caseroot), cesmroot, compset, res,
                        machine_name="cheyenne", driver="nuopc",
                        run_unsupported=True, answer="r",walltime=wall_time, 
                        project=project, queue="regular")
            # make sure that changing the casename will not affect these variables                                           
            case.set_value("EXEROOT",case.get_value("EXEROOT", resolved=True))
            case.set_value("RUNDIR",case.get_value("RUNDIR",resolved=True)+".00")

        case.set_value("RUN_TYPE",run_type)
        case.set_value("GET_REFCASE",ref_case_get)
        case.set_value("RUN_REFCASE",ref_case_name)
        case.set_value("RUN_REFDIR",ref_case_path)
        case.set_value("RUN_REFDATE",ref_case_date)
        case.set_value("STOP_OPTION","nyears")
        case.set_value("STOP_N",stop_yrs)
        case.set_value("REST_OPTION","nyears")
        case.set_value("REST_N",rest_yrs)
        case.set_value("RUN_STARTDATE",start_date)

        case.set_value("NTASKS",ntasks)
        case.set_value("NTASKS_GLC",1)
        case.set_value("NTASKS_WAV",1)
        case.set_value("NTASKS_ESP",1)

        case.set_value("CAM_CONFIG_OPTS", 
                       case.get_value("CAM_CONFIG_OPTS",resolved=True)+' '+cam_conopts)
        if specify_sstice_file:
            case.set_value("SSTICE_DATA_FILENAME",sstice_filename)
            case.set_value("SSTICE_YEAR_ALIGN",sstice_yr_align)
            case.set_value("SSTICE_YEAR_START",sstice_yr_start)
            case.set_value("SSTICE_YEAR_END",sstice_yr_end)

#            case.set_value("QUEUE","economy")

        rundir = case.get_value("RUNDIR")
        caseroot = case.get_value("CASEROOT")
        
        print(">> base case_setup...")
        case.case_setup()

        print(">> base SourceMods copy...")
        camSrcPath = os.path.join(caseroot,"SourceMods","src.cam")
        subprocess.check_call("cp "+camSourceModsPath+" "+camSrcPath, shell=True)
        clmSrcPath = os.path.join(caseroot,"SourceMods","src.clm")
        subprocess.check_call("cp "+clmSourceModsPath+" "+clmSrcPath, shell=True)
        
        print(">> base case write user_nl_cam...")
        usernlfile = os.path.join(caseroot,"user_nl_cam")
        funl = open(usernlfile, "a")
        funl.write(user_nl_cam_string)
        funl.close()

        print(">> base case write user_nl_clm...")
        usernlfile = os.path.join(caseroot,"user_nl_clm")
        funl = open(usernlfile, "a")
        funl.write(user_nl_clm_string)
        funl.close()
        
        print(">> base case_build...")
        #build.case_build(caseroot, case=case)
        subprocess.check_call("./preview_namelists", shell=True)
        subprocess.check_call("./case.build", shell=True)

        return caseroot

def clone_base_case(caseroot, ensemble, overwrite, paramdict, ensemble_num):
    print(">>>>> CLONING BASE CASE...")
    print(ensemble)
    startval = ensemble_startval
    nint = len(ensemble_startval)
    cloneroot = caseroot
    for i in range(int(startval), int(startval)+ensemble):
        ensemble_idx = '{{0:0{0:d}d}}'.format(nint).format(i)
        member_string = ensemble_num[i-1]
        print("member_string="+member_string)
        if ensemble > 1:
            caseroot = caseroot[:-nint] + member_string
        if overwrite and os.path.isdir(caseroot):
            shutil.rmtree(caseroot)
        if not os.path.isdir(caseroot):
            with Case(cloneroot, read_only=False) as clone:
                clone.create_clone(caseroot, keepexe=True)
        with Case(caseroot, read_only=False) as case:
            per_run_case_updates(case, member_string, nint, paramdict, ensemble_idx)


def take(n, iterable):
    "Return first n items of the iterable as a list"
    return list(islice(iterable, n))


def _main_func(description):

    print ("Starting SCAM PPE case creation, building, and submission script")
    print ("Base case name is {}".format(basecasename))
    print ("Parameter file is "+paramfile)

    overwrite = True

    # read in NetCDF parameter file
    inptrs = Dataset(paramfile,'r')
    print ("Variables in paramfile:")
    print (inptrs.variables.keys())
    print ("Dimensions in paramfile:")
    print (inptrs.dimensions.keys())
    num_sims = inptrs.dimensions['nmb_sim'].size
    num_vars = len(inptrs.variables.keys())-1
    ensemble_num = inptrs['Sample_nmb']
    print('ensemble_num')
    print(ensemble_num)
    print(ensemble_num[0])
    if specify_user_num_sims:
        num_sims = user_num_sims

    print ("Number of sims = {}".format(num_sims))
    print ("Number of params = {}".format(num_vars))


    # Save a pointer to the netcdf variables
    paramdict = inptrs.variables
    del paramdict['Sample_nmb']

    # Create and build the base case that all PPE cases are cloned from
    caseroot = build_base_case(baseroot, basecasename, res,
                                   compset, overwrite)

    # Pass in a dictionary with all of the parameters and their values
    # for each PPE simulation 
    # This code clones the base case, using the same build as the base case,
    # Adds the namelist parameters from the paramfile to the user_nl_cam 
    # file of each new case, does a case.set up, builds namelists, and 
    # submits the runs.

    clone_base_case(caseroot, num_sims, overwrite, paramdict, ensemble_num)

    inptrs.close()

if __name__ == "__main__":
    _main_func(__doc__)

