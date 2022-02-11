#!/bin/tcsh
../cime/scripts/create_newcase --compset FADIAB -res ne30_ne30_mg16 --pecount 108 --case /glade/scratch/pel/SEne30.L32.fadiab.acid --run-unsupported --project P93300642
cd /glade/scratch/pel/SEne30.L32.fadiab.acid
./xmlchange DOUT_S=FALSE
./xmlchange STOP_OPTION=ndays,STOP_N=4
./xmlquery CAM_CONFIG_OPTS
./xmlchange --file env_build.xml --id CAM_CONFIG_OPTS --val "--phys adiabatic --analytic_ic"
./xmlquery CAM_CONFIG_OPTS
./xmlchange JOB_WALLCLOCK_TIME=00:15:00
./pelayout
./case.setup
#cp /glade/u/home/cjablono/CAM_Jan22_cases/project_2_source_code/acid_test.F90 /glade/scratch/pel/CAM_Jan22.SEne30.L30_400m.fadiab.acid/SourceMods/src.cam/ic_baro_dry_jw06.F90

#Use the following namelist settings for user_nl_cam:
echo "empty_htapes     = .TRUE.   " >> user_nl_cam
echo "avgflag_pertape  = 'I'      " >> user_nl_cam
echo "MFILT            = 180      " >> user_nl_cam
echo "NDENS            = 2        " >> user_nl_cam
echo "fincl1 = 'PS','T','U','V','OMEGA','T850','U850','V850','OMEGA850','T700','T500','OMEGA500','U500','V500','PHIS','PSL','Z3','Z700','Z500','Z300'" >> user_nl_cam
echo "NHTFRQ           = -6       " >> user_nl_cam
echo "analytic_ic_type = 'dry_baroclinic_wave_jw2006'" >> user_nl_cam
echo "omega = 0.                  " >> user_nl_cam
echo "interpolate_output = .true. " >> user_nl_cam
echo "se_statefreq           = 1  " >> user_nl_cam


#NCDATA  	= '/glade/u/home/cjablono/CESM_vertical_grids/cam_vcoords_L30_dz400m_top_12km.nc"
./preview_namelists
qcmd -A P93300642 -- ./case.build
./case.submit
