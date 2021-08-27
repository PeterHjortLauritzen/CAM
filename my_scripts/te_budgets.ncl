load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

begin
  Ttilde  = "T~H-15V6F35~D~FV-6H3~~B~surf~N~"
  Ktilde  = "K~H-15V6F35~D~FV-6H3~~B~surf~N~"


  diags_fname = fname+".te_budgets.txt"
  fname =  systemfunc("ls "+fname)
  fall = addfile(fname,"r")       ; Open netCDF files.  

  TS       = fall->TS(:,:)
  PRECT    = fall->PRECT(:,:)
  SST      = fall->TS(:,:)
  tmp      = fall->T(:,:,:)
  nlev = dimsizes(tmp(0,:,0))
  print("nlev="+nlev)
  delete(tmp)
  TSatm      = fall->T(:,nlev-1,:)

  FTURB    = fall->FTURB(:,:)
  FLH      = fall->FLAT(:,:)
  FLH_T     = fall->FLAT_T(:,:)
  FLHE     = fall->FLATE(:,:)
  FLHP     = fall->FLATP(:,:)
  FH2O     = fall->FNH2O(:,:)
  FH2O     = fall->FNH2O(:,:)
  FH2O     = fall->FNH2O(:,:)
  FWV      = fall->FNWV(:,:)
  FLIQ     = fall->FNLIQ(:,:)
  FICE     = fall->FNICE(:,:)
  FNET_TS    = fall->FNET_TS(:,:)
  FNET_TBOT   = fall->FNET_TBOT(:,:)
  FKE      = fall->FKE(:,:)
  FPHIS    = fall->FPHIS(:,:)
  FTAU     = fall->FTAU(:,:)

  ;
  ; constants from CAM
  ;
  cpair    = 1004.64D0 
  cpliq    = 4188.00D0
  cpwv     = 1810.00D0
  cpice    = 2117.27D0
  LV = 2501000.0D0
  LS = 2834700.0D0
  LF = 333700.0D0

  print("File used for energy diagnostics "+fname)
  print("done input summary")
  print(" ")  
  print("Summary of input to ncl script:")
  print("===============================")
  print(" ")
  if (dycore.eq."SE") then
    print("Your data is for the SE dycore")
    print(" ")
    print("nsplit                   ="+nsplit)
    print("rsplit                   ="+rsplit)
    print("hypervis_subcycle        ="+hypervis_subcycle)
    print("hypervis_subcycle_sponge ="+hypervis_subcycle_sponge)
    print("ftype                    ="+ftype)
    print("qsize_condensate_loading ="+qsize_condensate_loading)
    print("lcp_moist                ="+lcp_moist)
    latlon = False

    ;
    ; check if physics-dynamics variables are present
    ;
    if (isfilevardim(fall,"SE_dED","ncol_d")) then
      print("File appears to have separate grids for physics and dynamics - setting pgx=True")
      pgx = True
    else
      print("File appears to be have the same grid for physics and dynamics - setting pgx=False")
      pgx = False
    end if

    tmp       = fall->SE_dED(:,:)
    nt        = dimsizes(tmp(:,0))
    print("Number of time-levels = "+nt)
    ncol      = dimsizes(tmp(0,:))
    lat1d         = fall->lat          ; Pull off lat/lon variable and
    lon1d        = fall->lon          ; convert both to 1D.
    if (pgx) then
      lat1d_gll     = fall->lat_d          ; Pull off lat/lon variable and
      lon1d_gll     = fall->lon_d          ; convert both to 1D.
    else
      lat1d_gll     = lat1d
      lon1d_gll     = lon1d
    end if


  end if
  if (dycore.eq."FV") then
    print("Your data is for the FV dycore")  
    latlon = True
    lcp_moist = "false"
    tmp       = fall->SE_pAP(:,:,:)
    nt        = dimsizes(tmp(:,0,0))
    ncol      = -1
  end if
  if (dycore.eq."FV3") then
    print("Your data is for the FV3 dycore")
    print(" ")
    print("nsplit                   ="+nsplit)
    print("rsplit                   ="+rsplit)
    print("hypervis_subcycle        ="+hypervis_subcycle)
    print("ftype                    ="+ftype)
    print("qsize_condensate_loading ="+qsize_condensate_loading)
    print("lcp_moist                ="+lcp_moist)
    latlon = False
    pgx = False

    tmp       = fall->SE_dED(:,:)
    nt        = dimsizes(tmp(:,0))
    ncol      = dimsizes(tmp(0,:))
    lat1d_gll     = fall->lat          ; Pull off lat/lon variable and
    lon1d_gll     = fall->lon          ; convert both to 1D.
    lat1d         = fall->lat          ; Pull off lat/lon variable and
    lon1d         = fall->lon          ; convert both to 1D.

  end if


  print("dtime                    ="+dtime)

;
; check
;
  if (lcp_moist.ne."true".and.lcp_moist.ne."false") then
    print("ERROR: lcp_moist should either be true or false; lcp_moist= "+lcp_moist)
    exit()
  end if
  ;
  ; initialization
  ;
  g         = 9.80665D0
  PI        = get_pi("double")
  deg2rad   = PI/180.0
  ncol_p    = ncol
  eps       = 1.0D-13
  area_sphere_p = eps
  if (latlon) then
    lat = fall->lat
    lat = deg2rad*lat
    lon = fall->lon
    lon = deg2rad*lon
    nlat=dimsizes(lat)
    nlon=dimsizes(lon)
    area  = new((/nlat,nlon/),typeof(tmp(0,0,0)))    
    dlat  = lat(3)-lat(2)
    dlat2 = 0.5*dlat    
    dlon  = lon(2)-lon(1)
    j=0
    do i=0,nlon-1
      area(j,i) = (sin(lat(j)+dlat2)-sin(lat(j)))*dlon
    end do    
    do j=1,nlat-2
      dsinlat =  (sin(lat(j)+dlat2)-sin(lat(j)-dlat2))      
      do i=0,nlon-1
        area(j,i) = dsinlat*dlon
      end do
    end do
    j=nlat-1
    do i=0,nlon-1
      area(j,i) = (sin(lat(j))-sin(lat(j)-dlat2))*dlon
    end do
    area_sphere = sum(area)
  else
    if (pgx) then
      tmp_p         = fall->SE_pAM(:,:)
      ncol_p        = dimsizes(tmp_p(0,:))
      area_p        = fall->area_fvm(:  )
      area_sphere_p = sum(area_p(0:ncol_p-1))          
      area          = fall->area_d (:  )
      area_sphere = sum(area(0:ncol-1))                
      print("Area of sphere (GLL)         : "+area_sphere)
      error = area_sphere-4.0*PI
      print("Error spanning the sphere    : "+error)
      print("Area of sphere (physgrid)    : "+area_sphere_p)
      error = area_sphere_p-4.0*PI
      print("Error spanning the sphere    : "+error)
      print("ncol_p="+ncol_p)
    else
      area      = fall->area   (:  )
      area_p    = area
      area_sphere = sum(area(0:ncol-1))          
      area_sphere_p = sum(area_p(0:ncol_p-1))          
    end if  
  end if
  print("xxx"+ncol+" "+ncol_p)
  

  ;
  ; check if physics-dynamics variables are present
  ;
  if (isfilevar(fall,"WI_pBF").and.isfilevar(fall,"WL_pBF")) then
    print("File appears to contain variables from a non-dry run - setting wet=True")
    wet = True
  else
    print("File appears to be from a dry run (no WL and no WI) - setting wet=False")
    wet = False
  end if

  if (isfilevar(fall,"WV_PDC")) then
    vars_pdc  = new((/4,nt,ncol/),typeof(tmp(0,0)))
    int_vars_pdc  = new((/4,nt/) ,typeof(tmp(0,0)))
    pdc_vars_present = True
    print("PDC vars are present on file")
  else
    print("No PDC vars present on file")
    pdc_vars_present = False
  end if
  ;
  ; check if physics-dynamics variables are present
  ;
  if (isfilevar(fall,"TT_pBF")) then
    TT_vars_present = True
    print("TT vars are present on file")
  else
    print("No TT vars present on file")
    TT_vars_present = False
  end if

  print(" ")
  print(" ")


    print("suffix (parameterization side)")
    print("pBF: state passed to parameterizations, before energy fixer")
    print("pBP: after energy fixer, before parameterizations")
    print("pAP: after last phys_update in parameterizations and state saved for energy fixer")
    print("pAM: after dry mass correction calculate but not used in SE")
    print(" ")
    print("history files saved off here")

  print(" ")
  print("suffix (dynamics)")
  print("dED: state from end of previous dynamics (= pBF + time sampling)")
  print("   loop over vertical remapping and physics dribbling -------- (nsplit) -------")
  print("            (dribbling and remapping always done together)                    |")
  print("          dAF: state from previous remapping                                  |")
  print("          dBD: state after physics dribble, before dynamics                   |")
  print("          loop over vertical Lagrangian dynamics --------rsplit-------------  |")
  print("              dynamics here                                                |  |")
  print("              loop over hyperviscosity ----------hypervis_sub------------  |  |")
  print("                 dBH   state before hyperviscosity                      |  |  |")
  print("                 dCH   state after hyperviscosity                       |  |  |")
  print("                 dAH   state after hyperviscosity momentum heating      |  |  |")
  print("              end hyperviscosity loop -----------------------------------  |  |")
  print("              dBS   state before del2 sponge                            |  |  |")
  print("              dAS   state after del2+mom heating sponge                 |  |  |")
  print("          end of vertical Lagrangian dynamics loop -------------------------  |")
  print("      dAD  state after dynamics, before vertical remapping                    |")
  print("      dAR     state after vertical remapping                                  |")
  print("   end of remapping loop ------------------------------------------------------")
  print("dBF  state passed to parameterizations = state after last remapping            ")
  print(" ")
  print(" ")
  print(" ")

;      " end of previous dynamics                           ",& !dED
;      " from previous remapping or state passed to dynamics",& !dAF - state in beginning of nsplit loop
;      " state after applying CAM forcing                   ",& !dBD - state after applyCAMforcing
;      " before vertical remapping                          ",& !dAD - state before vertical remapping
;      " after vertical remapping                           ",& !dAR - state at end of nsplit loop
;      " state passed to parameterizations                  ",& !dBF
;      " state before hypervis                              ",& !dBH
;      " state after hypervis but before adding heating term",& !dCH
;      " state after hypervis                               ",& !dAH
;      " state before sponge layer diffusion                ",& !dBS - state before sponge del2
;      " state after sponge layer diffusion                 ",& !dAS - state after sponge del2
;      " state before two dz filter                         ",& !dBZ
;      " state after two dz filter                          ",& !dAZ
;      " phys2dyn mapping errors (requires ftype-1)         " & !p2d - for asses

  if (wet) then
    qsize = 3
    str_vars      = (/"SE","KE","WV","WL","WI"/);,"TT"/)
    if (isfilevar(fall,"SE_dED")) then
      print("Dycore diagnostics seem to be on file")
      if (dycore.eq."FV3") then
        str_stage     = (/"pBF","pBP","pAP","pAM","dED","dAP","dBD","dAT","dAF","dAD"  ,"dAR" ,"dBF"/)
        subcycle      = (/1.0D0,1.0D0,1.0D0,1.0D0,1.0D0,1.0D0,1.0D0,1.0D0,nsplit,nsplit,nsplit,1.0D0/)
;        subcycle      = (/1D0  ,1D0  ,1D0  ,1D0  ,1.0D0,1.0D0,1.0D0,1D0,1D0,1D0,1D0/)	;averaging is not being done by CAM history
      else
        str_stage     = (/"pBF","pBP","pAP","pAM","dED","dAF","dBD"  ,"dAD" ,"dAR" ,"dBF","dBH","dCH","dAH","dBS","dAS","AP1","BP1","AP2","BP2","AP3","BP3","AP4","BP4","AP5","BP5","AP6","BP6","AP7","BP7","AP8","BP8","AP9","BP9"/);,"dX1","dX2","dX3"/);,"dBZ","dAZ"/)
        subcycle    =     (/1D0,1D0  ,1D0  ,1D0  ,1.0D0,nsplit,nsplit,nsplit,nsplit,1.0D0,nsplit*rsplit*hypervis_subcycle,nsplit*rsplit*hypervis_subcycle,nsplit*rsplit*hypervis_subcycle,nsplit*rsplit,nsplit*rsplit,1D0  ,1D0,1D0  ,1D0,1D0  ,1D0,1D0,1D0,1D0,1D0,1D0,1D0,1D0,1D0,1D0,1D0,1D0,1D0/)
      end if	
      dycoreDiags = True
    else
      dycoreDiags = False    
      print("No dycore diagnostics seem to be on file")
      str_stage     = (/"pBF","pBP","pAP","pAM"/)
      subcycle    = (/1D0,1D0,1D0,1D0/)          
    end if
  else
    ;
    ; dycore only diags
    ;
    qsize = 0
    str_vars      = (/"SE","KE"/)
    str_stage     = (/"pBF","pBP","pAP","pAM","dED","dAF","dBD"  ,"dAD" ,"dAR" ,"dBF","dBH","dCH","dAH","dBS","dAS"/);,"dX1","dX2","dX3"/);,"dBZ","dAZ"/)
    subcycle    =     (/1D0,1D0  ,1D0  ,1D0  ,1.0D0,nsplit,nsplit,nsplit,nsplit,1.0D0,nsplit*rsplit*hypervis_subcycle,nsplit*rsplit*hypervis_subcycle,nsplit*rsplit*hypervis_subcycle,nsplit*rsplit,nsplit*rsplit/);,nsplit*rsplit,nsplit*rsplit,nsplit*rsplit/);


;    str_stage     = (/"dED","dAF","dBD","dAD","dAR","dBF","dBH","dCH","dAH"/)
;    subcycle    = (/1.0D0,nsplit,nsplit,nsplit,nsplit,1.0D0,nsplit*rsplit*hypervis_subcycle,nsplit*rsplit*hypervis_subcycle,nsplit*rsplit*hypervis_subcycle/)    
    dycoreDiags = True
  end if

  nvars       = dimsizes(str_vars)
  nstage      = dimsizes(str_stage)
  if (latlon) then
    vars        = new((/nvars,nstage,nt,nlat,nlon/),typeof(tmp(0,0,0)))
    tot_energy  = new((/nstage/),typeof(tmp(0,0,0)))
    tracer_mass = new((/nt,3,nstage/),typeof(tmp(0,0,0)))
    int_vars  = new((/nvars,nstage,nt/)     ,typeof(tmp(0,0,0)))        
  else
    if (ncol.ge.ncol_p) then
      vars          = new((/nvars,nstage,nt,ncol/),typeof(tmp(0,0)))
      tot_energy_2d = new((/3,nt,nstage,ncol/),typeof(tmp(0,0)))
    else
      vars        = new((/nvars,nstage,nt,ncol_p/),typeof(tmp(0,0)))
      tot_energy_2d = new((/3,nt,nstage,ncol_p/),typeof(tmp(0,0)))
    end if
    vars = 0.0D0  
    tot_energy_2d =  0.0D0

    tot_energy  = new((/3,nt,nstage/),typeof(tmp(0,0)))


    tracer_mass = new((/nt,3,nstage/),typeof(tmp(0,0)))
    int_vars    = new((/nvars,nstage,nt/)     ,typeof(tmp(0,0)))    
  end if

  int_vars  = 0.0D0
  ;
  ; check that variables are on file and read them and compute global integrals
  ;
  inv_area_sphere   = 1.0D0/area_sphere
  inv_area_sphere_p = 1.0D0/area_sphere_p
  print("nvars="+nvars)
  do variable=0,nvars-1
    do stage=0,nstage-1
      str = str_vars(variable)+"_"+str_stage(stage)
      print("xxx "+str)
      if (isfilevar(fall,str)) then
        if (latlon) then
          vars(variable,stage,:,:,:) = fall->$str$(:,:,:)
          do t=0,nt-1
            int_vars(variable,stage,t) = sum(area(:,:)*vars(variable,stage,t,:,:))*inv_area_sphere		  
	  end do
	else
          fld = fall->$str$(:,:)
	  if (dimsizes(fld(0,:)).eq.ncol_p) then
              vars(variable,stage,:,0:ncol_p-1) = fall->$str$(:,0:ncol_p-1)
            do t=0,nt-1	      
              int_vars(variable,stage,t) = sum(area_p(0:ncol_p-1)*fld(t,0:ncol_p-1))*inv_area_sphere_p
	    end do
          else
            do t=0,nt-1	      	  
              int_vars(variable,stage,t) = sum(area(0:ncol-1)*fld(t,0:ncol-1))*inv_area_sphere
              vars(variable,stage,:,0:ncol-1) = fall->$str$(:,:)	  
	    end do
	  end if
	  delete(fld)
   	end if
	
      else
        print("Variable "+str+" not found")
        exit
      end if
    end do
  end do
  print("Done integrals")
  ;
  ; compute global integrals
  ;

;  do t=0,nt-1
;    do variable = 0,nvars-1
;      do stages=0,nstage-1
;        if (latlon) then
;          int_vars(variable,stages,t) = sum(area(:,:)*vars(variable,stages,t,:,:))*inv_area_sphere	
;	else
;          int_vars(variable,stages,t) = sum(area(0:ncol-1)*vars(variable,stages,t,0:ncol-1))*inv_area_sphere	
;	end if
;      end do
;    end do
;  end do

  ;
  ; For code clarity define indices corresponding to the different stages
  ;
  if (isfilevar(fall,"SE_pBF")) then
    ipBF = str_match_ind(str_stage,"pBF")
  else
    iPBF = -1
  end if
  if (isfilevar(fall,"SE_pBP")) then
    ipBP = str_match_ind(str_stage,"pBP")
  else
    iPBP = -1
  end if
  if (isfilevar(fall,"SE_pAP")) then
    ipAP = str_match_ind(str_stage,"pAP")
  else
    iPAP = -1
  end if

  if (isfilevar(fall,"SE_AP1")) then
    iAP1 = str_match_ind(str_stage,"AP1")
  else
    iAP1 = -1
  end if
  if (isfilevar(fall,"SE_AP2")) then
    iAP2 = str_match_ind(str_stage,"AP2")
  else
    iAP2 = -1
  end if
  if (isfilevar(fall,"SE_AP3")) then
    iAP3 = str_match_ind(str_stage,"AP3")
  else
    iAP3 = -1
  end if
  if (isfilevar(fall,"SE_AP4")) then
    iAP4 = str_match_ind(str_stage,"AP4")
  else
    iAP4 = -1
  end if
  if (isfilevar(fall,"SE_AP5")) then
    iAP5 = str_match_ind(str_stage,"AP5")
  else
    iAP5 = -1
  end if
  if (isfilevar(fall,"SE_AP6")) then
    iAP6 = str_match_ind(str_stage,"AP6")
  else
    iAP6 = -1
  end if
  if (isfilevar(fall,"SE_AP7")) then
    iAP7 = str_match_ind(str_stage,"AP7")
  else
    iAP7 = -1
  end if
  if (isfilevar(fall,"SE_AP8")) then
    iAP8 = str_match_ind(str_stage,"AP8")
  else
    iAP8 = -1
  end if
  if (isfilevar(fall,"SE_AP9")) then
    iAP9 = str_match_ind(str_stage,"AP9")
  else
    iAP9 = -1
  end if




  if (isfilevar(fall,"SE_BP1")) then
    iBP1 = str_match_ind(str_stage,"BP1")
  else
    iBP1 = -1
  end if
  if (isfilevar(fall,"SE_BP2")) then
    iBP2 = str_match_ind(str_stage,"BP2")
  else
    iBP2 = -1
  end if
  if (isfilevar(fall,"SE_BP3")) then
    iBP3 = str_match_ind(str_stage,"BP3")
  else
    iBP3 = -1
  end if
  if (isfilevar(fall,"SE_BP4")) then
    iBP4 = str_match_ind(str_stage,"BP4")
  else
    iBP4 = -1
  end if
  if (isfilevar(fall,"SE_BP5")) then
    iBP5 = str_match_ind(str_stage,"BP5")
  else
    iBP5 = -1
  end if
  if (isfilevar(fall,"SE_BP6")) then
    iBP6 = str_match_ind(str_stage,"BP6")
  else
    iBP6 = -1
  end if
  if (isfilevar(fall,"SE_BP7")) then
    iBP7 = str_match_ind(str_stage,"BP7")
  else
    iBP7 = -1
  end if
  if (isfilevar(fall,"SE_BP8")) then
    iBP8 = str_match_ind(str_stage,"BP8")
  else
    iBP8 = -1
  end if
  if (isfilevar(fall,"SE_BP9")) then
    iBP9 = str_match_ind(str_stage,"BP9")
  else
    iBP9 = -1
  end if


  if (isfilevar(fall,"SE_pAM")) then
    ipAM = str_match_ind(str_stage,"pAM")
  else
    iPAM = -1
  end if
  if (isfilevar(fall,"SE_dED")) then
    idED = str_match_ind(str_stage,"dED")
  else
    idED = -1
  end if
  if (isfilevar(fall,"SE_dBD")) then
    idBD = str_match_ind(str_stage,"dBD")
  else
    idBD = -1
  end if
  if (isfilevar(fall,"SE_dAF")) then
    idAF = str_match_ind(str_stage,"dAF")
  else
    idAF = -1
  end if
  if (isfilevar(fall,"SE_dAD")) then
    idAD = str_match_ind(str_stage,"dAD")
  else
    idAD = -1
  end if
  if (isfilevar(fall,"SE_dAR")) then
    idAR = str_match_ind(str_stage,"dAR")
  else
    idAR = -1
  end if
  if (isfilevar(fall,"SE_dBF")) then
    idBF = str_match_ind(str_stage,"dBF")
  else
    idBF = -1
  end if
  if (isfilevar(fall,"SE_dBH")) then
    idBH = str_match_ind(str_stage,"dBH")
  else
    idBH = -1
  end if
  if (isfilevar(fall,"SE_dCH")) then
    idCH = str_match_ind(str_stage,"dCH")
  else
    idCH = -1
  end if                  
  if (isfilevar(fall,"SE_dAH")) then
    idAH = str_match_ind(str_stage,"dAH")
  else
    idAH = -1
  end if
  if (isfilevar(fall,"SE_dAP")) then
    idAP = str_match_ind(str_stage,"dAP")
  else
    idAP = -1
  end if

  if (isfilevar(fall,"SE_dAT")) then
    idAT = str_match_ind(str_stage,"dAT")
  else
    idAT = -1
  end if                      


  if (isfilevar(fall,"SE_dBS")) then
    idBS = str_match_ind(str_stage,"dBS")
  else
    idBS = -1
  end if                      
  if (isfilevar(fall,"SE_dAS")) then
    idAS = str_match_ind(str_stage,"dAS")
  else
    idAS = -1
  end if                      
  if (isfilevar(fall,"SE_dBZ")) then
    idBZ = str_match_ind(str_stage,"dBZ")
  else
    idBZ = -1
  end if                      
  if (isfilevar(fall,"SE_dAZ")) then
    idAZ = str_match_ind(str_stage,"dAZ")
  else
    idAZ = -1
  end if                      
  if (isfilevar(fall,"SE_dX1")) then
    idX1 = str_match_ind(str_stage,"dX1")
  else
    idX1 = -1
  end if                      
  if (isfilevar(fall,"SE_dX2")) then
    idX2 = str_match_ind(str_stage,"dX2")
  else
    idX2 = -1
  end if                      
  if (isfilevar(fall,"SE_dX3")) then
    idX3 = str_match_ind(str_stage,"dX3")
  else
    idX3 = -1
  end if                      


  if (ismissing(ipBF).or.ismissing(ipBF).or.ismissing(ipAM).or.ismissing(ipAP)) then
    PhysDiags = False
    print("Physics diagnostics are not all on file")
  else
    PhysDiags = True
    print("Physics diagnostics are on file")
  end if

  do t=0,nt-1
    do j=0,nstage-1
      tot_energy(0,t,j)     = subcycle(j)*(int_vars(0,j,t)+int_vars(1,j,t))
      ;
      ;
      tot_energy(1,t,j)  = subcycle(j)*int_vars(0,j,t) ; thermal energy only
      tot_energy(2,t,j)  = subcycle(j)*int_vars(1,j,t) ; K energy only	       

      tot_energy_2d(0,t,j,:) = subcycle(j)*(vars(0,j,t,:)+vars(1,j,t,:))
      tot_energy_2d(1,t,j,:) = subcycle(j)*(vars(0,j,t,:))
      tot_energy_2d(2,t,j,:) = subcycle(j)*(vars(1,j,t,:))
      if (wet) then
        idx = str_match_ind(str_vars,"WV")
        tracer_mass(t,0,j) = subcycle(j)*(int_vars(idx,j,t))
        idx = str_match_ind(str_vars,"WL")
        tracer_mass(t,1,j) = subcycle(j)*(int_vars(idx,j,t))
		print("iii "+j+" "+tracer_mass(t,2,j)+" "+int_vars(idx,j,t))		
        idx = str_match_ind(str_vars,"WI")
        tracer_mass(t,2,j) = subcycle(j)*(int_vars(idx,j,t))
	print("qqq "+str_stage(j)+" "+tracer_mass(t,2,j)+" "+int_vars(idx,j,t))
      end if
    end do
  end do
    
  print("  ")
  print("  ")
  print("  ")
  print("ANALYSIS_START_ANALYSIS_START_ANALYSIS_START_ANALYSIS_START_ANALYSIS_START_ANALYSIS_START")
  print("  ")
  print("  ")
  print("  ")

  te_consistent_version=False
  if (qsize_condensate_loading.eq."1") then
    if (lcp_moist.eq."false") then
      print("Using total energy consistent version: qsize_condensate_loading=1 and cp=cpdry")
      te_consistent_version=True
    else
      print("WARNING: Total energy formulaes for dynamics and physics are different:")
      print("   Dynamics (cp includes water vapor; condensates not thermodynamically active).")
      print("   Physics (cp=cp_dry in internal energy).")
    end if
  else
    print("WARNING: Total energy formulaes for dynamics and physics are different")
    print("in dynamics (cp and dp includes all water variables) and physics (cp=cp_dry in internal energy).")
  end if

  energy_names = new(3, string)
  energy_names = (/"TOTAL ENERGY","THERMAL ENERGY","KINETIC ENERGY"/)

  do t=0,nt-1
    if (PhysDiags) then
      pEFIX  = new((/3/),typeof(tot_energy(0,0,0)))
      pDMEA  = new((/3/),typeof(tot_energy(0,0,0)))
      param  = new((/3/),typeof(tot_energy(0,0,0)))
      phys_total  = new((/3/),typeof(tot_energy(0,0,0)))
      pDEFAULT  = new((/3,10/),typeof(tot_energy(0,0,0)))
      do idx = 0,2   
        pEFIX(idx) = (tot_energy(idx,t,ipBP)-tot_energy(idx,t,ipBF))/dtime
        pDMEA(idx) = (tot_energy(idx,t,ipAM)-tot_energy(idx,t,ipAP))/dtime
        param(idx) = (tot_energy(idx,t,ipAP)-tot_energy(idx,t,ipBP))/dtime
        phys_total(idx) = (tot_energy(idx,t,ipAM)-tot_energy(idx,t,ipBF))/dtime
      end do
      diff_DME = pDMEA(0)      

      pDEFAULT(1,1) = (tot_energy(1,t,iAP1)-tot_energy(1,t,iBP1))/dtime
      pDEFAULT(2,1) = (tot_energy(2,t,iAP1)-tot_energy(2,t,iBP1))/dtime
      

      pDEFAULT(1,2) = (tot_energy(1,t,iAP2)-tot_energy(1,t,iBP2))/dtime
      pDEFAULT(2,2) = (tot_energy(2,t,iAP2)-tot_energy(2,t,iBP2))/dtime

      CAM_LHS_ave = pDEFAULT(1,2)
      CAM_RHS_ave = pDEFAULT(2,2)

      cp_dflt = (tot_energy_2d(1,t,iAP2,:)-tot_energy_2d(1,t,iBP2,:)-(tot_energy_2d(2,t,iAP2,:)-tot_energy_2d(2,t,iBP2,:)))/dtime

      CAM_LHS = (tot_energy_2d(1,t,iAP2,:)-tot_energy_2d(1,t,iBP2,:))/dtime
      CAM_RHS = (tot_energy_2d(2,t,iAP2,:)-tot_energy_2d(2,t,iBP2,:))/dtime	

      pDEFAULT(1,3) = (tot_energy(1,t,iAP3)-tot_energy(1,t,iBP3))/dtime
      pDEFAULT(2,3) = (tot_energy(2,t,iAP3)-tot_energy(2,t,iBP3))/dtime

      cpwv_2d = (tot_energy_2d(1,t,iAP3,:)-tot_energy_2d(1,t,iBP3,:)-(tot_energy_2d(2,t,iAP3,:)-tot_energy_2d(2,t,iBP3,:)))/dtime

      pDEFAULT(1,4) = (tot_energy(1,t,iAP4)-tot_energy(1,t,iBP4))/dtime
      pDEFAULT(2,4) = (tot_energy(2,t,iAP4)-tot_energy(2,t,iBP4))/dtime

      all_water = (tot_energy_2d(1,t,iAP4,:)-tot_energy_2d(1,t,iBP4,:)-(tot_energy_2d(2,t,iAP4,:)-tot_energy_2d(2,t,iBP4,:)))/dtime

      pDEFAULT(1,5) = (tot_energy(1,t,iAP5)-tot_energy(1,t,iBP5))/dtime
      pDEFAULT(2,5) = (tot_energy(2,t,iAP5)-tot_energy(2,t,iBP5))/dtime

      K_2d = (tot_energy_2d(1,t,iAP5,:)-tot_energy_2d(1,t,iBP5,:)-(tot_energy_2d(2,t,iAP5,:)-tot_energy_2d(2,t,iBP5,:)))/dtime

      pDEFAULT(1,6) = (tot_energy(1,t,iAP6)-tot_energy(1,t,iBP6))/dtime
      pDEFAULT(2,6) = (tot_energy(2,t,iAP6)-tot_energy(2,t,iBP6))/dtime

      cpT_2d = (tot_energy_2d(1,t,iAP6,:)-tot_energy_2d(1,t,iBP6,:)-(tot_energy_2d(2,t,iAP6,:)-tot_energy_2d(2,t,iBP6,:)))/dtime

      pDEFAULT(1,7) = (tot_energy(1,t,iAP7)-tot_energy(1,t,iBP7))/dtime
      pDEFAULT(2,7) = (tot_energy(2,t,iAP7)-tot_energy(2,t,iBP7))/dtime

      PHIS_2d = (tot_energy_2d(1,t,iAP7,:)-tot_energy_2d(1,t,iBP7,:))/dtime - (tot_energy_2d(2,t,iAP7,:)-tot_energy_2d(2,t,iBP7,:))/dtime

      pDEFAULT(1,8) = (tot_energy(1,t,iAP8)-tot_energy(1,t,iBP8))/dtime
      pDEFAULT(2,8) = (tot_energy(2,t,iAP8)-tot_energy(2,t,iBP8))/dtime

      lt_2d = (tot_energy_2d(1,t,iAP8,:)-tot_energy_2d(1,t,iBP8,:))/dtime - (tot_energy_2d(2,t,iAP8,:)-tot_energy_2d(2,t,iBP8,:))/dtime

      pDEFAULT(1,9) = (tot_energy(1,t,iAP9)-tot_energy(1,t,iBP9))/dtime
      pDEFAULT(2,9) = (tot_energy(2,t,iAP9)-tot_energy(2,t,iBP9))/dtime

      lvariableL_2d = (tot_energy_2d(1,t,iAP9,:)-tot_energy_2d(1,t,iBP9,:))/dtime - ((tot_energy_2d(2,t,iAP9,:)-tot_energy_2d(2,t,iBP9,:))/dtime )


      pwork = (tot_energy_2d(0,t,ipAM,:)-tot_energy_2d(0,t,ipAP,:))/dtime

      dparam = (tot_energy_2d(0,t,ipAP,:)-tot_energy_2d(0,t,ipBP,:))/dtime  

    end if

    if (dycoreDiags) then
      dyn_total       = (tot_energy(0,t,idBF)-tot_energy(0,t,idED))/dtime
      dyn_phys_total       = (tot_energy(0,t,idAD)-tot_energy(0,t,idAF))/dtime
    end if
   

    str_data  = new( 1, "string")
    if (PhysDiags) then
      print(" ")
      print("=================================================================================")
      print("|                                                                               |")
      print("| NEW ANALYSIS                                                                  |")
      print("|                                                                               |")
      print("=================================================================================")
      print("1 Default CAM")
      print("dE/dt parameterizations (AP1-BP1)"+sprintf("%4g",pDEFAULT(1,1))+" W/M^2")
      print("Flux terms              (AP1-BP1)"+sprintf("%4g",pDEFAULT(2,1))+" W/M^2")
      print(" ")
      print("2 Default CAM using equation (57)")
      print(" ")
      print("dE/dt parameterizations (AP2-BP2)"+sprintf("%4g",pDEFAULT(1,2))+" W/M^2")
      print("Flux terms              (AP2-BP2)"+sprintf("%4g",pDEFAULT(2,2))+" W/M^2")
      diff = pDEFAULT(1,2)-pDEFAULT(2,2)
      print("Diff                             "+sprintf("%4g",diff         )+" W/M^2")
      print(" ")
      print("3 Default CAM but using cpwv and not cpair")
      print(" ")
      print("dE/dt parameterizations (AP3-BP3)"+sprintf("%4g",pDEFAULT(1,3))+" W/M^2")
      print("Flux terms              (AP3-BP3)"+sprintf("%4g",pDEFAULT(2,3))+" W/M^2")
      diff_cpwv = pDEFAULT(1,3)-pDEFAULT(2,3)
      print("Diff                             "+sprintf("%4g",diff_cpwv    )+" W/M^2")
      print(" ")
      print("4 Default CAM but incl. all forms of water in pdel")
      print(" ")
      print("dE/dt parameterizations (AP4-BP4)"+sprintf("%4g",pDEFAULT(1,4))+" W/M^2")
      print("Flux terms              (AP4-BP4)"+sprintf("%4g",pDEFAULT(2,4))+" W/M^2")
      diff_all_water_pel = pDEFAULT(1,4)-pDEFAULT(2,4)
      print("Diff                             "+sprintf("%4g",diff_all_water_pel        )+" W/M^2")
      print(" ")
      print("5 Default CAM using equation (57) BUT ALL FORMS OF WATER IN PDEL AND UPDATE WATER TRACERS IN K AND PHIS TERM ")
      print(" ")
      print("dE/dt parameterizations (AP5-BP5)"+sprintf("%4g",pDEFAULT(1,5))+" W/M^2")
      print("Flux terms              (AP5-BP5)"+sprintf("%4g",pDEFAULT(2,5))+" W/M^2")
      diff_KPHIS = pDEFAULT(1,5)-pDEFAULT(2,5)
      print("Diff                             "+sprintf("%4g",diff_KPHIS         )+" W/M^2")

      print(" ")
      print("6 Default CAM using equation (57) BUT ALL FORMS OF WATER IN PDEL AND UPDATE WATER TRACERS IN CP*T TERM")
      print(" ")
      print("dE/dt parameterizations (AP6-BP6)"+sprintf("%4g",pDEFAULT(1,6))+" W/M^2")
      print("Flux terms              (AP6-BP6)"+sprintf("%4g",pDEFAULT(2,6))+" W/M^2")
      diff_cpT = pDEFAULT(1,6)-pDEFAULT(2,6)
      print("Diff                             "+sprintf("%4g",diff_cpT         )+" W/M^2")

      print(" ")
      print("7 Default CAM using equation (57)")
      print(" ")
      print("dE/dt parameterizations (AP7-BP7)"+sprintf("%4g",pDEFAULT(1,7))+" W/M^2")
      print("Flux terms              (AP7-BP7)"+sprintf("%4g",pDEFAULT(2,7))+" W/M^2")
      diff = pDEFAULT(1,7)-pDEFAULT(2,7)
      print("Diff                             "+sprintf("%4g",diff         )+" W/M^2")

      print(" ")
      print("8 Default CAM using equation (57)")
      print(" ")
      print("dE/dt parameterizations (AP8-BP8)"+sprintf("%4g",pDEFAULT(1,8))+" W/M^2")
      print("Flux terms              (AP8-BP8)"+sprintf("%4g",pDEFAULT(2,8))+" W/M^2")
      diff = pDEFAULT(1,8)-pDEFAULT(2,8)
      print("Diff                             "+sprintf("%4g",diff         )+" W/M^2")

      print(" ")
      print("9 Default CAM using equation (57)")
      print(" ")
      print("dE/dt parameterizations (AP9-BP9)"+sprintf("%4g",pDEFAULT(1,9))+" W/M^2")
      print("Flux terms              (AP9-BP9)"+sprintf("%4g",pDEFAULT(2,9))+" W/M^2")
      diff = pDEFAULT(1,9)-pDEFAULT(2,9)
      print("Diff                             "+sprintf("%4g",diff         )+" W/M^2")


;      do idx=0,2
;        print("dE/dt parameterizations + source sink terms  (AP1-BP1) "+sprintf("%4g",pDEFAULT(idx))+" W/M^2")
;      end do
      print(" ")
      print("=================================================================================")
      print("|                                                                               |")
      print("| ANALYSIS OF ENERGY DIAGNOSTICS IN PHYSCS                                      |")
      print("|                                                                               |")
      print("=================================================================================")
      do idx=0,2
        print(" ")
        print(energy_names(idx))
        print("                                                  "+"total E  " +"  thermal E  "+"  K E")
        print("dE/dt parameterizations (param)         (pAP-pBP) "+sprintf("%4g",param(idx))+" W/M^2")
        print("dE/dt dry mass adjustment (pwork)       (pAM-pAP) "+sprintf("%4g",pDMEA(idx))+" W/M^2")
        print("dE/dt energy fixer (efix)               (pBP-pBF) "+sprintf("%4g",pEFIX(idx))+" W/M^2")
        print("dE/dt physics total (phys)              (pAM-pBF) "+sprintf("%4g",phys_total(idx)) +" W/M^2")
      end do
      
    end if
    
    if (dycore.eq."FV".or.dycore.eq."FV3") then
      print(" ")    
      fv_dycore = -phys_total(0)-pDMEA(0)
      print("Dycore TE dissipation estimated from physics TE terms (assuming ftype=1 and no energy discretization errors is")
      print(sprintf("%6.3g",fv_dycore)     +" W/M^2")
      print(" ")      
    end if

    if (dycore.eq."FV3") then
      print(" ")
      print("=================================================================================")
      print("|                                                                               |")
      print("| ANALYSIS OF ENERGY DIAGNOSTICS IN DYNAMICS - specific for FV3 dycore          |")
      print("|                                                                               |")
      print("=================================================================================")
      print(" ")
      diff1 = (tot_energy(0,t,ipAM)-tot_energy(0,t,ipBF))/dtime
      diff2 = (tot_energy(0,t,idAP)-tot_energy(0,t,idED))/dtime
      
      diff  = diff2-diff1
      
      param_dyn= (tot_energy(0,t,idAP)-tot_energy(0,t,idED))/dtime	
      print("dE/dt parameterizations (param) in dyn on A-grid (dAP-dED) "+sprintf("%6.3g",param_dyn)     +" W/M^2")
      param_dyn_D= (tot_energy(0,t,idBD)-tot_energy(0,t,idED))/dtime	
      print("dE/dt parameterizations (param) in dyn on D-grid (dBD-dED) "+sprintf("%6.3g",param_dyn_D)     +" W/M^2")
      param_dyn_DA = (tot_energy(0,t,idBD)-tot_energy(0,t,idAP))/dtime
      print("dE/dt mapping wind tendencies from A to D grid (and back to A-grid for TE computations (dBD-dAP) "+sprintf("%6.3g",param_dyn_DA)     +" W/M^2")


      print(" Physics-dynamics coupling error (should be zero) : "+diff)
      print(" (dE/dt physics increment minus dE/dt physics increment in dynamics)")
      print(" ")


      dDycore = (tot_energy(0,t,idBF)-tot_energy(0,t,idBD))/dtime      
      dadiab  = (tot_energy_2d(0,t,idBF,:)-tot_energy_2d(0,t,idBD,:))/dtime      

      print(" dE/dt dynamical core total (dBF-dBD) = "+dDycore)
      print(" ")
      print(" Breakdown of dycore dE/dt:")
      print(" ")      

      Rayleigh = (tot_energy(0,t,idAT)-tot_energy(0,t,idBD))/dtime      
      print(" dE/dt Rayleigh friction (dAT-dBD) = "+Rayleigh)


      dDycore2d = (tot_energy(0,t,idAD)-tot_energy(0,t,idAF))/dtime
      print(" dE/dt 2D dynamical core (dAD-dAF) = "+dDycore2d)      

      dRemap = (tot_energy(0,t,idAR)-tot_energy(0,t,idAD))/dtime
      remap  = (tot_energy_2d(0,t,idAR,:)-tot_energy_2d(0,t,idAD,:))/dtime
      print(" dE/dt vertical remapping(dAR-dAD) = "+dRemap)      
      print(" ")
      diff = dDycore-(dDycore2d+dRemap+Rayleigh)
      print("Check that dE/dt 2D dynamics + vertical remapping + Rayleigh friction - dE/dt dynamical core total (should be zero)="+diff)

      hvis = remap
      hvis = 0.0

      fheat = remap
      fheat = 0.0
      
    end if

    if (dycore.eq."SE".and.dycoreDiags) then
      print(" ")
      print("=================================================================================")
      print("|                                                                               |")
      print("| ANALYSIS OF ENERGY DIAGNOSTICS IN DYNAMICS - specific for SE dycore           |")
      print("|                                                                               |")
      print("=================================================================================")
      print(" ")


      rate_of_change_2D_dyn = (tot_energy(0,t,idAD)-tot_energy(0,t,idBD))/dtime

      rate_of_change_vertical_remapping = (tot_energy(0,t,idAR)-tot_energy(0,t,idAD))/dtime

      remap = (tot_energy_2d(0,t,idAR,:)-tot_energy_2d(0,t,idAD,:))/dtime

      dADIA = rate_of_change_2D_dyn+rate_of_change_vertical_remapping

      dadiab = (tot_energy_2d(0,t,idAR,:)-tot_energy_2d(0,t,idBD,:))/dtime


;      x1       = (tot_energy(0,t,idX2)-tot_energy(0,t,idX1))/dtime
;      x2       = (tot_energy(0,t,idX3)-tot_energy(0,t,idX2))/dtime
;      print("dE/dt x1 "+sprintf("%6.3g",x1)+" W/M^2")
;      print("dE/dt x2 "+sprintf("%6.3g",x2)+" W/M^2")


      diffusion_del4       = (tot_energy(0,t,idCH)-tot_energy(0,t,idBH))/dtime
      diffusion_fric       = (tot_energy(0,t,idAH)-tot_energy(0,t,idCH))/dtime
      diffusion_del4_tot   = (tot_energy(0,t,idAH)-tot_energy(0,t,idBH))/dtime
      diffusion_sponge     = (tot_energy(0,t,idAS)-tot_energy(0,t,idBS))/dtime
      diffusion_total      = diffusion_del4_tot+diffusion_sponge

      print("dE/dt dyn total (dycore+phys tendency   (dBF-dED) "+sprintf("%6.3g",dyn_total)+" W/M^2")

      print("dE/dt total adiabatic dynamics (adiab)            "+sprintf("%6.3g",dADIA)+" W/M^2")    
      print("dE/dt 2D dynamics (2D)                  (dAD-dBD) "+sprintf("%6.3g",rate_of_change_2D_dyn)+" W/M^2")
      print("dE/dt vertical remapping (remap)        (dAR-dAD) "+sprintf("%6.3g",rate_of_change_vertical_remapping)+" W/M^2")
      rate_of_change_physics = (tot_energy(0,t,idBD)-tot_energy(0,t,idAF))/dtime
      print("dE/dt physics tendency in dynamics      (dBD-dAF) "+sprintf("%6.3g",rate_of_change_physics)+" W/M^2")
      print(" ")
      print("Breakdown of 2D dynamics:")
      print(" ")
      rate_of_change_heating_term_put_back_in = (tot_energy(0,t,idAH)-tot_energy(0,t,idCH))/dtime
      fheat = (tot_energy_2d(0,t,idAH,:)-tot_energy_2d(0,t,idCH,:))/dtime

      hvis = (tot_energy_2d(0,t,idAH,:)-tot_energy_2d(0,t,idBH,:))/dtime
      rate_of_change_hvis_sponge = (tot_energy(0,t,idAS)-tot_energy(0,t,idBS))/dtime

      print("      dE/dt hypervis del4 (hvis)             (dCH-dBH) "+sprintf("%6.3g",diffusion_del4)+" W/M^2")
      print("      dE/dt hypervis frictional heating del4 (dAH-dCH) "+sprintf("%6.3g",diffusion_fric)+" W/M^2")
      print("      dE/dt hypervis del4 total (hvis)       (dAH-dBH) "+sprintf("%6.3g",diffusion_del4_tot)+" W/M^2")
      print("      dE/dt hypervis sponge total            (dAS-dBS) "+sprintf("%6.3g",diffusion_sponge)+" W/M^2")
      print("                                                       ")
      print("      dE/dt explicit diffusion total                   "+sprintf("%6.3g",diffusion_total)+" W/M^2")


;      consistency = diffusion_del4_tot-(diffusion_del4+diffusion_fric)
;      print("      dE/dt hypervis consistency (dAH-dBH)=(dCH-dBH)+(dAH-dCH)="+sprintf("%6.3g",consistency)+" W/M^2")



      twoDresidual = rate_of_change_2D_dyn-diffusion_total
      print("      dE/dt residual (res)                        "+sprintf("%6.3g",twoDresidual)+" W/M^2")
    end if
    print(" ")
    print(" ")

    if (dycore.eq."SE".and.PhysDiags.and.dycoreDiags) then
      print(" ")
      print("=================================================================================")
      print("|                                                                               |")
      print("| ANALYSIS OF ENERGY DIAGNOSTICS IN DYNAMICS-PHYSICS COMBINED                   |")
      print("|                                                                               |")
      print("=================================================================================")
      print(" ")
      value_pdc = phys_total(0)-rate_of_change_physics
      if (te_consistent_version.eq."True") then
        print("Your model is energy consistent (qsize_condensate_loading=1 and cpdry)")
        if (ftype .eq. "1") then
          print(" ")
          print("You are using ftype=1 so PDC errors should be zero:")
          print(" ")
          print("    dE/dt physics tendency in dynamics (dBD-dAF) should exactly match dE/dt physics total (pAM-pBF): "+sprintf("%6.3g",value_pdc))
          print(" ")
        else
          print(" ")
          print("You are using ftype=0 or 2 so there are PDC errors (dribbling errors):")
          print(" ")
          print("   Dribbling errors (pAM-pBF-(dBD-dAF))/dt: "+sprintf("%6.3g",value_pdc))
        end if
        discr = "0       "
        str_pdc = sprintf("%6.3g",10*value_pdc)
      else
        print("Your model is energy inconsistent (qsize_condensate_loading<>1 and/or cp<>cpdry)")
        print(" ")
        print("PDC errors can not be assesed trhough ")
        print(" ")
        print("   dE/dt physics tendency in dynamics (dBD-dAF) does not match dE/dt physics total (pAM-pBF) due to energy discrepancy:"+sprintf("%6.3g",value_pdc))
	print(phys_total+"  "+rate_of_change_physics)
        str_pdc = "undef"
      end if
      print(" ")
      print("Some more consisitency/budget terms")
      print("===================================")
      print(" ")
      print("Energy fixer fixes dme_adjust (pDMEA), lack of energy conservation in adiabatic")
      print("dynamical core (dADIA), energy discrepancy (EDIFF) and energy lost/gained in physics-dynamics coupling")
      print(" ")
      print("dADIA                     "+sprintf("%6.3g",dADIA    )+" W/M^2")    
      print("pDMEA                     "+sprintf("%6.3g",pDMEA(0) )+" W/M^2")    
      print("physics-dynamics coupling "+sprintf("%6.3g",value_pdc)+" W/M^2")
      print(" ")
      str="dPDC+EDIFF"
      print(" ")
      print("		-energy fixer = DME_adjust+adaib dycore+phys-dyn errors+discr")
      print("            ")
      ttt = -pDMEA(0)-dADIA-value_pdc
      discr = -99.0
      print("          DME_adjust+adaib dycore+phys-dyn errors+discr = "+sprintf("%6.3g",ttt))
      print("          Energy fixer                                  = "+sprintf("%6.3g",pEFIX(0)))
      print(" ")
      fff = ttt-pEFIX(0)
      print("          Difference                                    = "+sprintf("%6.3g",fff))     


      qneg  = (tot_energy(0,t,ipBF)-tot_energy(0,t,idED))/dtime
      print(" ")
      print(" qneg: "+qneg)
      print(" ")

      ;
      ; write diagnostics to file
      ;
      str1 = qsize_condensate_loading+"                       &"
      str1 = str1 + " "+lcp_moist+"   &"      
      str1 = str1 + " "+ftype+"   &"      
      str1 = str1 + " "+sprintf("%6.3g",param(0))+"&"
      str1 = str1 + " "+sprintf("%6.3g",pDMEA(0))+"&"
      str1 = str1 + " "+sprintf("%6.3g",phys_total(0))+"&"
      str1 = str1 + " "+discr+"&"
      str1 = str1 + " "+sprintf("%6.3g",rate_of_change_2D_dyn)+"&"
      str1 = str1 + " "+sprintf("%6.3g",diffusion_total)+"&"
      str1 = str1 + " "+sprintf("%6.3g",rate_of_change_heating_term_put_back_in)+"&"
 ;      str1 = str1 + " "+sprintf("%6.3g",twoDresidual)+"&"
      str1 = str1 + " "+sprintf("%6.3g",twoDresidual)+"&"
      str1 = str1 + " "+sprintf("%6.3g",rate_of_change_vertical_remapping)+"&"
      str1 = str1 + " "+sprintf("%6.3g",dADIA)+"&"
      str1 = str1 + " "+sprintf("%6.3g",value_pdc)
      str1 = str1 + " "+str_pdc

      if (fileexists(diags_fname)) then
        c1 = (/str1/)
        alist   = [/c1/]
        write_table(diags_fname, "a", alist, "%s")
      else
        c1 = (/"#qsize_condensate_loading cpmoist  ftype param     pwork     efix      discr     2D         hvis       fheat     res       remap      adiab      PDC",str1/)
        alist   = [/c1/]
        write_table(diags_fname, "w", alist, "%s")
      end if
    end if
    if (qsize.gt.0) then
      print("    ")
      print("    ")
      print(" ")
      print("=================================================================================")
      print("|                                                                               |")
      print("| ANALYSIS OF WATER VAPOR, CLOUD LIQUID AND CLOUD ICE BUDGETS                   |")
      print("|                                                                               |")
      print("=================================================================================")
      print(" ")
      print(" ")
    end if

    do m_cnst=0,qsize-1
      if (m_cnst.eq.0) then 
        print("Water vapor")
        print("-----------")
      end if
      if (m_cnst.eq.1) then
        print("Cloud liquid")
        print("------------")
      end if
      if (m_cnst.eq.2) then 
        print("Cloud ice")
        print("---------")
      end if
      print(" ")
      phys_total(0)    = (tracer_mass(t,m_cnst,ipBP)-tracer_mass(t,m_cnst,ipBF))
      pDMEA(0) = (tracer_mass(t,m_cnst,ipAM)-tracer_mass(t,m_cnst,ipAP))
      param(0)           = (tracer_mass(t,m_cnst,ipAP)-tracer_mass(t,m_cnst,ipBP))
      dyn_total       = (tracer_mass(t,m_cnst,idBF)-tracer_mass(t,m_cnst,idED))
      phys_total      = (tracer_mass(t,m_cnst,ipAM)-tracer_mass(t,m_cnst,ipBF))

      print("dMASS/dt energy fixer                      (pBP-pBF) "+sprintf("%6.3g",phys_total(0))+" Pa")
      print("dMASS/dt parameterizations                 (pAP-pBP) "+sprintf("%6.3g",param(0))+" Pa")
      print("dMASS/dt dry mass adjustment               (pAM-pAP) "+sprintf("%6.3g",pDMEA(0))+" Pa")
      print(" ")
      value = phys_total(0)+pDMEA(0)
      print("=> dMASS/dt dynamical core (estimated from physics)  ")
      print("   dMASS/dt energy fixer + dMASS/dt dry mass adjustment "+sprintf("%6.3g",value)+" Pa")

      print("=> dMASS/dt physics total                   (pAM-pBF)"+sprintf("%6.3g",phys_total(0))+" Pa")
 

      print("  ")
      print("  ")
      print("  ")

      if (idAD.ge.0.and.idBD.ge.0.and.idAR.ge.0.and.idCH.ge.0) then
        mass_change__2D_dyn = (tracer_mass(t,m_cnst,idAD)-tracer_mass(t,m_cnst,idBD))
        mass_change__vertical_remapping = (tracer_mass(t,m_cnst,idAR)-tracer_mass(t,m_cnst,idAD))
        dADIA = mass_change__2D_dyn+mass_change__vertical_remapping
        print("dE/dt total adiabatic dynamics                    "+sprintf("%6.3g",dADIA)+" Pa")
        print("dE/dt 2D dynamics                       (dAD-dBD) "+sprintf("%6.3g",mass_change__2D_dyn)+" Pa")
        print(" ")
        print("Breakdown of 2D dynamics:")
        print(" ")
        mass_change__heating_term_put_back_in = (tracer_mass(t,m_cnst,idAH)-tracer_mass(t,m_cnst,idCH))
        mass_change__hypervis_total           = (tracer_mass(t,m_cnst,idAH)-tracer_mass(t,m_cnst,idBH))
        print("      dE/dt hypervis                    (dAH-dBH) "+sprintf("%6.3g",mass_change__hypervis_total)+" Pa")
        print("      dE/dt frictional heating          (dAH-dCH) "+sprintf("%6.3g",mass_change__heating_term_put_back_in)+" Pa")
        error = mass_change__2D_dyn-mass_change__hypervis_total
        print("      dE/dt residual (time truncation errors)     "+sprintf("%6.3g",error)+" Pa")	
      end if
      if (idAR.ge.0.and.idAD.ge.0) then
        mass_change__vertical_remapping = (tracer_mass(t,m_cnst,idAR)-tracer_mass(t,m_cnst,idAD))
        print("dE/dt vertical remapping                (dAR-dAD) "+sprintf("%6.3g",mass_change__vertical_remapping)+" Pa")
      end if
      print(" ")
      print(" ")

      if (idBD.ge.0.and.idAF.ge.0) then
        mass_change__physics = (tracer_mass(t,m_cnst,idBD)-tracer_mass(t,m_cnst,idAF))
        print("dE/dt physics tendency in dynamics      (dBD-dAF) "+sprintf("%6.3g",mass_change__physics)+" Pa")
        value = phys_total(0)-mass_change__physics
      end if
      if (idBD.ge.0.and.idAF.ge.0) then
        if (ftype .eq. "1") then
          print(" ")
          print("      Consistency check:")
          print(" ")
          print("      dE/dt physics tendency in dynamics (dBD-dAF) should exactly match dE/dt physics total (pAM-pBF):"+sprintf("%6.3g",value))
          print(" ")
        else
          print("Dribbling errors (pAM-pBF-(dBD-dAF))"+sprintf("%6.3g",value))
        end if
      end if

 
      if (dycore.eq."FV3") then
        print("Checking physics dynamics coupling")
        mass_change_phys = (tracer_mass(t,m_cnst,ipAM)-tracer_mass(t,m_cnst,ipBF))
        print("Mass change in physics:               "+mass_change_phys)
        mass_change_p_d_coupling  = (tracer_mass(t,m_cnst,idBD)-tracer_mass(t,m_cnst,idED))
        print("Mass change in d_p_coupling: "+mass_change_p_d_coupling)
	error = mass_change_p_d_coupling-mass_change_phys
	print("Error  "+error)
        mass_change_dyn  = (tracer_mass(t,m_cnst,idBF)-tracer_mass(t,m_cnst,idED))
        print("Mass change in dynamics from physics: "+mass_change_dyn)

        mass_error_in_coupling = (tracer_mass(t,m_cnst,ipAM)-tracer_mass(t,m_cnst,idBF))
        print("Mass error in coulpling "+mass_error_in_coupling)
       end if
       print(" ")
       print("---------------------------------------------------------------------------------------------------")
       print(" ")
     end do
  ;
  ;----------------------------------------------------------------------------------------------------------------
  ;
  ; plotting stuff
  ;
  ;----------------------------------------------------------------------------------------------------------------
  ;
  if (t.eq.nt-1) then
  output_format = "pdf"
  wks = gsn_open_wks(output_format,panelTitle)
  gsn_define_colormap(wks,"NCV_jaisnd");3gauss")      
; gsn_define_colormap(wks,"BlAqGrYeOrReVi200")   ; choose colormap
; gsn_define_colormap(wks,"GMT_polar")
    



  res                     = True         ; plot modifications desired
  res@gsnMaximize         = False         ; Maximize size of plot in frame
  res@gsnSpreadColors     = True         ; Use full colormap, but start
  res@gsnDraw             = False           ; don't draw
  res@gsnFrame            = False           ; don't advance frame
  res@cnFillOn            = True         ; Turn on contour fill
;  res@cnFillMode          = "AreaFill"   ; Style of fill. You can also
                                         ; use "CellFill" and "RasterFill"
  res@cnLinesOn           = False        ; Turn off contour lines
  res@cnLineLabelsOn      = False        ; Turn off contour line labels
  res@lbLabelAutoStride   = True         ; Clean up labelbar labels.
  res@cnFillMode="RasterFill"


;  res@cnLevelSelectionMode = "ManualLevels"    ; manually set the contour levels with the following 3 resources
;  res@cnMinLevelValF  = -30.0                        ; set the minimum contour level
;  res@cnMaxLevelValF  = 30.0                   ; set the maximum contour level
;  res@cnLevelSpacingF = 1.0                  ; set the interval between contours

;  res@mpMinLatF = -30.                   ; set the minimum latitude = -70.
;  res@mpMaxLatF = 30.                    ; set the maximum latitude = 70.

; res@mpMinLonF = -40.                    ; set the minimum latitude = -70.
;  res@mpMaxLonF = 40.

  res@gsnStringFontHeightF = 0.02
  res@tiMainFontHeightF = 0.025
  res@lbLabelFontHeightF = 0.02
  res@tmXBLabelFontHeightF = 0.015
  res@tmYLLabelFontHeightF = 0.015
  res@pmLabelBarOrthogonalPosF = 0.2

  res@lbBoxLinesOn = False

  res_p=res

  plot = new(6,graphic) 



  res_p@cnLevelSelectionMode = "ManualLevels"    ; manually set the contour levels with the following 3 resources
  fmin = -150.
      fmax = 150.
  res_p@cnMinLevelValF  = fmin                        ; set the minimum contour level
  res_p@cnMaxLevelValF  = fmax                   ; set the maximum contour level
  res_p@cnLevelSpacingF = (fmax-fmin)/50.0                  ; set the interval between contours
  res_p@gsnLeftString        = "Total physics TE tendency"      

;  res_p@tiMainString        = " [W/m^2]"
      ;-148.2925226269828   max=1770.210548306306


  min_f = min(dparam(0:ncol_p-1))
  max_f = max(dparam(0:ncol_p-1))
  print("dparam "+min_f+" "+max_f)
  res_p@lbTitleString    = "global min = "+sprintf("%6.4g",min_f)+"  global  max="+sprintf("%6.4g",max_f)
  res_p@sfXArray            = lon1d
  res_p@sfYArray            = lat1d
  plot(0) = gsn_csm_contour_map_ce(wks,dparam(0:ncol_p-1),res_p)

  res_p@cnLevelSelectionMode = "ManualLevels"    ; manually set the contour levels with the following 3 resources
  res_p@cnMinLevelValF  =-20.                        ; set the minimum contour level
  res_p@cnMaxLevelValF  = 20.                   ; set the maximum contour level
  res_p@cnLevelSpacingF = 0.50                  ; set the interval between contours
;  res_p@gsnRightString       = "W/m~S~2~N"
  res_p@gsnLeftString        = "Pressure work error TE tendency (pwork)"
  res@gsnRightString        = "mean: "+sprintf("%3.2g",diff_DME         )+"W/m^2"
;  res@tiMainString        = "Pressure work error (pwork) [W/m^2]"


  min_f = min(pwork(0:ncol_p-1))
  max_f = max(pwork(0:ncol_p-1))
  res_p@lbTitleString    = "global min = "+sprintf("%6.4g",min_f)+"  global  max="+sprintf("%6.4g",max_f)
  res_p@sfXArray            = lon1d
  res_p@sfYArray            = lat1d
  plot(1) = gsn_csm_contour_map_ce(wks,pwork(0:ncol_p-1),res_p)

;  res@cnLevelSelectionMode = "AutomaticLevels"    ; manually set the contour levels with the following 3 resources
  res@cnLevelSelectionMode = "ManualLevels"    ; manually set the contour levels with the following 3 resources
  fmin   = -0.1
      fmax = 0.1
  res@cnMinLevelValF  = fmin                        ; set the minimum contour level
  res@cnMaxLevelValF  = fmax                   ; set the maximum contour level
  res@cnLevelSpacingF = (fmax-fmin)/50.0                  ; set the interval between contours
  res@gsnLeftString       = "All water in dp"        
  res@gsnRightString        = "mean: "+sprintf("%3.2g",diff_all_water_pel         )+"W/m^2"
;  min_f = min(all_water)
;  max_f = max(all_water); max((tot_energy_2d(0,t,idAR,0:ncol-1)-tot_energy_2d(0,t,idBD,0:ncol-1))/dtime)

;  res@lbTitleString    = "global min = "+sprintf("%6.4g",min_f)+"  global  max="+sprintf("%6.4g",max_f)
  res@sfXArray            = lon1d_gll
  res@sfYArray            = lat1d_gll
  plot(2) = gsn_csm_contour_map_ce(wks,all_water(0:ncol-1),res)
      printMinMax((tot_energy_2d(0,t,idAR,:)-tot_energy_2d(0,t,idBD,:))/dtime,False)
      ;min=-1489.836497456233   max=122.200959737036

  fmin   = -2.0
  fmax = 2.0
  res@cnMinLevelValF  = fmin                        ; set the minimum contour level
  res@cnMaxLevelValF  = fmax                   ; set the maximum contour level
  res@cnLevelSpacingF = (fmax-fmin)/50.0                  ; set the interval between contours


  res@gsnLeftString        = "cpwv"        
  res@gsnRightString        = "mean: "+sprintf("%3.2g",diff_cpwv         )+"W/m^2"


;  min_f = min(hvis(0:ncol-1))
;  max_f = max(hvis(0:ncol-1))
;  res@lbTitleString    = "global min = "+sprintf("%6.4g",min_f)+"  global  max="+sprintf("%6.4g",max_f)

  plot(3) = gsn_csm_contour_map_ce(wks,cpwv_2d(0:ncol-1),res)	

;      average_hvis = sum(area(:)*(tot_energy_2d(0,t,idAH,:)-tot_energy_2d(0,t,idBH,:))/dtime)/(4.0*PI)
;      print("ave_hvis = "+average_hvis)
  res@cnLevelSelectionMode = "ManualLevels"    ; manually set the contour levels with the following 3 resources
 if (dycore.eq."SE") then
  res@cnMinLevelValF  = -2                        ; set the minimum contour level
  res@cnMaxLevelValF  = 2.0                   ; set the maximum contour level
  res@cnLevelSpacingF = 0.05                  ; set the interval between contours 
else
  res@cnMinLevelValF  = -10.0                        ; set the minimum contour level
  res@cnMaxLevelValF  = 10.0                   ; set the maximum contour level
  res@cnLevelSpacingF = 0.5                  ; set the interval between contours 
end if
  res@gsnLeftString        = "K TERM"      
  res@gsnRightString        = "mean: "+sprintf("%4g",diff_KPHIS         )+"W/m^2"

  min_f = min(K_2d(0:ncol-1))
  max_f = max(K_2d(0:ncol-1))
  res@lbTitleString    = "global min = "+sprintf("%6.4g",min_f)+"  global  max="+sprintf("%6.4g",max_f)

  plot(4) = gsn_csm_contour_map_ce(wks,K_2d(0:ncol-1),res)


  res@cnLevelSelectionMode = "ManualLevels"    ; manually set the contour levels with the following 3 resources
;  res@cnMinLevelValF  = -10.0                        ; set the minimum contour level
;  res@cnMaxLevelValF  = 10.0                   ; set the maximum contour level
;  res@cnLevelSpacingF = 0.5                  ; set the interval between contours
  res@cnMinLevelValF  = -20.0                        ; set the minimum contour level
  res@cnMaxLevelValF  = 20.0                   ; set the maximum contour level
  res@cnLevelSpacingF = 0.5                  ; set the interval between contours
  res@gsnLeftString        = "cpT terms"  
  res@gsnRightString        = "mean: "+sprintf("%4g",diff_cpT         )+"W/m^2"
  min_f = min(cpT_2d(0:ncol-1))
  max_f = max(cpT_2d(0:ncol-1))
  res@lbTitleString    = "global min = "+sprintf("%6.4g",min_f)+"  global  max="+sprintf("%6.4g",max_f)

  cppp = cpT_2d(0:ncol-1)

  plot(5) = gsn_csm_contour_map_ce(wks,cpT_2d(0:ncol-1),res)
 

  resP                  = True                   ; modify the panel plot
  resP@gsnFrame         = False                  ; don't advance panel plot
  resP@gsnPanelFigureStrings = (/"a","b","c","d","e","f"/)
  resP@gsnPanelYWhiteSpacePercent = 4.0
  resP@gsnPanelXWhiteSpacePercent = 4.0

;  resP@gsnPanelLabelBar = True                   ; add common colorbar
;  resP@gsnPanelBottom   = 0.05
  if (dycore.eq."FV3") then    
     resP@txString   = "TE tendencies for the default CAM-FV3 configuration (AMIP)"
  end if
  if (dycore.eq."SE") then    
     resP@txString   = "TE tendencies for the default CAM-SE configuration (AMIP)"
  end if
  resP@txString   = panelTitle
     gsn_panel(wks,plot,(/2,3/),resP)               ; now draw as one plot



   txres               = True
;  txres@txFontHeightF = 0.015
;  gsn_text_ndc(wks,"Figure 1: A nifty panel plot",0.5,0.02,txres)
   frame(wks)

;
;
; NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT
; NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT
; NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT
; NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT
; NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT
;
;
   delete(wks)
   delete(res)
   delete(plot)
   delete(res_p)



  filePlot = panelTitle+"_CAM_ENERGY_EQN"
  wks = gsn_open_wks(output_format,filePlot)
  gsn_define_colormap(wks,"NCV_jaisnd");3gauss")      
    

  res                     = True         ; plot modifications desired
  res@gsnMaximize         = False         ; Maximize size of plot in frame
  res@gsnSpreadColors     = True         ; Use full colormap, but start
  res@gsnDraw             = False           ; don't draw
  res@gsnFrame            = False           ; don't advance frame
  res@cnFillOn            = True         ; Turn on contour fill
;  res@cnFillMode          = "AreaFill"   ; Style of fill. You can also
                                         ; use "CellFill" and "RasterFill"
  res@cnLinesOn           = False        ; Turn off contour lines
  res@cnLineLabelsOn      = False        ; Turn off contour line labels
  res@lbLabelAutoStride   = True         ; Clean up labelbar labels.
  res@cnFillMode="RasterFill"

  res@gsnStringFontHeightF = 0.02
  res@tiMainFontHeightF = 0.025
  res@lbLabelFontHeightF = 0.02
  res@tmXBLabelFontHeightF = 0.015
  res@tmYLLabelFontHeightF = 0.015
  res@pmLabelBarOrthogonalPosF = 0.2

  res@lbBoxLinesOn = False

  res_p=res

  plot = new(7,graphic) 

; ***********************************************
; Left-hand side of CAM energy equation
; ***********************************************

  fld = CAM_LHS(0:ncol_p-1)
  str = "dE/dt parameterizations (left-hand side)"
  fmin = -150.
  fmax = 150.

  res_p@cnLevelSelectionMode = "ManualLevels"    ; manually set the contour levels with the following 3 resources
  res_p@cnMinLevelValF  = fmin                        ; set the minimum contour level
  res_p@cnMaxLevelValF  = fmax                   ; set the maximum contour level
  res_p@cnLevelSpacingF = (fmax-fmin)/50.0                  ; set the interval between contours

  min_f = min(fld(0:ncol_p-1))
  max_f = max(fld(0:ncol_p-1))
  res_p@lbTitleString    = "global min = "+sprintf("%6.4g",min_f)+"  global  max="+sprintf("%6.4g",max_f)
  print(str+" min/max: "+min_f+"  "+max_f)

  res_p@sfXArray            = lon1d
  res_p@sfYArray            = lat1d

  res_p@gsnLeftString        = str
  global_ave = sum(area_p(0:ncol_p-1)*fld(0:ncol_p-1))*inv_area_sphere_p
  print(str+" global ave: "+global_ave)
  res_p@gsnRightString        = "mean: "+sprintf("%3.2g", global_ave        )+" W/m^2"
  plot(0) = gsn_csm_contour_map_ce(wks,fld,res_p)

; ***********************************************
; Rhs-hand side of CAM energy equation
; ***********************************************

  fld = CAM_RHS(0:ncol_p-1)
  str = "Flux into atmosphere (right-hand side)"
  fmin = -150.
  fmax = 150.

  res_p@cnLevelSelectionMode = "ManualLevels"    ; manually set the contour levels with the following 3 resources
  res_p@cnMinLevelValF  = fmin                        ; set the minimum contour level
  res_p@cnMaxLevelValF  = fmax                   ; set the maximum contour level
  res_p@cnLevelSpacingF = (fmax-fmin)/50.0                  ; set the interval between contours

  min_f = min(fld(0:ncol_p-1))
  max_f = max(fld(0:ncol_p-1))
  res_p@lbTitleString    = "global min = "+sprintf("%6.4g",min_f)+"  global  max="+sprintf("%6.4g",max_f)
  print(str+" min/max: "+min_f+"  "+max_f)

  res_p@sfXArray            = lon1d
  res_p@sfYArray            = lat1d

  res_p@gsnLeftString        = str
  global_ave = sum(area_p(0:ncol_p-1)*fld(0:ncol_p-1))*inv_area_sphere_p
  print(str+" global ave: "+global_ave)
  res_p@gsnRightString        = "mean: "+sprintf("%3.2g", global_ave        )+" W/m^2"
  plot(1) = gsn_csm_contour_map_ce(wks,fld,res_p)

; ***********************************************
; Residual (right-hand side minus left-hand side
; ***********************************************

  fld = CAM_RHS(0:ncol_p-1)-CAM_LHS(0:ncol_p-1)
  str = "Residual: plot (a) minus (b)"
  fmin = -150.
  fmax = 150.

  res_p@cnLevelSelectionMode = "ManualLevels"    ; manually set the contour levels with the following 3 resources
  res_p@cnMinLevelValF  = fmin                        ; set the minimum contour level
  res_p@cnMaxLevelValF  = fmax                   ; set the maximum contour level
  res_p@cnLevelSpacingF = (fmax-fmin)/50.0                  ; set the interval between contours

  min_f = min(fld(0:ncol_p-1))
  max_f = max(fld(0:ncol_p-1))
  res_p@lbTitleString    = "global min = "+sprintf("%6.4g",min_f)+"  global  max="+sprintf("%6.4g",max_f)
  print(str+" min/max: "+min_f+"  "+max_f)

  res_p@sfXArray            = lon1d
  res_p@sfYArray            = lat1d

  res_p@gsnLeftString        = str
  global_ave = sum(area_p(0:ncol_p-1)*fld(0:ncol_p-1))*inv_area_sphere_p
  print(str+" global ave: "+global_ave)
  res_p@gsnRightString        = "mean: "+sprintf("%3.2g", global_ave        )+" W/m^2"
  plot(2) = gsn_csm_contour_map_ce(wks,fld,res_p)

; ***********************************************

; ***********************************************

  fld = LS*FWV(t,0:ncol_p-1)
  str = "Latent heat flux: Fnet^wv*LS00"
  fmin = -150.
  fmax = 150.

  res_p@cnLevelSelectionMode = "ManualLevels"    ; manually set the contour levels with the following 3 resources
  res_p@cnMinLevelValF  = fmin                        ; set the minimum contour level
  res_p@cnMaxLevelValF  = fmax                   ; set the maximum contour level
  res_p@cnLevelSpacingF = (fmax-fmin)/50.0                  ; set the interval between contours

  min_f = min(fld(0:ncol_p-1))
  max_f = max(fld(0:ncol_p-1))
  res_p@lbTitleString    = "global min = "+sprintf("%6.4g",min_f)+"  global  max="+sprintf("%6.4g",max_f)
  print(str+" min/max: "+min_f+"  "+max_f)

  res_p@sfXArray            = lon1d
  res_p@sfYArray            = lat1d

  res_p@gsnLeftString        = str
  global_ave = sum(area_p(0:ncol_p-1)*fld(0:ncol_p-1))*inv_area_sphere_p
  print(str+" global ave: "+global_ave)
  res_p@gsnRightString        = "mean: "+sprintf("%5.4g", global_ave        )+" W/m^2"
  plot(3) = gsn_csm_contour_map_ce(wks,fld,res_p)

; ***********************************************

; ***********************************************

  fld = LF*FLIQ(t,0:ncol_p-1)
;  fld = FTURB(t,0:ncol_p-1)
  str = "Latent heat flux: Fnet^liq*LF00"
  fmin = -150.
  fmax = 150.

  res_p@cnLevelSelectionMode = "ManualLevels"    ; manually set the contour levels with the following 3 resources
  res_p@cnMinLevelValF  = fmin                        ; set the minimum contour level
  res_p@cnMaxLevelValF  = fmax                   ; set the maximum contour level
  res_p@cnLevelSpacingF = (fmax-fmin)/50.0                  ; set the interval between contours

  min_f = min(fld(0:ncol_p-1))
  max_f = max(fld(0:ncol_p-1))
  res_p@lbTitleString    = "global min = "+sprintf("%6.4g",min_f)+"  global  max="+sprintf("%6.4g",max_f)
  print(str+" min/max: "+min_f+"  "+max_f)

  res_p@sfXArray            = lon1d
  res_p@sfYArray            = lat1d

  res_p@gsnLeftString        = str
  global_ave = sum(area_p(0:ncol_p-1)*fld(0:ncol_p-1))*inv_area_sphere_p
  print(str+" global ave: "+global_ave)
  res_p@gsnRightString        = "mean: "+sprintf("%6.4g", global_ave        )+" W/m^2"
  plot(4) = gsn_csm_contour_map_ce(wks,fld,res_p)

  delete(fld)

  fld = FTURB(t,0:ncol_p-1)
  str = "Turb+rad flux into atmosphere"
  fmin = -150.
  fmax = 150.

  res_p@cnLevelSelectionMode = "ManualLevels"    ; manually set the contour levels with the following 3 resources
  res_p@cnMinLevelValF  = fmin                        ; set the minimum contour level
  res_p@cnMaxLevelValF  = fmax                   ; set the maximum contour level
  res_p@cnLevelSpacingF = (fmax-fmin)/50.0                  ; set the interval between contours

  min_f = min(fld(0:ncol_p-1))
  max_f = max(fld(0:ncol_p-1))
  res_p@lbTitleString    = "global min = "+sprintf("%6.4g",min_f)+"  global  max="+sprintf("%6.4g",max_f)
  print(str+" min/max: "+min_f+"  "+max_f)

  res_p@sfXArray            = lon1d
  res_p@sfYArray            = lat1d

  res_p@gsnLeftString        = str
  global_ave = sum(area_p(0:ncol_p-1)*fld(0:ncol_p-1))*inv_area_sphere_p
  print(str+" global ave: "+global_ave)
  res_p@gsnRightString        = "mean: "+sprintf("%3.2g", global_ave )+" (kg/kg)/(s m^2)"
  plot(5) = gsn_csm_contour_map_ce(wks,fld,res_p)

; ***********************************************
; DME adjust
; ***********************************************

  delete(fld)
  fld = pwork(0:ncol_p-1)
  str = "Dry-mass adjustment"
  fmin = -30.
  fmax = 30.

  res_p@cnLevelSelectionMode = "ManualLevels"    ; manually set the contour levels with the following 3 resources
  res_p@cnMinLevelValF  = fmin                        ; set the minimum contour level
  res_p@cnMaxLevelValF  = fmax                   ; set the maximum contour level
  res_p@cnLevelSpacingF = (fmax-fmin)/50.0                  ; set the interval between contours

  min_f = min(fld(0:ncol_p-1))
  max_f = max(fld(0:ncol_p-1))
  res_p@lbTitleString    = "global min = "+sprintf("%6.4g",min_f)+"  global  max="+sprintf("%6.4g",max_f)
  print(str+" min/max: "+min_f+"  "+max_f)

  res_p@sfXArray            = lon1d
  res_p@sfYArray            = lat1d

  res_p@gsnLeftString        = str
  global_ave = sum(area_p(0:ncol_p-1)*fld(0:ncol_p-1))*inv_area_sphere_p
  print(str+" global ave: "+global_ave)
  res_p@gsnRightString        = "mean: "+sprintf("%3.2g", global_ave        )+" W/m^2"
  plot(6) = gsn_csm_contour_map_ce(wks,fld,res_p)


  delete(resP)

  resP                  = True                   ; modify the panel plot
  resP@gsnFrame         = False                  ; don't advance panel plot
  resP@gsnPanelFigureStrings = (/"a","b","c","d","e","f"/)
  resP@gsnPanelYWhiteSpacePercent = 4.0
  resP@gsnPanelXWhiteSpacePercent = 4.0

  resP@txString   = "CAM total energy equation terms"
  gsn_panel(wks,plot,(/3,3/),resP)               ; now draw as one plot

  txres               = True
  frame(wks)
;
; NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT
; NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT
; NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT
; NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT
; NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT
;
;
   delete(wks)
   delete(res)
   delete(plot)
   delete(res_p)



  filePlot = panelTitle+"water_flux"
  wks = gsn_open_wks(output_format,filePlot)
  gsn_define_colormap(wks,"NCV_jaisnd");3gauss")      
    

  res                     = True         ; plot modifications desired
  res@gsnMaximize         = False         ; Maximize size of plot in frame
  res@gsnSpreadColors     = True         ; Use full colormap, but start
  res@gsnDraw             = False           ; don't draw
  res@gsnFrame            = False           ; don't advance frame
  res@cnFillOn            = True         ; Turn on contour fill
;  res@cnFillMode          = "AreaFill"   ; Style of fill. You can also
                                         ; use "CellFill" and "RasterFill"
  res@cnLinesOn           = False        ; Turn off contour lines
  res@cnLineLabelsOn      = False        ; Turn off contour line labels
  res@lbLabelAutoStride   = True         ; Clean up labelbar labels.
  res@cnFillMode="RasterFill"

  res@gsnStringFontHeightF = 0.02
  res@tiMainFontHeightF = 0.025
  res@lbLabelFontHeightF = 0.02
  res@tmXBLabelFontHeightF = 0.015
  res@tmYLLabelFontHeightF = 0.015
  res@pmLabelBarOrthogonalPosF = 0.2

  res@lbBoxLinesOn = False

  res_p=res

  plot = new(4,graphic) 

; ***********************************************
; total water flux
; ***********************************************

  fld = FH2O(t,0:ncol_p-1)
  fld = 1000*24*3600*fld/dtime ;convert from m/s to mm/day
  str = "Flux of total water into atmosphere"
  fmin = -10.0
  fmax = 10.0

  res_p@cnLevelSelectionMode = "ManualLevels"    ; manually set the contour levels with the following 3 resources
  res_p@cnMinLevelValF  = fmin                        ; set the minimum contour level
  res_p@cnMaxLevelValF  = fmax                   ; set the maximum contour level
  res_p@cnLevelSpacingF = (fmax-fmin)/50.0                  ; set the interval between contours

  min_f = min(fld(0:ncol_p-1))
  max_f = max(fld(0:ncol_p-1))
  res_p@lbTitleString    = "global min = "+sprintf("%6.4g",min_f)+"  global  max="+sprintf("%6.4g",max_f)
  print(str+" min/max: "+min_f+"  "+max_f)

  res_p@sfXArray            = lon1d
  res_p@sfYArray            = lat1d

  res_p@gsnLeftString        = str
  global_ave = sum(area_p(0:ncol_p-1)*fld(0:ncol_p-1))*inv_area_sphere_p
  print(str+" global ave: "+global_ave)
  res_p@gsnRightString        = "mean: "+sprintf("%6.4g", global_ave        )+" mm/day"
  plot(0) = gsn_csm_contour_map_ce(wks,fld,res_p)
; ***********************************************
; water vapor flux
; ***********************************************
  fld = FWV(t,0:ncol_p-1)
  fld = 1000*24*3600*fld/dtime ;convert from m/s to mm/day
  str = "Flux of water vapor into atmosphere"
;  fmin = -0.0001
;  fmax = 0.0001

  res_p@cnLevelSelectionMode = "ManualLevels"    ; manually set the contour levels with the following 3 resources
  res_p@cnMinLevelValF  = fmin                        ; set the minimum contour level
  res_p@cnMaxLevelValF  = fmax                   ; set the maximum contour level
  res_p@cnLevelSpacingF = (fmax-fmin)/50.0                  ; set the interval between contours

  min_f = min(fld(0:ncol_p-1))
  max_f = max(fld(0:ncol_p-1))
  res_p@lbTitleString    = "global min = "+sprintf("%6.4g",min_f)+"  global  max="+sprintf("%6.4g",max_f)
  print(str+" min/max: "+min_f+"  "+max_f)

  res_p@sfXArray            = lon1d
  res_p@sfYArray            = lat1d

  res_p@gsnLeftString        = str
  global_ave = sum(area_p(0:ncol_p-1)*fld(0:ncol_p-1))*inv_area_sphere_p
  print(str+" global ave: "+global_ave)
  res_p@gsnRightString        = "mean: "+sprintf("%6.4g", global_ave        )+" mm/day"
  plot(1) = gsn_csm_contour_map_ce(wks,fld,res_p)
; ***********************************************
; liq flux
; ***********************************************
  fld = FLIQ(t,0:ncol_p-1)
  fld = 1000*24*3600*fld/dtime ;convert from m/s to mm/day
  str = "Flux of liquid into atmosphere"
;  fmin = -0.0001
;  fmax = 0.0001

  res_p@cnLevelSelectionMode = "ManualLevels"    ; manually set the contour levels with the following 3 resources
  res_p@cnMinLevelValF  = fmin                        ; set the minimum contour level
  res_p@cnMaxLevelValF  = fmax                   ; set the maximum contour level
  res_p@cnLevelSpacingF = (fmax-fmin)/50.0                  ; set the interval between contours

  min_f = min(fld(0:ncol_p-1))
  max_f = max(fld(0:ncol_p-1))
  res_p@lbTitleString    = "global min = "+sprintf("%6.4g",min_f)+"  global  max="+sprintf("%6.4g",max_f)
  print(str+" min/max: "+min_f+"  "+max_f)

  res_p@sfXArray            = lon1d
  res_p@sfYArray            = lat1d

  res_p@gsnLeftString        = str
  global_ave = sum(area_p(0:ncol_p-1)*fld(0:ncol_p-1))*inv_area_sphere_p
  print(str+" global ave: "+global_ave)
  res_p@gsnRightString        = "mean: "+sprintf("%6.4g", global_ave        )+" mm/day"
  plot(2) = gsn_csm_contour_map_ce(wks,fld,res_p)
; ***********************************************
; ice flux
; ***********************************************
  fld = FICE(t,0:ncol_p-1)
  fld = 1000*24*3600*fld/dtime ;convert from m/s to mm/day
  str = "Flux of ice into atmosphere"

;  fmin = -0.0001
;  fmax = 0.0001

  res_p@cnLevelSelectionMode = "ManualLevels"    ; manually set the contour levels with the following 3 resources
  res_p@cnMinLevelValF  = fmin                        ; set the minimum contour level
  res_p@cnMaxLevelValF  = fmax                   ; set the maximum contour level
  res_p@cnLevelSpacingF = (fmax-fmin)/50.0                  ; set the interval between contours

  min_f = min(fld(0:ncol_p-1))
  max_f = max(fld(0:ncol_p-1))
  res_p@lbTitleString    = "global min = "+sprintf("%6.4g",min_f)+"  global  max="+sprintf("%6.4g",max_f)
  print(str+" min/max: "+min_f+"  "+max_f)

  res_p@sfXArray            = lon1d
  res_p@sfYArray            = lat1d

  res_p@gsnLeftString        = str
  global_ave = sum(area_p(0:ncol_p-1)*fld(0:ncol_p-1))*inv_area_sphere_p
  print(str+" global ave: "+global_ave)
  res_p@gsnRightString        = "mean: "+sprintf("%6.4g", global_ave        )+" mm/day"
  plot(3) = gsn_csm_contour_map_ce(wks,fld,res_p)

  delete(resP)

  resP                  = True                   ; modify the panel plot
  resP@gsnFrame         = False                  ; don't advance panel plot
  resP@gsnPanelFigureStrings = (/"a","b","c","d"/)
  resP@gsnPanelYWhiteSpacePercent = 4.0
  resP@gsnPanelXWhiteSpacePercent = 4.0

  resP@txString   = "Water fluxes"
  gsn_panel(wks,plot,(/2,2/),resP)               ; now draw as one plot

  txres               = True
  frame(wks)

;
; NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT
; NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT
; NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT
; NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT
; NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT
;
;
   delete(wks)
   delete(res)
   delete(plot)
   delete(res_p)



  filePlot = panelTitle+"all_water_terms"
  wks = gsn_open_wks(output_format,filePlot)
  gsn_define_colormap(wks,"NCV_jaisnd");3gauss")      
    

  res                     = True         ; plot modifications desired
  res@gsnMaximize         = False         ; Maximize size of plot in frame
  res@gsnSpreadColors     = True         ; Use full colormap, but start
  res@gsnDraw             = False           ; don't draw
  res@gsnFrame            = False           ; don't advance frame
  res@cnFillOn            = True         ; Turn on contour fill
;  res@cnFillMode          = "AreaFill"   ; Style of fill. You can also
                                         ; use "CellFill" and "RasterFill"
  res@cnLinesOn           = False        ; Turn off contour lines
  res@cnLineLabelsOn      = False        ; Turn off contour line labels
  res@lbLabelAutoStride   = True         ; Clean up labelbar labels.
  res@cnFillMode="RasterFill"

  res@gsnStringFontHeightF = 0.02
  res@tiMainFontHeightF = 0.025
  res@lbLabelFontHeightF = 0.02
  res@tmXBLabelFontHeightF = 0.015
  res@tmYLLabelFontHeightF = 0.015
  res@pmLabelBarOrthogonalPosF = 0.2

  res@lbBoxLinesOn = False

  res_p=res

  plot = new(6,graphic) 
; ***********************************************
  ; DME adjust
; ***********************************************

  delete(fld)
  fld = pwork(0:ncol_p-1)
  str = "Dry-mass adjustment"
  fmin = -30.
  fmax = 30.

  res_p@cnLevelSelectionMode = "ManualLevels"    ; manually set the contour levels with the following 3 resources
  res_p@cnMinLevelValF  = fmin                        ; set the minimum contour level
  res_p@cnMaxLevelValF  = fmax                   ; set the maximum contour level
  res_p@cnLevelSpacingF = (fmax-fmin)/50.0                  ; set the interval between contours

  min_f = min(fld(0:ncol_p-1))
  max_f = max(fld(0:ncol_p-1))
  res_p@lbTitleString    = "global min = "+sprintf("%6.4g",min_f)+"  global  max="+sprintf("%6.4g",max_f)
  print(str+" min/max: "+min_f+"  "+max_f)

  res_p@sfXArray            = lon1d
  res_p@sfYArray            = lat1d

  res_p@gsnLeftString        = str
  global_ave = sum(area_p(0:ncol_p-1)*fld(0:ncol_p-1))*inv_area_sphere_p
  print(str+" global ave: "+global_ave)
  res_p@gsnRightString        = "mean: "+sprintf("%3.2g", global_ave        )+" W/m^2"
  plot(0) = gsn_csm_contour_map_ce(wks,fld,res_p)


; ***********************************************
; L(T)
; ***********************************************

  delete(fld)
  fld = lvariableL_2d(0:ncol_p-1)
  str = "L(T) imbalance"
  fmin = -30.
  fmax = 30.

  res_p@cnLevelSelectionMode = "ManualLevels"    ; manually set the contour levels with the following 3 resources
  res_p@cnMinLevelValF  = fmin                        ; set the minimum contour level
  res_p@cnMaxLevelValF  = fmax                   ; set the maximum contour level
  res_p@cnLevelSpacingF = (fmax-fmin)/50.0                  ; set the interval between contours

  min_f = min(fld(0:ncol_p-1))
  max_f = max(fld(0:ncol_p-1))
  res_p@lbTitleString    = "global min = "+sprintf("%6.4g",min_f)+"  global  max="+sprintf("%6.4g",max_f)
  print(str+" min/max: "+min_f+"  "+max_f)

  res_p@sfXArray            = lon1d
  res_p@sfYArray            = lat1d

  res_p@gsnLeftString        = str
  global_ave = sum(area_p(0:ncol_p-1)*fld(0:ncol_p-1))*inv_area_sphere_p
  print(str+" global ave: "+global_ave)
  res_p@gsnRightString        = "mean: "+sprintf("%3.2g", global_ave        )+" W/m^2"
  plot(1) = gsn_csm_contour_map_ce(wks,fld,res_p)
; ***********************************************
; 
; ***********************************************

  delete(fld)
    fld = FLH_T(t,0:ncol_p-1)
  str = "enthalpy flux with variable cp"
  fmin = -30.
  fmax = 30.

  res_p@cnLevelSelectionMode = "ManualLevels"    ; manually set the contour levels with the following 3 resources
  res_p@cnMinLevelValF  = fmin                        ; set the minimum contour level
  res_p@cnMaxLevelValF  = fmax                   ; set the maximum contour level
  res_p@cnLevelSpacingF = (fmax-fmin)/50.0                  ; set the interval between contours

  min_f = min(fld(0:ncol_p-1))
  max_f = max(fld(0:ncol_p-1))
  res_p@lbTitleString    = "global min = "+sprintf("%6.4g",min_f)+"  global  max="+sprintf("%6.4g",max_f)
  print(str+" min/max: "+min_f+"  "+max_f)

  res_p@sfXArray            = lon1d
  res_p@sfYArray            = lat1d

  res_p@gsnLeftString        = str
  global_ave = sum(area_p(0:ncol_p-1)*fld(0:ncol_p-1))*inv_area_sphere_p
  print(str+" global ave: "+global_ave)
  res_p@gsnRightString        = "mean: "+sprintf("%3.2g", global_ave        )+" W/m^2"
  plot(2) = gsn_csm_contour_map_ce(wks,fld,res_p)



  delete(resP)
  resP                  = True                   ; modify the panel plot
  resP@gsnFrame         = False                  ; don't advance panel plot
  resP@gsnPanelFigureStrings = (/"a","b","c","d","e","f","g","h","i","j"/)
  resP@gsnPanelYWhiteSpacePercent = 4.0
  resP@gsnPanelXWhiteSpacePercent = 4.0

  resP@txString   = "Epsilon terms"
  gsn_panel(wks,plot,(/1,3/),resP)               ; now draw as one plot

  txres               = True
  frame(wks)


;
;
; NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT
; NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT
; NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT
; NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT
; NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT
;
;
   delete(wks)
   delete(res)
   delete(plot)
   delete(res_p)



  filePlot = panelTitle+"3"
  wks = gsn_open_wks(output_format,filePlot)
  gsn_define_colormap(wks,"NCV_jaisnd");3gauss")      
    

  res                     = True         ; plot modifications desired
  res@gsnMaximize         = False         ; Maximize size of plot in frame
  res@gsnSpreadColors     = True         ; Use full colormap, but start
  res@gsnDraw             = False           ; don't draw
  res@gsnFrame            = False           ; don't advance frame
  res@cnFillOn            = True         ; Turn on contour fill
;  res@cnFillMode          = "AreaFill"   ; Style of fill. You can also
                                         ; use "CellFill" and "RasterFill"
  res@cnLinesOn           = False        ; Turn off contour lines
  res@cnLineLabelsOn      = False        ; Turn off contour line labels
  res@lbLabelAutoStride   = True         ; Clean up labelbar labels.
  res@cnFillMode="RasterFill"

  res@gsnStringFontHeightF = 0.02
  res@tiMainFontHeightF = 0.025
  res@lbLabelFontHeightF = 0.02
  res@tmXBLabelFontHeightF = 0.015
  res@tmYLLabelFontHeightF = 0.015
  res@pmLabelBarOrthogonalPosF = 0.2

  res@lbBoxLinesOn = False

  res_p=res

  plot = new(12,graphic) 

; ***********************************************
; cpwv
; ***********************************************

  delete(fld)
  fld = cpwv_2d(0:ncol_p-1)
  str = "dE/dt using cpwv"
  fmin = -1.
  fmax = 1.

  res_p@cnLevelSelectionMode = "ManualLevels"    ; manually set the contour levels with the following 3 resources
  res_p@cnMinLevelValF  = fmin                        ; set the minimum contour level
  res_p@cnMaxLevelValF  = fmax                   ; set the maximum contour level
  res_p@cnLevelSpacingF = (fmax-fmin)/50.0                  ; set the interval between contours

  min_f = min(fld(0:ncol_p-1))
  max_f = max(fld(0:ncol_p-1))
  res_p@lbTitleString    = "global min = "+sprintf("%6.4g",min_f)+"  global  max="+sprintf("%6.4g",max_f)
  print(str+" min/max: "+min_f+"  "+max_f)

  res_p@sfXArray            = lon1d
  res_p@sfYArray            = lat1d

  res_p@gsnLeftString        = str
  global_ave = sum(area_p(0:ncol_p-1)*fld(0:ncol_p-1))*inv_area_sphere_p
  print(str+" global ave: "+global_ave)
  res_p@gsnRightString        = "mean: "+sprintf("%3.2g", global_ave        )+" W/m^2"
  plot(0) = gsn_csm_contour_map_ce(wks,fld,res_p)

; ***********************************************
; all water in pdel
; ***********************************************

  fld = all_water(0:ncol_p-1)
  str = "dE/dt total H2O (total H2O const) "
  fmin = -0.05
  fmax = 0.05

  res_p@cnLevelSelectionMode = "ManualLevels"    ; manually set the contour levels with the following 3 resources
  res_p@cnMinLevelValF  = fmin                        ; set the minimum contour level
  res_p@cnMaxLevelValF  = fmax                   ; set the maximum contour level
  res_p@cnLevelSpacingF = (fmax-fmin)/50.0                  ; set the interval between contours

  min_f = min(fld(0:ncol_p-1))
  max_f = max(fld(0:ncol_p-1))
  res_p@lbTitleString    = "global min = "+sprintf("%6.4g",min_f)+"  global  max="+sprintf("%6.4g",max_f)
  print(str+" min/max: "+min_f+"  "+max_f)

  res_p@sfXArray            = lon1d
  res_p@sfYArray            = lat1d

  res_p@gsnLeftString        = str
  global_ave = sum(area_p(0:ncol_p-1)*fld(0:ncol_p-1))*inv_area_sphere_p
  print(str+" global ave: "+global_ave)
  res_p@gsnRightString        = "mean: "+sprintf("%3.2g", global_ave        )+" W/m^2"
  plot(1) = gsn_csm_contour_map_ce(wks,fld,res_p)

; ***********************************************
; d/dt(KE update all water)
; ***********************************************

  fld = K_2d(0:ncol_p-1)
  str = "d/dt(KE total H2O varying)"
  fmin = -0.01
  fmax =  0.01

  res_p@cnLevelSelectionMode = "ManualLevels"    ; manually set the contour levels with the following 3 resources
  res_p@cnMinLevelValF  = fmin                        ; set the minimum contour level
  res_p@cnMaxLevelValF  = fmax                   ; set the maximum contour level
  res_p@cnLevelSpacingF = (fmax-fmin)/50.0                  ; set the interval between contours

  min_f = min(fld(0:ncol_p-1))
  max_f = max(fld(0:ncol_p-1))
  res_p@lbTitleString    = "global min = "+sprintf("%6.4g",min_f)+"  global  max="+sprintf("%6.4g",max_f)
  print(str+" min/max: "+min_f+"  "+max_f)

  res_p@sfXArray            = lon1d
  res_p@sfYArray            = lat1d

  res_p@gsnLeftString        = str
  global_ave = sum(area_p(0:ncol_p-1)*fld(0:ncol_p-1))*inv_area_sphere_p
  print(str+" global ave: "+global_ave)
  res_p@gsnRightString        = "mean: "+sprintf("%3.2g", global_ave        )+" W/m^2"
  plot(2) = gsn_csm_contour_map_ce(wks,fld,res_p)

; ***********************************************
; d/dt(cp*T term total H2O varying
; ***********************************************

  fld = cpT_2d(0:ncol_p-1);*cpair
  str = "d/dt(cp*T term total H2O varying)"
  fmin = -30.
  fmax = 30.

  res_p@cnLevelSelectionMode = "ManualLevels"    ; manually set the contour levels with the following 3 resources
  res_p@cnMinLevelValF  = fmin                        ; set the minimum contour level
  res_p@cnMaxLevelValF  = fmax                   ; set the maximum contour level
  res_p@cnLevelSpacingF = (fmax-fmin)/50.0                  ; set the interval between contours

  min_f = min(fld(0:ncol_p-1))
  max_f = max(fld(0:ncol_p-1))
  res_p@lbTitleString    = "global min = "+sprintf("%6.4g",min_f)+"  global  max="+sprintf("%6.4g",max_f)
  print(str+" min/max: "+min_f+"  "+max_f)

  res_p@sfXArray            = lon1d
  res_p@sfYArray            = lat1d

  res_p@gsnLeftString        = str
  global_ave = sum(area_p(0:ncol_p-1)*fld(0:ncol_p-1))*inv_area_sphere_p
  print(str+" global ave: "+global_ave)
  res_p@gsnRightString        = "mean: "+sprintf("%3.2g", global_ave        )+" W/m^2"
  plot(3) = gsn_csm_contour_map_ce(wks,fld,res_p)


; ***********************************************
; update all H2O in PHIS
; ***********************************************

  fld = PHIS_2d(0:ncol_p-1)
  str = "d/dt(PHIS term total H2O vayring)"
  fmin = -1.
  fmax = 1.

;  res_p@cnLevelSelectionMode = "AutomaticLevels"    ; manually set the contour levels with the following 3 resources
  res_p@cnLevelSelectionMode = "ManualLevels"    ; manually set the contour levels with the following 3 resources
  res_p@cnMinLevelValF  = fmin                        ; set the minimum contour level
  res_p@cnMaxLevelValF  = fmax                   ; set the maximum contour level
  res_p@cnLevelSpacingF = (fmax-fmin)/50.0                  ; set the interval between contours

  min_f = min(fld(0:ncol_p-1))
  max_f = max(fld(0:ncol_p-1))
  res_p@lbTitleString    = "global min = "+sprintf("%6.4g",min_f)+"  global  max="+sprintf("%6.4g",max_f)
  print(str+" min/max: "+min_f+"  "+max_f)

  res_p@sfXArray            = lon1d
  res_p@sfYArray            = lat1d

  res_p@gsnLeftString        = str
  global_ave = sum(area_p(0:ncol_p-1)*fld(0:ncol_p-1))*inv_area_sphere_p
  print(str+" global ave: "+global_ave)
  res_p@gsnRightString        = "mean: "+sprintf("%3.2g", global_ave        )+" W/m^2"
  plot(4) = gsn_csm_contour_map_ce(wks,fld,res_p)

  delete(fld)
; ***********************************************
; Missing enthalpy flux at the surface
; ***********************************************

  fld = FNET_TS(t,0:ncol_p-1)*cpair
  str = "Enthalpy flux  using TS"
  fmin = -30.
  fmax = 30.

;  res_p@cnLevelSelectionMode = "AutomaticLevels"    ; manually set the contour levels with the following 3 resources
  res_p@cnLevelSelectionMode = "ManualLevels"    ; manually set the contour levels with the following 3 resources
  res_p@cnMinLevelValF  = fmin                        ; set the minimum contour level
  res_p@cnMaxLevelValF  = fmax                   ; set the maximum contour level
  res_p@cnLevelSpacingF = (fmax-fmin)/50.0                  ; set the interval between contours

  min_f = min(fld(0:ncol_p-1))
  max_f = max(fld(0:ncol_p-1))
  res_p@lbTitleString    = "global min = "+sprintf("%6.4g",min_f)+"  global  max="+sprintf("%6.4g",max_f)
  print(str+" min/max: "+min_f+"  "+max_f)

  res_p@sfXArray            = lon1d
  res_p@sfYArray            = lat1d

  res_p@gsnLeftString        = str
  global_ave2 = sum(area_p(0:ncol_p-1)*fld(0:ncol_p-1))*inv_area_sphere_p
  print(str+" global ave: "+global_ave2)
  res_p@gsnRightString        = "mean: "+sprintf("%3.2g", global_ave2        )+" W/m^2"
  plot(5) = gsn_csm_contour_map_ce(wks,fld,res_p)




; ***********************************************
; difference
; ***********************************************

  delete(fld)
;  fld = cpT_2d(0:ncol_p-1)+PHIS_2d(0:ncol_p-1)-FNET_TS(t,0:ncol_p-1)*cpair
  fld = cpT_2d(0:ncol_p-1)-FNET_TS(t,0:ncol_p-1)*cpair
  str = "(d)+(e)-(f)"
  fmin = -2.
  fmax = 2.

;  res_p@cnLevelSelectionMode = "AutomaticLevels"    ; manually set the contour levels with the following 3 resources
  res_p@cnLevelSelectionMode = "ManualLevels"    ; manually set the contour levels with the following 3 resources
  res_p@cnMinLevelValF  = fmin                        ; set the minimum contour level
  res_p@cnMaxLevelValF  = fmax                   ; set the maximum contour level
  res_p@cnLevelSpacingF = (fmax-fmin)/50.0                  ; set the interval between contours

  min_f = min(fld(0:ncol_p-1))
  max_f = max(fld(0:ncol_p-1))
  res_p@lbTitleString    = "global min = "+sprintf("%6.4g",min_f)+"  global  max="+sprintf("%6.4g",max_f)
  print(str+" min/max: "+min_f+"  "+max_f)

  res_p@sfXArray            = lon1d
  res_p@sfYArray            = lat1d

  res_p@gsnLeftString        = str
  global_ave = sum(area_p(0:ncol_p-1)*fld(0:ncol_p-1))*inv_area_sphere_p
  print(str+" global ave: "+global_ave)

  res_p@gsnRightString        = "mean: "+sprintf("%3.2g", global_ave        )
  plot(6) = gsn_csm_contour_map_ce(wks,fld,res_p)
; ***********************************************
; 
; ***********************************************

  delete(fld)
  fld = TSatm(t,:)-TS(t,:)
  str = "T(lowest model level) - TS"
  fmin = -5.
  fmax = 5.

;  res_p@cnLevelSelectionMode = "AutomaticLevels"    ; manually set the contour levels with the following 3 resources
  res_p@cnLevelSelectionMode = "ManualLevels"    ; manually set the contour levels with the following 3 resources
  res_p@cnMinLevelValF  = fmin                        ; set the minimum contour level
  res_p@cnMaxLevelValF  = fmax                   ; set the maximum contour level
  res_p@cnLevelSpacingF = (fmax-fmin)/50.0                  ; set the interval between contours

  min_f = min(fld(0:ncol_p-1))
  max_f = max(fld(0:ncol_p-1))
  res_p@lbTitleString    = "global min = "+sprintf("%6.4g",min_f)+"  global  max="+sprintf("%6.4g",max_f)
  print(str+" min/max: "+min_f+"  "+max_f)

  res_p@sfXArray            = lon1d
  res_p@sfYArray            = lat1d

  res_p@gsnLeftString        = str
  global_ave = sum(area_p(0:ncol_p-1)*fld(0:ncol_p-1))*inv_area_sphere_p
  print(str+" global ave: "+global_ave)

  res_p@gsnRightString        = "mean: "+sprintf("%3.2g", global_ave        )
  plot(7) = gsn_csm_contour_map_ce(wks,fld,res_p)


; ***********************************************
; difference
; ***********************************************

  delete(fld)
;  fld = cpT_2d(0:ncol_p-1)-FNET_TBOT(t,0:ncol_p-1)*cpair
;  str = "(e)-(f)"
  fld = (FNET_TS(t,0:ncol_p-1)-FNET_TBOT(t,0:ncol_p-1))*cpair
  str = "Enthalpy flux using T(lowest model level) - (f)"
  fmin = -0.2
  fmax = 0.2

;  res_p@cnLevelSelectionMode = "AutomaticLevels"    ; manually set the contour levels with the following 3 resources
  res_p@cnLevelSelectionMode = "ManualLevels"    ; manually set the contour levels with the following 3 resources
  res_p@cnMinLevelValF  = fmin                        ; set the minimum contour level
  res_p@cnMaxLevelValF  = fmax                   ; set the maximum contour level
  res_p@cnLevelSpacingF = (fmax-fmin)/50.0                  ; set the interval between contours

  min_f = min(fld(0:ncol_p-1))
  max_f = max(fld(0:ncol_p-1))
  res_p@lbTitleString    = "global min = "+sprintf("%6.4g",min_f)+"  global  max="+sprintf("%6.4g",max_f)
  print(str+" min/max: "+min_f+"  "+max_f)

  res_p@sfXArray            = lon1d
  res_p@sfYArray            = lat1d

  res_p@gsnLeftString        = str
  global_ave = sum(area_p(0:ncol_p-1)*fld(0:ncol_p-1))*inv_area_sphere_p
  print(str+" global ave: "+global_ave)

  res_p@gsnRightString        = "mean: "+sprintf("%3.2g", global_ave        )
  plot(8) = gsn_csm_contour_map_ce(wks,fld,res_p)




  delete(resP)
  resP                  = True                   ; modify the panel plot
  resP@gsnFrame         = False                  ; don't advance panel plot
  resP@gsnPanelFigureStrings = (/"a","b","c","d","e","f","g","h","i","j"/)
  resP@gsnPanelYWhiteSpacePercent = 4.0
  resP@gsnPanelXWhiteSpacePercent = 4.0

  resP@txString   = "Assess missing terms / more accurate thermodynamics"
  gsn_panel(wks,plot,(/4,3/),resP)               ; now draw as one plot

  txres               = True
  frame(wks)

;
;
; NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT
; NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT
;
; 4444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444
; 4444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444
;
;
   delete(wks)
   delete(res)
   delete(plot)
   delete(res_p)



  filePlot = panelTitle+"4"
  wks = gsn_open_wks(output_format,filePlot)
  gsn_define_colormap(wks,"NCV_jaisnd");3gauss")      
    

  res                     = True         ; plot modifications desired
  res@gsnMaximize         = False         ; Maximize size of plot in frame
  res@gsnSpreadColors     = True         ; Use full colormap, but start
  res@gsnDraw             = False           ; don't draw
  res@gsnFrame            = False           ; don't advance frame
  res@cnFillOn            = True         ; Turn on contour fill
;  res@cnFillMode          = "AreaFill"   ; Style of fill. You can also
                                         ; use "CellFill" and "RasterFill"
  res@cnLinesOn           = False        ; Turn off contour lines
  res@cnLineLabelsOn      = False        ; Turn off contour line labels
  res@lbLabelAutoStride   = True         ; Clean up labelbar labels.
  res@cnFillMode="RasterFill"

  res@gsnStringFontHeightF = 0.02
  res@tiMainFontHeightF = 0.025
  res@lbLabelFontHeightF = 0.02
  res@tmXBLabelFontHeightF = 0.015
  res@tmYLLabelFontHeightF = 0.015
  res@pmLabelBarOrthogonalPosF = 0.2

  res@lbBoxLinesOn = False

  res_p=res

  plot = new(15,graphic) 


; ***********************************************
; d/dt(KE water flux)
; ***********************************************

  fld = K_2d(0:ncol_p-1)
  str = "d/dt(KE total H2O varying)"
  fmin = -0.01
  fmax =  0.01

  res_p@cnLevelSelectionMode = "ManualLevels"    ; manually set the contour levels with the following 3 resources
  res_p@cnMinLevelValF  = fmin                        ; set the minimum contour level
  res_p@cnMaxLevelValF  = fmax                   ; set the maximum contour level
  res_p@cnLevelSpacingF = (fmax-fmin)/50.0                  ; set the interval between contours

  min_f = min(fld(0:ncol_p-1))
  max_f = max(fld(0:ncol_p-1))
  res_p@lbTitleString    = "global min = "+sprintf("%6.4g",min_f)+"  global  max="+sprintf("%6.4g",max_f)
  print(str+" min/max: "+min_f+"  "+max_f)

  res_p@sfXArray            = lon1d
  res_p@sfYArray            = lat1d

  res_p@gsnLeftString        = str
  global_ave = sum(area_p(0:ncol_p-1)*fld(0:ncol_p-1))*inv_area_sphere_p
  print(str+" global ave: "+global_ave)
  res_p@gsnRightString        = "mean: "+sprintf("%3.2g", global_ave        )+" W/m^2"
  plot(0) = gsn_csm_contour_map_ce(wks,fld,res_p)

; ***********************************************
; d/dt(I+PHI water varying)
; ***********************************************

  fld = cpT_2d(0:ncol_p-1)
  str = "d/dt(enthalpy total H2O varying)"
  fmin = -30.
  fmax = 30.

  res_p@cnLevelSelectionMode = "ManualLevels"    ; manually set the contour levels with the following 3 resources
  res_p@cnMinLevelValF  = fmin                        ; set the minimum contour level
  res_p@cnMaxLevelValF  = fmax                   ; set the maximum contour level
  res_p@cnLevelSpacingF = (fmax-fmin)/50.0                  ; set the interval between contours

  min_f = min(fld(0:ncol_p-1))
  max_f = max(fld(0:ncol_p-1))
  res_p@lbTitleString    = "global min = "+sprintf("%6.4g",min_f)+"  global  max="+sprintf("%6.4g",max_f)
  print(str+" min/max: "+min_f+"  "+max_f)

  res_p@sfXArray            = lon1d
  res_p@sfYArray            = lat1d

  res_p@gsnLeftString        = str
  global_ave = sum(area_p(0:ncol_p-1)*fld(0:ncol_p-1))*inv_area_sphere_p
  print(str+" global ave: "+global_ave)
  res_p@gsnRightString        = "mean: "+sprintf("%3.2g", global_ave        )+" W/m^2"
  plot(1) = gsn_csm_contour_map_ce(wks,fld,res_p)


; ***********************************************
; update all H2O (dme_adjust)
; ***********************************************

  fld = PHIS_2d(0:ncol_p-1)
  str = "d/dt(PHIS term; H2O vayring)"
  fmin = -20.
  fmax = 20.

;  res_p@cnLevelSelectionMode = "AutomaticLevels"    ; manually set the contour levels with the following 3 resources
  res_p@cnLevelSelectionMode = "ManualLevels"    ; manually set the contour levels with the following 3 resources
  res_p@cnMinLevelValF  = fmin                        ; set the minimum contour level
  res_p@cnMaxLevelValF  = fmax                   ; set the maximum contour level
  res_p@cnLevelSpacingF = (fmax-fmin)/50.0                  ; set the interval between contours

  min_f = min(fld(0:ncol_p-1))
  max_f = max(fld(0:ncol_p-1))
  res_p@lbTitleString    = "global min = "+sprintf("%6.4g",min_f)+"  global  max="+sprintf("%6.4g",max_f)
  print(str+" min/max: "+min_f+"  "+max_f)

  res_p@sfXArray            = lon1d
  res_p@sfYArray            = lat1d

  res_p@gsnLeftString        = str
  global_ave = sum(area_p(0:ncol_p-1)*fld(0:ncol_p-1))*inv_area_sphere_p
  print(str+" global ave: "+global_ave)
  res_p@gsnRightString        = "mean: "+sprintf("%3.2g", global_ave        )+" W/m^2"
  plot(2) = gsn_csm_contour_map_ce(wks,fld,res_p)

  delete(fld)
; ***********************************************
; Missing enthalpy flux 
; ***********************************************

  fld = FNET_TS(t,0:ncol_p-1)*cpair
  str = "Enthalpy flux  (using TS)"
  fmin = -30.
  fmax = 30.

;  res_p@cnLevelSelectionMode = "AutomaticLevels"    ; manually set the contour levels with the following 3 resources
  res_p@cnLevelSelectionMode = "ManualLevels"    ; manually set the contour levels with the following 3 resources
  res_p@cnMinLevelValF  = fmin                        ; set the minimum contour level
  res_p@cnMaxLevelValF  = fmax                   ; set the maximum contour level
  res_p@cnLevelSpacingF = (fmax-fmin)/50.0                  ; set the interval between contours

  min_f = min(fld(0:ncol_p-1))
  max_f = max(fld(0:ncol_p-1))
  res_p@lbTitleString    = "global min = "+sprintf("%6.4g",min_f)+"  global  max="+sprintf("%6.4g",max_f)
  print(str+" min/max: "+min_f+"  "+max_f)

  res_p@sfXArray            = lon1d
  res_p@sfYArray            = lat1d

  res_p@gsnLeftString        = str
  global_ave2 = sum(area_p(0:ncol_p-1)*fld(0:ncol_p-1))*inv_area_sphere_p
  print(str+" global ave: "+global_ave2)
  res_p@gsnRightString        = "mean: "+sprintf("%3.2g", global_ave2        )+" W/m^2"
  plot(3) = gsn_csm_contour_map_ce(wks,fld,res_p)


; ***********************************************
; Missing enthalpy flux 
; ***********************************************

  fld = FNET_TBOT(t,0:ncol_p-1)*cpair
  str = "Enthalpy flux  xxx"
  fmin = -30.
  fmax = 30.

;  res_p@cnLevelSelectionMode = "AutomaticLevels"    ; manually set the contour levels with the following 3 resources
  res_p@cnLevelSelectionMode = "ManualLevels"    ; manually set the contour levels with the following 3 resources
  res_p@cnMinLevelValF  = fmin                        ; set the minimum contour level
  res_p@cnMaxLevelValF  = fmax                   ; set the maximum contour level
  res_p@cnLevelSpacingF = (fmax-fmin)/50.0                  ; set the interval between contours

  min_f = min(fld(0:ncol_p-1))
  max_f = max(fld(0:ncol_p-1))
  res_p@lbTitleString    = "global min = "+sprintf("%6.4g",min_f)+"  global  max="+sprintf("%6.4g",max_f)
  print(str+" min/max: "+min_f+"  "+max_f)

  res_p@sfXArray            = lon1d
  res_p@sfYArray            = lat1d

  res_p@gsnLeftString        = str
  global_ave2 = sum(area_p(0:ncol_p-1)*fld(0:ncol_p-1))*inv_area_sphere_p
  print(str+" global ave: "+global_ave2)
  res_p@gsnRightString        = "mean: "+sprintf("%3.2g", global_ave2        )+" W/m^2"
  plot(4) = gsn_csm_contour_map_ce(wks,fld,res_p)


; ***********************************************
; difference
; ***********************************************

  delete(fld)
  fld = cpT_2d(0:ncol_p-1)-FNET_TS(t,0:ncol_p-1)*cpair
  str = "(d)-(e)"
  fmin = -2.
  fmax = 2.

;  res_p@cnLevelSelectionMode = "AutomaticLevels"    ; manually set the contour levels with the following 3 resources
  res_p@cnLevelSelectionMode = "ManualLevels"    ; manually set the contour levels with the following 3 resources
  res_p@cnMinLevelValF  = fmin                        ; set the minimum contour level
  res_p@cnMaxLevelValF  = fmax                   ; set the maximum contour level
  res_p@cnLevelSpacingF = (fmax-fmin)/50.0                  ; set the interval between contours

  min_f = min(fld(0:ncol_p-1))
  max_f = max(fld(0:ncol_p-1))
  res_p@lbTitleString    = "global min = "+sprintf("%6.4g",min_f)+"  global  max="+sprintf("%6.4g",max_f)
  print(str+" min/max: "+min_f+"  "+max_f)

  res_p@sfXArray            = lon1d
  res_p@sfYArray            = lat1d

  res_p@gsnLeftString        = str
  global_ave = sum(area_p(0:ncol_p-1)*fld(0:ncol_p-1))*inv_area_sphere_p
  print(str+" global ave: "+global_ave)

  res_p@gsnRightString        = "mean: "+sprintf("%3.2g", global_ave        )
  plot(5) = gsn_csm_contour_map_ce(wks,fld,res_p)


; ***********************************************
; difference
; ***********************************************

  delete(fld)
  fld = cpT_2d(0:ncol_p-1)-FNET_TS(t,0:ncol_p-1)*cpair
  str = "(b)-(d)"
  fmin = -0.7
  fmax = 0.7

;  res_p@cnLevelSelectionMode = "AutomaticLevels"    ; manually set the contour levels with the following 3 resources
  res_p@cnLevelSelectionMode = "ManualLevels"    ; manually set the contour levels with the following 3 resources
  res_p@cnMinLevelValF  = fmin                        ; set the minimum contour level
  res_p@cnMaxLevelValF  = fmax                   ; set the maximum contour level
  res_p@cnLevelSpacingF = (fmax-fmin)/50.0                  ; set the interval between contours

  min_f = min(fld(0:ncol_p-1))
  max_f = max(fld(0:ncol_p-1))
  res_p@lbTitleString    = "global min = "+sprintf("%6.4g",min_f)+"  global  max="+sprintf("%6.4g",max_f)
  print(str+" min/max: "+min_f+"  "+max_f)

  res_p@sfXArray            = lon1d
  res_p@sfYArray            = lat1d

  res_p@gsnLeftString        = str
  global_ave = sum(area_p(0:ncol_p-1)*fld(0:ncol_p-1))*inv_area_sphere_p
  print(str+" global ave: "+global_ave)

  res_p@gsnRightString        = "mean: "+sprintf("%3.2g", global_ave        )
  plot(6) = gsn_csm_contour_map_ce(wks,fld,res_p)
; ***********************************************
; difference
; ***********************************************

  delete(fld)
  fld = cpT_2d(0:ncol_p-1)-FNET_TBOT(t,0:ncol_p-1)*cpair
  str = "(b)-(e)"
  fmin = -0.7
  fmax = 0.7

;  res_p@cnLevelSelectionMode = "AutomaticLevels"    ; manually set the contour levels with the following 3 resources
  res_p@cnLevelSelectionMode = "ManualLevels"    ; manually set the contour levels with the following 3 resources
  res_p@cnMinLevelValF  = fmin                        ; set the minimum contour level
  res_p@cnMaxLevelValF  = fmax                   ; set the maximum contour level
  res_p@cnLevelSpacingF = (fmax-fmin)/50.0                  ; set the interval between contours

  min_f = min(fld(0:ncol_p-1))
  max_f = max(fld(0:ncol_p-1))
  res_p@lbTitleString    = "global min = "+sprintf("%6.4g",min_f)+"  global  max="+sprintf("%6.4g",max_f)
  print(str+" min/max: "+min_f+"  "+max_f)

  res_p@sfXArray            = lon1d
  res_p@sfYArray            = lat1d

  res_p@gsnLeftString        = str
  global_ave = sum(area_p(0:ncol_p-1)*fld(0:ncol_p-1))*inv_area_sphere_p
  print(str+" global ave: "+global_ave)

  res_p@gsnRightString        = "mean: "+sprintf("%3.2g", global_ave        )
  plot(7) = gsn_csm_contour_map_ce(wks,fld,res_p)

; ***********************************************
; difference
; ***********************************************

  delete(fld)
;  fld = cpT_2d(0:ncol_p-1)-FNET_TBOT(t,0:ncol_p-1)*cpair
;  str = "(e)-(f)"
  fld = (FNET_TS(t,0:ncol_p-1)-FNET_TBOT(t,0:ncol_p-1))*cpair
  str = "FNET_TS-FNET_TBOT"
  fmin = -0.7
  fmax = 0.7

;  res_p@cnLevelSelectionMode = "AutomaticLevels"    ; manually set the contour levels with the following 3 resources
  res_p@cnLevelSelectionMode = "ManualLevels"    ; manually set the contour levels with the following 3 resources
  res_p@cnMinLevelValF  = fmin                        ; set the minimum contour level
  res_p@cnMaxLevelValF  = fmax                   ; set the maximum contour level
  res_p@cnLevelSpacingF = (fmax-fmin)/50.0                  ; set the interval between contours

  min_f = min(fld(0:ncol_p-1))
  max_f = max(fld(0:ncol_p-1))
  res_p@lbTitleString    = "global min = "+sprintf("%6.4g",min_f)+"  global  max="+sprintf("%6.4g",max_f)
  print(str+" min/max: "+min_f+"  "+max_f)

  res_p@sfXArray            = lon1d
  res_p@sfYArray            = lat1d

  res_p@gsnLeftString        = str
  global_ave = sum(area_p(0:ncol_p-1)*fld(0:ncol_p-1))*inv_area_sphere_p
  print(str+" global ave: "+global_ave)

  res_p@gsnRightString        = "mean: "+sprintf("%3.2g", global_ave        )
  plot(8) = gsn_csm_contour_map_ce(wks,fld,res_p)


; ***********************************************
; 
; ***********************************************

  delete(fld)
  fld = TS(t,:)-TSatm(t,:)
  str = "TS (skin)-T(k=nlev) atm "
  fmin = -4.
  fmax = 4.

;  res_p@cnLevelSelectionMode = "AutomaticLevels"    ; manually set the contour levels with the following 3 resources
  res_p@cnLevelSelectionMode = "ManualLevels"    ; manually set the contour levels with the following 3 resources
  res_p@cnMinLevelValF  = fmin                        ; set the minimum contour level
  res_p@cnMaxLevelValF  = fmax                   ; set the maximum contour level
  res_p@cnLevelSpacingF = (fmax-fmin)/50.0                  ; set the interval between contours

  min_f = min(fld(0:ncol_p-1))
  max_f = max(fld(0:ncol_p-1))
  res_p@lbTitleString    = "global min = "+sprintf("%6.4g",min_f)+"  global  max="+sprintf("%6.4g",max_f)
  print(str+" min/max: "+min_f+"  "+max_f)

  res_p@sfXArray            = lon1d
  res_p@sfYArray            = lat1d

  res_p@gsnLeftString        = str
  global_ave = sum(area_p(0:ncol_p-1)*fld(0:ncol_p-1))*inv_area_sphere_p
  print(str+" global ave: "+global_ave)

  res_p@gsnRightString        = "mean: "+sprintf("%3.2g", global_ave        )
  plot(9) = gsn_csm_contour_map_ce(wks,fld,res_p)


; ***********************************************
; 
; ***********************************************

  delete(fld)
  fld = TS(t,:); -SST(t,:)
  str = "TS (= SST over ocea)"
  fmin = -10.
  fmax = 10.

  res_p@cnLevelSelectionMode = "AutomaticLevels"    ; manually set the contour levels with the following 3 resources
;  res_p@cnLevelSelectionMode = "ManualLevels"    ; manually set the contour levels with the following 3 resources
;  res_p@cnMinLevelValF  = fmin                        ; set the minimum contour level
;  res_p@cnMaxLevelValF  = fmax                   ; set the maximum contour level
;  res_p@cnLevelSpacingF = (fmax-fmin)/50.0                  ; set the interval between contours

  min_f = min(fld(0:ncol_p-1))
  max_f = max(fld(0:ncol_p-1))
  res_p@lbTitleString    = "global min = "+sprintf("%6.4g",min_f)+"  global  max="+sprintf("%6.4g",max_f)
  print(str+" min/max: "+min_f+"  "+max_f)

  res_p@sfXArray            = lon1d
  res_p@sfYArray            = lat1d

  res_p@gsnLeftString        = str
  global_ave = sum(area_p(0:ncol_p-1)*fld(0:ncol_p-1))*inv_area_sphere_p
  print(str+" global ave: "+global_ave)

  res_p@gsnRightString        = "mean: "+sprintf("%3.2g", global_ave        )
  plot(10) = gsn_csm_contour_map_ce(wks,fld,res_p)



; ***********************************************
; Water flux
; ***********************************************

  fld = FH2O(t,0:ncol_p-1)
  str = "Flux of total water"
  fmin = -0.0001
  fmax = 0.0001

  res_p@cnLevelSelectionMode = "ManualLevels"    ; manually set the contour levels with the following 3 resources
  res_p@cnMinLevelValF  = fmin                        ; set the minimum contour level
  res_p@cnMaxLevelValF  = fmax                   ; set the maximum contour level
  res_p@cnLevelSpacingF = (fmax-fmin)/50.0                  ; set the interval between contours

  min_f = min(fld(0:ncol_p-1))
  max_f = max(fld(0:ncol_p-1))
  res_p@lbTitleString    = "global min = "+sprintf("%6.4g",min_f)+"  global  max="+sprintf("%6.4g",max_f)
  print(str+" min/max: "+min_f+"  "+max_f)

  res_p@sfXArray            = lon1d
  res_p@sfYArray            = lat1d

  res_p@gsnLeftString        = str
  global_ave = sum(area_p(0:ncol_p-1)*fld(0:ncol_p-1))*inv_area_sphere_p
  print(str+" global ave: "+global_ave)
  res_p@gsnRightString        = "mean: "+sprintf("%3.2g", global_ave        )
  plot(11) = gsn_csm_contour_map_ce(wks,fld,res_p)


; ***********************************************
; 
; ***********************************************

  delete(fld)
  fld = cpT_2d(0:ncol_p-1)
  do i=0,ncol_p-1
;    if (abs(FH2O(t,i)).LE.1D-5) then
    if (abs(FH2O(t,i)).LE.1D-6) then
      fld(i) = 1000.0D0; TS(t,i)
    else
      fld(i) = cpT_2d(i)/(cpair*FH2O(t,i))
    end if
  end do
  testimate = fld
  str = "Estimate of TS to balance I+PHIS term"
  fmin = 220.
  fmax = 310.

;  res_p@cnLevelSelectionMode = "AutomaticLevels"    ; manually set the contour levels with the following 3 resources
  res_p@cnLevelSelectionMode = "ManualLevels"    ; manually set the contour levels with the following 3 resources
  res_p@cnMinLevelValF  = fmin                        ; set the minimum contour level
  res_p@cnMaxLevelValF  = fmax                   ; set the maximum contour level
  res_p@cnLevelSpacingF = (fmax-fmin)/50.0                  ; set the interval between contours

  min_f = min(fld(0:ncol_p-1))
  max_f = max(fld(0:ncol_p-1))
  res_p@lbTitleString    = "global min = "+sprintf("%6.4g",min_f)+"  global  max="+sprintf("%6.4g",max_f)
  print(str+" min/max: "+min_f+"  "+max_f)

  res_p@sfXArray            = lon1d
  res_p@sfYArray            = lat1d

  res_p@gsnLeftString        = str
  global_ave = sum(area_p(0:ncol_p-1)*fld(0:ncol_p-1))*inv_area_sphere_p
  print(str+" global ave: "+global_ave)

  res_p@gsnRightString        = "mean: "+sprintf("%3.2g", global_ave        )
  plot(12) = gsn_csm_contour_map_ce(wks,fld,res_p)

; ***********************************************
; 
; ***********************************************

  delete(fld)
  fld = testimate(0:ncol_p-1)-TS(t,0:ncol_p-1)
  str = "Estimate of TS to balance cp*T term minus TS"
  fmin = -30.
  fmax = 30.

;  res_p@cnLevelSelectionMode = "AutomaticLevels"    ; manually set the contour levels with the following 3 resources
  res_p@cnLevelSelectionMode = "ManualLevels"    ; manually set the contour levels with the following 3 resources
  res_p@cnMinLevelValF  = fmin                        ; set the minimum contour level
  res_p@cnMaxLevelValF  = fmax                   ; set the maximum contour level
  res_p@cnLevelSpacingF = (fmax-fmin)/50.0                  ; set the interval between contours

  min_f = min(fld(0:ncol_p-1))
  max_f = max(fld(0:ncol_p-1))
  res_p@lbTitleString    = "global min = "+sprintf("%6.4g",min_f)+"  global  max="+sprintf("%6.4g",max_f)
  print(str+" min/max: "+min_f+"  "+max_f)

  res_p@sfXArray            = lon1d
  res_p@sfYArray            = lat1d

  res_p@gsnLeftString        = str
  global_ave = sum(area_p(0:ncol_p-1)*fld(0:ncol_p-1))*inv_area_sphere_p
  print(str+" global ave: "+global_ave)

  res_p@gsnRightString        = "mean: "+sprintf("%3.2g", global_ave        )
  plot(13) = gsn_csm_contour_map_ce(wks,fld,res_p)

; ***********************************************
; 
; ***********************************************

  delete(fld)

  fld = PRECT(t,:)
  fld = 1000*24*3600*fld ;convert from m/s to mm/day
  str = "PRECT (mm/day)"
  fmin = -10.
  fmax = 10.

;  res_p@cnLevelSelectionMode = "AutomaticLevels"    ; manually set the contour levels with the following 3 resources
  res_p@cnLevelSelectionMode = "ManualLevels"    ; manually set the contour levels with the following 3 resources
  res_p@cnMinLevelValF  = fmin                        ; set the minimum contour level
  res_p@cnMaxLevelValF  = fmax                   ; set the maximum contour level
  res_p@cnLevelSpacingF = (fmax-fmin)/50.0                  ; set the interval between contours

  min_f = min(fld(0:ncol_p-1))
  max_f = max(fld(0:ncol_p-1))
  res_p@lbTitleString    = "global min = "+sprintf("%6.4g",min_f)+"  global  max="+sprintf("%6.4g",max_f)
  print(str+" min/max: "+min_f+"  "+max_f)

  res_p@sfXArray            = lon1d
  res_p@sfYArray            = lat1d

  res_p@gsnLeftString        = str
  global_ave = sum(area_p(0:ncol_p-1)*fld(0:ncol_p-1))*inv_area_sphere_p
  print(str+" global ave: "+global_ave)

  res_p@gsnRightString        = "mean: "+sprintf("%3.2g", global_ave        )
  plot(14) = gsn_csm_contour_map_ce(wks,fld,res_p)

 



 
  delete(resP)
  resP                  = True                   ; modify the panel plot
  resP@gsnFrame         = False                  ; don't advance panel plot
  resP@gsnPanelFigureStrings = (/"a","b","c","d","e","f","g","h","i","j"/)
  resP@gsnPanelYWhiteSpacePercent = 4.0
  resP@gsnPanelXWhiteSpacePercent = 4.0

  resP@txString   = "Assess missing terms"
  gsn_panel(wks,plot,(/5,3/),resP)               ; now draw as one plot

  txres               = True
  frame(wks)

;
;
; NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT
; NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT-NEWPLOT
;
; 4444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444
; 4444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444
;
;
   delete(wks)
   delete(res)
   delete(plot)
   delete(res_p)



  filePlot = panelTitle+"KE"
  wks = gsn_open_wks(output_format,filePlot)
  gsn_define_colormap(wks,"NCV_jaisnd");3gauss")      
    

  res                     = True         ; plot modifications desired
  res@gsnMaximize         = False         ; Maximize size of plot in frame
  res@gsnSpreadColors     = True         ; Use full colormap, but start
  res@gsnDraw             = False           ; don't draw
  res@gsnFrame            = False           ; don't advance frame
  res@cnFillOn            = True         ; Turn on contour fill
;  res@cnFillMode          = "AreaFill"   ; Style of fill. You can also
                                         ; use "CellFill" and "RasterFill"
  res@cnLinesOn           = False        ; Turn off contour lines
  res@cnLineLabelsOn      = False        ; Turn off contour line labels
  res@lbLabelAutoStride   = True         ; Clean up labelbar labels.
  res@cnFillMode="RasterFill"

  res@gsnStringFontHeightF = 0.02
  res@tiMainFontHeightF = 0.025
  res@lbLabelFontHeightF = 0.02
  res@tmXBLabelFontHeightF = 0.015
  res@tmYLLabelFontHeightF = 0.015
  res@pmLabelBarOrthogonalPosF = 0.2

  res@lbBoxLinesOn = False

  res_p=res

  plot = new(4,graphic) 


; ***********************************************
; d/dt(KE water flux)
; ***********************************************

  fld = K_2d(0:ncol_p-1)
  str = "d/dt(KE total H2O varying)"
  fmin = -0.01
  fmax =  0.01

  res_p@cnLevelSelectionMode = "ManualLevels"    ; manually set the contour levels with the following 3 resources
  res_p@cnMinLevelValF  = fmin                        ; set the minimum contour level
  res_p@cnMaxLevelValF  = fmax                   ; set the maximum contour level
  res_p@cnLevelSpacingF = (fmax-fmin)/50.0                  ; set the interval between contours

  min_f = min(fld(0:ncol_p-1))
  max_f = max(fld(0:ncol_p-1))
  res_p@lbTitleString    = "global min = "+sprintf("%6.4g",min_f)+"  global  max="+sprintf("%6.4g",max_f)
  print(str+" min/max: "+min_f+"  "+max_f)

  res_p@sfXArray            = lon1d
  res_p@sfYArray            = lat1d

  res_p@gsnLeftString        = str
  global_ave = sum(area_p(0:ncol_p-1)*fld(0:ncol_p-1))*inv_area_sphere_p
  print(str+" global ave: "+global_ave)
  res_p@gsnRightString        = "mean: "+sprintf("%3.2g", global_ave        )+" W/m^2"
  plot(0) = gsn_csm_contour_map_ce(wks,fld,res_p)

; ***********************************************
; Missing kinetic energy flux 
; ***********************************************

  fld = FKE(t,0:ncol_p-1)
  str = "FH2O*"+Ktilde+"=FH2O*K(nlev)"
  fmin = -0.01
  fmax =  0.01


;  res_p@cnLevelSelectionMode = "AutomaticLevels"    ; manually set the contour levels with the following 3 resources
  res_p@cnLevelSelectionMode = "ManualLevels"    ; manually set the contour levels with the following 3 resources
  res_p@cnMinLevelValF  = fmin                        ; set the minimum contour level
  res_p@cnMaxLevelValF  = fmax                   ; set the maximum contour level
  res_p@cnLevelSpacingF = (fmax-fmin)/50.0                  ; set the interval between contours

  min_f = min(fld(0:ncol_p-1))
  max_f = max(fld(0:ncol_p-1))
  res_p@lbTitleString    = "global min = "+sprintf("%6.4g",min_f)+"  global  max="+sprintf("%6.4g",max_f)
  print(str+" min/max: "+min_f+"  "+max_f)

  res_p@sfXArray            = lon1d
  res_p@sfYArray            = lat1d

  res_p@gsnLeftString        = str
  global_ave2 = sum(area_p(0:ncol_p-1)*fld(0:ncol_p-1))*inv_area_sphere_p
  print(str+" global ave: "+global_ave2)
  res_p@gsnRightString        = "mean: "+sprintf("%3.2g", global_ave2        )+" W/m^2"
  plot(1) = gsn_csm_contour_map_ce(wks,fld,res_p)


; ***********************************************
; Missing kinetic energy flux 
; ***********************************************

  fld = K_2d(0:ncol_p-1)-FKE(t,0:ncol_p-1)
  str = "(a)-(b)"
  fmin = -0.01
  fmax =  0.01


;  res_p@cnLevelSelectionMode = "AutomaticLevels"    ; manually set the contour levels with the following 3 resources
  res_p@cnLevelSelectionMode = "ManualLevels"    ; manually set the contour levels with the following 3 resources
  res_p@cnMinLevelValF  = fmin                        ; set the minimum contour level
  res_p@cnMaxLevelValF  = fmax                   ; set the maximum contour level
  res_p@cnLevelSpacingF = (fmax-fmin)/50.0                  ; set the interval between contours

  min_f = min(fld(0:ncol_p-1))
  max_f = max(fld(0:ncol_p-1))
  res_p@lbTitleString    = "global min = "+sprintf("%6.4g",min_f)+"  global  max="+sprintf("%6.4g",max_f)
  print(str+" min/max: "+min_f+"  "+max_f)

  res_p@sfXArray            = lon1d
  res_p@sfYArray            = lat1d

  res_p@gsnLeftString        = str
  global_ave2 = sum(area_p(0:ncol_p-1)*fld(0:ncol_p-1))*inv_area_sphere_p
  print(str+" global ave: "+global_ave2)
  res_p@gsnRightString        = "mean: "+sprintf("%3.2g", global_ave2        )+" W/m^2"
  plot(2) = gsn_csm_contour_map_ce(wks,fld,res_p)


; ***********************************************
; wind stress
; ***********************************************

  delete(fld)
  fld = FTAU(t,0:ncol_p-1)
  str = "Wind stress term"
  fmin = -4
  fmax = 4

;  res_p@cnLevelSelectionMode = "AutomaticLevels"    ; manually set the contour levels with the following 3 resources
  res_p@cnLevelSelectionMode = "ManualLevels"    ; manually set the contour levels with the following 3 resources
  res_p@cnMinLevelValF  = fmin                        ; set the minimum contour level
  res_p@cnMaxLevelValF  = fmax                   ; set the maximum contour level
  res_p@cnLevelSpacingF = (fmax-fmin)/50.0                  ; set the interval between contours

  min_f = min(fld(0:ncol_p-1))
  max_f = max(fld(0:ncol_p-1))
  res_p@lbTitleString    = "global min = "+sprintf("%6.4g",min_f)+"  global  max="+sprintf("%6.4g",max_f)
  print(str+" min/max: "+min_f+"  "+max_f)

  res_p@sfXArray            = lon1d
  res_p@sfYArray            = lat1d

  res_p@gsnLeftString        = str
  global_ave = sum(area_p(0:ncol_p-1)*fld(0:ncol_p-1))*inv_area_sphere_p
  print(str+" global ave: "+global_ave)

  res_p@gsnRightString        = "mean: "+sprintf("%3.2g", global_ave        )+" W/m^2"
  plot(3) = gsn_csm_contour_map_ce(wks,fld,res_p)

  delete(resP)
  resP                  = True                   ; modify the panel plot
  resP@gsnFrame         = False                  ; don't advance panel plot
  resP@gsnPanelFigureStrings = (/"a","b","c","d"/)
  resP@gsnPanelYWhiteSpacePercent = 4.0
  resP@gsnPanelXWhiteSpacePercent = 4.0

  resP@txString   = "Kinetic energy terms"
  gsn_panel(wks,plot,(/2,2/),resP)               ; now draw as one plot

  txres               = True
  frame(wks)


   delete(wks)
   delete(res)
   delete(plot)
   delete(res_p)

;EnthalpyEnthalpyEnthalpyEnthalpyEnthalpyEnthalpyEnthalpyEnthalpyEnthalpyEnthalpyEnthalpyEnthalpyEnthalpyEnthalpyEnthalpyEnthalpyEnthalpyEnthalpy
;EnthalpyEnthalpyEnthalpyEnthalpyEnthalpyEnthalpyEnthalpyEnthalpyEnthalpyEnthalpyEnthalpyEnthalpyEnthalpyEnthalpyEnthalpyEnthalpyEnthalpyEnthalpy
;EnthalpyEnthalpyEnthalpyEnthalpyEnthalpyEnthalpyEnthalpyEnthalpyEnthalpyEnthalpyEnthalpyEnthalpyEnthalpyEnthalpyEnthalpyEnthalpyEnthalpyEnthalpy

  filePlot = panelTitle+"H"
  wks = gsn_open_wks(output_format,filePlot)
  gsn_define_colormap(wks,"NCV_jaisnd");3gauss")      
    

  res                     = True         ; plot modifications desired
  res@gsnMaximize         = False         ; Maximize size of plot in frame
  res@gsnSpreadColors     = True         ; Use full colormap, but start
  res@gsnDraw             = False           ; don't draw
  res@gsnFrame            = False           ; don't advance frame
  res@cnFillOn            = True         ; Turn on contour fill
;  res@cnFillMode          = "AreaFill"   ; Style of fill. You can also
                                         ; use "CellFill" and "RasterFill"
  res@cnLinesOn           = False        ; Turn off contour lines
  res@cnLineLabelsOn      = False        ; Turn off contour line labels
  res@lbLabelAutoStride   = True         ; Clean up labelbar labels.
  res@cnFillMode="RasterFill"

  res@gsnStringFontHeightF = 0.02
  res@tiMainFontHeightF = 0.025
  res@lbLabelFontHeightF = 0.02
  res@tmXBLabelFontHeightF = 0.015
  res@tmYLLabelFontHeightF = 0.015
  res@pmLabelBarOrthogonalPosF = 0.2

  res@lbBoxLinesOn = False

  res_p=res

  plot = new(6,graphic) 

; ***********************************************
; d/dt(enthalpy water varying)
; ***********************************************

  fld = cpT_2d(0:ncol_p-1)
  str = "d/dt(enthalpy total H2O varying)"
  fmin = -30.
  fmax = 30.

  res_p@cnLevelSelectionMode = "ManualLevels"    ; manually set the contour levels with the following 3 resources
  res_p@cnMinLevelValF  = fmin                        ; set the minimum contour level
  res_p@cnMaxLevelValF  = fmax                   ; set the maximum contour level
  res_p@cnLevelSpacingF = (fmax-fmin)/50.0                  ; set the interval between contours

  min_f = min(fld(0:ncol_p-1))
  max_f = max(fld(0:ncol_p-1))
  res_p@lbTitleString    = "global min = "+sprintf("%6.4g",min_f)+"  global  max="+sprintf("%6.4g",max_f)
  print(str+" min/max: "+min_f+"  "+max_f)

  res_p@sfXArray            = lon1d
  res_p@sfYArray            = lat1d

  res_p@gsnLeftString        = str
  global_ave = sum(area_p(0:ncol_p-1)*fld(0:ncol_p-1))*inv_area_sphere_p
  print(str+" global ave: "+global_ave)
  res_p@gsnRightString        = "mean: "+sprintf("%3.2g", global_ave        )+" W/m^2"
  plot(0) = gsn_csm_contour_map_ce(wks,fld,res_p)

; ***********************************************
; Missing enthalpy flux 
; ***********************************************

  fld = FNET_TBOT(t,0:ncol_p-1)*cpair
  str = "Enthalpy flux using "+Ttilde+"=T~B~a,surf~N~=T(nlev)"
  fmin = -30.
  fmax = 30.

;  res_p@cnLevelSelectionMode = "AutomaticLevels"    ; manually set the contour levels with the following 3 resources
  res_p@cnLevelSelectionMode = "ManualLevels"    ; manually set the contour levels with the following 3 resources
  res_p@cnMinLevelValF  = fmin                        ; set the minimum contour level
  res_p@cnMaxLevelValF  = fmax                   ; set the maximum contour level
  res_p@cnLevelSpacingF = (fmax-fmin)/50.0                  ; set the interval between contours

  min_f = min(fld(0:ncol_p-1))
  max_f = max(fld(0:ncol_p-1))
  res_p@lbTitleString    = "global min = "+sprintf("%6.4g",min_f)+"  global  max="+sprintf("%6.4g",max_f)
  print(str+" min/max: "+min_f+"  "+max_f)

  res_p@sfXArray            = lon1d
  res_p@sfYArray            = lat1d

  res_p@gsnLeftString        = str
  global_ave2 = sum(area_p(0:ncol_p-1)*fld(0:ncol_p-1))*inv_area_sphere_p
  print(str+" global ave: "+global_ave2)
  res_p@gsnRightString        = "mean: "+sprintf("%3.2g", global_ave2        )+" W/m^2"
  plot(1) = gsn_csm_contour_map_ce(wks,fld,res_p)

; ***********************************************
; Missing enthalpy flux 
; ***********************************************

  fld = cpT_2d(0:ncol_p-1)-FNET_TBOT(t,0:ncol_p-1)*cpair
  str = "(d) - (e)"
  fmin = -1.
  fmax = 1.

;  res_p@cnLevelSelectionMode = "AutomaticLevels"    ; manually set the contour levels with the following 3 resources
  res_p@cnLevelSelectionMode = "ManualLevels"    ; manually set the contour levels with the following 3 resources
  res_p@cnMinLevelValF  = fmin                        ; set the minimum contour level
  res_p@cnMaxLevelValF  = fmax                   ; set the maximum contour level
  res_p@cnLevelSpacingF = (fmax-fmin)/50.0                  ; set the interval between contours

  min_f = min(fld(0:ncol_p-1))
  max_f = max(fld(0:ncol_p-1))
  res_p@lbTitleString    = "global min = "+sprintf("%6.4g",min_f)+"  global  max="+sprintf("%6.4g",max_f)
  print(str+" min/max: "+min_f+"  "+max_f)

  res_p@sfXArray            = lon1d
  res_p@sfYArray            = lat1d

  res_p@gsnLeftString        = str
  global_ave2 = sum(area_p(0:ncol_p-1)*fld(0:ncol_p-1))*inv_area_sphere_p
  print(str+" global ave: "+global_ave2)
  res_p@gsnRightString        = "mean: "+sprintf("%3.2g", global_ave2        )+" W/m^2"
  plot(2) = gsn_csm_contour_map_ce(wks,fld,res_p)

; ***********************************************
; equivalent T
; ***********************************************

  fld = cpT_2d(0:ncol_p-1)
  do i=0,ncol_p-1
    if (abs(FH2O(t,i)).LE.1D-6) then
      fld(i) = 0.0
    else
      fld(i) = (cpT_2d(i)/(cpair*FH2O(t,i)))-TSatm(t,i)
    end if
  end do
  str = Ttilde+" in (b) to balance (a) minus T(nlev)"
  fmin = -20.
  fmax = 20.


;  res_p@cnLevelSelectionMode = "AutomaticLevels"    ; manually set the contour levels with the following 3 resources
  res_p@cnLevelSelectionMode = "ManualLevels"    ; manually set the contour levels with the following 3 resources
  res_p@cnMinLevelValF  = fmin                        ; set the minimum contour level
  res_p@cnMaxLevelValF  = fmax                   ; set the maximum contour level
  res_p@cnLevelSpacingF = (fmax-fmin)/50.0                  ; set the interval between contours

  min_f = min(fld(0:ncol_p-1))
  max_f = max(fld(0:ncol_p-1))
  res_p@lbTitleString    = "global min = "+sprintf("%6.4g",min_f)+"  global  max="+sprintf("%6.4g",max_f)
  print(str+" min/max: "+min_f+"  "+max_f)

  res_p@sfXArray            = lon1d
  res_p@sfYArray            = lat1d

  res_p@gsnLeftString        = str
  global_ave2 = sum(area_p(0:ncol_p-1)*fld(0:ncol_p-1))*inv_area_sphere_p
  print(str+" global ave: "+global_ave2)
  res_p@gsnRightString        = "mean: "+sprintf("%3.2g", global_ave2        )+" K"
  plot(3) = gsn_csm_contour_map_ce(wks,fld,res_p)



; ***********************************************
; Fh2O
; ***********************************************

  delete(fld)
  fld = FH2O(t,0:ncol_p-1)
  str = "FH2O"
  fmin = -0.0001
  fmax = 0.0001

;  res_p@cnLevelSelectionMode = "AutomaticLevels"    ; manually set the contour levels with the following 3 resources
  res_p@cnLevelSelectionMode = "ManualLevels"    ; manually set the contour levels with the following 3 resources
  res_p@cnMinLevelValF  = fmin                        ; set the minimum contour level
  res_p@cnMaxLevelValF  = fmax                   ; set the maximum contour level
  res_p@cnLevelSpacingF = (fmax-fmin)/50.0                  ; set the interval between contours

  min_f = min(fld(0:ncol_p-1))
  max_f = max(fld(0:ncol_p-1))
  res_p@lbTitleString    = "global min = "+sprintf("%6.4g",min_f)+"  global  max="+sprintf("%6.4g",max_f)
  print(str+" min/max: "+min_f+"  "+max_f)

  res_p@sfXArray            = lon1d
  res_p@sfYArray            = lat1d

  res_p@gsnLeftString        = str
  global_ave2 = sum(area_p(0:ncol_p-1)*fld(0:ncol_p-1))*inv_area_sphere_p
  print(str+" global ave: "+global_ave2)
  res_p@gsnRightString        = "mean: "+sprintf("%3.2g", global_ave2        )+" kg/(s m^2)"
  plot(4) = gsn_csm_contour_map_ce(wks,fld,res_p)

; ***********************************************
; PRECT
; ***********************************************

  fld = PRECT(t,:)
  fld = 1000*24*3600*fld ;convert from m/s to mm/day
  str = "PRECT (mm/day)"

  fmin = -10.0
  fmax =  10.0

;  res_p@cnLevelSelectionMode = "AutomaticLevels"    ; manually set the contour levels with the following 3 resources
  res_p@cnLevelSelectionMode = "ManualLevels"    ; manually set the contour levels with the following 3 resources
  res_p@cnMinLevelValF  = fmin                        ; set the minimum contour level
  res_p@cnMaxLevelValF  = fmax                   ; set the maximum contour level
  res_p@cnLevelSpacingF = (fmax-fmin)/50.0                  ; set the interval between contours

  min_f = min(fld(0:ncol_p-1))
  max_f = max(fld(0:ncol_p-1))
  res_p@lbTitleString    = "global min = "+sprintf("%6.4g",min_f)+"  global  max="+sprintf("%6.4g",max_f)
  print(str+" min/max: "+min_f+"  "+max_f)

  res_p@sfXArray            = lon1d
  res_p@sfYArray            = lat1d

  res_p@gsnLeftString        = str
  global_ave2 = sum(area_p(0:ncol_p-1)*fld(0:ncol_p-1))*inv_area_sphere_p
  print(str+" global ave: "+global_ave2)
  res_p@gsnRightString        = "mean: "+sprintf("%3.2g", global_ave2        )+" mm/day"
  plot(5) = gsn_csm_contour_map_ce(wks,fld,res_p)



  delete(resP)
  resP                  = True                   ; modify the panel plot
  resP@gsnFrame         = False                  ; don't advance panel plot
  resP@gsnPanelFigureStrings = (/"a","b","c","d","e","f"/)
  resP@gsnPanelYWhiteSpacePercent = 4.0
  resP@gsnPanelXWhiteSpacePercent = 4.0

  resP@txString   = "Enthalpy terms"
  gsn_panel(wks,plot,(/2,3/),resP)               ; now draw as one plot

  txres               = True
  frame(wks)

;PHIPHISPHIPHISPHIPHISPHIPHISPHIPHISPHIPHISPHIPHISPHIPHISPHIPHISPHIPHISPHIPHISPHIPHISPHIPHISPHIPHISPHIPHISPHIPHIS
;PHIPHISPHIPHISPHIPHISPHIPHISPHIPHISPHIPHISPHIPHISPHIPHISPHIPHISPHIPHISPHIPHISPHIPHISPHIPHISPHIPHISPHIPHISPHIPHIS
;PHIPHISPHIPHISPHIPHISPHIPHISPHIPHISPHIPHISPHIPHISPHIPHISPHIPHISPHIPHISPHIPHISPHIPHISPHIPHISPHIPHISPHIPHISPHIPHIS


  filePlot = panelTitle+"PHIS"
  wks = gsn_open_wks(output_format,filePlot)
  gsn_define_colormap(wks,"NCV_jaisnd");3gauss")      
    

  res                     = True         ; plot modifications desired
  res@gsnMaximize         = False         ; Maximize size of plot in frame
  res@gsnSpreadColors     = True         ; Use full colormap, but start
  res@gsnDraw             = False           ; don't draw
  res@gsnFrame            = False           ; don't advance frame
  res@cnFillOn            = True         ; Turn on contour fill
;  res@cnFillMode          = "AreaFill"   ; Style of fill. You can also
                                         ; use "CellFill" and "RasterFill"
  res@cnLinesOn           = False        ; Turn off contour lines
  res@cnLineLabelsOn      = False        ; Turn off contour line labels
  res@lbLabelAutoStride   = True         ; Clean up labelbar labels.
  res@cnFillMode="RasterFill"

  res@gsnStringFontHeightF = 0.02
  res@tiMainFontHeightF = 0.025
  res@lbLabelFontHeightF = 0.02
  res@tmXBLabelFontHeightF = 0.015
  res@tmYLLabelFontHeightF = 0.015
  res@pmLabelBarOrthogonalPosF = 0.2

  res@lbBoxLinesOn = False

  res_p=res
  delete(plot)
  plot = new(3,graphic) 


; ***********************************************
; update all H2O (dme_adjust)
; ***********************************************

  fld = PHIS_2d(0:ncol_p-1)
  str = "d/dt(PHIS term H2O vayring)"
  fmin = -2.
  fmax = 2.

;  res_p@cnLevelSelectionMode = "AutomaticLevels"    ; manually set the contour levels with the following 3 resources
  res_p@cnLevelSelectionMode = "ManualLevels"    ; manually set the contour levels with the following 3 resources
  res_p@cnMinLevelValF  = fmin                        ; set the minimum contour level
  res_p@cnMaxLevelValF  = fmax                   ; set the maximum contour level
  res_p@cnLevelSpacingF = (fmax-fmin)/50.0                  ; set the interval between contours

  min_f = min(fld(0:ncol_p-1))
  max_f = max(fld(0:ncol_p-1))
  res_p@lbTitleString    = "global min = "+sprintf("%6.4g",min_f)+"  global  max="+sprintf("%6.4g",max_f)
  print(str+" min/max: "+min_f+"  "+max_f)

  res_p@sfXArray            = lon1d
  res_p@sfYArray            = lat1d

  res_p@gsnLeftString        = str
  global_ave = sum(area_p(0:ncol_p-1)*fld(0:ncol_p-1))*inv_area_sphere_p
  print(str+" global ave: "+global_ave)
  res_p@gsnRightString        = "mean: "+sprintf("%3.2g", global_ave        )+" W/m^2"
  plot(0) = gsn_csm_contour_map_ce(wks,fld,res_p)


; ***********************************************
; PHIS flux
; ***********************************************

  delete(fld)
  fld = FPHIS(t,0:ncol_p-1)
  str = "FPHIS=FH2O*PHIS"
  fmin = -2.
  fmax = 2.

;  res_p@cnLevelSelectionMode = "AutomaticLevels"    ; manually set the contour levels with the following 3 resources
  res_p@cnLevelSelectionMode = "ManualLevels"    ; manually set the contour levels with the following 3 resources
  res_p@cnMinLevelValF  = fmin                        ; set the minimum contour level
  res_p@cnMaxLevelValF  = fmax                   ; set the maximum contour level
  res_p@cnLevelSpacingF = (fmax-fmin)/50.0                  ; set the interval between contours

  min_f = min(fld(0:ncol_p-1))
  max_f = max(fld(0:ncol_p-1))
  res_p@lbTitleString    = "global min = "+sprintf("%6.4g",min_f)+"  global  max="+sprintf("%6.4g",max_f)
  print(str+" min/max: "+min_f+"  "+max_f)

  res_p@sfXArray            = lon1d
  res_p@sfYArray            = lat1d

  res_p@gsnLeftString        = str
  global_ave = sum(area_p(0:ncol_p-1)*fld(0:ncol_p-1))*inv_area_sphere_p
  print(str+" global ave: "+global_ave)

  res_p@gsnRightString        = "mean: "+sprintf("%3.2g", global_ave        )+"W/m^2"
  plot(1) = gsn_csm_contour_map_ce(wks,fld,res_p)
; ***********************************************
; PHIS flux
; ***********************************************

  delete(fld)
  fld = PHIS_2d(0:ncol_p-1)-FPHIS(t,0:ncol_p-1)
  str = "(g) - (h)"
  fmin = -2.
  fmax = 2.

;  res_p@cnLevelSelectionMode = "AutomaticLevels"    ; manually set the contour levels with the following 3 resources
  res_p@cnLevelSelectionMode = "ManualLevels"    ; manually set the contour levels with the following 3 resources
  res_p@cnMinLevelValF  = fmin                        ; set the minimum contour level
  res_p@cnMaxLevelValF  = fmax                   ; set the maximum contour level
  res_p@cnLevelSpacingF = (fmax-fmin)/50.0                  ; set the interval between contours

  min_f = min(fld(0:ncol_p-1))
  max_f = max(fld(0:ncol_p-1))
  res_p@lbTitleString    = "global min = "+sprintf("%6.4g",min_f)+"  global  max="+sprintf("%6.4g",max_f)
  print(str+" min/max: "+min_f+"  "+max_f)

  res_p@sfXArray            = lon1d
  res_p@sfYArray            = lat1d

  res_p@gsnLeftString        = str
  global_ave = sum(area_p(0:ncol_p-1)*fld(0:ncol_p-1))*inv_area_sphere_p
  print(str+" global ave: "+global_ave)

  res_p@gsnRightString        = "mean: "+sprintf("%3.2g", global_ave        )+"W/m^2"
  plot(2) = gsn_csm_contour_map_ce(wks,fld,res_p)


  delete(resP)
  resP                  = True                   ; modify the panel plot
  resP@gsnFrame         = False                  ; don't advance panel plot
  resP@gsnPanelFigureStrings = (/"a","b","c"/)
  resP@gsnPanelYWhiteSpacePercent = 4.0
  resP@gsnPanelXWhiteSpacePercent = 4.0

  resP@txString   = "Surface geopotential terms"
  gsn_panel(wks,plot,(/5,3/),resP)               ; now draw as one plot

  txres               = True
  frame(wks)


;otherotherotherotherotherotherotherotherotherotherotherotherotherotherotherotherotherotherotherotherotherotherotherotherotherotherotherotherotherotherotherother
;otherotherotherotherotherotherotherotherotherotherotherotherotherotherotherotherotherotherotherotherotherotherotherotherotherotherotherotherotherotherotherother
;otherotherotherotherotherotherotherotherotherotherotherotherotherotherotherotherotherotherotherotherotherotherotherotherotherotherotherotherotherotherotherother

 filePlot = panelTitle+"Other"
  wks = gsn_open_wks(output_format,filePlot)
  gsn_define_colormap(wks,"NCV_jaisnd");3gauss")      
    

  res                     = True         ; plot modifications desired
  res@gsnMaximize         = False         ; Maximize size of plot in frame
  res@gsnSpreadColors     = True         ; Use full colormap, but start
  res@gsnDraw             = False           ; don't draw
  res@gsnFrame            = False           ; don't advance frame
  res@cnFillOn            = True         ; Turn on contour fill
;  res@cnFillMode          = "AreaFill"   ; Style of fill. You can also
                                         ; use "CellFill" and "RasterFill"
  res@cnLinesOn           = False        ; Turn off contour lines
  res@cnLineLabelsOn      = False        ; Turn off contour line labels
  res@lbLabelAutoStride   = True         ; Clean up labelbar labels.
  res@cnFillMode="RasterFill"

  res@gsnStringFontHeightF = 0.02
  res@tiMainFontHeightF = 0.025
  res@lbLabelFontHeightF = 0.02
  res@tmXBLabelFontHeightF = 0.015
  res@tmYLLabelFontHeightF = 0.015
  res@pmLabelBarOrthogonalPosF = 0.2

  res@lbBoxLinesOn = False

  res_p=res
  delete(plot)
  plot = new(4,graphic) 
; ***********************************************
; difference
; ***********************************************

  delete(fld)
;  fld = cpT_2d(0:ncol_p-1)-FNET_TBOT(t,0:ncol_p-1)*cpair
;  str = "(e)-(f)"
  fld = (FNET_TS(t,0:ncol_p-1)-FNET_TBOT(t,0:ncol_p-1))*cpair
  str = "FNET_TS-FNET_TBOT"
  fmin = -0.7
  fmax = 0.7

;  res_p@cnLevelSelectionMode = "AutomaticLevels"    ; manually set the contour levels with the following 3 resources
  res_p@cnLevelSelectionMode = "ManualLevels"    ; manually set the contour levels with the following 3 resources
  res_p@cnMinLevelValF  = fmin                        ; set the minimum contour level
  res_p@cnMaxLevelValF  = fmax                   ; set the maximum contour level
  res_p@cnLevelSpacingF = (fmax-fmin)/50.0                  ; set the interval between contours

  min_f = min(fld(0:ncol_p-1))
  max_f = max(fld(0:ncol_p-1))
  res_p@lbTitleString    = "global min = "+sprintf("%6.4g",min_f)+"  global  max="+sprintf("%6.4g",max_f)
  print(str+" min/max: "+min_f+"  "+max_f)

  res_p@sfXArray            = lon1d
  res_p@sfYArray            = lat1d

  res_p@gsnLeftString        = str
  global_ave = sum(area_p(0:ncol_p-1)*fld(0:ncol_p-1))*inv_area_sphere_p
  print(str+" global ave: "+global_ave)

  res_p@gsnRightString        = "mean: "+sprintf("%3.2g", global_ave        )
  plot(0) = gsn_csm_contour_map_ce(wks,fld,res_p)


; ***********************************************
; 
; ***********************************************

  delete(fld)
  fld = TS(t,:)-TSatm(t,:)
  str = "TS (skin)-T(k=nlev) atm "
  fmin = -4.
  fmax = 4.

;  res_p@cnLevelSelectionMode = "AutomaticLevels"    ; manually set the contour levels with the following 3 resources
  res_p@cnLevelSelectionMode = "ManualLevels"    ; manually set the contour levels with the following 3 resources
  res_p@cnMinLevelValF  = fmin                        ; set the minimum contour level
  res_p@cnMaxLevelValF  = fmax                   ; set the maximum contour level
  res_p@cnLevelSpacingF = (fmax-fmin)/50.0                  ; set the interval between contours

  min_f = min(fld(0:ncol_p-1))
  max_f = max(fld(0:ncol_p-1))
  res_p@lbTitleString    = "global min = "+sprintf("%6.4g",min_f)+"  global  max="+sprintf("%6.4g",max_f)
  print(str+" min/max: "+min_f+"  "+max_f)

  res_p@sfXArray            = lon1d
  res_p@sfYArray            = lat1d

  res_p@gsnLeftString        = str
  global_ave = sum(area_p(0:ncol_p-1)*fld(0:ncol_p-1))*inv_area_sphere_p
  print(str+" global ave: "+global_ave)

  res_p@gsnRightString        = "mean: "+sprintf("%3.2g", global_ave        )
  plot(1) = gsn_csm_contour_map_ce(wks,fld,res_p)


; ***********************************************
; 
; ***********************************************



; ***********************************************
; Water flux
; ***********************************************

  fld = FH2O(t,0:ncol_p-1)
  str = "Flux of total water"
  fmin = -0.0001
  fmax = 0.0001

  res_p@cnLevelSelectionMode = "ManualLevels"    ; manually set the contour levels with the following 3 resources
  res_p@cnMinLevelValF  = fmin                        ; set the minimum contour level
  res_p@cnMaxLevelValF  = fmax                   ; set the maximum contour level
  res_p@cnLevelSpacingF = (fmax-fmin)/50.0                  ; set the interval between contours

  min_f = min(fld(0:ncol_p-1))
  max_f = max(fld(0:ncol_p-1))
  res_p@lbTitleString    = "global min = "+sprintf("%6.4g",min_f)+"  global  max="+sprintf("%6.4g",max_f)
  print(str+" min/max: "+min_f+"  "+max_f)

  res_p@sfXArray            = lon1d
  res_p@sfYArray            = lat1d

  res_p@gsnLeftString        = str
  global_ave = sum(area_p(0:ncol_p-1)*fld(0:ncol_p-1))*inv_area_sphere_p
  print(str+" global ave: "+global_ave)
  res_p@gsnRightString        = "mean: "+sprintf("%3.2g", global_ave        )
  plot(2) = gsn_csm_contour_map_ce(wks,fld,res_p)


; ***********************************************
; 
; ***********************************************

  delete(fld)

  fld = PRECT(t,:)
  fld = 1000*24*3600*fld ;convert from m/s to mm/day
  str = "PRECT (mm/day)"
  fmin = -10.
  fmax = 10.

;  res_p@cnLevelSelectionMode = "AutomaticLevels"    ; manually set the contour levels with the following 3 resources
  res_p@cnLevelSelectionMode = "ManualLevels"    ; manually set the contour levels with the following 3 resources
  res_p@cnMinLevelValF  = fmin                        ; set the minimum contour level
  res_p@cnMaxLevelValF  = fmax                   ; set the maximum contour level
  res_p@cnLevelSpacingF = (fmax-fmin)/50.0                  ; set the interval between contours

  min_f = min(fld(0:ncol_p-1))
  max_f = max(fld(0:ncol_p-1))
  res_p@lbTitleString    = "global min = "+sprintf("%6.4g",min_f)+"  global  max="+sprintf("%6.4g",max_f)
  print(str+" min/max: "+min_f+"  "+max_f)

  res_p@sfXArray            = lon1d
  res_p@sfYArray            = lat1d

  res_p@gsnLeftString        = str
  global_ave = sum(area_p(0:ncol_p-1)*fld(0:ncol_p-1))*inv_area_sphere_p
  print(str+" global ave: "+global_ave)

  res_p@gsnRightString        = "mean: "+sprintf("%3.2g", global_ave        )
  plot(3) = gsn_csm_contour_map_ce(wks,fld,res_p)

 



 
  delete(resP)
  resP                  = True                   ; modify the panel plot
  resP@gsnFrame         = False                  ; don't advance panel plot
  resP@gsnPanelFigureStrings = (/"a","b","c","d","e","f","g","h","i","j"/)
  resP@gsnPanelYWhiteSpacePercent = 4.0
  resP@gsnPanelXWhiteSpacePercent = 4.0

  resP@txString   = "Assess missing terms"
  gsn_panel(wks,plot,(/2,2/),resP)               ; now draw as one plot

  txres               = True
  frame(wks)



   end if
  end do


end

