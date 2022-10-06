load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

begin
  print("adsfasdf")
  fname  =  systemfunc("ls "+fname)
  fnamez =  systemfunc("ls "+fnamez)

  fall  = addfile(fname ,"r")       ; Open netCDF files.  
  fallz = addfile(fnamez,"r")       ; Open netCDF files.  

  print("File used for energy diagnostics "+fname)
  print("done input summary")
  print(" ")  
  print("Summary of input to ncl script:")
  print("===============================")
  print(" ")

  print("Your data is for the MPAS dycore")  
  lcp_moist = "false"
  tmp       = fall->SE_phAP(:,:)
  nt        = dimsizes(tmp(:,0))
  ncol      = dimsizes(tmp(0,:))
  ;
  ; initialization
  ;
  g         = 9.80665D0
  PI        = get_pi("double")
  deg2rad   = PI/180.0
  eps       = 1.0D-13

  area      = fall->area   (:  )
  area_sphere = sum(area(0:ncol-1))      
  accuracy_area_sphere = (area_sphere-4D0*PI)/(4D0*PI)
  print("Area of sphere in dataset "+area_sphere)
  print("Relative error            "+accuracy_area_sphere)
  inv_area_sphere   = 1.0D0/area_sphere
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

  qsize = 3;xxx


  print(" ")
  print(" ")

  print("suffix (parameterization side)")
  print("pBF: state passed to parameterizations, before energy fixer")
  print("pBP: after energy fixer, before parameterizations")
  print("pAP: after last phys_update in parameterizations and state saved for energy fixer")
  print("pAM: after dry mass correction calculate but not used in SE and MPAS")
  print(" ")
  print("history files saved off here")


  str_vars      = (/"PO","SE","KE","PT","PS","WV","WL","WI"/);,"TT"/)
  str_stage     = (/"phBF","phBP","phAP","phAM","dyBF","dyBP","dyAP","dyAM","dBF","dAP","dAM"/)  ; "dyAP","dyBD","dyAT","dyAF","dyAD"  ,"dyAR" ,"dyBF"/)
  subcycle      = (/1.0D0,1.0D0,1.0D0,1.0D0,1.0D0,1.0D0,1.0D0,1.0D0,1.0D0,1.0D0,1.0D0,1.0D0,1.0D0,1.0D0/)  ;,1.0D0,1.0D0,1.0D0,nsplit,nsplit,nsplit,1.0D0/)

  nvars       = dimsizes(str_vars)
  nstage      = dimsizes(str_stage)

  vars          = new((/nvars,nstage,nt,ncol/),typeof(tmp(0,0)))
  tot_energy    = new((/nt,nstage/),typeof(tmp(0,0)))
  tot_energy_t  = new((/nt,nstage/),typeof(tmp(0,0)))
  tot_energy_2d = new((/nt,6,nstage,ncol/),typeof(tmp(0,0)))

  tracer_mass = new((/nt,3,nstage/),typeof(tmp(0,0)))
  int_vars    = new((/nvars,nstage,nt/)     ,typeof(tmp(0,0)))    

  int_vars  = 0.0D0
  ;
  ; check that variables are on file and read them and compute global integrals
  ;

  print("nvars="+nvars)
  do variable=0,nvars-1
    do stage=0,nstage-1
      str = str_vars(variable)+"_"+str_stage(stage)
      print("xxx "+str)
      if (isfilevar(fall,str)) then
        fld = fall->$str$(:,:)
        do t=0,nt-1	      	  
          int_vars(variable,stage,t) = sum(area(0:ncol-1)*fld(t,0:ncol-1))*inv_area_sphere
          vars(variable,stage,t,:) = fld(t,0:ncol-1)
        end do
        delete(fld)	
      else
        print("Variable "+str+" not found")
;        exit
      end if
    end do
  end do
  print("Done integrals")
  ;
  ; compute global integrals
  ;

  ;
  ; For code clarity define indices corresponding to the different stages
  ;
  if (isfilevar(fall,"SE_phBF")) then
    ipBF = str_match_ind(str_stage,"phBF")
  else
    iPBF = -1
  end if
  if (isfilevar(fall,"SE_phBP")) then
    ipBP = str_match_ind(str_stage,"phBP")
  else
    iPBP = -1
  end if
  if (isfilevar(fall,"SE_phAP")) then
    ipAP = str_match_ind(str_stage,"phAP")
  else
    iPAP = -1
  end if
  if (isfilevar(fall,"SE_phAM")) then
    ipAM = str_match_ind(str_stage,"phAM")
  else
    iPAM = -1
  end if

  if (isfilevar(fall,"SE_dyBF")) then
    idyBF = str_match_ind(str_stage,"dyBF")
  else
    idyBF = -1
  end if
  if (isfilevar(fall,"SE_dyBP")) then
    idyBP = str_match_ind(str_stage,"dyBP")
  else
    idyBP = -1
  end if
  if (isfilevar(fall,"SE_dyAP")) then
    idyAP = str_match_ind(str_stage,"dyAP")
  else
    idyAP = -1
  end if
  if (isfilevar(fall,"SE_dyAM")) then
    idyAM = str_match_ind(str_stage,"dyAM")
  else
    idyAM = -1
  end if



  if (isfilevar(fall,"SE_dDP")) then
    idDP = str_match_ind(str_stage,"dDP")
  else
    idDP = -1
  end if   
  if (isfilevar(fall,"SE_dAM")) then
    idAM = str_match_ind(str_stage,"dAM")
  else
    idAM = -1
  end if                      
  if (isfilevar(fall,"SE_dPD")) then
    idPD = str_match_ind(str_stage,"dPD")
  else
    idPD = -1
  end if                      
  if (isfilevar(fall,"SE_dAP")) then
    idAP = str_match_ind(str_stage,"dAP")
  else
    idAP = -1
  end if                      
  if (isfilevar(fall,"SE_dBF")) then
    idBF = str_match_ind(str_stage,"dBF")
  else
    idBF = -1
  end if                      
  if (isfilevar(fall,"SE_dED")) then
    idED = str_match_ind(str_stage,"dED")
  else
    idED = -1
  end if                      

  do t=0,nt-1
    do j=0,nstage-1
      print("stage "+j)
      print("time  "+t)
      tot_energy(t,j)     = subcycle(j)*(int_vars(0,j,t)+int_vars(1,j,t)+int_vars(2,j,t))
;       tot_energy(t,j)  = subcycle(j)*int_vars(0,j,t) ; potential
;       tot_energy(t,j)  = subcycle(j)*int_vars(1,j,t) ; thermal energy only
;        tot_energy(t,j)  = subcycle(j)*int_vars(2,j,t); K energy only	       
      print("t="+t+": "+str_stage(j)+"   "+tot_energy(t,j))
      ;
      ; For spatial plot "PO","SE","KE"
      ;
      tot_energy_2d(t,0,j,:) = subcycle(j)*(vars(0,j,t,:)+vars(1,j,t,:)+vars(2,j,t,:))
      tot_energy_2d(t,1,j,:) = subcycle(j)*(vars(0,j,t,:)); themal energy only
      tot_energy_2d(t,2,j,:) = subcycle(j)*(vars(1,j,t,:)); themal energy only
      tot_energy_2d(t,3,j,:) = subcycle(j)*(vars(2,j,t,:)); K only
      tot_energy_2d(t,4,j,:) = subcycle(j)*(vars(3,j,t,:)); ptop
      tot_energy_2d(t,5,j,:) = subcycle(j)*(vars(4,j,t,:)); ps
      ;
      ;
      ;
      idx = str_match_ind(str_vars,"WV")
      tracer_mass(t,0,j) = subcycle(j)*(int_vars(idx,j,t))
      idx = str_match_ind(str_vars,"WL")
      tracer_mass(t,1,j) = subcycle(j)*(int_vars(idx,j,t))
;      print("iii "+j+" "+tracer_mass(t,2,j)+" "+int_vars(idx,j,t))		
      idx = str_match_ind(str_vars,"WI")
      tracer_mass(t,2,j) = subcycle(j)*(int_vars(idx,j,t))
;      print("qqq "+str_stage(j)+" "+tracer_mass(t,2,j)+" "+int_vars(idx,j,t))
    end do
  end do
    
  print("  ")
  print("  ")
  print("  ")
  print("ANALYSIS_START_ANALYSIS_START_ANALYSIS_START_ANALYSIS_START_ANALYSIS_START_ANALYSIS_START")
  print("  ")
  print("  ")
  print("  ")


  do t=0,nt-1
    ph_param  = (tot_energy(t,ipAP)-tot_energy(t,ipBP))/dtime
    ph_EFIX  = (tot_energy(t,ipBP)-tot_energy(t,ipBF))/dtime
    ph_DMEA  = (tot_energy(t,ipAM)-tot_energy(t,ipAP))/dtime
    ph_DMEA_2d  = (tot_energy_2d(t,0,ipAM,:)-tot_energy_2d(t,0,ipAP,:))/dtime
    ph_param_and_efix  = (tot_energy(t,ipAP)-tot_energy(t,ipBF))/dtime
    ph_param_and_efix_2d  = (tot_energy_2d(t,0,ipAP,:)-tot_energy_2d(t,0,ipBF,:))/dtime


    ph_phys_total      = (tot_energy(t,ipAM)-tot_energy(t,ipBF))/dtime

    dy_param  = (tot_energy(t,idyAP)-tot_energy(t,idyBP))/dtime
    dy_EFIX  = (tot_energy(t,idyBP)-tot_energy(t,idyBF))/dtime
    dy_DMEA  = (tot_energy(t,idyAM)-tot_energy(t,idyAP))/dtime
    dy_DMEA_2d  = (tot_energy_2d(t,0,idyAM,:)-tot_energy_2d(t,0,idyAP,:))/dtime
    dy_param_and_efix  = (tot_energy(t,idyAP)-tot_energy(t,idyBF))/dtime
    dy_phys_total      = (tot_energy(t,idyAM)-tot_energy(t,idyBF))/dtime

    pwork = (tot_energy_2d(t,0,ipAM,:)-tot_energy_2d(t,0,ipAP,:))/dtime






    mpas_param       = (tot_energy(t,idAP)-tot_energy(t,idBF))/dtime
    mpas_param_2d    = (tot_energy_2d(t,0,idAP,:)-tot_energy_2d(t,0,idBF,:))/dtime

    mpas_dmea        = (tot_energy(t,idAM)-tot_energy(t,idAP))/dtime
    mpas_dmea_2d     = (tot_energy_2d(t,0,idAM,:)-tot_energy_2d(t,0,idAP,:))/dtime
    mpas_dmea_Potential_2d     = (tot_energy_2d(t,1,idAM,:)-tot_energy_2d(t,1,idAP,:))/dtime
    mpas_dmea_Internal_2d     = (tot_energy_2d(t,2,idAM,:)-tot_energy_2d(t,2,idAP,:))/dtime
    mpas_dmea_Kinetic_2d     = (tot_energy_2d(t,3,idAM,:)-tot_energy_2d(t,3,idAP,:))/dtime

    mpas_phys_total  = (tot_energy(t,idAM)-tot_energy(t,idBF))/dtime

;  str_vars      = (/"PO","SE","KE","PT","PS","WV","WL","WI"/);,"TT"/)
    mpas_dmea_ps_2d     = (tot_energy_2d(t,4,idAM,:)-tot_energy_2d(t,4,idAP,:))/dtime	
;    mpas_dmea_ps_2d     = tot_energy_2d(t,4,idAM,:)

    print(" ")
    print("=================================================================================")
    print("|                                                                               |")
    print("| ANALYSIS OF ENERGY DIAGNOSTICS IN PHYSICS                                      |")
    print("|                                                                               |")
    print("=================================================================================")
    print(" ")
    print("-------------------------------------------------------")
    print(" CAM physics energy increments (in pressure coordinate)")
    print("-------------------------------------------------------")
    print(" ")
    print("dE/dt parameterizations no efix (param) (pAP-pBP) "+sprintf("%6.3g",ph_param)     +" W/M^2")
    print("dE/dt energy fixer (efix)               (pBP-pBF) "+sprintf("%6.3g",ph_EFIX)     +" W/M^2")
    print("NOTE: energy fixer uses energy formula consistent with dycore (so this is not p-based for MPAS) ")
    print("dE/dt parameterizations + efix          (pAP-pBF) "+sprintf("%6.3g",ph_param_and_efix)     +" W/M^2")
    print(" ")
    print("dE/dt dry mass adjustment (pwork)       (pAM-pAP) "+sprintf("%6.3g",ph_DMEA)     +" W/M^2")
    print("dE/dt physics total (phys)              (pAM-pBF) "+sprintf("%6.3g",ph_phys_total)+" W/M^2")
    print(" ")
    dycore = -ph_EFIX-ph_DMEA
    print("Dycore TE dissipation estimated from physics in pressure coordinate "+sprintf("%6.3g",dycore)    +" W/M^2")
    print("(assuming no physics-dynamics coupling errors)    ")
    print(" ")
    print("-----------------------------------------------------------------------------------")
    print(" CAM physics dynamical core consistent energy increments (for MPAS in z coordinate)")
    print("-----------------------------------------------------------------------------------")
    print(" ")
    print("dE/dt parameterizations no efix (param) (dyAP-dyBP) "+sprintf("%6.3g",dy_param)     +" W/M^2")
    print("dE/dt energy fixer (efix)               (dyBP-dyBF) "+sprintf("%6.3g",dy_EFIX)     +" W/M^2")
    print("dE/dt parameterizations + efix          (dyAP-dyBF) "+sprintf("%6.3g",dy_param_and_efix)     +" W/M^2")
    print(" ")
    print("dE/dt dry mass adjustment (pwork)       (dyAM-dyAP) "+sprintf("%6.3g",dy_DMEA)     +" W/M^2")
    print("dE/dt physics total (phys)              (dyAM-dyBF) "+sprintf("%6.3g",dy_phys_total)+" W/M^2")
    print(" ")
    dycore = -dy_EFIX-dy_DMEA
    print("Dycore TE dissipation estimated from physics with dycore energy "+sprintf("%6.3g",dycore)    +" W/M^2")
    print("(assuming no physics-dynamics coupling errors; -efix-dme_adjust)    ")
    print(" ")

    print(" ")
    print("=================================================================================")
    print("|                                                                               |")
    print("| ANALYSIS OF ENERGY DIAGNOSTICS IN PHYSICS dp_coupling (MPAS)                  |")
    print("|                                                                               |")
    print("=================================================================================")
    print(" ")
    print("  ")
    print("dE/dt parameterizations + efix (total physics increment) in MPAS   ")
    print("when adding as one increment - no dribbling (dAP-dBF) "+sprintf("%6.3g",mpas_param)+" W/M^2")
    err = ph_param_and_efix-mpas_param
    print("compare to same tendency in physics (MUST BE SMALL!)  "+sprintf("%6.3g",err)+" W/M^2")
    print(" ")
    print("dE/dt dry mass adjustment in MPAS           (dAM-dAP) "+sprintf("%6.3g",mpas_dmea) +" W/M^2")
    err = dy_DMEA-mpas_dmea
    print("compare to same tendency in physics (MUST BE SMALL!)  "+sprintf("%6.3g",err)+" W/M^2")
    print(" ")
    print(" ")
    print("=================================================================================")
    print("|                                                                               |")
    print("| ANALYSIS OF WATER CONERVATION                                                 |")
    print("|                                                                               |")
    print("=================================================================================")
    print(" ")

    do m_cnst=0,qsize-1

      if (m_cnst.eq.0) then 

        print("Water vapor budget")
        print("------------------")
      end if
      if (m_cnst.eq.1) then
        print("Cloud liquid budget")
        print("-------------------")
      end if
      if (m_cnst.eq.2) then 
        print("Cloud ice budget")
        print("----------------")
      end if
      print(" ")
      pEFIX    = (tracer_mass(t,m_cnst,ipBP)-tracer_mass(t,m_cnst,ipBF))
      pDMEA = (tracer_mass(t,m_cnst,ipAM)-tracer_mass(t,m_cnst,ipAP))
      param           = (tracer_mass(t,m_cnst,ipAP)-tracer_mass(t,m_cnst,ipBP))
      param_mpas      = (tracer_mass(t,m_cnst,idAM)-tracer_mass(t,m_cnst,idBF))
      
;xxx      dyn_total       = (tracer_mass(t,m_cnst,idBF)-tracer_mass(t,m_cnst,idED))
      phys_total      = (tracer_mass(t,m_cnst,ipAM)-tracer_mass(t,m_cnst,ipBF))
  
      print("dMASS/dt energy fixer                      (pBP-pBF) "+sprintf("%6.3g",pEFIX)+" Pa")
      print("dMASS/dt parameterizations                 (pAP-pBP) "+sprintf("%6.3g",param)+" Pa")
      print("dMASS/dt dry mass adjustment               (pAM-pAP) "+sprintf("%6.3g",pDMEA)+" Pa")
      print(" ")
      print(" ")
      print("dMass/dt physics total in MPAS             (dAM-dBF) "+sprintf("%6.3g",param_mpas)+" Pa")
      err = (param_mpas-param);/param
      print("Is mass budget closed?          (pAP-pBP)- (dAM-dBF) "+err)
      print("---------------------------------------------------------------------------------------------------")
      print(" ")
    end do
  end do

  print("start plotting")

  output_format = "pdf"
  filePlot = "mpas_energy"
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

  res@gsnStringFontHeightF = 0.02
  res@tiMainFontHeightF = 0.025
  res@lbLabelFontHeightF = 0.02
  res@tmXBLabelFontHeightF = 0.015
  res@tmYLLabelFontHeightF = 0.015
  res@pmLabelBarOrthogonalPosF = 0.2

  res@lbBoxLinesOn = False


  cellfill = True
  scripFile = "/glade/p/cesmdata/inputdata/atm/cam/coords/mpasa120_SCRIP_desc.200911.nc"
  sf =  addfile(scripFile,"r")
  if (cellfill) then
    lat1d   = sf->grid_center_lat
    lon1d   = sf->grid_center_lon
    latvert = sf->grid_corner_lat
    lonvert = sf->grid_corner_lon
	
    if (lonvert@units.eq."radians") then
      r2d  = get_r2d("double")
      latvert = r2d*latvert
      lonvert = r2d*lonvert
    end if

    res@cnFillMode = "CellFill"
    res@sfXCellBounds = lonvert
    res@sfYCellBounds = latvert

    delete([/sf,latvert,lonvert/])
  else
    lat1d = a->lat
    lon1d = a->lon
    res@cnFillMode          ="RasterFill"
  end if   


  res_p=res

  plot = new(9,graphic) 

  res_p@lbTitleFontHeightF= .02                ; make title smaller
;  res_p@lbTitleString    = "W/m^s"                ; title string
;  res_p@lbTitlePosition  = "Right"              ; title position
;  res_p@lbTitleFontHeightF= .015                ; make title smaller
;  res@lbTitleDirection = "Across"             ; title direction



; ***********************************************
; (1) DME adjust computed by MPAS dyagnostics
; ***********************************************

  fld = mpas_dmea_2d(:)
  str = "DMA (MPAS z-coordinates)"
  fmin = -20.
  fmax = 20.

  res_p@cnLevelSelectionMode = "ManualLevels"    ; manually set the contour levels with the following 3 resources
  res_p@cnMinLevelValF  = fmin                        ; set the minimum contour level
  res_p@cnMaxLevelValF  = fmax                   ; set the maximum contour level
  res_p@cnLevelSpacingF = (fmax-fmin)/50.0                  ; set the interval between contours

  min_f = min(fld)
  max_f = max(fld)
  res_p@lbTitleString    = "global min = "+sprintf("%6.4g",min_f)+" W/m^2"+"  global  max="+sprintf("%6.4g",max_f)+" W/m^2"
  print(str+" min/max: "+min_f+"  "+max_f)

  res_p@sfXArray            = lon1d
  res_p@sfYArray            = lat1d

  res_p@gsnLeftString        = str
  global_ave = sum(area*fld)*inv_area_sphere
  print(str+" global ave: "+global_ave)
  res_p@gsnRightString        = "mean: "+sprintf("%3.2g", global_ave        )+" W/m^2"
  plot(0) = gsn_csm_contour_map_ce(wks,fld,res_p)
  fld_mpas = fld
; ***********************************************
; (2) DME adjust computed by MPAS dyagnostics
; ***********************************************

  fld = ph_DMEA_2d(:)
  str = "DMA (CAM physics p-coordinates)"
  fmin = -20.
  fmax = 20.

  res_p@cnLevelSelectionMode = "ManualLevels"    ; manually set the contour levels with the following 3 resources
  res_p@cnMinLevelValF  = fmin                        ; set the minimum contour level
  res_p@cnMaxLevelValF  = fmax                   ; set the maximum contour level
  res_p@cnLevelSpacingF = (fmax-fmin)/50.0                  ; set the interval between contours

  min_f = min(fld)
  max_f = max(fld)
  res_p@lbTitleString    = "global min = "+sprintf("%6.4g",min_f)+" W/m^2"+"  global  max="+sprintf("%6.4g",max_f)+" W/m^2"
  print(str+" min/max: "+min_f+"  "+max_f)

  res_p@sfXArray            = lon1d
  res_p@sfYArray            = lat1d

  res_p@gsnLeftString        = str
  global_ave = sum(area*fld)*inv_area_sphere
  print(str+" global ave: "+global_ave)
  res_p@gsnRightString        = "mean: "+sprintf("%3.2g", global_ave        )+" W/m^2"
  plot(1) = gsn_csm_contour_map_ce(wks,fld,res_p)

  fld_CAM = fld


; ***********************************************
; (3) diff
; ***********************************************
    fld = fld_mpas-fld_CAM
;  fld = fld+mpas_dmea_Potential_2d(:)
  str = "(a) - (b)"
  fmin = -10.0
  fmax = 10.0

  res_p@cnLevelSelectionMode = "ManualLevels"    ; manually set the contour levels with the following 3 resources
  res_p@cnMinLevelValF  = fmin                        ; set the minimum contour level
  res_p@cnMaxLevelValF  = fmax                   ; set the maximum contour level
  res_p@cnLevelSpacingF = (fmax-fmin)/50.0                  ; set the interval between contours

  min_f = min(fld)
  max_f = max(fld)
  res_p@lbTitleString    = "global min = "+sprintf("%6.4g",min_f)+" W/m^2"+"  global  max="+sprintf("%6.4g",max_f)+" W/m^2"
  print(str+" min/max: "+min_f+"  "+max_f)

  res_p@sfXArray            = lon1d
  res_p@sfYArray            = lat1d

  res_p@gsnLeftString        = str
  global_ave = sum(area*fld)*inv_area_sphere
  print(str+" global ave: "+global_ave)
  res_p@gsnRightString        = "mean: "+sprintf("%3.2g", global_ave        )+" W/m^2"
  plot(2) = gsn_csm_contour_map_ce(wks,fld,res_p)


  fld = mpas_dmea_Internal_2d(:)
  str = "DMA internal energy tendency"
  fmin = -20.
  fmax = 20.

  res_p@cnLevelSelectionMode = "ManualLevels"    ; manually set the contour levels with the following 3 resources
  res_p@cnMinLevelValF  = fmin                        ; set the minimum contour level
  res_p@cnMaxLevelValF  = fmax                   ; set the maximum contour level
  res_p@cnLevelSpacingF = (fmax-fmin)/50.0                  ; set the interval between contours

  min_f = min(fld)
  max_f = max(fld)
  res_p@lbTitleString    = "global min = "+sprintf("%6.4g",min_f)+" W/m^2"+"  global  max="+sprintf("%6.4g",max_f)+" W/m^2"
  print(str+" min/max: "+min_f+"  "+max_f)

  res_p@sfXArray            = lon1d
  res_p@sfYArray            = lat1d

  res_p@gsnLeftString        = str
  global_ave = sum(area*fld)*inv_area_sphere
  print(str+" global ave: "+global_ave)
  res_p@gsnRightString        = "mean: "+sprintf("%3.2g", global_ave        )+" W/m^2"
  plot(3) = gsn_csm_contour_map_ce(wks,fld,res_p)

; ***********************************************
; (1) DME adjust computed by MPAS dyagnostics
; ***********************************************

  fld = mpas_dmea_Kinetic_2d(:)
  str = "DMA kinetic energy tendency"
  fmin = -0.01
  fmax = 0.01

  res_p@cnLevelSelectionMode = "ManualLevels"    ; manually set the contour levels with the following 3 resources
  res_p@cnMinLevelValF  = fmin                        ; set the minimum contour level
  res_p@cnMaxLevelValF  = fmax                   ; set the maximum contour level
  res_p@cnLevelSpacingF = (fmax-fmin)/50.0                  ; set the interval between contours

  min_f = min(fld)
  max_f = max(fld)
  res_p@lbTitleString    = "global min = "+sprintf("%6.4g",min_f)+" W/m^2"+"  global  max="+sprintf("%6.4g",max_f)+" W/m^2"
  print(str+" min/max: "+min_f+"  "+max_f)

  res_p@sfXArray            = lon1d
  res_p@sfYArray            = lat1d

  res_p@gsnLeftString        = str
  global_ave = sum(area*fld)*inv_area_sphere
  print(str+" global ave: "+global_ave)
  res_p@gsnRightString        = "mean: "+sprintf("%3.2g", global_ave        )+" W/m^2"
  plot(4) = gsn_csm_contour_map_ce(wks,fld,res_p)

; ***********************************************
; (1) DME adjust computed by MPAS dyagnostics
; ***********************************************

  fld = mpas_dmea_Potential_2d(:)
  str = "DMA potential energy tendency"
  fmin = -2.
  fmax = 2.

  res_p@cnLevelSelectionMode = "ManualLevels"    ; manually set the contour levels with the following 3 resources
  res_p@cnMinLevelValF  = fmin                        ; set the minimum contour level
  res_p@cnMaxLevelValF  = fmax                   ; set the maximum contour level
  res_p@cnLevelSpacingF = (fmax-fmin)/50.0                  ; set the interval between contours

  min_f = min(fld)
  max_f = max(fld)
  res_p@lbTitleString    = "global min = "+sprintf("%6.4g",min_f)+" W/m^2"+"  global  max="+sprintf("%6.4g",max_f)+" W/m^2"
  print(str+" min/max: "+min_f+"  "+max_f)

  res_p@sfXArray            = lon1d
  res_p@sfYArray            = lat1d

  res_p@gsnLeftString        = str
  global_ave = sum(area*fld)*inv_area_sphere
  print(str+" global ave: "+global_ave)
  res_p@gsnRightString        = "mean: "+sprintf("%3.2g", global_ave        )+" W/m^2"
  plot(5) = gsn_csm_contour_map_ce(wks,fld,res_p)

; ***********************************************
; (1) DME adjust computed by MPAS dyagnostics
; ***********************************************





; ***********************************************
; (3) work term
; ***********************************************
  fld = fall->ztop_dPtopdt(0,:)
;  fld = fall->ptop_dztopdt(0,:); ztop_dPtopdt(0,:)

  str = "Model top work term"
  fmin = -0.01
  fmax =  0.01

  res_p@cnLevelSelectionMode = "ManualLevels"    ; manually set the contour levels with the following 3 resources
  res_p@cnMinLevelValF  = fmin                        ; set the minimum contour level
  res_p@cnMaxLevelValF  = fmax                   ; set the maximum contour level
  res_p@cnLevelSpacingF = (fmax-fmin)/50.0                  ; set the interval between contours

  min_f = min(fld)
  max_f = max(fld)
  res_p@lbTitleString    = "global min = "+sprintf("%6.4g",min_f)+" W/m^2"+"  global  max="+sprintf("%6.4g",max_f)+" W/m^2"
  print(str+" min/max: "+min_f+"  "+max_f)

  res_p@sfXArray            = lon1d
  res_p@sfYArray            = lat1d

  res_p@gsnLeftString        = str
  global_ave = sum(area*fld)*inv_area_sphere
  print(str+" global ave: "+global_ave)
  res_p@gsnRightString        = "mean: "+sprintf("%3.2g", global_ave        )+" W/m^2"
  plot(6) = gsn_csm_contour_map_ce(wks,fld,res_p)

; ***********************************************
; (1) DME adjust computed by MPAS dyagnostics
; ***********************************************

  fld = fallz->ptop_mpas(0,:)
  str = "MPAS model top pressure"

  res_p@cnLevelSelectionMode = "AutomaticLevels";ManualLevels"    ; manually set the contour levels with the following 3 resources
;  fmin = -20.
;  fmax = 20.

;  res_p@cnLevelSelectionMode = "ManualLevels"    ; manually set the contour levels with the following 3 resources
;  res_p@cnMinLevelValF  = fmin                        ; set the minimum contour level
;  res_p@cnMaxLevelValF  = fmax                   ; set the maximum contour level
;  res_p@cnLevelSpacingF = (fmax-fmin)/50.0                  ; set the interval between contours

  min_f = min(fld)
  max_f = max(fld)
  res_p@lbTitleString    = "global min = "+sprintf("%6.4g",min_f)+" W/m^2"+"  global  max="+sprintf("%6.4g",max_f)+" W/m^2"
  print(str+" min/max: "+min_f+"  "+max_f)

  res_p@sfXArray            = lon1d
  res_p@sfYArray            = lat1d

  res_p@gsnLeftString        = str
  global_ave = sum(area*fld)*inv_area_sphere
  print(str+" global ave: "+global_ave)
  res_p@gsnRightString        = "mean: "+sprintf("%3.2g", global_ave        )+" W/m^2"
  plot(7) = gsn_csm_contour_map_ce(wks,fld,res_p)



; ***********************************************
; 
; ***********************************************
  
  delete(fld)
  fld = ph_param_and_efix_2d-mpas_param_2d
  fld = fld+(mpas_dmea_2d(:)-dy_DMEA_2d(:))
  str = "dE phys+dma - dE MPAS phys+dma"

  res_p@cnLevelSelectionMode = "AutomaticLevels";ManualLevels"    ; manually set the contour levels with the following 3 resources

  min_f = min(fld)
  max_f = max(fld)
  res_p@lbTitleString    = "global min = "+sprintf("%6.4g",min_f)+" W/m^2"+"  global  max="+sprintf("%6.4g",max_f)+" W/m^2"
  print(str+" min/max: "+min_f+"  "+max_f)

  res_p@sfXArray            = lon1d
  res_p@sfYArray            = lat1d

  res_p@gsnLeftString        = str
  global_ave = sum(area*fld)*inv_area_sphere
  print(str+" global ave: "+global_ave)
  res_p@gsnRightString        = "mean: "+sprintf("%3.2g", global_ave        )+" W/m^2"
  plot(8) = gsn_csm_contour_map_ce(wks,fld,res_p)



; ***********************************************
; Left-hand side of CAM energy equation
; ***********************************************
  
;  delete(fld)
;  lev = 25
;  zm_phBF = fallz->zm_phBF(0,lev,:)
;  zm_phAP = fallz->zm_phAP(0,lev,:)
;  fld = zm_phAP-zm_phBF
;
;  str = "dz float level "+lev
;  fmin = -300.
;  fmax = 300.
;
;  res_p@cnLevelSelectionMode = "AutomaticLevels";ManualLevels"    ; manually set the contour levels with the following 3 resources
;;  res_p@cnMinLevelValF  = fmin                        ; set the minimum contour level
;;  res_p@cnMaxLevelValF  = fmax                   ; set the maximum contour level
;;  res_p@cnLevelSpacingF = (fmax-fmin)/50.0                  ; set the interval between contours
;
;  min_f = min(fld)
;  max_f = max(fld)
;
;  res_p@lbTitleString    = "global min = "+sprintf("%6.4g",min_f)+" W/m^2"+"  global  max="+sprintf("%6.4g",max_f)+" W/m^2"
;  print(str+" min/max: "+min_f+"  "+max_f)
;
;  res_p@sfXArray            = lon1d
;  res_p@sfYArray            = lat1d
;
;  res_p@gsnLeftString        = str
;  global_ave = sum(area*fld)*inv_area_sphere
;  print(str+" global ave: "+global_ave)
;  res_p@gsnRightString        = "mean: "+sprintf("%3.2g", global_ave        )+" W/m^2"
;
;  plot(5) = gsn_csm_contour_map_ce(wks,fld,res_p)
;
;

; ***********************************************
; plot panel
; ***********************************************

  resP                  = True                   ; modify the panel plot
  resP@gsnFrame         = False                  ; don't advance panel plot
  resP@gsnPanelFigureStrings = (/"a","b","c","d","e","f","g","h","i","j"/)
  resP@gsnPanelYWhiteSpacePercent = 4.0
  resP@gsnPanelXWhiteSpacePercent = 4.0

  resP@txString   = "Energy diagnostics (1yr average)"
  gsn_panel(wks,plot,(/3,3/),resP)               ; now draw as one plot

  txres               = True
  frame(wks)


  filePlot = "mpas_energy_detail"
  wks = gsn_open_wks(output_format,filePlot)

; ***********************************************
; (1) DME adjust computed by MPAS dyagnostics
; ***********************************************

  fld = mpas_dmea_Potential_2d(:)
  str = "Potential energy tendency"
  fmin = -2.
  fmax = 2.

  res_p@cnLevelSelectionMode = "ManualLevels"    ; manually set the contour levels with the following 3 resources
  res_p@cnMinLevelValF  = fmin                        ; set the minimum contour level
  res_p@cnMaxLevelValF  = fmax                   ; set the maximum contour level
  res_p@cnLevelSpacingF = (fmax-fmin)/50.0                  ; set the interval between contours

  min_f = min(fld)
  max_f = max(fld)
  res_p@lbTitleString    = "global min = "+sprintf("%6.4g",min_f)+" W/m^2"+"  global  max="+sprintf("%6.4g",max_f)+" W/m^2"
  print(str+" min/max: "+min_f+"  "+max_f)

  res_p@sfXArray            = lon1d
  res_p@sfYArray            = lat1d

  res_p@gsnLeftString        = str
  global_ave = sum(area*fld)*inv_area_sphere
  print(str+" global ave: "+global_ave)
  res_p@gsnRightString        = "mean: "+sprintf("%3.2g", global_ave        )+" W/m^2"
  plot(0) = gsn_csm_contour_map_ce(wks,fld,res_p)

; ***********************************************
; (1) DME adjust computed by MPAS dyagnostics
; ***********************************************

  fld = mpas_dmea_Internal_2d(:)
  str = "Internal energy tendency"
  fmin = -15.
  fmax = 15.

  res_p@cnLevelSelectionMode = "ManualLevels"    ; manually set the contour levels with the following 3 resources
  res_p@cnMinLevelValF  = fmin                        ; set the minimum contour level
  res_p@cnMaxLevelValF  = fmax                   ; set the maximum contour level
  res_p@cnLevelSpacingF = (fmax-fmin)/50.0                  ; set the interval between contours

  min_f = min(fld)
  max_f = max(fld)
  res_p@lbTitleString    = "global min = "+sprintf("%6.4g",min_f)+" W/m^2"+"  global  max="+sprintf("%6.4g",max_f)+" W/m^2"
  print(str+" min/max: "+min_f+"  "+max_f)

  res_p@sfXArray            = lon1d
  res_p@sfYArray            = lat1d

  res_p@gsnLeftString        = str
  global_ave = sum(area*fld)*inv_area_sphere
  print(str+" global ave: "+global_ave)
  res_p@gsnRightString        = "mean: "+sprintf("%3.2g", global_ave        )+" W/m^2"
  plot(1) = gsn_csm_contour_map_ce(wks,fld,res_p)

; ***********************************************
; (1) DME adjust computed by MPAS dyagnostics
; ***********************************************

  fld = mpas_dmea_Kinetic_2d(:)
  str = "Kinetic energy tendency"
  fmin = -0.01
  fmax = 0.01

  res_p@cnLevelSelectionMode = "ManualLevels"    ; manually set the contour levels with the following 3 resources
  res_p@cnMinLevelValF  = fmin                        ; set the minimum contour level
  res_p@cnMaxLevelValF  = fmax                   ; set the maximum contour level
  res_p@cnLevelSpacingF = (fmax-fmin)/50.0                  ; set the interval between contours

  min_f = min(fld)
  max_f = max(fld)
  res_p@lbTitleString    = "global min = "+sprintf("%6.4g",min_f)+" W/m^2"+"  global  max="+sprintf("%6.4g",max_f)+" W/m^2"
  print(str+" min/max: "+min_f+"  "+max_f)

  res_p@sfXArray            = lon1d
  res_p@sfYArray            = lat1d

  res_p@gsnLeftString        = str
  global_ave = sum(area*fld)*inv_area_sphere
  print(str+" global ave: "+global_ave)
  res_p@gsnRightString        = "mean: "+sprintf("%3.2g", global_ave        )+" W/m^2"
  plot(2) = gsn_csm_contour_map_ce(wks,fld,res_p)

; ***********************************************
; plot panel
; ***********************************************

  resP                  = True                   ; modify the panel plot
  resP@gsnFrame         = False                  ; don't advance panel plot
  resP@gsnPanelFigureStrings = (/"a","b","c","d","e","f","g","h","i","j"/)
  resP@gsnPanelYWhiteSpacePercent = 4.0
  resP@gsnPanelXWhiteSpacePercent = 4.0

  resP@txString   = "Dry-mass adjustment: Potential, internal and kinetic energy tendency"
  gsn_panel(wks,plot,(/1,3/),resP)               ; now draw as one plot

  txres               = True
  frame(wks)



  filePlot = "mpas_misc"
  wks = gsn_open_wks(output_format,filePlot)

; ***********************************************
; (1) DME adjust computed by MPAS dyagnostics
; ***********************************************

  fld = fallz->ptop_mpas(0,:)
  str = "MPAS model top pressure"

  res_p@cnLevelSelectionMode = "AutomaticLevels";ManualLevels"    ; manually set the contour levels with the following 3 resources
;  fmin = -20.
;  fmax = 20.

;  res_p@cnLevelSelectionMode = "ManualLevels"    ; manually set the contour levels with the following 3 resources
;  res_p@cnMinLevelValF  = fmin                        ; set the minimum contour level
;  res_p@cnMaxLevelValF  = fmax                   ; set the maximum contour level
;  res_p@cnLevelSpacingF = (fmax-fmin)/50.0                  ; set the interval between contours

  min_f = min(fld)
  max_f = max(fld)
  res_p@lbTitleString    = "global min = "+sprintf("%6.4g",min_f)+" W/m^2"+"  global  max="+sprintf("%6.4g",max_f)+" W/m^2"
  print(str+" min/max: "+min_f+"  "+max_f)

  res_p@sfXArray            = lon1d
  res_p@sfYArray            = lat1d

  res_p@gsnLeftString        = str
  global_ave = sum(area*fld)*inv_area_sphere
  print(str+" global ave: "+global_ave)
  res_p@gsnRightString        = "mean: "+sprintf("%3.2g", global_ave        )+" W/m^2"
  plot(0) = gsn_csm_contour_map_ce(wks,fld,res_p)

; ***********************************************
; 
; ***********************************************

  fld = mpas_dmea_ps_2d 
  str = "dma dPSdt"

  res_p@cnLevelSelectionMode = "AutomaticLevels";ManualLevels"    ; manually set the contour levels with the following 3 resources
;  fmin = -20.
;  fmax = 20.

;  res_p@cnLevelSelectionMode = "ManualLevels"    ; manually set the contour levels with the following 3 resources
;  res_p@cnMinLevelValF  = fmin                        ; set the minimum contour level
;  res_p@cnMaxLevelValF  = fmax                   ; set the maximum contour level
;  res_p@cnLevelSpacingF = (fmax-fmin)/50.0                  ; set the interval between contours

  min_f = min(fld)
  max_f = max(fld)
  res_p@lbTitleString    = "global min = "+sprintf("%6.4g",min_f)+" W/m^2"+"  global  max="+sprintf("%6.4g",max_f)+" W/m^2"
  print(str+" min/max: "+min_f+"  "+max_f)

  res_p@sfXArray            = lon1d
  res_p@sfYArray            = lat1d

  res_p@gsnLeftString        = str
  global_ave = sum(area*fld)*inv_area_sphere
  print(str+" global ave: "+global_ave)
  res_p@gsnRightString        = "mean: "+sprintf("%3.2g", global_ave        )+" W/m^2"
  plot(1) = gsn_csm_contour_map_ce(wks,fld,res_p)

; ***********************************************
; plot panel
; ***********************************************

  resP                  = True                   ; modify the panel plot
  resP@gsnFrame         = False                  ; don't advance panel plot
  resP@gsnPanelFigureStrings = (/"a","b","c","d","e","f","g","h","i","j"/)
  resP@gsnPanelYWhiteSpacePercent = 4.0
  resP@gsnPanelXWhiteSpacePercent = 4.0

  resP@txString   = "Dry-mass adjustment: Potential, internal and kinetic energy tendency"
  gsn_panel(wks,plot,(/1,2/),resP)               ; now draw as one plot

  txres               = True
  frame(wks)

end