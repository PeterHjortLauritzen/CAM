load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
begin
  option = "nhtfrq1-B-case"

  output_format="pdf"

;++++++++++++++++++++++++++++++++++++++++++++++++++++
;
; 
;
;++++++++++++++++++++++++++++++++++++++++++++++++++++

;  casename = "b.e30_alpha03c_enthalpy.BLT1850.ne30_t232.enthalpy.09"
;;  casename = "b.e30_alpha03c_enthalpy.BLT1850.ne30_t232.enthalpy.06"
;  dir="/glade/derecho/scratch/gmarques/"+casename+"/run/"
;  fname = casename+".cam.h7i."
;  fnameall =  systemfunc("ls "+dir+fname+"*.nc")
;  fall = addfiles(fnameall,"r")       ; Open netCDF files.
;;;
;;; if energy diags on file
;;;
;  fname_energy = casename+".cam.h8i."
;  fnameall_energy =  systemfunc("ls "+dir+fname_energy+"*.nc")
;    print(fnameall_energy)
;  fall_energy = addfiles(fnameall_energy,"r")       ; Open netCDF files.
;
;  tidx = 191


;++++++++++++++++++++++++++++++++++++++++++++++++++++
;
; look at particular time-stap from nhtfrq=1 runs
;
;++++++++++++++++++++++++++++++++++++++++++++++++++++
  if (option .eq. "nhtfrq1-B-case") then
;    tt = "0001-01-03-84600"
    tt = "0001-01-04-00000"
    casename = "b.e30_enthalpy_branch.BLT1850.ne30_t232.enthalpy.12"
    dir="/glade/derecho/scratch/gmarques/"+casename+"/run/"

   ;  casename = "cam_enth_simple_debug"
   ;  dir="/glade/derecho/scratch/pel/"+casename+"/run/"

   fname = casename+".cam.h8i."+tt
   fnameall =  systemfunc("ls "+dir+fname+"*.nc")
   fall = addfiles(fnameall,"r")       ; Open netCDF files.
   print(fnameall)
   ;
   ;  if energy diags on file
   ;
   fname_energy = casename+".cam.h7i."+tt
   fnameall_energy =  systemfunc("ls "+dir+fname_energy+"*.nc")
   print(fnameall_energy)
   fall_energy = addfiles(fnameall_energy,"r")       ; Open netCDF files.
   tidx = 0
  end if

  OCNFRAC=fall[:]->OCNFRAC
  
  tmp=fall[:]->EFIX
  date=fall[:]->date
  time=fall[:]->time

  print(time)

  print("date="+date(tidx))
    print("time="+time(tidx))
  dt=1800.0
  firsttime = 0
  lasttime  = dimsizes(fall[:]->time)-1
  print("Total number of time samples "+lasttime)

  ncol        = dimsizes(tmp(0,:))
  ncol_p      = ncol
  area        = fall[0]->area
  area_sphere_ocn = sum(OCNFRAC(0,0:ncol-1)*area(0:ncol-1))
 
  area_sphere = sum(area(0:ncol-1))

  area_d    = fall_energy[0]->area_d
  area_sphere_d = sum(area_d(:))
    
  lat1d_d        = fall_energy[0]->lat_d          ; Pull off lat/lon variable and
  lon1d_d        = fall_energy[0]->lon_d  
  print("area_sphere_ocn="+area_sphere_ocn)
  print("ncol="+ncol)
  lat1d        = fall[0]->lat          ; Pull off lat/lon variable and
  lon1d        = fall[0]->lon  

  hsnow_liq_ref=fall[:]->hsnow_liq_ref
  hrain_liq_ref=fall[:]->hrain_liq_ref
  hevap_liq_ref=fall[:]->hevap_liq_ref

  enth_flux_to_not_ocn=fall[:]->enth_flux_to_not_ocn
  enth_flux_to_ocn=fall[:]->enth_flux_to_ocn

  EFIX=fall[:]->EFIX
  dEdt_cpdycore=fall[:]->dEdt_cpdycore
  TE_dEdt_param_efix_dynE=fall_energy[:]->TE_dEdt_param_efix_dynE


;  TE_dEdt_dycore=fall_energy[:]->TE_dEdt_dycore   - for some reason all zeros

  TE_dEdt_floating_dyn=fall_energy[:]->TE_dEdt_floating_dyn
  TE_dEdt_vert_remap=fall_energy[:]->TE_dEdt_vert_remap
;  TE_dEdt_dycore = TE_dEdt_floating_dyn+TE_dEdt_vert_remap
  TE_dEdt_dycore = TE_dEdt_vert_remap
  TE_dEdt_phys_tot_in_dyn=fall_energy[:]->TE_dEdt_phys_tot_in_dyn

  TE_dEdt_dme_adjust_dynE=fall_energy[:]->TE_dEdt_dme_adjust_dynE


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
; *NEWPAGENEWPAGENEWPAGENEWPAGENEWPAGENEWPAGE
; *NEWPAGENEWPAGENEWPAGENEWPAGENEWPAGENEWPAGE
; ***********************************************
; ***********************************************
  filePlot = "ocnfrac-flux-to-coupler-from-CAM."+tt
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
  plot = new(9,graphic)
  
; ***********************************************
; 
; ***********************************************

  fld = hsnow_liq_ref(tidx,0:ncol_p-1)
  str = "hsnow_liq_ref (OCN only)"
  fmin = -2.
  fmax = 2.

  res_p@cnLevelSelectionMode = "ManualLevels"    ; manually set the contour levels with the following 3 resources
  res_p@cnMinLevelValF  = fmin                        ; set the minimum contour level
  res_p@cnMaxLevelValF  = fmax                   ; set the maximum contour level
  res_p@cnLevelSpacingF = (fmax-fmin)/50.0                  ; set the interval between contours

  min_f = min(fld)
  max_f = max(fld)
  res_p@lbTitleString    = "global min = "+sprintf("%6.4g",min_f)+"  global  max="+sprintf("%6.4g",max_f)
  print(str+" min/max: "+min_f+"  "+max_f)

  res_p@sfXArray            = lon1d
  res_p@sfYArray            = lat1d

  res_p@gsnLeftString        = str
  global_ave = sum(area*fld)/area_sphere_ocn
  print(str+" global ave: "+global_ave)
  res_p@gsnRightString        = "mean: "+sprintf("%6.4g", global_ave        )+" W/m^2"
  plot(0) = gsn_csm_contour_map_ce(wks,fld,res_p)

; ***********************************************
; 
; ***********************************************

  fld = hrain_liq_ref(tidx,0:ncol_p-1)
  str = "hrain_liq_ref (OCN only)"
  fmin = -25.
  fmax = 25.

  res_p@cnLevelSelectionMode = "ManualLevels"    ; manually set the contour levels with the following 3 resources
  res_p@cnMinLevelValF  = fmin                        ; set the minimum contour level
  res_p@cnMaxLevelValF  = fmax                   ; set the maximum contour level
  res_p@cnLevelSpacingF = (fmax-fmin)/50.0                  ; set the interval between contours

  min_f = min(fld)
  max_f = max(fld)
  res_p@lbTitleString    = "global min = "+sprintf("%6.4g",min_f)+"  global  max="+sprintf("%6.4g",max_f)
  print(str+" min/max: "+min_f+"  "+max_f)

  res_p@sfXArray            = lon1d
  res_p@sfYArray            = lat1d

  res_p@gsnLeftString        = str
  global_ave = sum(area*fld)/area_sphere_ocn
  print(str+" global ave: "+global_ave)
  res_p@gsnRightString        = "mean: "+sprintf("%6.4g", global_ave        )+" W/m^2"
  plot(1) = gsn_csm_contour_map_ce(wks,fld,res_p)

; ***********************************************
; 
; ***********************************************

  fld = hevap_liq_ref(tidx,0:ncol_p-1)
  str = "hevap_liq_ref (OCN only)"
  fmin = -25.
  fmax = 25.

  res_p@cnLevelSelectionMode = "ManualLevels"    ; manually set the contour levels with the following 3 resources
  res_p@cnMinLevelValF  = fmin                        ; set the minimum contour level
  res_p@cnMaxLevelValF  = fmax                   ; set the maximum contour level
  res_p@cnLevelSpacingF = (fmax-fmin)/50.0                  ; set the interval between contours

  min_f = min(fld)
  max_f = max(fld)
  res_p@lbTitleString    = "global min = "+sprintf("%6.4g",min_f)+"  global  max="+sprintf("%6.4g",max_f)
  print(str+" min/max: "+min_f+"  "+max_f)

  res_p@sfXArray            = lon1d
  res_p@sfYArray            = lat1d

  res_p@gsnLeftString        = str
  global_ave = sum(area*fld)/area_sphere_ocn
  print(str+" global ave: "+global_ave)
  res_p@gsnRightString        = "mean: "+sprintf("%6.4g", global_ave        )+" W/m^2"
  plot(2) = gsn_csm_contour_map_ce(wks,fld,res_p)

; ***********************************************
; 
; ***********************************************

  fld = hevap_liq_ref(tidx,0:ncol_p-1)+hrain_liq_ref(tidx,0:ncol_p-1)+hsnow_liq_ref(tidx,0:ncol_p-1)
  str = "(a)+(b)+(c) (OCN only)"
  fmin = -25.
  fmax = 25.

  res_p@cnLevelSelectionMode = "ManualLevels"    ; manually set the contour levels with the following 3 resources
  res_p@cnMinLevelValF  = fmin                        ; set the minimum contour level
  res_p@cnMaxLevelValF  = fmax                   ; set the maximum contour level
  res_p@cnLevelSpacingF = (fmax-fmin)/50.0                  ; set the interval between contours

  min_f = min(fld)
  max_f = max(fld)
  res_p@lbTitleString    = "global min = "+sprintf("%6.4g",min_f)+"  global  max="+sprintf("%6.4g",max_f)
  print(str+" min/max: "+min_f+"  "+max_f)

  res_p@sfXArray            = lon1d
  res_p@sfYArray            = lat1d

  res_p@gsnLeftString        = str
  global_ave = sum(area*fld)/area_sphere_ocn
  print(str+" global ave: "+global_ave)
  res_p@gsnRightString        = "mean: "+sprintf("%6.4g", global_ave        )+" W/m^2"
  plot(3) = gsn_csm_contour_map_ce(wks,fld,res_p)



; ***********************************************
; Plot everything
; ***********************************************
  resP                  = True                   ; modify the panel plot
  resP@gsnFrame         = False                  ; don't advance panel plot
  resP@gsnPanelFigureStrings = (/"a","b","c","d","e","f","g","h","i","j"/)
  resP@gsnPanelYWhiteSpacePercent = 4.0
  resP@gsnPanelXWhiteSpacePercent = 4.0

  resP@txString   = "h-flux to coupler from CAM (means using ocn area only): "+tt


  gsn_panel(wks,plot,(/3,2/),resP)               ; now draw as one plot

  txres               = True
  frame(wks)



  delete(plot)
  plot = new(9,graphic)
; ***********************************************
; *NEWPAGENEWPAGENEWPAGENEWPAGENEWPAGENEWPAGE
; *NEWPAGENEWPAGENEWPAGENEWPAGENEWPAGENEWPAGE
; ***********************************************
; ***********************************************
  filePlot = "enthalpy-flux-end-of-CAM-physics."+tt
  wks = gsn_open_wks(output_format,filePlot)  
; ***********************************************
; 
; ***********************************************

  fld = enth_flux_to_ocn(tidx,0:ncol_p-1)
  str = "Enthalpy flux over ocean end of CAM"
  fmin = -25.
  fmax = 25.

  res_p@cnLevelSelectionMode = "ManualLevels"    ; manually set the contour levels with the following 3 resources
  res_p@cnMinLevelValF  = fmin                        ; set the minimum contour level
  res_p@cnMaxLevelValF  = fmax                   ; set the maximum contour level
  res_p@cnLevelSpacingF = (fmax-fmin)/50.0                  ; set the interval between contours

  min_f = min(fld)
  max_f = max(fld)
  res_p@lbTitleString    = "global min = "+sprintf("%6.4g",min_f)+"  global  max="+sprintf("%6.4g",max_f)
  print(str+" min/max: "+min_f+"  "+max_f)

  res_p@sfXArray            = lon1d
  res_p@sfYArray            = lat1d

  res_p@gsnLeftString        = str
  global_ave = sum(area*fld)/area_sphere
  global_ave_enth_flux_to_ocn = global_ave
  print(str+" global ave: "+global_ave)
  res_p@gsnRightString        = "mean: "+sprintf("%6.4g", global_ave        )+" W/m^2"
  plot(0) = gsn_csm_contour_map_ce(wks,fld,res_p)

; ***********************************************
; 
; ***********************************************

  fld = enth_flux_to_not_ocn(tidx,0:ncol_p-1)
  str = "Enthalpy flux over not ocean end of CAM"
  fmin = -25.
  fmax = 25.

  res_p@cnLevelSelectionMode = "ManualLevels"    ; manually set the contour levels with the following 3 resources
  res_p@cnMinLevelValF  = fmin                        ; set the minimum contour level
  res_p@cnMaxLevelValF  = fmax                   ; set the maximum contour level
  res_p@cnLevelSpacingF = (fmax-fmin)/50.0                  ; set the interval between contours

  min_f = min(fld)
  max_f = max(fld)
  res_p@lbTitleString    = "global min = "+sprintf("%6.4g",min_f)+"  global  max="+sprintf("%6.4g",max_f)
  print(str+" min/max: "+min_f+"  "+max_f)

  res_p@sfXArray            = lon1d
  res_p@sfYArray            = lat1d

  res_p@gsnLeftString        = str
  global_ave = sum(area*fld)/area_sphere
  global_ave_enth_flux_to_ocn = global_ave
  print(str+" global ave: "+global_ave)
  res_p@gsnRightString        = "mean: "+sprintf("%6.4g", global_ave        )+" W/m^2"
  plot(1) = gsn_csm_contour_map_ce(wks,fld,res_p)


; ***********************************************
; 
; ***********************************************

  delete(fld)
  fld = TE_dEdt_dme_adjust_dynE(tidx,:)/1800.0D
  str = "TE_dEdt_dme_adjust_dynE"
  fmin = -25.
  fmax = 25.

  res_p@cnLevelSelectionMode = "ManualLevels"    ; manually set the contour levels with the following 3 resources
  res_p@cnMinLevelValF  = fmin                        ; set the minimum contour level
  res_p@cnMaxLevelValF  = fmax                   ; set the maximum contour level
  res_p@cnLevelSpacingF = (fmax-fmin)/50.0                  ; set the interval between contours

  delete(min_f)
  min_f = min(fld)
  delete(max_f)
  max_f = max(fld)
  res_p@lbTitleString    = "global min = "+sprintf("%0.2f",min_f)+"  global  max="+sprintf("%6.4g",max_f)
  print(str+" min/max: "+min_f+"  "+max_f)

  delete(res_p@sfXArray)
  res_p@sfXArray            = lon1d
  delete(res_p@sfYArray)
  res_p@sfYArray            = lat1d

  res_p@gsnLeftString        = str
  global_ave = sum(area*fld)/area_sphere
  global_ave_TE_dEdt_dme_adjust_dynE=global_ave
  print(str+" global ave: "+global_ave)
  res_p@gsnRightString        = "mean: "+sprintf("%4.3f", global_ave        )+" W/m^2"
  plot(2) = gsn_csm_contour_map_ce(wks,fld,res_p)

 

; ***********************************************
; 
; ***********************************************

  fldE = TE_dEdt_dycore(tidx,:)/1800.0
  str = "TE_dEdt_dycore"
  fmin = -100.
  fmax = 100.

  res_p@cnLevelSelectionMode = "ManualLevels"    ; manually set the contour levels with the following 3 resources
  res_p@cnMinLevelValF  = fmin                        ; set the minimum contour level
  res_p@cnMaxLevelValF  = fmax                   ; set the maximum contour level
  res_p@cnLevelSpacingF = (fmax-fmin)/50.0                  ; set the interval between contours

  delete(min_f)
  delete(max_f)

  min_f = min(fldE)
  max_f = max(fldE)

  res_p@lbTitleString    = "global min = "+sprintf("%6.4g",min_f)+"  global  max="+sprintf("%6.4g",max_f)
  print(str+" min/max: "+min_f+"  "+max_f)

  delete(res_p@sfXArray)
  delete(res_p@sfYArray)

  res_p@sfXArray            = lon1d_d
  res_p@sfYArray            = lat1d_d

  res_p@gsnLeftString        = str
  global_ave = sum(area_d*fldE)/area_sphere_d
  global_ave_dycore = global_ave
  print(str+" global ave: "+global_ave)
  res_p@gsnRightString        = "mean: "+sprintf("%6.4g", global_ave        )+" W/m^2"
  plot(3) = gsn_csm_contour_map_ce(wks,fldE,res_p)

; ***********************************************
; 
; ***********************************************

  fld = EFIX(tidx,:)
  str = "Energy fixer"
  fmin = -25.
  fmax = 25.

  res_p@cnLevelSelectionMode = "ManualLevels"    ; manually set the contour levels with the following 3 resources
  res_p@cnMinLevelValF  = fmin                        ; set the minimum contour level
  res_p@cnMaxLevelValF  = fmax                   ; set the maximum contour level
  res_p@cnLevelSpacingF = (fmax-fmin)/50.0                  ; set the interval between contours

  min_f = min(fld)
  max_f = max(fld)
  res_p@lbTitleString    = "global min = "+sprintf("%6.4g",min_f)+"  global  max="+sprintf("%6.4g",max_f)
  print(str+" min/max: "+min_f+"  "+max_f)

  delete(res_p@sfXArray)
  delete(res_p@sfYArray)
  res_p@sfXArray            = lon1d
  res_p@sfYArray            = lat1d

  res_p@gsnLeftString        = str
  global_ave = sum(area*fld)/area_sphere
  global_ave_efix = global_ave
  print(str+" global ave: "+global_ave)
  res_p@gsnRightString        = "mean: "+sprintf("%6.4g", global_ave        )+" W/m^2"
  plot(4) = gsn_csm_contour_map_ce(wks,fld,res_p)




; ***********************************************
; Plot everything
; ***********************************************
  resP                  = True                   ; modify the panel plot
  resP@gsnFrame         = False                  ; don't advance panel plot
  resP@gsnPanelFigureStrings = (/"a","b","c","d","e","f","g","h","i","j"/)
  resP@gsnPanelYWhiteSpacePercent = 4.0
  resP@gsnPanelXWhiteSpacePercent = 4.0

  resP@txString   = "dEdt budget CAM: instantaneous output"
  gsn_panel(wks,plot,(/3,2/),resP)               ; now draw as one plot
  txres@txFontHeightF = 0.015
  str = "- dE/dt EFIX = - flux passed to MOM6 (a) + dE/dt DME (b) + dE/dt dycore (c) + dE/dt coupling errors"
  gsn_text_ndc(wks, str, 0.5, 0.3, txres)
  dEdt_coupling_errors = -global_ave_enth_flux_to_ocn+global_ave_TE_dEdt_dme_adjust_dynE+global_ave_dycore+global_ave_efix
  str = "dE/dt coupling errors ="+dEdt_coupling_errors+" W/m2"
  gsn_text_ndc(wks, str, 0.5, 0.25, txres)
  txres               = True
  frame(wks)
  delete(plot)
  plot = new(12,graphic) 
; ***********************************************
; *NEWPAGENEWPAGENEWPAGENEWPAGENEWPAGENEWPAGE
; *NEWPAGENEWPAGENEWPAGENEWPAGENEWPAGENEWPAGE
; ***********************************************
; ***********************************************
  filePlot = "mom6-flux-to-coupler-from-CAM."+tt
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
  plot = new(9,graphic)
  
; ***********************************************
; 
; ***********************************************

  fld = hsnow_liq_ref(tidx,0:ncol_p-1)
  str = "hsnow_liq_ref (OCN only)"
  fmin = -2.
  fmax = 2.

  res_p@cnLevelSelectionMode = "ManualLevels"    ; manually set the contour levels with the following 3 resources
  res_p@cnMinLevelValF  = fmin                        ; set the minimum contour level
  res_p@cnMaxLevelValF  = fmax                   ; set the maximum contour level
  res_p@cnLevelSpacingF = (fmax-fmin)/50.0                  ; set the interval between contours

  min_f = min(fld)
  max_f = max(fld)
  res_p@lbTitleString    = "global min = "+sprintf("%6.4g",min_f)+"  global  max="+sprintf("%6.4g",max_f)
  print(str+" min/max: "+min_f+"  "+max_f)

  res_p@sfXArray            = lon1d
  res_p@sfYArray            = lat1d

  res_p@gsnLeftString        = str
  global_ave = sum(area*fld)/area_sphere_ocn
  print(str+" global ave: "+global_ave)
  res_p@gsnRightString        = "mean: "+sprintf("%6.4g", global_ave        )+" W/m^2"
  plot(0) = gsn_csm_contour_map_ce(wks,fld,res_p)

; ***********************************************
; 
; ***********************************************

  fld = hrain_liq_ref(tidx,0:ncol_p-1)
  str = "hrain_liq_ref (OCN only)"
  fmin = -25.
  fmax = 25.

  res_p@cnLevelSelectionMode = "ManualLevels"    ; manually set the contour levels with the following 3 resources
  res_p@cnMinLevelValF  = fmin                        ; set the minimum contour level
  res_p@cnMaxLevelValF  = fmax                   ; set the maximum contour level
  res_p@cnLevelSpacingF = (fmax-fmin)/50.0                  ; set the interval between contours

  min_f = min(fld)
  max_f = max(fld)
  res_p@lbTitleString    = "global min = "+sprintf("%6.4g",min_f)+"  global  max="+sprintf("%6.4g",max_f)
  print(str+" min/max: "+min_f+"  "+max_f)

  res_p@sfXArray            = lon1d
  res_p@sfYArray            = lat1d

  res_p@gsnLeftString        = str
  global_ave = sum(area*fld)/area_sphere_ocn
  print(str+" global ave: "+global_ave)
  res_p@gsnRightString        = "mean: "+sprintf("%6.4g", global_ave        )+" W/m^2"
  plot(1) = gsn_csm_contour_map_ce(wks,fld,res_p)

; ***********************************************
; 
; ***********************************************

  fld = hevap_liq_ref(tidx,0:ncol_p-1)
  str = "hevap_liq_ref (OCN only)"
  fmin = -25.
  fmax = 25.

  res_p@cnLevelSelectionMode = "ManualLevels"    ; manually set the contour levels with the following 3 resources
  res_p@cnMinLevelValF  = fmin                        ; set the minimum contour level
  res_p@cnMaxLevelValF  = fmax                   ; set the maximum contour level
  res_p@cnLevelSpacingF = (fmax-fmin)/50.0                  ; set the interval between contours

  min_f = min(fld)
  max_f = max(fld)
  res_p@lbTitleString    = "global min = "+sprintf("%6.4g",min_f)+"  global  max="+sprintf("%6.4g",max_f)
  print(str+" min/max: "+min_f+"  "+max_f)

  res_p@sfXArray            = lon1d
  res_p@sfYArray            = lat1d

  res_p@gsnLeftString        = str
  global_ave = sum(area*fld)/area_sphere_ocn
  print(str+" global ave: "+global_ave)
  res_p@gsnRightString        = "mean: "+sprintf("%6.4g", global_ave        )+" W/m^2"
  plot(2) = gsn_csm_contour_map_ce(wks,fld,res_p)

; ***********************************************
; 
; ***********************************************

  fld = hevap_liq_ref(tidx,0:ncol_p-1)+hrain_liq_ref(tidx,0:ncol_p-1)+hsnow_liq_ref(tidx,0:ncol_p-1)
  str = "(a)+(b)+(c) (OCN only)"
  fmin = -25.
  fmax = 25.

  res_p@cnLevelSelectionMode = "ManualLevels"    ; manually set the contour levels with the following 3 resources
  res_p@cnMinLevelValF  = fmin                        ; set the minimum contour level
  res_p@cnMaxLevelValF  = fmax                   ; set the maximum contour level
  res_p@cnLevelSpacingF = (fmax-fmin)/50.0                  ; set the interval between contours

  min_f = min(fld)
  max_f = max(fld)
  res_p@lbTitleString    = "global min = "+sprintf("%6.4g",min_f)+"  global  max="+sprintf("%6.4g",max_f)
  print(str+" min/max: "+min_f+"  "+max_f)

  res_p@sfXArray            = lon1d
  res_p@sfYArray            = lat1d

  res_p@gsnLeftString        = str
  global_ave = sum(area*fld)/area_sphere_ocn
  print(str+" global ave: "+global_ave)
  res_p@gsnRightString        = "mean: "+sprintf("%6.4g", global_ave        )+" W/m^2"
  plot(3) = gsn_csm_contour_map_ce(wks,fld,res_p)



; ***********************************************
; Plot everything
; ***********************************************
  resP                  = True                   ; modify the panel plot
  resP@gsnFrame         = False                  ; don't advance panel plot
  resP@gsnPanelFigureStrings = (/"a","b","c","d","e","f","g","h","i","j"/)
  resP@gsnPanelYWhiteSpacePercent = 4.0
  resP@gsnPanelXWhiteSpacePercent = 4.0

  resP@txString   = "MOM6 enthalpy flux to coupler from CAM (means using ocn area only)"


  gsn_panel(wks,plot,(/3,2/),resP)               ; now draw as one plot

  txres               = True
  frame(wks)


  delete(plot)
  plot = new(9,graphic)
; ***********************************************
; *NEWPAGENEWPAGENEWPAGENEWPAGENEWPAGENEWPAGE
; *NEWPAGENEWPAGENEWPAGENEWPAGENEWPAGENEWPAGE
; ***********************************************
; ***********************************************
  filePlot = "ocnfrac."+tt
  wks = gsn_open_wks(output_format,filePlot)  
; ***********************************************
; 
; ***********************************************
  delete(fld)
  fld=OCNFRAC(tidx,:)
  str = "OCNFRAC"
  fmin = 0.
  fmax = 1.

  res_p@cnLevelSelectionMode = "ManualLevels"    ; manually set the contour levels with the following 3 resources
  res_p@cnMinLevelValF  = fmin                        ; set the minimum contour level
  res_p@cnMaxLevelValF  = fmax                   ; set the maximum contour level
  res_p@cnLevelSpacingF = (fmax-fmin)/50.0                  ; set the interval between contours

  min_f = min(fld)
  max_f = max(fld)
  res_p@lbTitleString    = "global min = "+sprintf("%6.4g",min_f)+"  global  max="+sprintf("%6.4g",max_f)
  print(str+" min/max: "+min_f+"  "+max_f)

  res_p@sfXArray            = lon1d
  res_p@sfYArray            = lat1d

  res_p@gsnLeftString        = str
  global_ave = sum(area*fld)/area_sphere
  global_ave_enth_flux_to_ocn = global_ave
  print(str+" global ave: "+global_ave)
  res_p@gsnRightString        = "mean: "+sprintf("%6.4g", global_ave        )+" W/m^2"
  plot(0) = gsn_csm_contour_map_ce(wks,fld,res_p)

; ***********************************************
; Plot everything
; ***********************************************
  resP                  = True                   ; modify the panel plot
  resP@gsnFrame         = False                  ; don't advance panel plot
  resP@gsnPanelFigureStrings = (/"a","b","c","d","e","f","g","h","i","j"/)
  resP@gsnPanelYWhiteSpacePercent = 4.0
  resP@gsnPanelXWhiteSpacePercent = 4.0

  resP@txString   = "Ocean fraction"
  gsn_panel(wks,plot,(/1,1/),resP)               ; now draw as one plot

  txres               = True
  frame(wks)

  delete(plot)
  plot = new(9,graphic)
; ***********************************************
; *NEWPAGENEWPAGENEWPAGENEWPAGENEWPAGENEWPAGE
; *NEWPAGENEWPAGENEWPAGENEWPAGENEWPAGENEWPAGE
; ***********************************************
; ***********************************************
  filePlot = "dycore-budget."+tt
  wks = gsn_open_wks(output_format,filePlot)  
; ***********************************************
; 
; ***********************************************
  delete(fld)
  TE_dyAM= fall_energy[:]->TE_dyAM
  TE_dyBF= fall_energy[:]->TE_dyBF
  fld=(TE_dyAM(tidx,:)-TE_dyBF(tidx,:))/1800.0D
  str = "Physics tendency (efix+param+dme) in physics"
  fmin = -25.
  fmax = 25.

  res_p@cnLevelSelectionMode = "ManualLevels"    ; manually set the contour levels with the following 3 resources
  res_p@cnMinLevelValF  = fmin                        ; set the minimum contour level
  res_p@cnMaxLevelValF  = fmax                   ; set the maximum contour level
  res_p@cnLevelSpacingF = (fmax-fmin)/50.0                  ; set the interval between contours

  min_f = min(fld)
  max_f = max(fld)
  res_p@lbTitleString    = "global min = "+sprintf("%6.4g",min_f)+"  global  max="+sprintf("%6.4g",max_f)
  print(str+" min/max: "+min_f+"  "+max_f)

  res_p@sfXArray            = lon1d
  res_p@sfYArray            = lat1d

  res_p@gsnLeftString        = str
  global_ave = sum(area*fld)/area_sphere
  global_ave_enth_flux_to_ocn = global_ave
  print(str+" global ave: "+global_ave)
  res_p@gsnRightString        = "mean: "+sprintf("%6.4g", global_ave        )+" W/m^2"
  plot(0) = gsn_csm_contour_map_ce(wks,fld,res_p)

; ***********************************************
; 
; ***********************************************

  delete(fld)
;  fld = TE_dEdt_phys_tot_in_dyn(tidx,:)/1800.0D
  TE_dBD= fall_energy[:]->TE_dBD
  TE_dAF= fall_energy[:]->TE_dAF
  fld=(TE_dBD(tidx,:)-TE_dAF(tidx,:))/1800.0D
  str = "Physics tendency in dynamics"
  fmin = -25.
  fmax = 25.

  res_p@cnLevelSelectionMode = "ManualLevels"    ; manually set the contour levels with the following 3 resources
  res_p@cnMinLevelValF  = fmin                        ; set the minimum contour level
  res_p@cnMaxLevelValF  = fmax                   ; set the maximum contour level
  res_p@cnLevelSpacingF = (fmax-fmin)/50.0                  ; set the interval between contours

  min_f = min(fld)
  max_f = max(fld)
  res_p@lbTitleString    = "global min = "+sprintf("%6.4g",min_f)+"  global  max="+sprintf("%6.4g",max_f)
  print(str+" min/max: "+min_f+"  "+max_f)

  delete(res_p@sfXArray)
  res_p@sfXArray            = lon1d_d
  delete(res_p@sfYArray)
  res_p@sfYArray            = lat1d_d

  res_p@gsnLeftString        = str
  global_ave = sum(area_d*fld)/area_sphere
  global_ave_enth_flux_to_ocn = global_ave
  print(str+" global ave: "+global_ave)
  res_p@gsnRightString        = "mean: "+sprintf("%6.4g", global_ave        )+" W/m^2"
  plot(1) = gsn_csm_contour_map_ce(wks,fld,res_p)


; ***********************************************
; 
; ***********************************************

  delete(fld)
  fld = TE_dEdt_dme_adjust_dynE(tidx,:)/1800.0D
  str = "TE_dEdt_dme_adjust_dynE"
  fmin = -25.
  fmax = 25.

  res_p@cnLevelSelectionMode = "ManualLevels"    ; manually set the contour levels with the following 3 resources
  res_p@cnMinLevelValF  = fmin                        ; set the minimum contour level
  res_p@cnMaxLevelValF  = fmax                   ; set the maximum contour level
  res_p@cnLevelSpacingF = (fmax-fmin)/50.0                  ; set the interval between contours

  delete(min_f)
  min_f = min(fld)
  delete(max_f)
  max_f = max(fld)
  res_p@lbTitleString    = "global min = "+sprintf("%0.2f",min_f)+"  global  max="+sprintf("%6.4g",max_f)
  print(str+" min/max: "+min_f+"  "+max_f)

  delete(res_p@sfXArray)
  res_p@sfXArray            = lon1d
  delete(res_p@sfYArray)
  res_p@sfYArray            = lat1d

  res_p@gsnLeftString        = str
  global_ave = sum(area*fld)/area_sphere
  global_ave_TE_dEdt_dme_adjust_dynE=global_ave
  print(str+" global ave: "+global_ave)
  res_p@gsnRightString        = "mean: "+sprintf("%4.3f", global_ave        )+" W/m^2"
  plot(2) = gsn_csm_contour_map_ce(wks,fld,res_p)

 

; ***********************************************
; 
; ***********************************************

  fldE = TE_dEdt_dycore(tidx,:)/1800.0
  str = "TE_dEdt_dycore"
  fmin = -100.
  fmax = 100.

  res_p@cnLevelSelectionMode = "ManualLevels"    ; manually set the contour levels with the following 3 resources
  res_p@cnMinLevelValF  = fmin                        ; set the minimum contour level
  res_p@cnMaxLevelValF  = fmax                   ; set the maximum contour level
  res_p@cnLevelSpacingF = (fmax-fmin)/50.0                  ; set the interval between contours

  delete(min_f)
  delete(max_f)

  min_f = min(fldE)
  max_f = max(fldE)

  res_p@lbTitleString    = "global min = "+sprintf("%6.4g",min_f)+"  global  max="+sprintf("%6.4g",max_f)
  print(str+" min/max: "+min_f+"  "+max_f)

  delete(res_p@sfXArray)
  delete(res_p@sfYArray)

  res_p@sfXArray            = lon1d_d
  res_p@sfYArray            = lat1d_d

  res_p@gsnLeftString        = str
  global_ave = sum(area_d*fldE)/area_sphere
  global_ave_dycore = global_ave
  print(str+" global ave: "+global_ave)
  res_p@gsnRightString        = "mean: "+sprintf("%6.4g", global_ave        )+" W/m^2"
  plot(3) = gsn_csm_contour_map_ce(wks,fldE,res_p)

; ***********************************************
; 
; ***********************************************

  fld = EFIX(tidx,0:ncol_p-1)
  str = "Energy fixer"
  fmin = -25.
  fmax = 25.

  res_p@cnLevelSelectionMode = "ManualLevels"    ; manually set the contour levels with the following 3 resources
  res_p@cnMinLevelValF  = fmin                        ; set the minimum contour level
  res_p@cnMaxLevelValF  = fmax                   ; set the maximum contour level
  res_p@cnLevelSpacingF = (fmax-fmin)/50.0                  ; set the interval between contours

  min_f = min(fld)
  max_f = max(fld)
  res_p@lbTitleString    = "global min = "+sprintf("%6.4g",min_f)+"  global  max="+sprintf("%6.4g",max_f)
  print(str+" min/max: "+min_f+"  "+max_f)

  delete(res_p@sfXArray)
  delete(res_p@sfYArray)
  res_p@sfXArray            = lon1d
  res_p@sfYArray            = lat1d

  res_p@gsnLeftString        = str
  global_ave = sum(area*fld)/area_sphere
  global_ave_efix = global_ave
  print(str+" global ave: "+global_ave)
  res_p@gsnRightString        = "mean: "+sprintf("%6.4g", global_ave        )+" W/m^2"
  plot(4) = gsn_csm_contour_map_ce(wks,fld,res_p)




; ***********************************************
; Plot everything
; ***********************************************
  resP                  = True                   ; modify the panel plot
  resP@gsnFrame         = False                  ; don't advance panel plot
  resP@gsnPanelFigureStrings = (/"a","b","c","d","e","f","g","h","i","j"/)
  resP@gsnPanelYWhiteSpacePercent = 4.0
  resP@gsnPanelXWhiteSpacePercent = 4.0

  resP@txString   = "Dycore budget"
  gsn_panel(wks,plot,(/3,2/),resP)               ; now draw as one plot
  txres@txFontHeightF = 0.015
  str = "- dE/dt EFIX = - flux passed to MOM6 (a) + dE/dt DME (b) + dE/dt dycore (c) + dE/dt coupling errors"
  gsn_text_ndc(wks, str, 0.5, 0.3, txres)
  dEdt_coupling_errors = -global_ave_enth_flux_to_ocn+global_ave_TE_dEdt_dme_adjust_dynE+global_ave_dycore+global_ave_efix
  str = "dE/dt coupling errors ="+dEdt_coupling_errors+" W/m2"
  gsn_text_ndc(wks, str, 0.5, 0.25, txres)
  txres               = True
  frame(wks)


  delete(plot)
  plot = new(9,graphic)
; ***********************************************
; *NEWPAGENEWPAGENEWPAGENEWPAGENEWPAGENEWPAGE
; *NEWPAGENEWPAGENEWPAGENEWPAGENEWPAGENEWPAGE
; ***********************************************
; ***********************************************
  filePlot = "physics-budget."+tt
  wks = gsn_open_wks(output_format,filePlot)  


; ***********************************************
; 
; ***********************************************
  delete(fld)
  fld = EFIX(tidx,:)
  str = "Energy fixer (dyBP-dyBF)"
  fmin = -25.
  fmax = 25.

  res_p@cnLevelSelectionMode = "ManualLevels"    ; manually set the contour levels with the following 3 resources
  res_p@cnMinLevelValF  = fmin                        ; set the minimum contour level
  res_p@cnMaxLevelValF  = fmax                   ; set the maximum contour level
  res_p@cnLevelSpacingF = (fmax-fmin)/50.0                  ; set the interval between contours

  min_f = min(fld)
  max_f = max(fld)
  res_p@lbTitleString    = "global min = "+sprintf("%6.4g",min_f)+"  global  max="+sprintf("%6.4g",max_f)
  print(str+" min/max: "+min_f+"  "+max_f)

  delete(res_p@sfXArray)
  delete(res_p@sfYArray)
  res_p@sfXArray            = lon1d
  res_p@sfYArray            = lat1d

  res_p@gsnLeftString        = str
  global_ave = sum(area*fld)/area_sphere
  global_ave_efix = global_ave
  print(str+" global ave: "+global_ave)
  res_p@gsnRightString        = "mean: "+sprintf("%6.4g", global_ave        )+" W/m^2"
  plot(0) = gsn_csm_contour_map_ce(wks,fld,res_p)

; ***********************************************
; 
; ***********************************************
  delete(fld)
  TE_dyAP= fall_energy[:]->TE_dyAP
  TE_dyBP= fall_energy[:]->TE_dyBP
  fld=(TE_dyAP(tidx,:)-TE_dyBP(tidx,:))/1800.0D
  str = "Parameterization: dyAP-dyBP"
  fmin = -25.
  fmax = 25.

  res_p@cnLevelSelectionMode = "ManualLevels"    ; manually set the contour levels with the following 3 resources
  res_p@cnMinLevelValF  = fmin                        ; set the minimum contour level
  res_p@cnMaxLevelValF  = fmax                   ; set the maximum contour level
  res_p@cnLevelSpacingF = (fmax-fmin)/50.0                  ; set the interval between contours

  min_f = min(fld)
  max_f = max(fld)
  res_p@lbTitleString    = "global min = "+sprintf("%6.4g",min_f)+"  global  max="+sprintf("%6.4g",max_f)
  print(str+" min/max: "+min_f+"  "+max_f)

  delete(res_p@sfXArray)
  delete(res_p@sfYArray)
  res_p@sfXArray            = lon1d
  res_p@sfYArray            = lat1d

  res_p@gsnLeftString        = str
  global_ave = sum(area*fld)/area_sphere
  global_ave_param = global_ave
  print(str+" global ave: "+global_ave)
  res_p@gsnRightString        = "mean: "+sprintf("%6.4g", global_ave        )+" W/m^2"
  plot(1) = gsn_csm_contour_map_ce(wks,fld,res_p)


; ***********************************************
; 
; ***********************************************
  delete(fld)
  TE_dyAM= fall_energy[:]->TE_dyAM
  TE_dyAP= fall_energy[:]->TE_dyAP
  fld=(TE_dyAM(tidx,:)-TE_dyAP(tidx,:))/1800.0D

;  SE_dyAM= fall_energy[:]->SE_dyAM
;  SE_dyAP= fall_energy[:]->SE_dyAP
;  fld=(SE_dyAM(tidx,:)-SE_dyAP(tidx,:))/1800.0D
  
  str = "dme+cp (dyAM-dyAP)"
  fmin = -25.
  fmax = 25.

  res_p@cnLevelSelectionMode = "ManualLevels"    ; manually set the contour levels with the following 3 resources
  res_p@cnMinLevelValF  = fmin                        ; set the minimum contour level
  res_p@cnMaxLevelValF  = fmax                   ; set the maximum contour level
  res_p@cnLevelSpacingF = (fmax-fmin)/50.0                  ; set the interval between contours

  min_f = min(fld)
  max_f = max(fld)
  res_p@lbTitleString    = "global min = "+sprintf("%6.4g",min_f)+"  global  max="+sprintf("%6.4g",max_f)
  print(str+" min/max: "+min_f+"  "+max_f)

  delete(res_p@sfXArray)
  delete(res_p@sfYArray)
  res_p@sfXArray            = lon1d
  res_p@sfYArray            = lat1d

  res_p@gsnLeftString        = str
  global_ave = sum(area*fld)/area_sphere
  global_ave_dme_cp = global_ave
  print(str+" global ave: "+global_ave)
  res_p@gsnRightString        = "mean: "+sprintf("%6.4g", global_ave        )+" W/m^2"
  plot(2) = gsn_csm_contour_map_ce(wks,fld,res_p)

; ***********************************************
; 
; ***********************************************
  delete(fld)
  TE_dyAM= fall_energy[:]->TE_dyAM
  TE_dyBF= fall_energy[:]->TE_dyBF
  fld=(TE_dyAM(tidx,:)-TE_dyBF(tidx,:))/1800.0D
  
  str = "(a)+(b)+(c) (dyAM-dyBF)"
  fmin = -25.
  fmax = 25.

  res_p@cnLevelSelectionMode = "ManualLevels"    ; manually set the contour levels with the following 3 resources
  res_p@cnMinLevelValF  = fmin                        ; set the minimum contour level
  res_p@cnMaxLevelValF  = fmax                   ; set the maximum contour level
  res_p@cnLevelSpacingF = (fmax-fmin)/50.0                  ; set the interval between contours

  min_f = min(fld)
  max_f = max(fld)
  res_p@lbTitleString    = "global min = "+sprintf("%6.4g",min_f)+"  global  max="+sprintf("%6.4g",max_f)
  print(str+" min/max: "+min_f+"  "+max_f)

  delete(res_p@sfXArray)
  delete(res_p@sfYArray)
  res_p@sfXArray            = lon1d
  res_p@sfYArray            = lat1d

  res_p@gsnLeftString        = str
  global_ave = sum(area*fld)/area_sphere
  global_ave_am_bf = global_ave
  print(str+" global ave: "+global_ave)
  res_p@gsnRightString        = "mean: "+sprintf("%6.4g", global_ave        )+" W/m^2"
  plot(3) = gsn_csm_contour_map_ce(wks,fld,res_p)
; ***********************************************
; 
; ***********************************************
;  fld = enth_flux_to_ocn(tidx,0:ncol_p-1)
;  str = "Enthalpy flux to ocean"
;  fmin = -25.
;  fmax = 25.
;
;  res_p@cnLevelSelectionMode = "ManualLevels"    ; manually set the contour levels with the following 3 resources
;  res_p@cnMinLevelValF  = fmin                        ; set the minimum contour level
;  res_p@cnMaxLevelValF  = fmax                   ; set the maximum contour level
;  res_p@cnLevelSpacingF = (fmax-fmin)/50.0                  ; set the interval between contours
;
;  min_f = min(fld)
;  max_f = max(fld)
;  res_p@lbTitleString    = "global min = "+sprintf("%6.4g",min_f)+"  global  max="+sprintf("%6.4g",max_f)
;  print(str+" min/max: "+min_f+"  "+max_f);
;
;  res_p@sfXArray            = lon1d
;  res_p@sfYArray            = lat1d;
;
;  res_p@gsnLeftString        = str
;  global_ave = sum(area*fld)/area_sphere
;  global_ave_enth_flux_to_ocn = global_ave
;  print(str+" global ave: "+global_ave)
;  res_p@gsnRightString        = "mean: "+sprintf("%6.4g", global_ave        )+" W/m^2"
;  plot(4) = gsn_csm_contour_map_ce(wks,fld,res_p)

; ***********************************************
; 
; ***********************************************
;  delete(fld)
;  fld=dEdt_cpdycore(tidx,:)
;  str = "cp_dycore: enthalpy energy only"
;  fmin = -25.
;  fmax = 25.
;
;  res_p@cnLevelSelectionMode = "ManualLevels"    ; manually set the contour levels with the following 3 resources
;  res_p@cnMinLevelValF  = fmin                        ; set the minimum contour level
;  res_p@cnMaxLevelValF  = fmax                   ; set the maximum contour level
;  res_p@cnLevelSpacingF = (fmax-fmin)/50.0                  ; set the interval between contours
;
;  min_f = min(fld)
;  max_f = max(fld)
;  res_p@lbTitleString    = "global min = "+sprintf("%6.4g",min_f)+"  global  max="+sprintf("%6.4g",max_f)
;  print(str+" min/max: "+min_f+"  "+max_f)
;
;  delete(res_p@sfXArray)
;  delete(res_p@sfYArray)
;  res_p@sfXArray            = lon1d
;  res_p@sfYArray            = lat1d
;
;  res_p@gsnLeftString        = str
;  global_ave = sum(area*fld)/area_sphere
;  global_ave_efix = global_ave
;  print(str+" global ave: "+global_ave)
;  res_p@gsnRightString        = "mean: "+sprintf("%6.4g", global_ave        )+" W/m^2"
;  plot(3) = gsn_csm_contour_map_ce(wks,fld,res_p)


; ***********************************************
; Plot everything
; ***********************************************
  resP                  = True                   ; modify the panel plot
  resP@gsnFrame         = False                  ; don't advance panel plot
  resP@gsnPanelFigureStrings = (/"a","b","c","d","e","f","g","h","i","j"/)
  resP@gsnPanelYWhiteSpacePercent = 4.0
  resP@gsnPanelXWhiteSpacePercent = 4.0

  resP@txString   = "Physics budget"
  gsn_panel(wks,plot,(/3,2/),resP)               ; now draw as one plot
  txres@txFontHeightF = 0.015
  str = "Check if (a)+(b)+(c) is the same as (E(dyAM)-E(dyBF)/dt: "
  gsn_text_ndc(wks, str, 0.5, 0.3, txres)
;  dEdt_coupling_errors = -global_ave_enth_flux_to_ocn+global_ave_TE_dEdt_dme_adjust_dynE+global_ave_dycore+global_ave_efix
;  str = "dE/dt coupling errors ="+dEdt_coupling_errors+" W/m2"

  str=global_ave_efix+global_ave_param +global_ave_dme_cp-global_ave_am_bf

  gsn_text_ndc(wks, str, 0.5, 0.25, txres)
  txres               = True
  frame(wks)

end