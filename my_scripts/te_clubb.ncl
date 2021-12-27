load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

begin

  fname_ref =  systemfunc("ls "+fname_ref)
  fall_ref = addfile(fname_ref,"r")       ; Open netCDF files.  

;  te_leak_ref = fall_ref->ELEAK_CLUBB(:,:)
;  ke_leak_ref = fall_ref->KLEAK_CLUBB(:,:)
;  se_leak_ref = fall_ref->SLEAK_CLUBB(:,:)

  te_leak_ref = fall_ref->PRECC(:,:)+fall_ref->PRECL(:,:)
  ke_leak_ref = fall_ref->PRECC(:,:)+fall_ref->PRECL(:,:)
  se_leak_ref = fall_ref->PRECC(:,:)+fall_ref->PRECL(:,:) 
  te_leak_ref = 1000*24*3600*te_leak_ref ;convert from m/s to mm/day  




  fname_flx =  systemfunc("ls "+fname_flx)
  fall_flx = addfile(fname_flx,"r")       ; Open netCDF files.  

  te_leak_flx = fall_flx->PRECC(:,:)+fall_flx->PRECL(:,:)
  ke_leak_flx = fall_flx->PRECC(:,:)+fall_flx->PRECL(:,:)
  se_leak_flx = fall_flx->PRECC(:,:)+fall_flx->PRECL(:,:)

  te_leak_flx = 1000*24*3600*te_leak_flx ;convert from m/s to mm/day  

;  te_leak_flx = fall_flx->ELEAK_CLUBB(:,:)
;  ke_leak_flx = fall_flx->KLEAK_CLUBB(:,:)
;  se_leak_flx = fall_flx->SLEAK_CLUBB(:,:)



 
  lat1d         = fall_ref->lat          ; Pull off lat/lon variable and
  lon1d         = fall_ref->lon          ; convert both to 1D.
  area          = fall_ref->area(:)
  area_sphere   = sum(area)


  output_format = "pdf"
  filePlot = "clubb"
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
; 
; ***********************************************

  fld = te_leak_ref(0,:)
  str = "TE imbalance /w frictional heat"
  fmin = -5.
  fmax = 5.

;  res_p@cnLevelSelectionMode = "ManualLevels"    ; manually set the contour levels with the following 3 resources
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
  print(str+" global ave: "+global_ave)
  res_p@gsnRightString        = "mean: "+sprintf("%3.2g", global_ave        )+" W/m^2"
  plot(0) = gsn_csm_contour_map_ce(wks,fld,res_p)


; ***********************************************
; 
; ***********************************************

  fld = te_leak_flx(0,:)
  str = "TE imbalance /w frictional heat & flux scaling"
  fmin = -5.
  fmax = 5.

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
  print(str+" global ave: "+global_ave)
  res_p@gsnRightString        = "mean: "+sprintf("%3.2g", global_ave        )+" W/m^2"
  plot(1) = gsn_csm_contour_map_ce(wks,fld,res_p)

; ***********************************************
; 
; ***********************************************



; ***********************************************
; 
; ***********************************************

  fld = te_leak_flx(0,:)-te_leak_ref(0,:)
  str = "Difference: (b)-(a)"
  fmin = -5.
  fmax = 5.

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
  print(str+" global ave: "+global_ave)
  res_p@gsnRightString        = "mean: "+sprintf("%3.2g", global_ave        )+" W/m^2"
  plot(2) = gsn_csm_contour_map_ce(wks,fld,res_p)





  resP                  = True                   ; modify the panel plot
  resP@gsnFrame         = False                  ; don't advance panel plot
  resP@gsnPanelFigureStrings = (/"a","b","c","d","e","f","g"/)
  resP@gsnPanelYWhiteSpacePercent = 4.0
  resP@gsnPanelXWhiteSpacePercent = 4.0

  resP@txString   = "CLUBB sensible heat flux consistency experiments"
  gsn_panel(wks,plot,(/1,3/),resP)               ; now draw as one plot

  txres               = True
  frame(wks)

end