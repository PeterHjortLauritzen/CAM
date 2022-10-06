load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

begin
  fnamez =  systemfunc("ls "+fnamez)

  fallz = addfile(fnamez,"r")       ; Open netCDF files.  


  print("File used for energy diagnostics "+fnamez)
  print("done input summary")
  print(" ")  
  print("Summary of input to ncl script:")
  print("===============================")
  print(" ")

  print("Your data is for the MPAS dycore")  
  tmp       = fallz->zmid_mpas(:,0,:)
  nt        = dimsizes(tmp(:,0))
  ncol      = dimsizes(tmp(0,:))

  ;
  ; initialization
  ;
  g         = 9.80665D0
  PI        = get_pi("double")
  deg2rad   = PI/180.0
  eps       = 1.0D-13

  area      = fallz->area   (:  )
  area_sphere = sum(area(0:ncol-1))      
  accuracy_area_sphere = (area_sphere-4D0*PI)/(4D0*PI)
  print("Area of sphere in dataset "+area_sphere)
  print("Relative error            "+accuracy_area_sphere)
  inv_area_sphere   = 1.0D0/area_sphere
 

 
  print("start plotting")

  output_format = "pdf"
  filePlot = "z2d"
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


  PHIS   = fallz->PHIS(0,:)    ; Read some data; convert to 1D.        
  zs     = PHIS/9.80665
  delete(PHIS)

  ; 
  ; get heights over ocean (MPAS)
  ;
  lon=fallz->lon
  lat=fallz->lat
  dist=2.0D0
  lon_point = 277.0D0
  lat_point = 7.0D0
  do i=0,ncol-1
    sq = (lon(i)-lon_point)^2+(lat(i)-lat_point)^2
    if (sq.lt.dist) then
      icol = i
      dist = sq
    end if
  end do
  print("Point is at index "+i)
  print("Point location (lon,lat): "+lon(icol)+" "+lat(icol))

  zm_mpas = fallz->zmid_mpas(time_idx,:,icol)
  nlev = dimsizes(zm_mpas)
  do k=0,nlev-1  
    zm_mpas(k)   = zm_mpas(k);-zs(icol)
    print("zm_mpas("+k+")="+zm_mpas(k))
  end do

; ***********************************************
; 
; ***********************************************
  
  lev = 30
  zm_phBF = fallz->zm_dyBF(time_idx,lev,:)
  zmid_mpas = fallz->zmid_mpas(time_idx,lev,:)
  fld = zm_phBF-(zmid_mpas)

  str = "z~B4~k=31~S~(CAM)~N~(begin phys) - z~B4~k=31~S~(MPAS)~N~; z~B4~k=31~S~(MPAS)~N~="+zm_mpas(lev)+"m"
  fmin = -300.
  fmax = 300.

  res_p@cnLevelSelectionMode = "AutomaticLevels";ManualLevels"    ; manually set the contour levels with the following 3 resources
;  res_p@cnMinLevelValF  = fmin                        ; set the minimum contour level
;;  res_p@cnMaxLevelValF  = fmax                   ; set the maximum contour level
;;  res_p@cnLevelSpacingF = (fmax-fmin)/50.0                  ; set the interval between contours
;
  min_f = min(fld)
  max_f = max(fld)

  res_p@lbTitleString    = "global min = "+sprintf("%6.4g",min_f)+" m"+"  global  max="+sprintf("%6.4g",max_f)+" m"
  print(str+" min/max: "+min_f+"  "+max_f)

  res_p@sfXArray            = lon1d
  res_p@sfYArray            = lat1d

  res_p@gsnLeftString        = str
  global_ave = sum(area*fld)*inv_area_sphere
  print(str+" global ave: "+global_ave)
  res_p@gsnRightString        = "mean: "+sprintf("%3.2g", global_ave        )+" m"
;
  plot(0) = gsn_csm_contour_map_ce(wks,fld,res_p)

  lev = 20
  zm_phBF = fallz->zm_dyBF(time_idx,lev,:)
  zmid_mpas = fallz->zmid_mpas(time_idx,lev,:)
  fld = zm_phBF-(zmid_mpas)

  str = "z~B4~k=21~S~(CAM)~N~(begin phys) - z~B4~k=21~S~(MPAS)~N~; z~B4~k=21~S~(MPAS)~N~="+zm_mpas(lev)+"m"
  fmin = -300.
  fmax = 300.

  res_p@cnLevelSelectionMode = "AutomaticLevels";ManualLevels"    ; manually set the contour levels with the following 3 resources
;  res_p@cnMinLevelValF  = fmin                        ; set the minimum contour level
;;  res_p@cnMaxLevelValF  = fmax                   ; set the maximum contour level
;;  res_p@cnLevelSpacingF = (fmax-fmin)/50.0                  ; set the interval between contours
;
  min_f = min(fld)
  max_f = max(fld)

  res_p@lbTitleString    = "global min = "+sprintf("%6.4g",min_f)+" m"+"  global  max="+sprintf("%6.4g",max_f)+" m"
  print(str+" min/max: "+min_f+"  "+max_f)

  res_p@sfXArray            = lon1d
  res_p@sfYArray            = lat1d

  res_p@gsnLeftString        = str
  global_ave = sum(area*fld)*inv_area_sphere
  print(str+" global ave: "+global_ave)
  res_p@gsnRightString        = "mean: "+sprintf("%3.2g", global_ave        )+" m"
;
  plot(1) = gsn_csm_contour_map_ce(wks,fld,res_p)


  lev = 10
  zm_phBF = fallz->zm_phBF(time_idx,lev,:)
  zmid_mpas = fallz->zmid_mpas(time_idx,lev,:)
  fld = zm_phBF-(zmid_mpas)

  str = "z~B4~k=11~S~(CAM)~N~(begin phys) - z~B4~k=11~S~(MPAS)~N~; z~B4~k=11~S~(MPAS)~N~="+zm_mpas(lev)+"m"
  fmin = -300.
  fmax = 300.

  res_p@cnLevelSelectionMode = "AutomaticLevels";ManualLevels"    ; manually set the contour levels with the following 3 resources
;  res_p@cnMinLevelValF  = fmin                        ; set the minimum contour level
;;  res_p@cnMaxLevelValF  = fmax                   ; set the maximum contour level
;;  res_p@cnLevelSpacingF = (fmax-fmin)/50.0                  ; set the interval between contours
;
  min_f = min(fld)
  max_f = max(fld)

  res_p@lbTitleString    = "global min = "+sprintf("%6.4g",min_f)+" m"+"  global  max="+sprintf("%6.4g",max_f)+" m"
  print(str+" min/max: "+min_f+"  "+max_f)

  res_p@sfXArray            = lon1d
  res_p@sfYArray            = lat1d

  res_p@gsnLeftString        = str
  global_ave = sum(area*fld)*inv_area_sphere
  print(str+" global ave: "+global_ave)
  res_p@gsnRightString        = "mean: "+sprintf("%3.2g", global_ave        )+" m"
;
  plot(2) = gsn_csm_contour_map_ce(wks,fld,res_p)
;

; ***********************************************
; 
; ***********************************************
  
  lev = 30
  zm_phBF = fallz->zm_phBF(time_idx,lev,:)
  zm_phAP = fallz->zm_phAP(time_idx,lev,:)
  fld = zm_phAP-zm_phBF

  str = "z~B4~k=31~S~(CAM)~N~(end phys before dma) - z~B4~k=31~S~(CAM)~N~(begin phys)"
  fmin = -300.
  fmax = 300.

  res_p@cnLevelSelectionMode = "AutomaticLevels";ManualLevels"    ; manually set the contour levels with the following 3 resources
;  res_p@cnMinLevelValF  = fmin                        ; set the minimum contour level
;;  res_p@cnMaxLevelValF  = fmax                   ; set the maximum contour level
;;  res_p@cnLevelSpacingF = (fmax-fmin)/50.0                  ; set the interval between contours
;
  min_f = min(fld)
  max_f = max(fld)

  res_p@lbTitleString    = "global min = "+sprintf("%6.4g",min_f)+" m"+"  global  max="+sprintf("%6.4g",max_f)+" m"
  print(str+" min/max: "+min_f+"  "+max_f)

  res_p@sfXArray            = lon1d
  res_p@sfYArray            = lat1d

  res_p@gsnLeftString        = str
  global_ave = sum(area*fld)*inv_area_sphere
  print(str+" global ave: "+global_ave)
  res_p@gsnRightString        = "mean: "+sprintf("%3.2g", global_ave        )+" m"
;
  plot(3) = gsn_csm_contour_map_ce(wks,fld,res_p)

  lev = 20
  zm_phBF = fallz->zm_phBF(time_idx,lev,:)
  zm_phAP = fallz->zm_phAP(time_idx,lev,:)
  fld = zm_phAP-zm_phBF

  str = "z~B4~k=21~S~(CAM)~N~(end phys before dma) - z~B4~k=21~S~(CAM)~N~(begin phys)"
  fmin = -300.
  fmax = 300.

  res_p@cnLevelSelectionMode = "AutomaticLevels";ManualLevels"    ; manually set the contour levels with the following 3 resources
;  res_p@cnMinLevelValF  = fmin                        ; set the minimum contour level
;;  res_p@cnMaxLevelValF  = fmax                   ; set the maximum contour level
;;  res_p@cnLevelSpacingF = (fmax-fmin)/50.0                  ; set the interval between contours
;
  min_f = min(fld)
  max_f = max(fld)

  res_p@lbTitleString    = "global min = "+sprintf("%6.4g",min_f)+" m"+"  global  max="+sprintf("%6.4g",max_f)+" m"
  print(str+" min/max: "+min_f+"  "+max_f)

  res_p@sfXArray            = lon1d
  res_p@sfYArray            = lat1d

  res_p@gsnLeftString        = str
  global_ave = sum(area*fld)*inv_area_sphere
  print(str+" global ave: "+global_ave)
  res_p@gsnRightString        = "mean: "+sprintf("%3.2g", global_ave        )+" m"
;
  plot(4) = gsn_csm_contour_map_ce(wks,fld,res_p)


  lev = 10
  zm_phBF = fallz->zm_phBF(time_idx,lev,:)
  zm_phAP = fallz->zm_phAP(time_idx,lev,:)
  fld = zm_phAP-zm_phBF

  str = "z~B4~k=11~S~(CAM)~N~(end phys before dma) - z~B4~k=11~S~(CAM)~N~(begin phys)"
  fmin = -300.
  fmax = 300.

  res_p@cnLevelSelectionMode = "AutomaticLevels";ManualLevels"    ; manually set the contour levels with the following 3 resources
;  res_p@cnMinLevelValF  = fmin                        ; set the minimum contour level
;;  res_p@cnMaxLevelValF  = fmax                   ; set the maximum contour level
;;  res_p@cnLevelSpacingF = (fmax-fmin)/50.0                  ; set the interval between contours
;
  min_f = min(fld)
  max_f = max(fld)

  res_p@lbTitleString    = "global min = "+sprintf("%6.4g",min_f)+" m"+"  global  max="+sprintf("%6.4g",max_f)+" m"
  print(str+" min/max: "+min_f+"  "+max_f)

  res_p@sfXArray            = lon1d
  res_p@sfYArray            = lat1d

  res_p@gsnLeftString        = str
  global_ave = sum(area*fld)*inv_area_sphere
  print(str+" global ave: "+global_ave)
  res_p@gsnRightString        = "mean: "+sprintf("%3.2g", global_ave        )+" m"
;
  plot(5) = gsn_csm_contour_map_ce(wks,fld,res_p)
;
;;
; ***********************************************
; 
; ***********************************************
  
  lev = 30
  zm_phAM = fallz->zm_phAM(time_idx,lev,:)
  zm_phAP = fallz->zm_phAP(time_idx,lev,:)
  fld = zm_phAM-zm_phAP

  str = "z~B4~k=31~S~(CAM)~N~(after dma) - z~B4~k=31~S~(CAM)~N~(before dma)"
  fmin = -300.
  fmax = 300.

  res_p@cnLevelSelectionMode = "AutomaticLevels";ManualLevels"    ; manually set the contour levels with the following 3 resources
;  res_p@cnMinLevelValF  = fmin                        ; set the minimum contour level
;;  res_p@cnMaxLevelValF  = fmax                   ; set the maximum contour level
;;  res_p@cnLevelSpacingF = (fmax-fmin)/50.0                  ; set the interval between contours
;
  min_f = min(fld)
  max_f = max(fld)

  res_p@lbTitleString    = "global min = "+sprintf("%6.4g",min_f)+" m"+"  global  max="+sprintf("%6.4g",max_f)+" m"
  print(str+" min/max: "+min_f+"  "+max_f)

  res_p@sfXArray            = lon1d
  res_p@sfYArray            = lat1d

  res_p@gsnLeftString        = str
  global_ave = sum(area*fld)*inv_area_sphere
  print(str+" global ave: "+global_ave)
  res_p@gsnRightString        = "mean: "+sprintf("%3.2g", global_ave        )+" m"
;
  plot(6) = gsn_csm_contour_map_ce(wks,fld,res_p)

  lev = 20
  zm_phAM = fallz->zm_phAM(time_idx,lev,:)
  zm_phAP = fallz->zm_phAP(time_idx,lev,:)
  fld = zm_phAM-zm_phAP

  str = "z~B4~k=21~S~(CAM)~N~(end phys before dma) - z~B4~k=21~S~(CAM)~N~(after dma)"
  fmin = -300.
  fmax = 300.

  res_p@cnLevelSelectionMode = "AutomaticLevels";ManualLevels"    ; manually set the contour levels with the following 3 resources
;  res_p@cnMinLevelValF  = fmin                        ; set the minimum contour level
;;  res_p@cnMaxLevelValF  = fmax                   ; set the maximum contour level
;;  res_p@cnLevelSpacingF = (fmax-fmin)/50.0                  ; set the interval between contours
;
  min_f = min(fld)
  max_f = max(fld)

  res_p@lbTitleString    = "global min = "+sprintf("%6.4g",min_f)+" m"+"  global  max="+sprintf("%6.4g",max_f)+" m"
  print(str+" min/max: "+min_f+"  "+max_f)

  res_p@sfXArray            = lon1d
  res_p@sfYArray            = lat1d

  res_p@gsnLeftString        = str
  global_ave = sum(area*fld)*inv_area_sphere
  print(str+" global ave: "+global_ave)
  res_p@gsnRightString        = "mean: "+sprintf("%3.2g", global_ave        )+" m"
;
  plot(7) = gsn_csm_contour_map_ce(wks,fld,res_p)


  lev = 0
  ptop_mpas = fallz->ptop_mpas(time_idx,:)
  zi_phAM = fallz->zi_phAM(time_idx,lev,:)
  zi_phAP = fallz->zi_phAP(time_idx,lev,:)
;  fld = ptop_mpas*(zi_phAM-zi_phAP)/1800.0
  fld = (zi_phAM-zi_phAP)
;  fld = 150.0*(zm_phAM-zm_phAP)


;  fld = ptop_mpas


  str = "z~B4~k=1~S~(CAM)~N~(end phys before dma) - z~B4~k=1~S~(CAM)~N~(after dma)"
  fmin = -10.0
  fmax = 10.0

;  res_p@cnLevelSelectionMode = "AutomaticLevels";ManualLevels"    ; manually set the contour levels with the following 3 resources
  res_p@cnLevelSelectionMode = "ManualLevels"    ; manually set the contour levels with the following 3 resources
  res_p@cnMinLevelValF  = fmin                        ; set the minimum contour level
  res_p@cnMaxLevelValF  = fmax                   ; set the maximum contour level
  res_p@cnLevelSpacingF = (fmax-fmin)/50.0                  ; set the interval between contours
;
  min_f = min(fld)
  max_f = max(fld)

  res_p@lbTitleString    = "global min = "+sprintf("%6.4g",min_f)+" m"+"  global  max="+sprintf("%6.4g",max_f)+" m"
  print(str+" min/max: "+min_f+"  "+max_f)

  res_p@sfXArray            = lon1d
  res_p@sfYArray            = lat1d

  res_p@gsnLeftString        = str
  global_ave = sum(area*fld)*inv_area_sphere
  print(str+" global ave: "+global_ave)
  res_p@gsnRightString        = "mean: "+sprintf("%3.2g", global_ave        )+" m"
;
  plot(8) = gsn_csm_contour_map_ce(wks,fld,res_p)
;
;;

; ***********************************************
; plot panel
; ***********************************************

  resP                  = True                   ; modify the panel plot
  resP@gsnFrame         = False                  ; don't advance panel plot
  resP@gsnPanelFigureStrings = (/"a","b","c","d","e","f","g","h","i","j"/)
  resP@gsnPanelYWhiteSpacePercent = 4.0
  resP@gsnPanelXWhiteSpacePercent = 4.0

  resP@txString   = "Height (z) discrepancies and differences"
  gsn_panel(wks,plot,(/3,3/),resP)               ; now draw as one plot

  txres               = True
  frame(wks)

end