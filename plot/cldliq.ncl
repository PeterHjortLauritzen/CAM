load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
begin
  output_format="pdf"
  dir="/glade/derecho/scratch/pel/cam_enthalpy2_384_FLTHIST_ne30pg3_ne30pg3_mg17/run/"
  fname = "cam_enthalpy2_384_FLTHIST_ne30pg3_ne30pg3_mg17.cam.h3i."
  fnameall =  systemfunc("ls "+dir+fname+"*.nc")
  fall = addfiles(fnameall,"r")       ; Open netCDF files.
    
  tidx = 20

  T00 = 290D0
  cpair    = 1004.64D0 
  cpliq    = 4188.00D0
  cpwv     = 1810.00D0
  cpice    = 2117.27D0
  LS = 2834700.0
  LF =  333700.0
  dt=1800.0

;  T=fall[:]->T
;  TS=fall[:]->TS

  DCLDLIQ=fall[:]->DCLDLIQ
  DCLDICE=fall[:]->DCLDICE


  dt=1800.0
  firsttime = 0
  lasttime  = dimsizes(fall[:]->time)-1
  print("Total number of time samples "+lasttime)

  ncol        = dimsizes(DCLDLIQ(0,1,:))
  ncol_p      = ncol
  area        = fall[0]->area
  area_sphere = sum(area(0:ncol-1))
  print("area_sphere="+area_sphere)
  print("ncol="+ncol)
  lat1d        = fall[0]->lat          ; Pull off lat/lon variable and
  lon1d        = fall[0]->lon  


  filePlot = "cldliq"
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
; *NEWPAGENEWPAGENEWPAGENEWPAGENEWPAGENEWPAGE
; *NEWPAGENEWPAGENEWPAGENEWPAGENEWPAGENEWPAGE
; ***********************************************
; ***********************************************
  filePlot = "cldliq"
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
  plot = new(38,graphic)

  iplot = 0
  do k=20,57
    fld = DCLDLIQ(tidx,k,0:ncol_p-1)
    str = "DCLDLIQ level "+k
    fmin = -400.
    fmax = 400.

;    res_p@cnLevelSelectionMode = "ManualLevels"    ; manually set the contour levels with the following 3 resources
;    res_p@cnMinLevelValF  = fmin                        ; set the minimum contour level
;    res_p@cnMaxLevelValF  = fmax                   ; set the maximum contour level
;    res_p@cnLevelSpacingF = (fmax-fmin)/50.0                  ; set the interval between contours

    min_f = min(fld(0:ncol_p-1))
    max_f = max(fld(0:ncol_p-1))
    res_p@lbTitleString    = "global min = "+sprintf("%6.4g",min_f)+"  global  max="+sprintf("%6.4g",max_f)
    print(str+" min/max: "+min_f+"  "+max_f)

    res_p@sfXArray            = lon1d
    res_p@sfYArray            = lat1d

    res_p@gsnLeftString        = str
    global_ave = sum(area(0:ncol_p-1)*fld(0:ncol_p-1))/area_sphere
    print(str+" global ave: "+global_ave)
    res_p@gsnRightString        = "mean: "+sprintf("%3.2g", global_ave        )+" W/m^2"
    plot(iplot) = gsn_csm_contour_map_ce(wks,fld,res_p)
    iplot = iplot+1
  end do


  
; ***********************************************
; Plot everything
; ***********************************************
  resP                  = True                   ; modify the panel plot
  resP@gsnFrame         = False                  ; don't advance panel plot
  resP@gsnPanelFigureStrings = (/"a","b","c","d","e","f","g","h","i","j"/)
  resP@gsnPanelYWhiteSpacePercent = 4.0
  resP@gsnPanelXWhiteSpacePercent = 4.0

  resP@txString   = "Terms in energy equation"
  gsn_panel(wks,plot,(/6,7/),resP)               ; now draw as one plot

  txres               = True
  frame(wks)



end