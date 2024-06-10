load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
begin
  output_format="pdf"

  dir="/glade/derecho/scratch/pel/cam_enth_384_FLTHIST_ne30pg3_ne30pg3_mg17/run/"
  fname = "cam_enth_384_FLTHIST_ne30pg3_ne30pg3_mg17.cam.h3i."
  fnameall =  systemfunc("ls "+dir+fname+"*.nc")
  fall = addfiles(fnameall,"r")       ; Open netCDF files.
    
  fname_ediags = "cam_enth_384_FLTHIST_ne30pg3_ne30pg3_mg17.cam.h1i."
  fnameall_ediags =  systemfunc("ls "+dir+fname_ediags+"*.nc")
  print(fnameall_ediags)
  fall_ediags = addfiles(fnameall_ediags,"r")       ; Open netCDF files.

  tidx = 20

  T00 = 290D0
  cpair    = 1004.64D0 
  cpliq    = 4188.00D0
  cpwv     = 1810.00D0
  cpice    = 2117.27D0
  LS = 2834700.0
  LF =  333700.0
  dt=1800.0

  TE_dEdt_dme_adjust_physE=fall_ediags[:]->TE_dEdt_dme_adjust_physE
  TE_dEdt_dme_adjust_dynE =fall_ediags[:]->TE_dEdt_dme_adjust_dynE

  SE_dEdt_dme_adjust_physE=fall_ediags[:]->SE_dEdt_dme_adjust_physE
  SE_dEdt_dme_adjust_dynE =fall_ediags[:]->SE_dEdt_dme_adjust_dynE

  KE_dEdt_dme_adjust_physE=fall_ediags[:]->KE_dEdt_dme_adjust_physE
  KE_dEdt_dme_adjust_dynE =fall_ediags[:]->KE_dEdt_dme_adjust_dynE

  TE_phAM  =fall_ediags[:]->TE_phAM
  TE_phAP  =fall_ediags[:]->TE_phAP
  TE_phBP  =fall_ediags[:]->TE_phBP
  TE_phBF  =fall_ediags[:]->TE_phBF

  TE_dyAM  =fall_ediags[:]->TE_dyAM
    TE_dyAF  =fall_ediags[:]->TE_dyAF
  TE_dyAP  =fall_ediags[:]->TE_dyAP
  TE_dyBP  =fall_ediags[:]->TE_dyBP
  TE_dyBF  =fall_ediags[:]->TE_dyBF
  TE_dyAC  =fall_ediags[:]->TE_dyAC


  KE_phAP  =fall_ediags[:]->KE_phAP
  KE_phBP  =fall_ediags[:]->KE_phBP

  TE_dEdt_param_physE=fall_ediags[:]->TE_dEdt_param_physE
  TE_dEdt_param_efix_physE=fall_ediags[:]->TE_dEdt_param_efix_physE

  TE_dEdt_param_dynE=fall_ediags[:]->TE_dEdt_param_dynE
  TE_dEdt_param_efix_dynE=fall_ediags[:]->TE_dEdt_param_efix_dynE	

  TE_dEdt_cpupdate_dynE=fall_ediags[:]->TE_dEdt_cpupdate_dynE
  TE_dEdt_Q_LT_dynE=fall_ediags[:]->TE_dEdt_Q_LT_dynE

  T=fall[:]->T
  TS=fall[:]->TS
  enth_flx=fall[:]->enth_flx
  heating=fall[:]->heating
  imbalance=fall[:]->imbalance
  te_sen=fall[:]->te_sen
  te_lat=fall[:]->te_lat
  radiation=fall[:]->radiation
  te_tnd=fall[:]->te_tnd
  FRHS_FLX=fall[:]->FRHS_FLX
  FRHS_FLXA=fall[:]->FRHS_FLXA
  FRHS_FLXB=fall[:]->FRHS_FLXB
  FRHS_FLXC=fall[:]->FRHS_FLXC

  prect=fall[:]->PRECT
 
  FRAIN_coupler= fall[:]->FRAIN_coupler
  FSNOW_coupler= fall[:]->FSNOW_coupler
  FEVAP_coupler= fall[:]->FEVAP_coupler

  HRAIN= fall[:]->HRAIN
  HSNOW= fall[:]->HSNOW
  HEVAP= fall[:]->HEVAP

  rliqbc= fall[:]->rliqbc

  FRAIN_AC= fall[:]->FRAIN_AC
  FSNOW_AC= fall[:]->FSNOW_AC
  FEVAP_AC = fall[:]->FEVAP_AC

  FRAIN_BC= fall[:]->FRAIN_BC
  FSNOW_BC= fall[:]->FSNOW_BC
  FEVAP_BC= fall[:]->FEVAP_BC

  HRAIN_AC= fall[:]->HRAIN_AC
  HSNOW_AC= fall[:]->HSNOW_AC
  HEVAP_AC = fall[:]->HEVAP_AC

  HRAIN_BC= fall[:]->HRAIN_BC
  HSNOW_BC= fall[:]->HSNOW_BC
  HEVAP_BC= fall[:]->HEVAP_BC

  FRAIN_tot= fall[:]->FRAIN_tot
  FSNOW_tot= fall[:]->FSNOW_tot
  FEVAP_tot= fall[:]->FEVAP_tot

  dtot_wv= fall[:]->dtot_wv
  dtot_ice= fall[:]->dtot_ice
  dtot_liq= fall[:]->dtot_liq

  dtot_wv_coupler= fall[:]->dtot_wv_coupler
  dtot_ice_coupler= fall[:]->dtot_ice_coupler
  dtot_liq_coupler= fall[:]->dtot_liq_coupler


  dt=1800.0
  firsttime = 0
  lasttime  = dimsizes(fall[:]->time)-1
  print("Total number of time samples "+lasttime)

  ncol        = dimsizes(FEVAP_AC(0,:))
  ncol_p      = ncol
  area        = fall[0]->area
  area_sphere = sum(area(0:ncol-1))
  print("area_sphere="+area_sphere)
  print("ncol="+ncol)
  lat1d        = fall[0]->lat          ; Pull off lat/lon variable and
  lon1d        = fall[0]->lon  


  filePlot = "budget2"
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
  filePlot = "budget2"
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
  plot = new(15,graphic)
  
; ***********************************************
; dEdt_param_efix_physE
; ***********************************************

  fld = (TE_dyAF(tidx,0:ncol_p-1)-TE_dyBP(tidx,0:ncol_p-1))/dt
  str = "dE/dt (dyAF-dyBP): tot energy change phys no efix"
  fmin = -400.
  fmax = 400.

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
  global_ave = sum(area(0:ncol_p-1)*fld(0:ncol_p-1))/area_sphere
  print(str+" global ave: "+global_ave)
  res_p@gsnRightString        = "mean: "+sprintf("%3.2g", global_ave        )+" W/m^2"
  plot(0) = gsn_csm_contour_map_ce(wks,fld,res_p)



; ***********************************************
; 
; ***********************************************

  fld = te_tnd(tidx,0:ncol_p-1)
  str = "te_tnd"
  fmin = -400.
  fmax = 400.

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
  global_ave = sum(area(0:ncol_p-1)*fld(0:ncol_p-1))/area_sphere
  print(str+" global ave: "+global_ave)
  res_p@gsnRightString        = "mean: "+sprintf("%3.2g", global_ave        )+" W/m^2"
  plot(1) = gsn_csm_contour_map_ce(wks,fld,res_p)




; ***********************************************
; 
; ***********************************************
  enthalpy_flux = enth_flx(tidx,0:ncol_p-1)
  fld = enthalpy_flux
  str = "enth_flx"

  fmin = -400.0
  fmax = 400.0

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
  global_ave = sum(area(0:ncol_p-1)*fld(0:ncol_p-1))/area_sphere
  print(str+" global ave: "+global_ave)
  res_p@gsnRightString        = "mean: "+sprintf("%3.2g", global_ave        )+" W/m^2"
  plot(2) = gsn_csm_contour_map_ce(wks,fld,res_p)


; ***********************************************
; 
; ***********************************************

  fld = heating(tidx,0:ncol_p-1)/dt
  str = "heating"
;  fld = (TE_phAM(tidx,0:ncol_p-1)-TE_phBF(tidx,0:ncol_p-1))/dt
;  str = "dE/dt (phAM-phBF)"

  fmin = -400.
  fmax = 400.

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
  global_ave = sum(area(0:ncol_p-1)*fld(0:ncol_p-1))/area_sphere
  print(str+" global ave: "+global_ave)
  res_p@gsnRightString        = "mean: "+sprintf("%3.2g", global_ave        )+" W/m^2"
  plot(3) = gsn_csm_contour_map_ce(wks,fld,res_p)



; ***********************************************
; 
; ***********************************************

  
  fld = HEVAP_AC(tidx,0:ncol_p-1)+HEVAP_BC(tidx,0:ncol_p-1)
  str = "hevap_ac+hevap_bc"

  fmin = -400.
  fmax = 400.

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
  global_ave = sum(area(0:ncol_p-1)*fld(0:ncol_p-1))/area_sphere
  print(str+" global ave: "+global_ave)
  res_p@gsnRightString        = "mean: "+sprintf("%3.2g", global_ave        )+" W/m^2"
  plot(4) = gsn_csm_contour_map_ce(wks,fld,res_p)


; ***********************************************
; 
; ***********************************************

  
  fld = HEVAP_AC(tidx,0:ncol_p-1)+HEVAP_BC(tidx,0:ncol_p-1)+HRAIN_AC(tidx,0:ncol_p-1)+HRAIN_BC(tidx,0:ncol_p-1)+HSNOW_AC(tidx,0:ncol_p-1)+HSNOW_BC(tidx,0:ncol_p-1)
  str = "(d)+(e)"

  fmin = -400.
  fmax = 400.

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
  global_ave = sum(area(0:ncol_p-1)*fld(0:ncol_p-1))/area_sphere
  print(str+" global ave: "+global_ave)
  res_p@gsnRightString        = "mean: "+sprintf("%3.2g", global_ave        )+" W/m^2"
  plot(5) = gsn_csm_contour_map_ce(wks,fld,res_p)

; ***********************************************
; 
; ***********************************************


  fld = HRAIN_AC(tidx,0:ncol_p-1)+HRAIN_BC(tidx,0:ncol_p-1)
  str = "HRAIN_AC/BC"

  fmin = -400.
  fmax = 400.

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
  global_ave = sum(area(0:ncol_p-1)*fld(0:ncol_p-1))/area_sphere
  print(str+" global ave: "+global_ave)
  res_p@gsnRightString        = "mean: "+sprintf("%3.2g", global_ave        )+" W/m^2"
  plot(6) = gsn_csm_contour_map_ce(wks,fld,res_p)

; ***********************************************
; 
; ***********************************************

  TT = T(tidx,57,ncol_p-1)
  fld = -TT*FRAIN_BC(tidx,0:ncol_p-1)*cpliq-TT*FRAIN_AC(tidx,0:ncol_p-1)*cpliq ;+TT*FSNOW_BC(tidx,0:ncol_p-1)*cpice+TT*FSNOW_AC(tidx,0:ncol_p-1)*cpice
  str = "(g) diagnosed from F*TT*cp"

  fmin = -400.
  fmax = 400.

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
  global_ave = sum(area(0:ncol_p-1)*fld(0:ncol_p-1))/area_sphere
  print(str+" global ave: "+global_ave)
  res_p@gsnRightString        = "mean: "+sprintf("%3.2g", global_ave        )+" W/m^2"
  plot(7) = gsn_csm_contour_map_ce(wks,fld,res_p)

; ***********************************************
; 
; ***********************************************

; TT = T(tidx,47,ncol_p-1)
;  fld = -TT*FRAIN_BC(tidx,0:ncol_p-1)*cpliq-TT*FRAIN_AC(tidx,0:ncol_p-1)*cpliq ;+TT*FSNOW_BC(tidx,0:ncol_p-1)*cpice+TT*FSNOW_AC(tidx,0:ncol_p-1)*cpice
  str = "(g) diagnosed from F*TT*cp, k=47"

  fld = (TE_dyAC(tidx,0:ncol_p-1)-TE_dyAP(tidx,0:ncol_p-1))/dt
  str = "dme_adjust diagnosted from phAC-phAP"

  fmin = -400.
  fmax = 400.

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
  global_ave = sum(area(0:ncol_p-1)*fld(0:ncol_p-1))/area_sphere
  print(str+" global ave: "+global_ave)
  res_p@gsnRightString        = "mean: "+sprintf("%3.2g", global_ave        )+" W/m^2"
  plot(8) = gsn_csm_contour_map_ce(wks,fld,res_p); ***********************************************
; 
; ***********************************************

  fld = TE_dEdt_cpupdate_dynE(tidx,0:ncol_p-1)/dt
  str = "cpupdate"
  
  fmin = -4.
  fmax = 4.

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
  global_ave = sum(area(0:ncol_p-1)*fld(0:ncol_p-1))/area_sphere
  print(str+" global ave: "+global_ave)
  res_p@gsnRightString        = "mean: "+sprintf("%3.2g", global_ave        )+" W/m^2"
  plot(9) = gsn_csm_contour_map_ce(wks,fld,res_p)

; ***********************************************
; 
; ***********************************************
  fld =   TE_dEdt_Q_LT_dynE(tidx,0:ncol_p-1)/dt
  str = "TE_dEdt_Q_LT"
  fmin = -400.0
  fmax = 400.0

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
  global_ave = sum(area(0:ncol_p-1)*fld(0:ncol_p-1))/area_sphere
  print(str+" global ave: "+global_ave)
  res_p@gsnRightString        = "mean: "+sprintf("%3.2g", global_ave        )+" W/m^2"
  plot(10) = gsn_csm_contour_map_ce(wks,fld,res_p)

; ***********************************************
; 
; ***********************************************
  fld = imbalance(tidx,0:ncol_p-1)
  str = "imbalance"
  fmin = -400.0
  fmax = 400.0

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
  global_ave = sum(area(0:ncol_p-1)*fld(0:ncol_p-1))/area_sphere
  print(str+" global ave: "+global_ave)
  res_p@gsnRightString        = "mean: "+sprintf("%3.2g", global_ave        )+" W/m^2"
  plot(11) = gsn_csm_contour_map_ce(wks,fld,res_p)


; ***********************************************
; 
; ***********************************************
  TT = T(tidx,40,ncol_p-1)
  fld = 1000.0*prect(tidx,0:ncol_p-1)*cpliq*TT
  str = "PRECT"
  fmin = -400.0
  fmax = 400.0

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
  global_ave = sum(area(0:ncol_p-1)*fld(0:ncol_p-1))/area_sphere
  print(str+" global ave: "+global_ave)
  res_p@gsnRightString        = "mean: "+sprintf("%3.2g", global_ave        )+" W/m^2"
  plot(12) = gsn_csm_contour_map_ce(wks,fld,res_p)





  
; ***********************************************
; Plot everything
; ***********************************************
  resP                  = True                   ; modify the panel plot
  resP@gsnFrame         = False                  ; don't advance panel plot
  resP@gsnPanelFigureStrings = (/"a","b","c","d","e","f","g","h","i","j"/)
  resP@gsnPanelYWhiteSpacePercent = 4.0
  resP@gsnPanelXWhiteSpacePercent = 4.0

  resP@txString   = "Total energy equation terms/budget - cp dyn energy"
  gsn_panel(wks,plot,(/5,3/),resP)               ; now draw as one plot

  txres               = True
  frame(wks)



end