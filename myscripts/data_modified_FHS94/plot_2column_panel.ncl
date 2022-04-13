load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

begin
  ;
  ; user set variables
  ;
  
  print(vname)
;  vname      = (/var_left_panel,var_right_panel/)
;  vname                 = "PS"      ; variable you wish to output
  horizontal_resolution = "ne30_ne30" ; -res setting
  compset="FKESSLER"                ; compset used
;  iday  = 10                        ; plot day
;  iday  = 3                        ; plot day
  iday = tointeger(time_idx)
  output_format = "pdf"             ; output format
  

;******************************************************
; create plot
;******************************************************

  fnames = (/"modified-FHS94_cam_development_f09_f09_mg17_a.ave.h0.nc",\
             "modified-FHS94_topo-mods_ne30_ne30_mg17_b.ave.h0.nc",\
             "modified-FHS94_topo-mods_ne30_ne30_mg17_c.ave.h0.nc",\
	     "modified-FHS94_topo-mods_ne30_ne30_mg17_d.ave.h0.nc",\
             "modified-FHS94_topo-mods_ne30pg3_ne30pg3_mg17_e.ave.h0.nc",\
	     "modified-FHS94_topo-mods_ne30pg3_ne30pg3_mg17_f.ave.h0.nc"/)

  leftString = (/"FV","SE: old PGF (CESM2.0)","SE: Exner PGF","SE: Exner PGF + ref profiles in del4 (CESM2.2)","SE-CSLAM: Exner PGF","SE-CSLAM: Exner PGF + ref profiles in del4 (CESM2.2)"/);,"SE: Exner PGF + ref profiles in del4 + lap_p correction","SE: Exner PGF + lap_p correction"/)
  wks = gsn_open_wks(output_format,outfname)
  print("output file "+outfname)
  fall = addfiles(fnames,"r")       ; Open netCDF files.	
  fcnt = dimsizes(fnames)
  print(fnames)


  if (vname.eq."PRECT") then
    gsn_define_colormap(wks,"WhBlGrYeRe")
  else
    gsn_define_colormap(wks,"BlAqGrYeOrReVi200")   ; choose colormap
  end if

  res                     = True         ; plot modifications desired
  res@gsnMaximize         = False         ; Maximize size of plot in frame
  res@gsnSpreadColors     = True         ; Use full colormap, but start
  res@gsnDraw             = False           ; don't draw
  res@gsnFrame            = False           ; don't advance frame
  res@cnFillOn            = True         ; Turn on contour fill

  res@cnLinesOn           = False        ; Turn off contour lines
  res@cnLineLabelsOn      = False        ; Turn off contour line labels
  res@lbLabelAutoStride   = True         ; Clean up labelbar labels.
;  res@mpOutlineOn     = False       ; No continents
;  res@cnFillMode="RasterFill"

;  res@mpMinLatF = 20.                   ; set the minimum latitude = -70.
;  res@mpMaxLatF = 70.                    ; set the maximum latitude = 70.

;  res@mpMinLonF = 30.                    ; set the minimum latitude = -70.
;  res@mpMaxLonF = 210.

  res@gsnStringFontHeightF = 0.02
  res@tiMainFontHeightF = 0.025
  res@lbLabelFontHeightF = 0.02
  res@tmXBLabelFontHeightF = 0.015
  res@tmYLLabelFontHeightF = 0.015
  res@pmLabelBarOrthogonalPosF = 0.2



  plot = new(2*fcnt,graphic) 
  plotk=0
  do k=0,fcnt-1
;    time= fall[k]->time(iday)
;    print("time "+time)
  
    do nvar=0,0
      print("plot "+k)
      if (vname.eq."PRECT") then
          vartmp = fall[k]->PRECL;+PRECC
      else
         vartmp = fall[k]->$vname(nvar)$
      end if
      latlon=all(isdim(vartmp,(/"time","lat","lon"/)))
      delete(vartmp)
      
      if (latlon) then
        print("lat lon data")
        if (vname.eq."PRECT") then
           var1 = fall[k]->PRECL(iday,:,:)
           var2 = fall[k]->PRECC(iday,:,:)
           var = var1
           var  = var1+var2
        else
           var = fall[k]->$vname(nvar)$(iday,:,:)
        end if


        lat = fall[k]->lat          ; Pull off lat/lon variable and
        lon = fall[k]->lon          ; convert both to 1D.    	
      else
        print("1D data")      
        if (vname.eq."PRECT") then         
           var1 = fall[k]->PRECL(iday,:)
           var2 = fall[k]->PRECC(iday,:)
           var  = var1
           var = var1+var2
        else
           var = fall[k]->$vname(nvar)$(iday,:)
        end if
        lat1d = fall[k]->lat          ; Pull off lat/lon variable and
        lon1d = fall[k]->lon          ; convert both to 1D.    	
        res@sfXArray            = lon1d        ; Required to tell NCL where to
        res@sfYArray            = lat1d        ; overlay data on globe.	
      end if

      res@cnLevelSelectionMode = "AutomaticLevels"

    res@gsnLeftString = leftString(k)
      if (vname(nvar).eq."PRECT") then
;        res@cnFillPalette         = "WhiteBlueGreenYellowRed"
        var = 1000*24*3600*var ;convert from m/s to mm/day  
        res@cnLevelSelectionMode = "ExplicitLevels"   ; set explicit contour levels
        res@cnLevels    = (/0.2,0.5,1.0,2,3,4,5,6,7,8,9,10,12,14,17/)

;        res@cnLevelSelectionMode = "ManualLevels"
;        res@gsnRightString = "mm/day"
;        res@cnMinLevelValF  = 0.0
;        res@cnMaxLevelValF  = 22.0
;        res@cnLevelSpacingF = 1.0 
;        res@cnFillPalette         = "WhiteBlueGreenYellowRed"
      end if

      if (vname(nvar).eq."OMEGA500") then
        gsn_define_colormap(wks,"MPL_seismic");
        res@cnLevelSelectionMode = "ManualLevels"
        res@cnMinLevelValF  = -0.1
        res@cnMaxLevelValF  = 0.1              
        res@cnLevelSpacingF = 0.005              
;      res@gsnRightString = "mm/day"


      end if

      if (vname(nvar).eq."OMEGA850") then
        res@cnLevelSelectionMode = "ManualLevels"
        gsn_define_colormap(wks,"MPL_seismic");
        res@cnMinLevelValF  = -1.0
        res@cnMaxLevelValF  = 1.0              
        res@cnLevelSpacingF = 0.05              
;      res@gsnRightString = "mm/day"
      end if

      if (vname(nvar).eq."PS") then
        var = var/100.0 ;convert from Pa to hPa
        res@gsnRightString = "hPa"
        if (time.eq.10) then
          res@cnLevelSelectionMode = "ManualLevels"
          res@cnMinLevelValF  = 945.0
          res@cnMaxLevelValF  = 1015.0
          res@cnLevelSpacingF = 5.0
        end if
        if (time.eq.3) then
          res@cnLevelSelectionMode = "ManualLevels"
          res@cnMinLevelValF  = 990.0
          res@cnMaxLevelValF  = 1010.0
          res@cnLevelSpacingF = 2.0
        end if
      end if

      res@mpCenterLonF = 180.0           ; Rotate the projection.
      plot(plotk) = gsn_csm_contour_map_ce(wks,var,res)
      plotk=plotk+1
      delete(var)
;      delete(res@cnMinLevelValF)
;      delete(res@cnMaxLevelValF)
;      delete(res@cnLevelSpacingF)
      if (latlon) then
        delete(lat)
        delete(lon)
      else
        delete(lat1d)
        delete(lon1d)
        delete(res@sfXArray)
        delete(res@sfYArray)      
      end if      
;      delete(res@cnLevels)
;      delete(res@lbLabelStrings)

    end do
;    delete(time)    
  end do
;************************************************
; create panel
;************************************************
  resP                  = True                   ; modify the panel plot
  resP@gsnFrame         = False                  ; don't advance panel plot
;  resP@gsnPanelLabelBar = True                  ; add common colorbar
;  resP@gsnPanelBottom   = 0.05
   resP@gsnPanelMainFontHeightF = 0.01
  if (vname.eq."PRECT") then
    resP@txString   = "PRECT, 1 year F2000climo, 32 levels"
  else
    resP@txString   = "OMEGA500, 18 months average, FHS94 forcing, 32 levels"
  end if
  resP@gsnPanelYWhiteSpacePercent = 2.0
;    resP@gsnPanelYWhiteSpacePercent = 10.0
  resP@gsnPanelXWhiteSpacePercent = 4.0
  resP@gsnPanelFigureStrings = (/"a","b","c","d","e","f"/)
  resP@gsnPanelFigureStringsFontHeightF = 0.01
  gsn_panel(wks,plot,(/fcnt,2/),resP)            ; now draw as one plot

  txres               = True
  txres@txFontHeightF = 0.015  
;  txres@txFontHeightF = 0.015
;  gsn_text_ndc(wks,"Figure: OMEGA500, 1 year F2000climo, 32 levels",0.5,0.02,txres)
  frame(wks)
end
