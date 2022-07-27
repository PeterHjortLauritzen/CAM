begin
  frames = True
  fv = False

  f = addfile("h0.nc","r")

  lat = f->lat
  lon = f->lon 

  lon2d = new((/1,dimsizes(lon)/),typeof(lon(0)))
  lat2d = new((/1,dimsizes(lat)/),typeof(lat(0)))
  r2d  = get_r2d("double")
  lon2d(0,:)= lon(:)
  lat2d(0,:)= lat(:)


;  target_lat = 46.0
;  target_lon = 185.0
;  nm = getind_latlon2d (lat2d,lon2d, target_lat, target_lon) ; middle of Nortern Pacific
;  icol = nm(0,1)
;  str_point = "Pacific"


  target_lat = 31.0
  target_lon = 80.0
  nm = getind_latlon2d (lat2d,lon2d, target_lat, target_lon) ; Himalayas
  icol = nm(0,1)
  str_point = "Himalayas"

  print("Column location (lat,lon): "+lat(icol)+" "+lon(icol))


  if (fv) then
    N2O      = f->N2O(:,:,ilat,ilon)

    N2O_BC1 = f->N2O_BC1(:,:,ilat,ilon)
    N2O_BC2 = f->N2O_BC2(:,:,ilat,ilon)
 
    N2O_AC1 = f->N2O_AC1(:,:,ilat,ilon)
    N2O_AC4 = f->N2O_AC4(:,:,ilat,ilon)
    N2O_AC5 = f->N2O_AC5(:,:,ilat,ilon)
    N2O_AC6 = f->N2O_AC6(:,:,ilat,ilon)
    N2O_AC7 = f->N2O_AC7(:,:,ilat,ilon)
    N2O_AC8 = f->N2O_AC8(:,:,ilat,ilon)
    N2O_AC9 = f->N2O_AC9(:,:,ilat,ilon)
    N2O_AC10 = f->N2O_AC10(:,:,ilat,ilon)
    N2O_AC11 = f->N2O_AC11(:,:,ilat,ilon)
    N2O_AC11d = f->N2O_AC11d(:,:,ilat,ilon)
    N2O_AC11e = f->N2O_AC11e(:,:,ilat,ilon)
    N2O_AC12 = f->N2O_AC12(:,:,ilat,ilon)
  else
    N2O      = f->N2O(:,:,icol)

    N2O_dyn1 = f->N2O_dyn1(:,:,icol)
    N2O_dyn2 = f->N2O_dyn2(:,:,icol)
    N2O_DP   = f->N2O_DP(:,:,icol)
  
    N2O_BC1 = f->N2O_BC1(:,:,icol)
    N2O_BC2 = f->N2O_BC2(:,:,icol)
 
    N2O_AC1 = f->N2O_AC1(:,:,icol)
    N2O_AC4 = f->N2O_AC4(:,:,icol)
    N2O_AC5 = f->N2O_AC5(:,:,icol)
    N2O_AC6 = f->N2O_AC6(:,:,icol)
    N2O_AC7 = f->N2O_AC7(:,:,icol)
    N2O_AC8 = f->N2O_AC8(:,:,icol)
    N2O_AC9 = f->N2O_AC9(:,:,icol)
    N2O_AC10 = f->N2O_AC10(:,:,icol)
    N2O_AC11 = f->N2O_AC11(:,:,icol)
    N2O_AC11d = f->N2O_AC11d(:,:,icol)
    N2O_AC11e = f->N2O_AC11e(:,:,icol)
    N2O_AC12 = f->N2O_AC12(:,:,icol)

    N2O_dyn_remap1 = f->N2O_dyn_remap1(:,:,icol)
    N2O_dyn_remap2 = f->N2O_dyn_remap2(:,:,icol)
  end if

  nt   = dimsizes(N2O(:,0))
  nlev = dimsizes(N2O_AC1(0,:))


  if (fv) then
    str = "fv-"
  else
    str = "se-"
   end if

  if (frames) then
    str   = str + "vert-profile-frames."+str_point
    wks   = gsn_open_wks ("pdf",str)               ; send graphics to PNG file
  else
    str   = str + "vert-profile."+str_point
    wks   = gsn_open_wks ("pdf",str)               ; send graphics to PNG file
  end if

  plot = new(12,graphic)

  res          = True
  res@gsnDraw  = False                          ; don't draw
  res@gsnFrame = False                          ; don't advance frame
  res@trYReverse = True		; reverse the Y-axis

  ilev = ispan(0,nlev-1,1)



  printVarSummary(N2O_AC1)

  if (frames) then
    tott = nt-1
  else
    tott = 0
  end if

  do t=0,tott

    res@tiMainString     = "Incr. tphysac all"          ; add title
    if (frames) then
      fld=N2O_AC12(t,:)-N2O_AC1(t,:)
      plot(1)  = gsn_csm_xy (wks,1E9*fld,ilev,res) ; create plot
    else
      fld=N2O_AC12-N2O_AC1
      plot(1)  = gsn_csm_xy (wks,1E9*fld(:,:),ilev,res) ; create plot 
    end if


    res@tiMainString     = "Incr. tphysbc"          ; add title
    if (frames) then
      fld=N2O_BC2(t,:)-N2O_BC1(t,:)
      plot(2)  = gsn_csm_xy (wks,1E9*fld,ilev,res) ; create plot
    else
      fld=N2O_BC2-N2O_BC1
    end if
    plot(2)  = gsn_csm_xy (wks,1E9*fld,ilev,res) ; create plot
;;
;
    res@tiMainString     = "Incr. Carma+contrails (1-4)"          ; add title
    if (frames) then
      fld=N2O_AC4(t,:)-N2O_AC1(t,:)
    else
      fld=N2O_AC4-N2O_AC1
    end if
    plot(3)  = gsn_csm_xy (wks,1E9*fld,ilev,res) ; create plot

    res@tiMainString     = "Incr. CLUBB (5)"          ; add title
    if (frames) then
      fld=N2O_AC5(t,:)-N2O_AC4(t,:)
    else
      fld=N2O_AC5-N2O_AC4
    end if
    plot(4)  = gsn_csm_xy (wks,1E9*fld,ilev,res) ; create plot

    res@tiMainString     = "Incr. process (7-9)"          ; add title
    if (frames) then
      fld=N2O_AC9(t,:)-N2O_AC5(t,:)
    else
      fld=N2O_AC9-N2O_AC5
    end if
    plot(5)  = gsn_csm_xy (wks,1E9*fld,ilev,res) ; create plot

    res@tiMainString     = "Incr. co2 cycle+chemistry (10)"          ; add title
    if (frames) then
      fld=N2O_AC10(t,:)-N2O_AC9(t,:)
    else
      fld=N2O_AC10-N2O_AC9
    end if
    plot(6)  = gsn_csm_xy (wks,1E9*fld,ilev,res) ; create plot

    res@tiMainString     = "Incr. vertical diffusion+... (11-11d)"          ; add title
    if (frames) then
      fld=N2O_AC11d(t,:)-N2O_AC10(t,:)
    else
      fld=N2O_AC11d-N2O_AC10
    end if
    plot(7)  = gsn_csm_xy (wks,1E9*fld,ilev,res) ; create plot

    res@tiMainString     = "Incr. gw_drag"          ; add title
    if (frames) then
      fld=N2O_AC11e(t,:)-N2O_AC11d(t,:)
    else
      fld=N2O_AC11e-N2O_AC11d
    end if
    plot(8)  = gsn_csm_xy (wks,1E9*fld,ilev,res) ; create plot

    res@tiMainString     = "Incr. end tphysac"          ; add title
    if (frames) then
      fld=N2O_AC12(t,:)-N2O_AC11e(t,:)
    else
      fld=N2O_AC12-N2O_AC11e
    end if
    plot(9)  = gsn_csm_xy (wks,1E9*fld,ilev,res) ; create plot


    if (fv) then
    else
      res@tiMainString     = "Floating dyn incr. "          ; add title
      if (frames) then
        fld=N2O_dyn_remap1(t,:)-N2O_dyn2(t,:)
      else
        fld=N2O_dyn_remap1-N2O_dyn2
      end if
      plot(10)  = gsn_csm_xy (wks,1E9*fld,ilev,res) ; create plot

      res@tiMainString     = "vertical remap (dycore)"          ; add title
      if (frames) then
        fld=N2O_dyn_remap2(t,:)-N2O_dyn_remap1(t,:)
      else
        fld=N2O_dyn_remap2-N2O_dyn_remap1
      end if
      plot(11)  = gsn_csm_xy (wks,1E9*fld,ilev,res) ; create plot
    end if



;    res@tiMainString     = "phys incr in dynamics"          ; add title
;    if (frames) then
;      fld=N2O_dyn2(t,:)-N2O_dyn1(t,:)
;    else
;      fld=N2O_dyn2-N2O_dyn1
;    end if
;    plot(11)  = gsn_csm_xy (wks,1E9*fld,ilev,res) ; create plot

;    res@trXMinF          =  490.0    ; Set minimum Y axis value.
;    res@trXMaxF          =  492.5    ; Set maximum Y axis value.


    res@tiMainString     = "N2O_BC1"          ; add title
    if (frames) then
      fld=N2O_BC1(t,:)
    else
      fld=N2O_BC1
    end if
    plot(0)  = gsn_csm_xy (wks,1E9*fld,ilev,res) ; create plot


    resP = True
    resP@gsnPanelMainString = "CAM-SE; instantaneous every time-step; lat,lon="+target_lat+","+target_lon
    gsn_panel(wks,plot,(/3,4/),resP)             ; now draw as one plot
 
  end do
;************************************************
; create panel
;************************************************
end