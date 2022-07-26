begin
  frames = False
  fv = True

  str_TT = "TT_LW" ;1

;  var = "TT5"
;  str_TT = "TT_UN" ;5

;  str_TT = "TTRMD" ;4

  var = "Q"
  str_TT = "Q"


  f = addfile("h0.nc","r")

  lat = f->lat
  lon = f->lon 

  if (fv) then
    nlat = dimsizes(lat)
    nlon = dimsizes(lon)
 
    lon2d = new((/nlat,nlon/),typeof(lon(0)))
    lat2d = new((/nlat,nlon/),typeof(lat(0)))
    do j=0,nlat-1
      lon2d(j,:)= lon(:)
    end do
    do i=0,nlon-1
      lat2d(:,i)= lat(:)
    end do


;    target_lat = 46.0
;    target_lon = 185.0
;    nm = getind_latlon2d (lat2d,lon2d, target_lat, target_lon) ; middle of Nortern Pacific
;    print(nm)
;    ilat = nm(0,0)
;    ilon = nm(0,1)
;    print("Column location (lat,lon): "+lat(ilat)+" "+lon(ilon))
  ;  str_point = "Pacific"

    target_lat = 31.0
    target_lon = 80.0
    nm = getind_latlon2d (lat2d,lon2d, target_lat, target_lon) ; Himalayas
    ilat = nm(0,0)
    ilon = nm(0,1)
    print("Column location (lat,lon): "+lat(ilat)+" "+lon(ilon))
    str_point = "Himalayas"

  else 
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
  end if




  if (fv) then
    sr = var+"_BC1"
    TT_BC1 = f->$sr$(:,:,ilat,ilon)

    sr = var+"_AC4"
    TT_AC4 = f->$sr$(:,:,ilat,ilon)
    sr = var+"_AC5"
    TT_AC5 = f->$sr$(:,:,ilat,ilon)
    sr = var+"_AC11d"
    TT_AC11d = f->$sr$(:,:,ilat,ilon)
    sr = var+"_AC11e"
    TT_AC11e = f->$sr$(:,:,ilat,ilon)
    sr = var+"_AC12"
    TT_AC12 = f->$sr$(:,:,ilat,ilon)
  else
    sr = var+"_BC1"
    TT_BC1 = f->$sr$(:,:,icol)
    sr = var+"_AC4"
    TT_AC4   = f->$sr$(:,:,icol)
    sr = var+"_AC5"
    TT_AC5   = f->$sr$(:,:,icol)
    sr = var+"_AC11d"
    TT_AC11d = f->$sr$(:,:,icol)
    sr = var+"_AC11e"
    TT_AC11e = f->$sr$(:,:,icol)
    sr = var+"_AC12"
    TT_AC12  = f->$sr$(:,:,icol)

    sr = var+"_dyn1"
    TT_dyn1   = f->$sr$(:,:,icol)
    sr = var+"_dyn2"
    TT_dyn2   = f->$sr$(:,:,icol)
    sr = var+"_dyn_remap1"
    TT_remap1 = f->$sr$(:,:,icol)
    sr = var+"_dyn_remap2"
    TT_remap2 = f->$sr$(:,:,icol)
    sr = var+"_DP"
    TT_DP     = f->$sr$(:,:,icol)
  end if

  nt   = dimsizes(TT_BC1(:,0))
  nlev = dimsizes(TT_BC1(0,:))


  if (fv) then
    str = "fv-"
  else
    str = "se-"
   end if
  str = str+str_TT
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




  if (frames) then
    tott = nt-1
  else
    tott = 0
  end if

  do t=0,tott



    res@tiMainString     = "Incr. tphysac all"          ; add title
    if (frames) then
      fld=TT_AC12(t,:)-TT_BC1(t,:)
      plot(1)  = gsn_csm_xy (wks,fld,ilev,res) ; create plot
    else
      fld=TT_AC12-TT_BC1
      plot(1)  = gsn_csm_xy (wks,fld(:,:),ilev,res) ; create plot 
    end if


    res@tiMainString     = "Incr. CLUBB (5)"          ; add title
    if (frames) then
      fld=TT_AC5(t,:)-TT_AC4(t,:)
    else
      fld=TT_AC5-TT_AC4
    end if
    plot(2)  = gsn_csm_xy (wks,fld,ilev,res) ; create plot


    res@tiMainString     = "Incr. gw_drag"          ; add title
    if (frames) then
      fld=TT_AC11e(t,:)-TT_AC11d(t,:)
    else
      fld=TT_AC11e-TT_AC11d
    end if
    plot(3)  = gsn_csm_xy (wks,fld,ilev,res) ; create plot

    if (fv) then
    else
      res@tiMainString     = "Floating Lagrangian dycore"          ; add title
      if (frames) then
        fld=TT_remap1(t,:)-TT_dyn2(t,:)
      else
        fld=TT_remap1-TT_dyn2
      end if
      plot(4)  = gsn_csm_xy (wks,fld,ilev,res) ; create plot

      res@tiMainString     = "vertical remap dycore"          ; add title
      if (frames) then
        fld=TT_remap2(t,:)-TT_remap1(t,:)
      else
        fld=TT_remap2-TT_remap1
      end if
      plot(5)  = gsn_csm_xy (wks,fld,ilev,res) ; create plot
    end if

    res@tiMainString     = str_TT +" (tphysbc start)"         ; add title
    if (frames) then
      fld=TT_BC1(t,:)
    else
      fld=TT_BC1
    end if
    plot(0)  = gsn_csm_xy (wks,fld,ilev,res) ; create plot


    resP = True
    if (fv) then
      resP@gsnPanelMainString = "CAM-FV; instantaneous every time-step;"+str_TT+" lat,lon="+target_lat+","+target_lon
    else
      resP@gsnPanelMainString = "CAM-SE; instantaneous every time-step;"+str_TT+" lat,lon="+target_lat+","+target_lon
    end if
    gsn_panel(wks,plot,(/2,3/),resP)             ; now draw as one plot
 
  end do
;************************************************
; create panel
;************************************************
end