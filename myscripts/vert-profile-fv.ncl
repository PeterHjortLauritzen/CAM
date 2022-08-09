begin
  f = addfile("h0.nc","r")

  lat = f->lat
  lon = f->lon 
 
  nlat = dimsizes(lat)
  nlon = dimsizes(lon)

  lon2d = new((/nlat,nlon/),typeof(lon(0)))
  lat2d = new((/nlat,nlon/),typeof(lat(0)))
  r2d  = get_r2d("double")
  do j=0,nlat-1
    lon2d(j,:)= lon(:)
  end do
  do i=0,nlon-1
    lat2d(:,i)= lat(:)
  end do

  N2O = f->N2O

  N2O_BC1 = f->N2O_BC1
  N2O_BC2 = f->N2O_BC2

  N2O_AC1 = f->N2O_AC1
;  N2O_AC2 = f->N2O_AC2
;  N2O_AC3 = f->N2O_AC3
  N2O_AC4 = f->N2O_AC4
  N2O_AC5 = f->N2O_AC5
  N2O_AC6 = f->N2O_AC6
  N2O_AC7 = f->N2O_AC7
  N2O_AC8 = f->N2O_AC8
  N2O_AC9 = f->N2O_AC9
  N2O_AC10 = f->N2O_AC10
  N2O_AC11 = f->N2O_AC11
;  N2O_AC11b = f->N2O_AC11b
;  N2O_AC11c = f->N2O_AC11c
  N2O_AC11d = f->N2O_AC11d
  N2O_AC11e = f->N2O_AC11e
;  N2O_AC11f = f->N2O_AC11f
;  N2O_AC11g = f->N2O_AC11g
;  N2O_AC11h = f->N2O_AC11h
  N2O_AC12 = f->N2O_AC12
;  N2O_AC13 = f->N2O_AC13

  nt = dimsizes(N2O(:,0,0,0))

  wks   = gsn_open_wks ("pdf","vert-profile-fv")               ; send graphics to PNG file

  plot = new(12,graphic)

  res          = True
  res@gsnDraw  = False                          ; don't draw
  res@gsnFrame = False                          ; don't advance frame

  nlev = dimsizes(N2O_AC1(0,:,0,0))
  ilev = ispan(0,nlev-1,1)

  target_lat = 46.0
  target_lon = 185.0
  nm = getind_latlon2d (lat2d,lon2d, target_lat, target_lon) ; middle of Nortern Pacific
  print(nm)
  ilat = nm(0,0)
  ilon = nm(0,1)
  print("Column location (lat,lon): "+lat(ilat)+" "+lon(ilon))

;  target_lat = 31.0
;  target_lon = 80.0
;  nm = getind_latlon2d (lat2d,lon2d, target_lat, target_lon) ; Himalayas
;  ilat = nm(0,0)
;  ilon = nm(0,1)
;  print("Column location (lat,lon): "+lat(ilat)+" "+lon(ilon))


  printVarSummary(N2O_AC1)

  do t=0,0
    res@tiMainString     = "Incr. tphysac all"          ; add title
    fld=N2O_AC12-N2O_AC1
    plot(1)  = gsn_csm_xy (wks,1E9*fld(:,:,ilat,ilon),ilev,res) ; create plot

    res@tiMainString     = "Incr. tphysbc"          ; add title
    fld=N2O_BC2-N2O_BC1
    plot(2)  = gsn_csm_xy (wks,1E9*fld(:,:,ilat,ilon),ilev,res) ; create plot


    res@tiMainString     = "Incr. Carma+contrails (1-4)"          ; add title
    fld=N2O_AC4-N2O_AC1
    plot(3)  = gsn_csm_xy (wks,1E9*fld(:,:,ilat,ilon),ilev,res) ; create plot

    res@tiMainString     = "Incr. CLUBB (5)"          ; add title
    fld=N2O_AC5-N2O_AC4
    plot(4)  = gsn_csm_xy (wks,1E9*fld(:,:,ilat,ilon),ilev,res) ; create plot

    res@tiMainString     = "Incr. process (7-9)"          ; add title
    fld=N2O_AC9-N2O_AC5
    plot(5)  = gsn_csm_xy (wks,1E9*fld(:,:,ilat,ilon),ilev,res) ; create plot

    res@tiMainString     = "Incr. co2 cycle+chemistry (10)"          ; add title
    fld=N2O_AC10-N2O_AC9
    plot(6)  = gsn_csm_xy (wks,1E9*fld(:,:,ilat,ilon),ilev,res) ; create plot

    res@tiMainString     = "Incr. vertical diffusion+... (11-11d)"          ; add title
    fld=N2O_AC11d-N2O_AC10
    plot(7)  = gsn_csm_xy (wks,1E9*fld(:,:,ilat,ilon),ilev,res) ; create plot

    res@tiMainString     = "Incr. gw_drag"          ; add title
    fld=N2O_AC11e-N2O_AC11d
    plot(8)  = gsn_csm_xy (wks,1E9*fld(:,:,ilat,ilon),ilev,res) ; create plot

    res@tiMainString     = "Incr. end tphysac"          ; add title
    fld=N2O_AC12-N2O_AC11e
    plot(9)  = gsn_csm_xy (wks,1E9*fld(:,:,ilat,ilon),ilev,res) ; create plot

    res@trXMinF          =  490.0    ; Set minimum Y axis value.
    res@trXMaxF          =  492.5    ; Set maximum Y axis value.

    res@tiMainString     = "N2O_BC1"          ; add title
    fld=N2O_BC1
    plot(0)  = gsn_csm_xy (wks,1E9*fld(:,:,ilat,ilon),ilev,res) ; create plot



  end do
;************************************************
; create panel
;************************************************
  resP = True
  resP@gsnPanelMainString = "CAM-FV; instantaneous every time-step; lat,lon="+target_lat+","+target_lon
  gsn_panel(wks,plot,(/3,4/),resP)             ; now draw as one plot
end