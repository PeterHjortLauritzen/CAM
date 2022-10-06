load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

begin
  fname = "cam-mpas-pdc-paper_mpasa120_mpasa120_192_short.cam.h1.*.nc"
  fnames =  systemfunc("ls "+fname)
  print(fnames)
  fall = addfiles(fnames,"r")       ; Open netCDF files.	

  ptop = fall[1]->ptop_mpas
  time = fall[1]->time

  lon=fall[0]->lon
  lat=fall[0]->lat
  ncol = dimsizes(lon(:))
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


  wks   = gsn_open_wks ("pdf","ptop")  
  res                  = True                     ; plot mods desired
  res@tiMainString     = "Basic XY plot"          ; add title

  plot  = gsn_csm_xy (wks,time,ptop(:,icol),res) 

end
