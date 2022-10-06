begin
;************************************************
; open netCDF file
;************************************************
  a = addfile(" cam-mpas-pdc-paper_mpasa120_mpasa120_192_short.cam.h0.0001-02.nc","r")
  zm_phBF   = a->zm_phBF(:,:,:)
  zm_phAP   = a->zm_phAP(:,:,:) 
  k=10
  t=0
  fld = zm_phAP(t,k,:)-zm_phBF(t,k,:)
  

  ;
  ; SCRIP file
  ;
  scripFile = "/glade/p/cesmdata/inputdata/atm/cam/coords/mpasa120_SCRIP_desc.200911.nc"
  sf =  addfile(scripFile,"r")
  cellfill = True

;************************************************
; create plot
;************************************************

  wks = gsn_open_wks("pdf","raster")         ; send graphics to PNG file
  gsn_define_colormap(wks,"BlAqGrYeOrReVi200")   ; choose colormap

  res                     = True         ; plot modifications desired
  res@gsnMaximize         = False         ; Maximize size of plot in frame
  res@cnFillOn            = True         ; Turn on contour fill
                                         ; use "CellFill" and "RasterFill"
  res@cnLinesOn           = False        ; Turn off contour lines
  res@cnLineLabelsOn      = False        ; Turn off contour line labels
  res@lbLabelAutoStride   = True         ; Clean up labelbar labels.

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

  res@sfXArray            = lon1d        ; Required to tell NCL where to
  res@sfYArray            = lat1d        ; overlay data on globe.

;  res@cnFillMode       = "RasterFill"       ; Raster Mode

  res@mpCenterLonF     = 180                ; set map center at 180
                                            ; must be in colormap  
  plot = gsn_csm_contour_map_ce(wks,fld, res) ; create plot
end