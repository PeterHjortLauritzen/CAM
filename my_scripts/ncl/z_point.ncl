load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

;function ave1d(fld,area)
;local ncol,i,sum_
;begin
;  ncol = dimsizes(area(:))    
;  sum_area = sum(area)
;  sum_fld  = sum(fld(:)*area(:))/sum(area)
;  ave1d    = sum_fld/sum_area
;  do i=0,i-1
;    sum_area = sum_area+area(i)
;    
;  end do
;end

begin
  latlon = False

;  threeD = False
;  vname = "dPSdt"

  threeD = True
  vname = "zi_phBF"

  fname = "cam-mpas-pdc-paper_mpasa120_mpasa120_192_short.cam.h1.*.nc";data/cam-mpas-pdc-paper_mpasa120_mpasa120_192_short.cam.h1.*.nc"
;  fname = "cam-mpas-pdc-paper_mpasa120_mpasa120_192_short.cam.h1.0001-07-15-00000.nc"
  fnames =  systemfunc("ls "+fname)
  output_format = "pdf"

  wks = gsn_open_wks(output_format,"2d_rasta_"+vname)
  fall = addfiles(fnames,"r")       ; Open netCDF files.	
  fcnt = dimsizes(fnames)
  print(fnames)

  time = fall[:]->time
  ;  
  ; get dimensions
  ;
  aa = fall[:]->zi_phBF(:,:,:) 
  nlev = dimsizes(aa(0,:,0))
  ncol = dimsizes(aa(0,0,:))
  nt   = dimsizes(aa(:,0,0))
  delete(aa)
  print("nt="+nt)

  area     = fall[0]->area
  inv_area = 1.0D0/sum(area(:))
  PHIS   = fall[:]->PHIS(0,:)    ; Read some data; convert to 1D.        
  zs     = PHIS/9.80665
  delete(PHIS)

  lon=fall[0]->lon
  lat=fall[0]->lat
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


  ave_zi_phBF = new((/nlev/),typeof(area(0)))
  ave_zi_phAP = new((/nlev/),typeof(area(0)))
  ave_zi_mpas = new((/nlev/),typeof(area(0)))
  
  ave_zm_phBF = new((/nlev/),typeof(area(0)))
  ave_zm_phAP = new((/nlev/),typeof(area(0)))
  ave_zm_mpas = new((/nlev/),typeof(area(0)))

  time_series_zm_phAP_minus_phBF = new((/3,nt,nlev/),typeof(area(0)))

  zi_phBF   = fall[:]->zi_phBF(:,:,icol)    ; Read some data; convert to 1D.        
;  zi_dyBF   = fall[:]->zi_dyBF(:,:,icol)    ; Read some data; convert to 1D.        
  zi_phAP   = fall[:]->zi_phAP(:,:,icol)    ; Read some data; convert to 1D.        

  zint_mpas = fall[:]->zint_mpas(:,:,icol)    ; Read some data; convert to 1D.      
 
  zm_phBF   = fall[:]->zm_phBF(:,:,icol)    ; Read some data; convert to 1D.        
  zm_dyBF   = fall[:]->zm_dyBF(:,:,icol)    ; Read some data; convert to 1D.        
  zm_phAP   = fall[:]->zm_phAP(:,:,icol)    ; Read some data; convert to 1D.        
  zm_phAM   = fall[:]->zm_phAM(:,:,icol)    ; Read some data; convert to 1D.        
  zmid_mpas = fall[:]->zmid_mpas(:,:,icol)    ; Read some data; convert to 1D.        


  do t=0,nt-1
    print("===============t "+t+"=================")


    ;
    ; MPAS heights area relative to sea-level
    ;
    zm_mpas = zmid_mpas(0,:)
    zi_mpas = zint_mpas(0,:)
    do k=0,nlev-1
      zm_mpas(k)   = zmid_mpas(t,k)-zs(icol)
      zi_mpas(k)   = zint_mpas(t,k)-zs(icol) 
    end do

    do k=0,nlev-1
    ;
    ; float diags
    ;
;      time_series_zm_phAP_minus_phBF(0,t,k) = zm_phAP(t,k)-zm_phBF(t,k)
      time_series_zm_phAP_minus_phBF(1,t,k) = zm_phAP(t,k)-zm_phBF(t,k)
      time_series_zm_phAP_minus_phBF(2,t,k) = zm_phAM(t,k)-zm_phAP(t,k)
    ;
    ; pmid bias
    ;
      time_series_zm_phAP_minus_phBF(0,t,k) = zm_phBF(t,k)-zmid_mpas(t,k)
;      time_series_zm_phAP_minus_phBF(1,t,k) = zm_dyBF(t,k)-zmid_mpas(t,k)
;      time_series_zm_phAP_minus_phBF(2,t,k) = zm_dyBF(t,k)-zmid_mpas(t,k)

    end do
  end do

  number_of_plots = 17
  offset =0 ;14
  
 
;  plot = new(nlev,graphic)
  plot = new(number_of_plots,graphic)
  poly = new(16,graphic)   ; For attaching grey-filled rectangles

  wtype          = "pdf"
;  wtype@wkWidth  = 2500
;  wtype@wkHeight = 2500
  wks = gsn_open_wks(wtype ,"z_point")                   ; ps,pdf,x11,png

  trxmin = time(0)
  trxmax = time(nt-1)

  res                       = True          ; plot mods desired
  res@gsnDraw               = False
  res@gsnFrame              = False

  res@vpWidthF              = 0.70       ; change aspect ratio of plot
  tmp = 1.0/tofloat(number_of_plots-1)
  res@vpHeightF             = tmp

  res@trXMinF               = trxmin
  res@trXMaxF               = trxmax

  res@xyLineColor           = "black"
  res@xyLineThicknessF      = 2.0


  res@tiYAxisFontColor      = "black" ;res@xyLineColor
  res@tiYAxisFontHeightF    = 0.03

  res@tfPolyDrawOrder       = "PreDraw"    ; Draw any attached primitives
                                           ; before rest of plot.

;  res@tmXBOn                = False
  res@tmYLLabelFontHeightF  = 0.01      ; resize tick labels
  res@tmYLLabelFontColor    = res@xyLineColor
  res@xyDashPatterns    = 0

  res@xyLineColors      = (/"blue","green","red"/)
;;;;; TEMP
  
;---Top right string
  amres=True
  amres@amParallelPosF   = -0.475    ; This is the right edge of the plot.
  amres@amOrthogonalPosF =  0.39    ; This is the bottom edge of the plo

  amres2=True
  amres2@amParallelPosF   =  0.428    ; This is the right edge of the plot.
  amres2@amOrthogonalPosF =  0.39    ; This is the bottom edge of the plo

  txres                       = True
  txres@txPerimOn             = True
  txres@txBackgroundFillColor = "Salmon"
  txres@txFontHeightF         = 0.006

  do k=0,number_of_plots-1
    idx = k+offset+1
    print("idx="+idx)
    if ((k % 2) .eq. 0) then ; only works for integers
      res@tmYRLabelsOn         = True
      res@tmYLLabelsOn         = False
    else
       res@tmYRLabelsOn         = False
       res@tmYLLabelsOn         = True
    end if
    txt = "level "+idx
    txid = gsn_create_text(wks, txt, txres)
    txt = "z mid MPAS = "+zmid_mpas(0,k+offset)+" m"
    txid2 = gsn_create_text(wks, txt, txres)
    if (k.eq.0) then
      res@xyExplicitLegendLabels = (/"beginning CAM physics mid-level z minus MPAS mid-level z",\
                                     "end of CAM physics (but before dry-mass adjustmet) z minus beginning CAM physics z",\
				     "CAM physics z after dry-mass adjustment (dma) minus z before dma"/) 
      res@pmLegendDisplayMode    = "Always"            ; turn on legend
      res@pmLegendSide           = "Top"               ; Change location of 
      res@pmLegendWidthF         = 0.15                ; Change width and
      res@pmLegendHeightF        = 0.1                ; height of legend.
      res@lgPerimOn              = False               ; turn off box around
      res@lgLabelFontHeightF     = .01                 ; label font height
    else
      res@pmLegendDisplayMode    = "Never"            ; turn on legend
    end if
;  res@pmLegendParallelPosF   = .1                  ; move units right
;  res@pmLegendOrthogonalPosF = -0.3                ; move units down

    plot(k) = gsn_csm_xy(wks,time,time_series_zm_phAP_minus_phBF(:,:,k+offset),res)
    annoid = gsn_add_annotation(plot(k), txid , amres)
    annoid = gsn_add_annotation(plot(k), txid2, amres2)
  end do

;;;;;; Attach plots along the X axes
  attachres1                     = True
  attachres1@gsnAttachPlotsXAxis = True  ;; attaches along x-axis
  attachres1@gsnAttachBorderOn   = True; False ;; No border please
  attachres2                     = True
  attachres2@gsnAttachPlotsXAxis = True  ;; attaches along x-axis
  attachres2@gsnAttachBorderOn   = True; False ;; No border please

  attachid1  = gsn_attach_plots(plot(0),plot(1:number_of_plots-1),attachres1,attachres2)



;;;;;;; Maximize output on the page and draw everything
  pres = False               ; No resources needed
  maximize_output(wks,pres)  ; Maximize plot on page  


end
