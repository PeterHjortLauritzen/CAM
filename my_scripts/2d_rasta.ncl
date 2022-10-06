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

;  fname = "data/cam-mpas-pdc-paper_mpasa120_mpasa120_192_short.cam.h1.0001-06-*.nc"
  fname = "data/cam-mpas-pdc-paper_mpasa120_mpasa120_192_short.cam.h1.*.nc"
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

  ave_zi_phBF = new((/nlev/),typeof(area(0)))
  ave_zi_phAP = new((/nlev/),typeof(area(0)))
  ave_zi_mpas = new((/nlev/),typeof(area(0)))
  
  ave_zm_phBF = new((/nlev/),typeof(area(0)))
  ave_zm_phAP = new((/nlev/),typeof(area(0)))
  ave_zm_mpas = new((/nlev/),typeof(area(0)))

  time_series_zm_phAP_minus_phBF = new((/3,nt,nlev/),typeof(area(0)))

  ;
  ; get time-space averages
  ;
  fname_ave  = "data/cam-mpas-pdc-paper_mpasa120_mpasa120_192_short.ave.ANN.h0.nc"
  fnames_ave =  systemfunc("ls "+fname_ave)
  fall_ave = addfiles(fnames_ave,"r")       ; Open netCDF files.	
  zm_phBF   = fall_ave[:]->zm_phBF(0,:,:)
  zm_phAP   = fall_ave[:]->zm_phAP(0,:,:) 
  ; 
  ; time-space average (ave_ave)
  ; 
  ave_ave_zm_phAP_minus_phBF    = new((/nlev/),typeof(area(0)))
  do k=0,nlev-1
    ave_ave_zm_phAP_minus_phBF(k) = sum((zm_phAP(k,:)-zm_phBF(k,:))*area(:))*inv_area
  end do
  delete(zm_phBF)
  delete(zm_phAP)

  zi_phBF   = fall[:]->zi_phBF(:,:,:)    ; Read some data; convert to 1D.        
;  zi_dyBF   = fall[:]->zi_dyBF(:,:,:)    ; Read some data; convert to 1D.        
  zi_phAP   = fall[:]->zi_phAP(:,:,:)    ; Read some data; convert to 1D.        

  zint_mpas = fall[:]->zint_mpas(:,:,:)    ; Read some data; convert to 1D.      
 
  zm_phBF   = fall[:]->zm_phBF(:,:,:)    ; Read some data; convert to 1D.        
;    zm_dyBF   = fall[:]->zm_dyBF(:,:,:)    ; Read some data; convert to 1D.        
  zm_phAP   = fall[:]->zm_phAP(:,:,:)    ; Read some data; convert to 1D.        
  zmid_mpas = fall[:]->zmid_mpas(:,:,:)    ; Read some data; convert to 1D.        


  do t=0,nt-1
    print("===============t "+t+"=================")


    ;
    ; MPAS heights area relative to sea-level
    ;
    zm_mpas = zmid_mpas(0,:,:)
    zi_mpas = zint_mpas(0,:,:)
    do k=0,nlev-1
      zm_mpas(k,:)   = zmid_mpas(t,k,:)-zs
      zi_mpas(k,:)   = zint_mpas(t,k,:)-zs 
    end do

    do k=0,nlev-1
;      ave_zi_phBF(k) = sum(zi_phBF(t,k,:)*area(:))*inv_area
;      ave_zi_phAP(k) = sum(zi_phAP(t,k,:)*area(:))*inv_area
;      ave_zi_mpas(k) = sum(zi_mpas(t,k,:)*area(:))*inv_area

      ave_zm_phBF(k) = sum(zm_phBF(t,k,:)*area(:))*inv_area
      ave_zm_phAP(k) = sum(zm_phAP(t,k,:)*area(:))*inv_area
;      ave_zm_mpas(k) = sum(zm_mpas(t,k,:)*area(:))*inv_area 
    end do

    ;
    ; print resultss
    ;
;  do k=0,0 
;    tmp = ave_zi_mpas(k)-ave_zi_phBF(k)
;    abs_tmp = tmp/ave_zi_mpas(k)
;    print("k="+k)
;    print("zi diff dp_coupling (mpas zi minus zi phBF absolute and normalized) = "+tmp+" "+abs_tmp)
;  end do

    do k=0,nlev-1
      print("Looking at k="+k+" at MPAS height "+ave_zm_mpas(k))
      tmp = ave_zm_mpas(k)-ave_zm_phBF(k)
      abs_tmp = 100D0*tmp/ave_zm_mpas(k)
      print("k="+k)
      print("zm diff dp_coupling (mpas zm minus zm phBF absolute and %) = "+tmp+" "+abs_tmp)
      tmp = ave_zm_phAP(k)-ave_zm_phBF(k)
      abs_tmp = 100D0*tmp/ave_zm_phBF(k)
      print(ave_zm_phAP(k)+"  "+ave_zm_phBF(k))
      print("zm float in physics (zm_AP-zm_BF absolute and %)           = "+tmp+" "+abs_tmp)
      print("zm_phAP(k,0) "+zm_phAP(t,k,0))

      time_series_zm_phAP_minus_phBF(0,t,k) = tmp;abs_tmp
      time_series_zm_phAP_minus_phBF(1,t,k) = min(zm_phAP(t,k,:)-zm_phBF(t,k,:)) ;ave_ave_zm_phAP_minus_phBF(k)
      time_series_zm_phAP_minus_phBF(2,t,k) = max(zm_phAP(t,k,:)-zm_phBF(t,k,:)) ;ave_ave_zm_phAP_minus_phBF(k)
;      delete(tmp)
;      tmp = zm_phAP(k,:)-zm_phBF(k,:)
;      time_series_zm_phAP_minus_phBF(2,t,k) = tmp;min(tmp);max_zm_phAP_minus_phBF(k) 
;      time_series_zm_phAP_minus_phBF(3,t,k) = tmp;max_zm_phAP_minus_phBF(k) 
      delete(tmp)
    end do
  end do

;  do k=0,nlev-1
;    filename = "time_series_zm_phAP_minus_phBF_percent_dev_lev"+k+".txt"
;    print("writing file "+filename)
;    asciiwrite("data/"+filename,time_series_zm_phAP_minus_phBF(:,k))
;  end do

 ;
 ; check that zint consistency
 ;
;  do k=0,nlev-1
;    do i=0,ncol-1
;      tmp = zi_dyBF(k,i)-zint_mpas(k,i)
;      if (abs(tmp).gt.1.0) then
;        print("k="+k+" i="+i+"  "+tmp+"  "+zi_phBF(k,i)+"   "+zint_mpas(k,i)+"  "+PHIS(i)/9.80665)
;      end if
;     end do
;  end do
 
;  do k=0,nlev-1
;    do i=0,ncol-1
;      tmp = zm_dyBF(k,i)-zm_mpas(k,i)
;      if (abs(tmp).gt.1.0E-1) then
;        print("k="+k+" i="+i+"  "+tmp)
;        print("zm_dyBF   "+zm_dyBF(k,i)+" zm_mpas "+zm_mpas)
;      end if
;     end do
;  end do
;  xindexc = (/0,50000,100000,150000,200000,250000,300000,350000,400000/)
;  xlabelc = (/"0","50","100","150","200","250","300","350","400"/)

  number_of_plots = 17
  offset =15
  
 
;  plot = new(nlev,graphic)
  plot = new(number_of_plots,graphic)
  poly = new(16,graphic)   ; For attaching grey-filled rectangles

  wtype          = "pdf"
;  wtype@wkWidth  = 2500
;  wtype@wkHeight = 2500
  wks = gsn_open_wks(wtype ,"unique")                   ; ps,pdf,x11,png

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

  txres                       = True
  txres@txPerimOn             = True
  txres@txBackgroundFillColor = "Salmon"
  txres@txFontHeightF         = 0.006

  do k=0,number_of_plots-1
    idx = k+offset+1
    if ((k % 2) .eq. 0) then ; only works for integers
       res@tmYRLabelsOn         = True
       res@tmYLLabelsOn         = False
     else
        res@tmYRLabelsOn         = False
        res@tmYLLabelsOn         = True
     end if
    txt = "level "+idx
   txid = gsn_create_text(wks, txt, txres)
    plot(k) = gsn_csm_xy(wks,time,time_series_zm_phAP_minus_phBF(:,:,k+offset),res)
    annoid = gsn_add_annotation(plot(k), txid, amres)
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