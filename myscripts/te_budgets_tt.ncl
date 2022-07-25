load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

begin

  f = addfile("h1.nc","r")

  TT_phB1 = f->WI_phB1
  TT_phA1 = f->WI_phA1

  area    = f->area(:)
  area_sphere = sum(area)

  nt = dimsizes(TT_phB1(:,0))

  int_vars  = new((/2,nt/)     ,typeof(TT_phB1(0,0)))        

  do t=0,nt-1
     print(" ")
     print(" ")
     print("t = "+t)

     int_vars(0,t) = sum(area(:)*TT_phB1(t,:))/area_sphere	  
     print("WI_phB1="+int_vars(0,t))
     int_vars(1,t) = sum(area(:)*TT_phA1(t,:))/area_sphere	  
     print("WI_phA1="+int_vars(1,t))
     print(" ")
     print("rel. mass change = "+(int_vars(1,t)-int_vars(0,t))/int_vars(0,t))
  end do


end
