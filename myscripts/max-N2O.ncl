load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

begin
  f = addfile("h0.nc","r")

  N2O = f->N2O
  nt = dimsizes(N2O(:,0,0))
  do i=0,nt-1
    print(i+"  "+max(N2O(i,:,:)))
  end do
  print("max in lowest level")
  nlev = dimsizes(N2O(0,:,0))
    do i=0,nt-1
    print(i+"  "+max(N2O(i,nlev-1,:)))
  end do

end
