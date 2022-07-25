load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

begin
  f = addfile("h0.nc","r")
  N2O = f->N2O


  N2O_BC1 = f->N2O_BC1
  N2O_BC2 = f->N2O_BC2


  N2O_AC1 = f->N2O_AC1
  N2O_AC2 = f->N2O_AC2
  N2O_AC3 = f->N2O_AC3
  N2O_AC4 = f->N2O_AC4
  N2O_AC5 = f->N2O_AC5
  N2O_AC6 = f->N2O_AC6
  N2O_AC7 = f->N2O_AC7
  N2O_AC8 = f->N2O_AC8
  N2O_AC9 = f->N2O_AC9
  N2O_AC10 = f->N2O_AC10
  N2O_AC11 = f->N2O_AC11
  N2O_AC11b = f->N2O_AC11b
  N2O_AC11c = f->N2O_AC11c
  N2O_AC11d = f->N2O_AC11d
  N2O_AC11e = f->N2O_AC11e
  N2O_AC11f = f->N2O_AC11f
  N2O_AC11g = f->N2O_AC11g
  N2O_AC11h = f->N2O_AC11h
  N2O_AC12 = f->N2O_AC12
  N2O_AC13 = f->N2O_AC13

  nt = dimsizes(N2O(:,0,0,0))

  print("max "+max(N2O))

  do i=0,nt-1
     print(" ")
     print("time stamp "+i)
     print(" ")


     print(" ")
     print("begin tphysbc                 "+max(N2O_BC1(i,:,:,:)))
     print("end   tphysbc                 "+max(N2O_BC2(i,:,:,:)))
     print("begin tphysac                 "+max(N2O_AC1(i,:,:,:)))
     print("Carma emission                "+max(N2O_AC2(i,:,:,:)))
     print("Carma                         "+max(N2O_AC3(i,:,:,:)))
     print("Contrails                     "+max(N2O_AC4(i,:,:,:)))
     print("CLUBB MG                      "+max(N2O_AC5(i,:,:,:)))
     print("Wet scavenging                "+max(N2O_AC6(i,:,:,:)))
     print("CArma wet depo                "+max(N2O_AC7(i,:,:,:)))
     print("ZM mixing                     "+max(N2O_AC8(i,:,:,:)))
     print("radiation+aoa tracer          "+max(N2O_AC9(i,:,:,:)))
     print("co2 cyelc+chemistry           "+max(N2O_AC10(i,:,:,:)))
     print("Vertical diffusion            "+max(N2O_AC11(i,:,:,:)))
     print("rayleigh_friction_tend        "+max(N2O_AC11b(i,:,:,:)))
     print("aero_model_drydep             "+max(N2O_AC11c(i,:,:,:)))
     print("carma_do_aerosol              "+max(N2O_AC11d(i,:,:,:)))
     print("gw_tend                       "+max(N2O_AC11e(i,:,:,:)))
     print("qbo_relax                     "+max(N2O_AC11f(i,:,:,:)))
     print("lunar_tides_tend              "+max(N2O_AC11g(i,:,:,:)))
     print("waccmx stuff                  "+max(N2O_AC11h(i,:,:,:)))
     print("12                            "+max(N2O_AC12(i,:,:,:)))
     print("DME_ADJUST and end of tphysac "+max(N2O_AC13(i,:,:,:)))


  end do
end
