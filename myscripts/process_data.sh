#!/bin/tcsh
if ( "$#argv" != 1) then
  echo "Wrong number of arguments specified:"
  echo "  -arg 1 caze dir on Cheyenne"
endif
set n = 1
unset res
set caze = "$argv[$n]" 
rm *.cam.h1.*.nc
rm *.cam.h0.0001-01.nc
source ~/git-scripts/ncl_scripts/ave.sh $caze h0 ANN
#ncks -v OMEGA500,OMEGA850,PSL $caze.ave.ANN.h0.nc $caze.ave.h0.nc
ncks -v OMEGA500,OMEGA850 $caze.ave.ANN.h0.nc $caze.ave.h0.nc
rm $caze.ave.ANN.h0.nc
