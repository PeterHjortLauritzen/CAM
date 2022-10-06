#!/bin/tcsh
if ( "$#argv" != 2) then
  echo "Wrong number of arguments specified:"
  echo "  -arg 1 file with energy diagnostics (averaged)"
  echo "  -arg 2 time index"
  exit
endif
set n = 1
set filez = "$argv[$n]" 
set n = 2
set timei = "$argv[$n]" 

echo $filez

ncl 'dir="'$PWD'"' 'fnamez="'$filez'"' 'time_idx='$timei'' z2d.ncl  

