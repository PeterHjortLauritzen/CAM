#!/bin/tcsh
if ( "$#argv" != 1) then
  echo "Wrong number of arguments specified:"
  echo "  -arg 1 file with energy diagnostics (averaged)"
  exit
endif
set n = 1
set filez = "$argv[$n]" 

echo $filez

ncl 'dir="'$PWD'"' 'fnamez="'$filez'"' z2d-tmp.ncl  

