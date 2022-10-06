#!/bin/tcsh
if ( "$#argv" != 2) then
  echo "Wrong number of arguments specified:"
  echo "  -arg 1 file with energy diagnostics (averaged)"
  echo "  -arg 2 file with z diagnostics (averaged)"
  exit
endif
set n = 1
set file = "$argv[$n]" 

set n = 2
set filez = "$argv[$n]" 

echo $file
echo $filez

if (! -e nuopc.runconfig) then
  echo "This script expects nuopc.runconfig to get physics time-step"
  exit
endif
grep "atm_cpl_dt[[:space:]]" nuopc.runconfig >> tmp_file
set dtime = `grep -o "[0-9][0-9][0-9][0-9]" tmp_file || echo "none"`
echo "dtime ="$dtime
rm tmp_file

ncl 'dir="'$PWD'"' 'fname="'$file'"' 'fnamez="'$filez'"' 'dtime='$dtime'' te_budgets.ncl  

