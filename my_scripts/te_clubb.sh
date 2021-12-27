#!/bin/tcsh
if ( "$#argv" != 2 ) then
  echo "Wrong number of arguments specified:"
  echo "  -arg 1 reference run"
  echo "  -arg 2 scale surface flux run"
  exit
endif
set n = 1
set fileA = "$argv[$n]"
set n = 2
set fileB = "$argv[$n]"

set dir='/glade/u/home/pel/src/cam_energy_analysis/my_scripts'
ncl 'dir="'$PWD'"' 'fname_ref="'$fileA'"' 'fname_flx="'$fileB'"' $dir/te_clubb.ncl  

