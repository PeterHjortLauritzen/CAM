#!/bin/tcsh
if ( "$#argv" != 1) then
  echo "Wrong number of arguments specified:"
  echo "  -arg 1 time index to plot"
endif  
set n = 1
set time_idx = "$argv[$n]"
set vname = "(/"\"OMEGA500\""/)"
echo $vname
set hnum = "h0.ave"
#source set_datafile_location.sh
ncl 'time_idx="'$time_idx'"' 'vname='$vname''  'hnum="'$hnum'"' 'outfname="'modified-FHS94.OMEGA.pdf'"'  plot_2column_panel.ncl
#cd work/
#pdfcrop --margins 25 OMEGA.pdf
#ls OMEGA-crop.pdf
#mv OMEGA-crop.pdf ../plots/OMEGA.pdf
#cd ..

