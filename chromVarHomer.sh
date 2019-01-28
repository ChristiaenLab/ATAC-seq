#!/bin/bash

#./rscript.s chromVarHomer.R --peaks 2018-06-12/peakome.bed --out peakome
#./rscript.s chromVarHomer.R --peaks 2018-06-12/TSC.bed --out TSC
#./rscript.s chromVarHomer.R --peaks 2018-06-12/TSCpromoter.bed --out TSCpromoter
#./rscript.s chromVarHomer.R --peaks 2018-06-12/TSCpromoter.bed --out TSCcontrol --condition control --tissue "mesenchyme,B7.5" --noKO
#./rscript.s chromVarHomer.R --peaks 2018-06-12/TSCpromoter.bed --out laczTSC --condition control --tissue "B7.5" --noKO
#./rscript.s chromVarHomer.R --peaks 2018-06-12/TSCpromoter.bed --out gfpTSC --condition control --tissue mesenchyme
#./rscript.s chromVarHomer.R --peaks 2018-06-12/TSCpromoter.bed --out "6.10.TSC" --condition control --tissue "mesenchyme,B7.5" --time "6hpf,10hpf" --noKO
#./rscript.s chromVarHomer.R --peaks 2018-06-12/TSCpromoter.bed --out "6.10.18.TSC" --condition control --tissue "mesenchyme,B7.5" --time "6hpf,10hpf,18hpf" --noKO
#./rscript.s chromVarHomer.R --peaks 2018-06-12/TSCpromoter.bed --out "18.TSC" --condition control --tissue "mesenchyme,B7.5" --time "18hpf"
./rscript.s chromVarHomer.R --peaks 2019-01-18/denovoCardiacPeaks.bed --out denovoCardiac  --tissue "B7.5"
./rscript.s chromVarHomer.R --peaks 2019-01-18/denovoASMpeaks.bed --out denovoASM  --tissue "B7.5"

./rscript.s chromVarHomer.R --peaks 2018-10-23/daTime.bed --out time --condition "control" --tissue "B7.5" --noKO
./rscript.s chromVarHomer.R --peaks 2018-10-23/daMesp.bed --out mesp --tissue "B7.5" --time "6hpf,10hpf"
./rscript.s chromVarHomer.R --peaks 2018-10-23/daHandr.bed --out handr --tissue "B7.5" --time "18hpf"
./rscript.s chromVarHomer.R --peaks 2018-10-23/daTissue.bed --out tissue --tissue "mesenchyme,B7.5" --time "10hpf,18hpf" --condition "control" --noKO

./rscript.s chromVarHomer.R --peaks 2018-10-23/daTime.bed --out time_all  --tissue "B7.5"
./rscript.s chromVarHomer.R --peaks 2018-10-23/daMesp.bed --out mesp_all --tissue "B7.5"
./rscript.s chromVarHomer.R --peaks 2018-10-23/daHandr.bed --out handr_all --tissue "B7.5"
./rscript.s chromVarHomer.R --peaks 2018-10-23/daTissue.bed --out tissue_all --tissue "mesenchyme,B7.5"
