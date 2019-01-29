#!/bin/bash

# data for Fig. 2A
./rscript.s chromVarHomer.R --peaks daMesp.bed --out mesp --tissue "B7.5" --time "6hpf,10hpf"

# data for Fig. S20
./rscript.s chromVarHomer.R --peaks denovoCardiacPeaks.bed --out denovoCardiac  --tissue "B7.5"
./rscript.s chromVarHomer.R --peaks denovoASMpeaks.bed --out denovoASM  --tissue "B7.5"
