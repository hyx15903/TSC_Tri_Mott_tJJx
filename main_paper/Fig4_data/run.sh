#!/bin/sh

arg='SC1'

gnuplot Gnuplot_SC_density_correlation_${arg}.sh
latex SC_density_correlation_${arg}.tex
dvips SC_density_correlation_${arg}.dvi
ps2pdf SC_density_correlation_${arg}.ps
