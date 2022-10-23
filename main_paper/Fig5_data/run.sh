#!/bin/sh

arg='_v4'

gnuplot Gnuplot_phase_transition${arg}.sh
latex Phase_transition${arg}.tex
dvips Phase_transition${arg}.dvi
ps2pdf Phase_transition${arg}.ps
