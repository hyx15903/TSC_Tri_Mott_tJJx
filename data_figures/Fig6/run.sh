#!/bin/sh


gnuplot Gnuplot_Fig6.sh
latex Fig6.tex
dvips Fig6.dvi
ps2pdf Fig6.ps
