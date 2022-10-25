#!/bin/sh

gnuplot Gnuplot_Fig4.sh
latex Fig4.tex
dvips Fig4.dvi
ps2pdf Fig4.ps
