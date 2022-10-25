#!/bin/sh

gnuplot Gnuplot_Fig5.sh
latex Fig5.tex
dvips Fig5.dvi
ps2pdf Fig5.ps
