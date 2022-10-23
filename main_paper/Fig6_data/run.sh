#!/bin/sh

arg=''

gnuplot Gnuplot_Nature_Phase_transition${arg}.sh
latex Nature_Phase_tranistion${arg}.tex
dvips Nature_Phase_tranistion${arg}.dvi
ps2pdf Nature_Phase_tranistion${arg}.ps
