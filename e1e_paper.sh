#!/bin/tcsh -f
rm e1e_paper.pdf
rm e1e_paper.aux
rm e1e_paper.bbl
#rm text/aliases/*.aux text/conclusions/*.aux text/event_selection/*.aux text/motivation/*.aux text/setups/*.aux text/title_page/*.aux text/toc/*.aux text/appendix/*.aux text/rad_corr/rad_corr.aux text/sys_err/*.aux text/bin_corr/*.aux text/cross_sction/*.aux
pdflatex e1e_paper
bibtex e1e_paper
pdflatex e1e_paper
pdflatex e1e_paper
acroread e1e_paper.pdf >& /dev/null
rm e1e_paper.log e1e_paper.aux e1e_paper.toc e1e_paper.out e1e_paper.blg e1e_paper.bbl e1e_paperNotes.bib
