#!/bin/tcsh -f

foreach file (*.eps)

setenv file_pdf `echo $file | sed -e 's/.eps/.pdf/g'`

echo $file_pdf

gs -o $file_pdf -sDEVICE=pdfwrite -dEPSCrop $file


end
