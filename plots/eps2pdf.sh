#!/bin/tcsh -f
foreach file_eps (*.eps)

setenv file_pdf `echo ${file_eps} | sed -e s/.eps/.pdf/g` 

#echo ${file_eps}' '${file_pdf}

gs -o ${file_pdf} -sDEVICE=pdfwrite -dEPSCrop ${file_eps}

end
