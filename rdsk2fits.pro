;+
;
;  NAME: 
;     rdsk2fits
;
;  PURPOSE: To restore an RDSK file and save it as a FITS file.
;   
;
;  CATEGORY:
;      UTILITIES
;
;  CALLING SEQUENCE:
;
;      rdsk2fits
;
;  INPUTS:
;
;  OPTIONAL INPUTS:
;
;  OUTPUTS:
;
;  OPTIONAL OUTPUTS:
;
;  KEYWORD PARAMETERS:
;    
;  EXAMPLE:
;      rdsk2fits
;
;  MODIFICATION HISTORY:
;        c. Matt Giguere 2013.08.13 10:59:43
;
;-
pro rdsk2fits, $
filename = filename, $
rdskfname = rdskfname, $
data = data, $
header = header

if ~keyword_set(header) then header=''
if keyword_set(rdskfname) then rdsk, data, rdskfname

writefits, filename, data, header

end;rdsk2fits.pro