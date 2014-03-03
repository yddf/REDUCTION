;***********************************************************************
;+
; NAME:  NEXTNAMEEPS.PRO
;																	   
; PURPOSE: This procedure will check to see if a name already exists. 
;	If the name already exists, it will add the next available integer
;	to the end of the name. 
;																	   
; CATEGORY: UTILITIES							   
;																	   
; MODIFICATION HISTORY:												   
;     c. Matt Giguere, Thursday, June 05, 2008		
;-
;***********************************************************************
function nextnameeps, innm, suf, nosuf=nosuf
if size(suf, /typ) ne 7 then suf = ''

int = 2
if ~(file_test(innm+suf+'.eps')) then begin
fnlnm = innm+suf
xst = 0 
endif else xst = 1
while xst eq 1 do begin
fnlnmeps = innm+strt(int)+suf+'.eps'
print, 'fnlnmeps is: ', fnlnmeps
xst = file_test(fnlnmeps)
fnlnm=innm+strt(int++)+suf
endwhile
print, 'fnlnm is: ', fnlnm
return, fnlnm
end; fnlnm.pro