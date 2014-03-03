;***********************************************************************
; NAME:  NEXTNAME.PRO
;																	   
; PURPOSE: This procedure will check to see if a name already exists. 
;	If the name already exists, it will add the next available integer
;	to the end of the name. 
;																	   
; CATEGORY: EXOPLANETS							   
;																	   
; CALLING SEQUENCE:													   
;																	   
;																	   
; INPUTS:															   
;
; OPTIONAL INPUTS:													   
;																	   
; KEYWORD PARAMETERS:	
;
; OUTPUTS:															   
;																	   
; OPTIONAL OUTPUTS:													   
;																	   
; SIDE EFFECTS:														   
;																	   
; RESTRICTIONS:														   
;																	   
; PROCEDURE:														   
;																	   
; EXAMPLE:			
;																	   
; MODIFICATION HISTORY:												   
;     c. Matt Giguere, Thursday, June 05, 2008		
;***********************************************************************
function nextname, innm, suf
if size(suf, /typ) ne 7 then suf = '.dat'

int = 2
if ~(file_test(innm+suf)) then begin
fnlnm = innm+suf
xst = 0 
endif else xst = 1
while xst eq 1 do begin
fnlnm = innm+'_'+strt(int++)+suf
xst = file_test(fnlnm)
endwhile
return, fnlnm
end; fnlnm.pro