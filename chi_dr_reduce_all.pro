;+
;
;  NAME: 
;     chi_dr_reduce
;
;  PURPOSE: 
;   This routine will cycle through a bunch of consecutive nights 
;		reducing all of the data. 
;
;  CATEGORY:
;      CHIRON
;
;  CALLING SEQUENCE:
;
;      chi_dr_reduce
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
;		To reduce everything, convert everything to fits, but skip
;		THID creation, and then check to make sure everything has
;		been processed at the end, type:
;
;      chi_dr_reduce, start_date='110410', nnights=1, run='rqa33', $
;							/allofit, /skipthid
;
;  MODIFICATION HISTORY:
;        c. Matt Giguere 2012.06.12
;
;-
pro chi_dr_reduce_all

datearr = [$
'120310', $
'120312', $
'120314', $
'120315', $
'120318', $
'120320', $
'120322', $
'120324', $
;'120325', $ <<-- NO DATA
'120326', $
'120327', $
'120328', $
'120329', $
'120330', $
'120331']

for i=0, n_elements(datearr)-1 do begin
print, '******************************************'
print, '******************************************'
print, '******************************************'
print, 'CHI_DR_REDUCE_ALL.PRO'
PRINT, 'NOW ON DATE: ', datearr[i]
print, i/n_elements(i)+' % done...'
print, '******************************************'
print, '******************************************'
print, '******************************************'
print, '******************************************'
  chi_reduce_all, date = datearr[i]
endfor

end;chi_dr_reduce_all.pro