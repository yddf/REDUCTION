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
;        c. Matt Giguere 2011.03.30 07:29:39 PM
;
;-
pro chi_dr_reduce, $
start_date = start_date, $
nnights = nnights, $
run = run, $
reduceonly = reduceonly, $
fitsonly = fitsonly, $
skipthid = skipthid, $
allofit = allofit

;110309 - 110325: rqa31
;110325 - 110407: rqa32
;110408 - :rqa33

if ~keyword_set(run) then run = 'rqa33'
if ~keyword_set(start_date) then start_date = 110410L
if ~keyword_set(nnights) then nnights = 2

datearr = lindgen(nnights) + long(start_date)
print, 'the date arr is: '
print, datearr
stop
for i=0, nnights - 1 do begin
print, '********************************************************'
print, ' NOW DOING ', strt(datearr[i])
print, '********************************************************'

chi_reduce, date = strt(datearr[i]), $
run = run, reduceonly = reduceonly, $
skipreduce = keyword_set(fitsonly), $
skipthid = keyword_set(skipthid), $
allofit = allofit
endfor




end;chi_dr_reduce.pro