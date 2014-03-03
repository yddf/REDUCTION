;+
;  NAME: 
;     chi_redistrib
;  PURPOSE: 
;   
;
;  CATEGORY:
;      CHIRON
;
;  CALLING SEQUENCE:
;
;      chi_redistrib
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
;      chi_redistrib, month = '1111'
;
;  MODIFICATION HISTORY:
;        c. Matt Giguere 2011.10.25 10:09:48 AM
;
;-
pro chi_redistrib, $
year = year, $
month = month, $
date = date, $
help = help
print, 'Now starting chi_redistrib @: ', systime()
if keyword_set(help) then begin
	print, '*************************************************'
	print, '*******HELP INFORMATION FOR CHI_REDISTRIB*******'
	print, '*************************************************'
	print, 'KEYWORDS: '
	print, ''
	print, 'HELP: Use this keyword to print all available arguments'
	print, ''
	print, 'DATE: A date of array of dates in the yymmdd format that'
	print, '      chi_redistrib will redistribute the data for.'
	print, ''
	print, 'MONTH: The entire month, in the yymm STRING format, that you would'
	print, '       like chi_redistrib to redistribute the data for'
	print, ''
	print, 'YEAR: The entire year, in the yy format, that you would'
	print, '      like chi_redistrib to redistribute the data for.'
	print, '*************************************************'
endif

if keyword_set(year) then begin
	spawn, 'ls -1d /mir7/raw/'+strt(long(yy), f='(I02)')+'*', date 
endif ;kw:year

if keyword_set(month) then begin
   print, 'MONTH ENTERED: ', month
	if strlen(month) lt 4 then month=calform(systime(/julian), f='yy')+strt(long(month), f='(I02)')
	spawn, 'ls -1d /mir7/raw/'+month+'*', date 
	print, 'MONTH NOW: ', month
endif ;kw:month

if ~keyword_set(date) and ~keyword_set(month) then begin
month='1110'
date = [ $
'111014', $
'111018', $
'111020']
endif;kw:date

if month eq '1109' then begin
date = [ $
'110902', $
'110904', $
'110906', $
'110908', $
'110910', $
'110912', $
'110914', $
'110918', $
'110919', $
'110920', $
'110922', $
'110924', $
'110926', $
'110928', $
'110930']
endif

if month eq '1111' then begin
date = [ $
'111102', $
'111104', $
'111106', $
'111107', $
'111108', $
'111109', $
;'111110', $ <--NO DATA
'111111', $
'111112', $
'111114', $
'111117', $
'111118', $
'111120', $
'111122', $
'111124', $
'111126', $
'111127', $
'111128']
endif

if month eq '1203' then begin
date = [ $
'120302', $
'120305', $
'120306', $
'120307', $
'120309']
endif

print, 'The date array is: '
print, transpose(date)
print, '*************************************************'
for i=0, n_elements(date)-1 do begin
	chi_que_distrib, date=date[i]
endfor

stop, '~~END OF CHI_REDISTRIB~~'
end;chi_redistrib.pro
