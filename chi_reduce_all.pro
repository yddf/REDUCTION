;+
;
;  NAME: 
;     chi_reduce_all
;
;  PURPOSE: 
;   A procedure designed to run Andrei's sorting hat routine for all modes
;
;  CATEGORY:
;      CHIRON
;
;  CALLING SEQUENCE:
;
;      chi_reduce_all
;
;  KEYWORD PARAMETERS:
;    
;  EXAMPLE:
;      chi_reduce_all
;
;  MODIFICATION HISTORY:
;        c. Matt Giguere 2011.11.18 11:05:47 AM
;
;-
pro chi_reduce_all, $
help = help, $
date = date, $
doppler = doppler, $
skipbary = skipbary, $
skipdistrib = skipdistrib, $
skipqc = skipqc


if keyword_set(help) then begin
	print, '*************************************************'
	print, '*************************************************'
	print, '        HELP INFORMATION FOR chi_reduce_all'
	print, 'KEYWORDS: '
	print, ''
	print, 'HELP: Use this keyword to print all available arguments'
	print, ''
	print, ''
	print, ''
	print, '*************************************************'
	print, '                     EXAMPLE                     '
	print, "IDL>"
	print, 'IDL> '
	print, '*************************************************'
	stop
endif

print, 'Started @: ', systime()
spawn, 'hostname', host

if host eq 'ctimac1.ctio.noao.edu' then begin
  rawdir = '/mir7/raw/'
endif else begin
  rawdir = '/raw/mir7/'
endelse
lfn = '/tous/mir7/logsheets/20'+strmid(date, 0, 2)+'/'+strt(date)+'.log'

;This part gets the image prefix:
spawn, 'ls -1 '+rawdir+date+'/', filearr
nel = n_elements(filearr)
nfa = strarr(nel)
for i=0, nel-1 do nfa[i] = strmid(filearr[i], 0, strlen(filearr[i])-9)
uniqprefs =  nfa(uniq(nfa))
pref = 'junk'
ii=-1
repeat begin
  ii++
  pref = uniqprefs[ii]
  print, 'pref is: ', uniqprefs[ii]
  print, 'pref 1st 2 are: ', strmid(uniqprefs[ii], 0,2)
endrep until ( ((strmid(uniqprefs[ii],0,2) eq 'qa') or $
        (strmid(uniqprefs[ii],0,2) eq 'ch')) and $
        (strmid(uniqprefs[ii], 0, 4) ne 'chir') and $
        (strmid(uniqprefs[ii], 0, 1, /reverse) eq '.'))

print, 'pref is: ', pref
;stop, 'pref is: ', pref

modearr = [$
'narrow', $ 
'slicer', $
'slit', $
'fiber']

for i=0, 3 do begin
  print, '*************************************************'
  print, ' NOW ON TO THE ', MODEARR[I], ' MODE...'
  print, '*************************************************'
  sorting_hat,date,run=pref,mode=modearr[i],/reduce,/getthid,/iod2fits
endfor

if ~keyword_set(skipbary) then begin
  print, '**************************************************'
  print, 'NOW ADDING OBSERVATIONS TO THE BARYLOG...'
  print, '**************************************************'

  ;No Dot Prefix (ndpref) - The prefix without the dot (e.g. 'chi120402')
  ndpref = strmid(pref, 0, strlen(pref)-1)
  qbarylog, lfn, prefix=ndpref
  chi_barystruct, asciifn='/tous/mir7/bary/qbcvel.ascii', $
  				structfn='/tous/mir7/bary/qbcvel.dat'
endif ;~KW(skipbary)


if keyword_set(doppler) then begin
  print, '**************************************************'
  print, 'NOW RUNNING THE DOPPLER CODE...'
  print, '**************************************************'

  print, '**********************'
  print, 'NARROW SLIT DOPPLER CODE'
  print, '**********************'
  sorting_hat, date, run=run, $
  mode='narrow', $
  /doppler, doptag='t'

  print, '**********************'
  print, 'SLICER DOPPLER CODE'
  print, '**********************'
  ;sorting_hat, date, run=run, $
  ;mode='slicer', $
  ;/doppler, doptag = 'ks'
endif ;KW(doppler)

if ~keyword_set(skipdistrib) then begin
  print, '******************************'
  print, 'NOW DISTRIBUTING THE DATA'
  print, '******************************'
  print, systime()
  chi_que_distrib, date=date
endif;KW(skipdistrib)

if ~keyword_set(skipqc) then begin
  print, '******************************'
  print, 'NOW RERUNNING QC'
  print, '******************************'
  print, systime()
  chi_quality, date=date
endif;KW(skipqc)

print, 'Finished @: ', systime()

end;chi_reduce_all.pro