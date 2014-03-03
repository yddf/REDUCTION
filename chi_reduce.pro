;+
;
;  NAME: 
;     chi_reduce
;
;  PURPOSE: 
;   To reduce the CHIRON data from each night
;
;  CATEGORY:
;      CHIRON
;
;  CALLING SEQUENCE:
;
;      chi_reduce, date, /keywords
;
;  INPUTS:
;
;		DATE: The date of the observing night. In the form 'yymmdd'. 
;		For example, chi_reduce, '110315'
;
;  KEYWORD PARAMETERS:
;
;		ALLOFIT: Use this keyword to run all 4 modes at the same
;		time. 
;
;		AUTOTHID: Use this keyword to read in the logsheets
;		and automatically determine what observations were
;		thar in each slit position. 
;
;		OBSNUM: The observations you are going to run THID.PRO
;		on. The first OBSNUM must be a file for which the .THID
;		file already exists. The subsequent OBSNUMs will run 
;		through thid producing the crude wavelength calibrations
;
;		CTHARS: The number of ThAr observations taken with the 
;		slicer
;
;		FTHARS: The number of ThAr observations taken with the 
;		fiber
;    
;		NTHARS: The number of ThAr observations taken with the
;		narrow slit
;
;		TTHARS: The number of ThAr observations taken with the 
;		normal slit
;    
;		REDUCEONLY: Use this keyword to ONLY reduce the data, 
;		and NOT run the rest of the procedure that solves for the
;		wavelength, calibrates, and checks. 
;	
;		SKIPFITS: Use this keyword to indicate NOT to convert
;		everything from iodspec to fits files at the end of 
;		the routine. 
;
;		SKIPREDUCE: Use this keyword to skip the reduction step
;
;		SKIPTHID: Use this keyword if you want to skip the 
;		wavelength solution calculations (THIDs)
;    
;  EXAMPLE:
;
;    To just convert the iodspecs to fits files:
;		chi_reduce, run='rqa31', /skipreduce, /skipthid, $
;					date='110312'
;
;		The normal way to run it at the end of the night. This 
;		reduces all of the data from each mode, runs thid on 
;		each ThAr exposure, and converts all of the iodspec
;		files to fits files and includes the wavelength
;		calibrations from thid in the fits files:
;
;		chi_reduce, /autothid, /allofit, run='rqa33', date='110415'		
;
;		
;
;
;  MODIFICATION HISTORY:
;        c. Matt Giguere 2011.03.24 06:28:32 PM
;		- 20110326: greatly increased the automation and
;		robustness of the code to handle all dates and all
;		slicer positions.  ~MJG
;		- 20110407 - MJG: started work on autothid keyword 
;		- 20110412 - MJG: added skipthid kw
;
;-
pro chi_reduce, $
allofit = allofit, $
autothid = autothid, $
bary = bary, $
cthars = cthars, $
date = date, $
doppler = doppler, $
fiber = fiber, $
fthars = fthars, $
narrow_slit = narrow_slit, $
nthars = nthars, $
reduceonly = reduceonly, $
run = run, $
skipfits = skipfits, $
skipreduce = skipreduce, $
skipthid = skipthid, $
slicer = slicer, $
slit = slit, $
tthars = tthars


if ~keyword_set(date) then date='110326'
if ~keyword_set(run) then run='rqa32'
print, 'The date is: ', date
print, 'The run is: ', run

if (keyword_set(allofit) OR keyword_set(reduceonly)) then begin
narrow_slit = 1
slit = 1
fiber = 1
slicer = 1
endif


if ~keyword_set(skipreduce) then begin
  if keyword_set(narrow_slit) then begin
     print, '**************************************************'
     print, 'Now reducing narrow slit spectra...'
     print, '**************************************************'
	 sorting_hat, date, run=run, /reduce, 	 mode = 'narrow'
  endif
  if keyword_set(slit) then begin
     print, '**************************************************'
     print, 'Now reducing slit spectra...'
     print, '**************************************************'
	 sorting_hat, date, run=run, /reduce, 	 mode = 'slit'
  endif
  if keyword_set(fiber) then begin
     print, '**************************************************'
     print, 'Now reducing fiber spectra...'
     print, '**************************************************'
	 sorting_hat, date, run=run, /reduce, 	 mode = 'fiber'
  endif
  if keyword_set(slicer) then begin
     print, '**************************************************'
     print, 'Now reducing slicer spectra...'
     print, '**************************************************'
	 sorting_hat, date, run=run, /reduce, 	 mode = 'slicer'
  endif
endif 

if ~keyword_set(reduceonly) then begin
if ~keyword_set(skipthid) then begin
if keyword_set(autothid) then begin
lfn = '/mir7/logsheets/'+strt(date)+'.log'
readcol, lfn, obsnmbr, objnm, i2, midtime, exptime, binsize, slitpos, $
delimiter=' ', f='A, A, A, A, D, A, A'

;find the thar exposures for each slit position:
nts = where(((strt(objnm) eq 'thar') and (strt(slitpos) eq 'narrow') and (strt(binsize) eq '3x1')))
tts = where(((strt(objnm) eq 'thar') and (strt(slitpos) eq 'slit') and (strt(binsize) eq '3x1')))
fts = where(((strt(objnm) eq 'thar') and (strt(slitpos) eq 'fiber') and (strt(binsize) eq '4x4')))
cts = where(((strt(objnm) eq 'thar') and (strt(slitpos) eq 'slicer') and (strt(binsize) eq '3x1')))

if strmid(run, 0, 1) eq 'q' then run = 'r'+run
;retrieve the thid files already processed:
spawn, 'ls -1 /mir7/thid/thidfile/'+run+'.*', res
;take the result, shorten it to just the observation number, then convert to long:
lress = long(strmid(res, 26, 4))
lastrunflag=0

;insert the most recently processed thid file (before these observations)
;as the one to use for reference:
initthar = where(lress lt long(obsnmbr[nts[0]]))

;stop

;this is for the case that we are just starting a new run and no thidfiles
;exist yet for the current run number. 
if res[0] eq '' then begin
runnum = long(strmid(run, 3, strlen(run) - 3))
runnum--
lastrun = strmid(run, 0, 3)+strt(runnum)

;retrieve the thid files already processed:
spawn, 'ls -1 /mir7/thid/thidfile/'+lastrun+'.*', res
;take the result, shorten it to just the observation nubmer, then convert to long:
lress = long(strmid(res[n_elements(res)-1], 26, 4))
initthar=0
lastrunflag=1
dumim = readfits('/mir7/fitspec/'+lastrun+'.'+strt(lress)+'.fits', hd, /silent)
endif;res=''


;**************************************************
;FIRST TO MAKE THE NARROW SLIT ARRAY:
;**************************************************
;This loop will cycle through the previous THARs looking for one taken
;with the same decker position as the array of interest:
if n_elements(lress) gt 1 then begin
dumvar = 1
ii=0L
while dumvar do begin
tharidx = initthar[n_elements(initthar) - 1L - ii]
dumim = readfits('/mir7/fitspec/'+run+'.'+strt(lress[tharidx], f='(I04)')+'.fits', hd, /silent)
ii++
if sxpar(hd, 'decker') eq 'narrow_slit' then dumvar = 0
endwhile
endif else tharidx = 1L
;the extra -1L compensates for the ii++ at the end of the DO..WHILE:
nthars = [string(lress[tharidx-1L]), obsnmbr[nts]]



print, '********************************************'
print, 'The narrow slit ThArs are: '
print, nthars
print, sxpar(hd, 'decker'), slitpos[nts]
print, '********************************************'

;stop

;**************************************************
;NEXT TO MAKE THE NORMAL SLIT ARRAY:
;**************************************************
;This loop will cycle through the previous THARs looking for one taken
;with the same decker position as the array of interest:
if n_elements(lress) gt 1 then begin
dumvar = 1
ii=0L
while dumvar do begin
tharidx = initthar[n_elements(initthar) - 1L - ii]
dumim = readfits('/mir7/fitspec/'+run+'.'+strt(lress[tharidx], f='(I04)')+'.fits', hd, /silent)
ii++
dpos = strt(sxpar(hd, 'decker'))
print, 'dpos is: ', dpos
if  dpos eq 'slit' then dumvar = 0
endwhile
endif else tharidx = 1L

;the extra -1L compensates for the ii++ at the end of the DO..WHILE:
tthars = [string(lress[tharidx-1L]), obsnmbr[tts]]

;tthars = [string(lress[initthar[n_elements(initthar)-1]]), obsnmbr[tts]]
print, '********************************************'
print, 'The normal slit ThArs are: '
print, tthars
print, sxpar(hd, 'decker'), slitpos[tts]
print, '********************************************'

;**************************************************
;THE FIBER ARRAY:
;**************************************************
;This loop will cycle through the previous THARs looking for one taken
;with the same decker position as the array of interest:
if n_elements(lress) gt 1 then begin
dumvar = 1
ii=0L
while dumvar do begin
tharidx = initthar[n_elements(initthar) - 1L - ii]
dumim = readfits('/mir7/fitspec/'+run+'.'+strt(lress[tharidx], f='(I04)')+'.fits', hd, /silent)
ii++
dpos = strt(sxpar(hd, 'decker'))
print, 'dpos is: ', dpos
if dpos eq 'fiber' then dumvar = 0
endwhile
endif else tharidx = 1L

;the extra -1L compensates for the ii++ at the end of the DO..WHILE:
fthars = [string(lress[tharidx-1L]), obsnmbr[fts]]

;fthars = [string(lress[initthar[n_elements(initthar)-1]]), obsnmbr[fts]]
print, '********************************************'
print, 'The fiber ThArs are: '
print, fthars
print, sxpar(hd, 'decker'), slitpos[fts]
print, '********************************************'


;**************************************************
;THE SLICER ARRAY:
;**************************************************
;This loop will cycle through the previous THARs looking for one taken
;with the same decker position as the array of interest:
if n_elements(lress) gt 1 then begin
dumvar = 1
ii=0L
while dumvar do begin
tharidx = initthar[n_elements(initthar) - 1L - ii]
dumim = readfits('/mir7/fitspec/'+run+'.'+strt(lress[tharidx], f='(I04)')+'.fits', hd, /silent)
ii++
dpos = strt(sxpar(hd, 'decker'))
print, 'dpos is: ', dpos
if dpos eq 'slicer' then dumvar = 0
endwhile
endif else tharidx = 1L

;the extra -1L compensates for the ii++ at the end of the DO..WHILE:
cthars = [string(lress[tharidx-1L]), obsnmbr[cts]]

;cthars = [string(lress[initthar[n_elements(initthar)-1]]), obsnmbr[cts]]
print, '********************************************'
print, 'The slicer ThArs are: '
print, cthars
print, sxpar(hd, 'decker'), slitpos[cts]
print, '********************************************'

endif;autothid
;stop


obsnum=0
for jm = 0, 3 do begin ;cycle through the # modes
nskip = 1
print, 'jm is: ', jm
;stop
if jm eq 0 then begin 
if keyword_set(nthars) then obsnum = nthars else nskip = 0
endif
if jm eq 1 then begin 
if keyword_set(tthars) then obsnum = tthars else nskip = 0
endif
if jm eq 2 then begin 
if keyword_set(fthars) then obsnum = fthars else nskip = 0
endif
if jm eq 3 then begin 
if keyword_set(cthars) then obsnum = cthars else nskip = 0
endif

print, 'jm is: ', jm
print, 'n skip is: ', nskip
print, 'obsnum is: ', obsnum
;stop
if nskip then begin
print, '**************************************************'
print, 'Now to run THID.PRO to get the wavelength solution'
print, ' using all of the ThAr observations with the input'
print, ' decker position...'
print, '**************************************************'

print, '*************************************************'
print, 'Now on ThAr SET: ',strt(jm+1),' of 4'
print, 'Obsnum is: ', obsnum
;print, 'The decker position for this set is: ', slitpos[nts]
print, '*************************************************'

for i=1, n_elements(obsnum) - 1 do begin


if ~lastrunflag then begin
  restore, '/mir7/thid/thidfile/'+run+'.'+strt(obsnum[i-1], f='(I04)')+'.thid'
endif
if (lastrunflag and (i eq 1)) then begin
  restore, '/mir7/thid/thidfile/'+lastrun+'.'+strt(obsnum[i-1], f='(I04)')+'.thid'
endif
  rdsk, thar, '/mir7/iodspec/'+run+'.'+strt(obsnum[i], f='(I04)'), 1
  
  print, 'NOW CALCULATING THE WAVELENGTH SOLUTION FOR: ', $
  				'/mir7/iodspec/'+run+'.'+strt(obsnum[i], f='(I04)')
  
  thid, thar, 84d, 84d * [6654d, 6734d], wvc, thid, $
	  init = thid.wvc, /orev
  ;a for auto
  ;15
  ;5
  ;9
  ;f
  ;6
  ;6
  ;o
  ;drag and drop above to remove outliers
  ;drag and drop below to remove outliers
  ;click margin twice
  ;repeat the above a two more times
  ;f
  ;6
  ;6
  ;o
  ;drag and drop above to remove outliers
  ;drag and drop below to remove outliers
  ;click margin twice
  ;q to quit
  ;(type cmd: h at any point for help)
  
  save, thid, file='/mir7/thid/thidfile/'+run+'.'+strt(obsnum[i], f='(I04)')+'.thid'
  
  mkwave, w, thid.wvc
  
  w = reverse(w, 2)
  
  save, w, file='/mir7/thid/wavfile/ctio_'+run+'.'+strt(obsnum[i], f='(I04)')+'.dat'

endfor ;run through thids
endif ;nskip
endfor ;run through modes

endif;KW: skipthid

;stop
if ~keyword_set(skipfits) then begin
print, '**************************************************'
print, 'Now to match wavelengths to the reduced spectra...'
print, '**************************************************'

print, '**********************'
print, 'NARROW SLIT IOD2FITS'
print, '**********************'
sorting_hat, date, run=run, mode='narrow', /iod2fits

print, '**********************'
print, 'NORMAL SLIT IOD2FITS'
print, '**********************'
sorting_hat, date, run=run,  /normslit,  /iod2fits

print, '**********************'
print, 'SLICER IOD2FITS'
print, '**********************'
sorting_hat, date, run=run,  mode='slicer', /iod2fits

print, '**********************'
print, 'FIBER IOD2FITS'
print, '**********************'
sorting_hat, date, run=run, mode='fiber',  /iod2fits
endif ;KW:skipfits
;stop

print, '**************************************************'
print, 'Finally to check to make sure everything has been run...'
print, '**************************************************'
sorting_hat, date, run=run, $
narrow=narrow_slit, $
normslit = slit, $
slicer = slicer, $
fiber = fiber, $
/end_check

endif;KW:reduceonly

if keyword_set(bary) then begin
print, '**************************************************'
print, 'NOW RUNNING THE DOPPLER CODE...'
print, '**************************************************'

qbarylog, '/mir7/logsheets/'+date+'.log', prefix=run
barystruct_dbl, observatory='ctio'
endif ;KW(bary)


if keyword_set(doppler) then begin
print, '**************************************************'
print, 'NOW RUNNING THE DOPPLER CODE...'
print, '**************************************************'

print, '**********************'
print, 'NARROW SLIT DOPPLER CODE'
print, '**********************'
sorting_hat, date, run=run,  mode='narrow', /doppler, doptag = 'kd'

print, '**********************'
print, 'SLICER DOPPLER CODE'
print, '**********************'
sorting_hat, date, run=run,  mode='slicer',  /doppler, doptag = 'ks'

endif ;KW(doppler)


;stop
end;chi_reduce.pro
