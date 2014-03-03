;+
; PURPOSE: To find the latest THID file for a given night
;
; INPUT: night='110820', parfile, [mode=mode]
; run is needed for the old format, otherwise run  'chi'+night
; if mode is specified, only THID files for this mode will be searched
;
; OUTPUT: found = name of the file
;
; CREATED:
; Oct 17, 2011 AT
;
; MODIFICATION HISTORY:
; *routine wasn't working since the 2012 upgrade. 
; Changed a number of things, removed hard coding 
; and got it up and running. 20120504 ~MJG
; *Modified the file search line to work with multiple years 20130111 ~MJG
; Made sure thar lamp was in for observation, got rid of 1 goto. 20130901 ~MJG
;-

pro findthid, night, redpar, thidfile, mode=mode, run=run

depth = 10 ; how many nights back to look?
thidfile='none'
logdir = redpar.rootdir+redpar.logdir
if ~keyword_set(mode) then mode = redpar.modes[redpar.mode]
if ~keyword_set(run) then run = redpar.prefix
if ~keyword_set(night) then night = redpar.date
pfxtg = redpar.prefix_tag

;The extra "../*/" takes care of the year directory:
logs = file_search(logdir+'../*/*.log', count=nlogs)
if nlogs eq 0 then return ; No logfiles
print, nlogs,' logsheets found by FINDTHID'
logs = logs[sort(logs)]

start = (where(logdir+night+'.log' eq logs))[0] 
if start lt 0 then stop, 'Night '+night+' is not found in the logs!'

found = 0
lookback = (start - depth) > 0
; start with the current night
for i=start, lookback ,-1 do begin ; try all nights backwards
   ;for the case where the automatic hourly log generator is currently 
   ;creating the logs, wait 180 seconds for it to finish before proceeding:
   spawn, 'tail '+logs[i], logout
   if logout[0] eq '' then wait, 180
   ; read logsheets
   readcol, logs[i], skip=9, obnm, objnm, i2, mdpt, exptm, bin, slit,  f='(a5, a13, a4, a14, a8, a3, a6)'
   if keyword_set(mode) then sel = where((objnm eq 'thar') and (slit eq mode)) else  sel = where(objnm eq 'thar')
   if sel[0] lt 0 then continue ; nothing for this night! 

   curnight = strmid(logs[i], strlen(logdir))
   curnight = strmid(curnight, 0, strpos(curnight,'.log'))
   crun = 'chi'+curnight

   if strpos(crun,'.') lt 0 then crun=crun+'.' ; add the point

   thidfile =''
   j=0L
   while j lt n_elements(sel)-1 and thidfile eq '' do begin ; search thids
	 fn = redpar.rootdir+redpar.thidfiledir+pfxtg+crun+obnm[sel[j]]+'.thid'
	 print, 'Looking for '+fn
	 thids = file_search(fn, count=found)
	 rawfn = redpar.rootdir+redpar.rawdir+curnight+'/'+crun+obnm[sel[j]]+'.fits'
	 header = headfits(rawfn)
	 if strt(fxpar(header, 'COMPLAMP')) ne 'TH-AR' then begin
	 	print, 'WARNING! NO LAMP IN FOR: '
	 	print, rawfn
	 	print, 'TYPE THE IDL COMMAND: '
	 	print, "chi_junk, date='"+curnight+"', seqnum='"+obnm[sel[j]]+"', reason = 'No ThAr Lamp.', /chi_q, /log"
	 	print, 'TO GET RID OF IT.'
	 	stop
	 endif
	 if found and strt(fxpar(header, 'COMPLAMP')) eq 'TH-AR' then begin 
		thidfile = fn
		print, 'FOUND thid file '+thidfile
		return
	 endif 
	 j++
   endwhile
endfor

end;findthid.pro
