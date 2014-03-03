;+
;
;  NAME: 
;     chi_que_distrib
;
;  PURPOSE: 
;   To distribute data to all of the CHIRON users. 
;
;  CATEGORY:
;      CHIRON
;
;  CALLING SEQUENCE:
;
;      chi_que_distrib, date=date
;
;  INPUTS:
;	       DATE: The date the data was taken. 
;
;  OPTIONAL INPUTS:
;          NOCOPY: Use this keyword if you don't want to
;			copy files to each plan directory
;
;			NOCALIB: Use this keyword if you don't
;			want to redo the calibrations.
;
;			NOLOG:
;
;  OUTPUTS:
;
;  OPTIONAL OUTPUTS:
;
;  KEYWORD PARAMETERS:
;    
;  EXAMPLE:
;      chi_que_distrib, date='110808'
;
;  MODIFICATION HISTORY:
;        c. Matt Giguere 2011.08.03 10:31:39 AM
;
;-
pro chi_que_distrib, date=date, $
nocopy=nocopy, nocalib = nocalib, $
nolog=nolog, noballs=noballs

spawn, 'unalias cp'

if ~keyword_set(date) then date='110802'
qodir = '/tous/mir7/queueobs/'+date+'/'

logpath='/tous/mir7/logsheets/'
logsheet=logpath+'20'+strmid(date, 0, 2)+'/'+date+'.log'
fits_path='/tous/mir7/fitspec/'

;FIRST TO PROCESS THE CALIBRATIONS FOR THE NIGHT:
spawn, 'mkdir '+qodir
caldir = qodir+date+'_cals/'
spawn, 'mkdir '+caldir


;read in the image prefix used for the night:
readcol,logsheet, obs, fn, ln, tel, ct, 1.5, pre, prefix, f='A,A,A,A,A,A,A,A', skip=3, num=1

;or just hardcode it now:
prefix = 'chi'+date

;read in the logsheet from the night:
readcol,logsheet, obnm, objnm, i2, mdpt, exptm, bin, slit, propid, $
	skip=9, delim=' ', f = 'A,A,A,A,A,A,A,A'

;****************************************************************************
;				COPY CALIBRATION DATA TO COMMON DIRECTORY:
;****************************************************************************
if ~keyword_set(nocalib) then begin
;first to copy all of the quartz files into the cal dir:
print, '****************************************************'
print, 'COPYING QUARTZ EXPOSURES...'
print, '****************************************************'
qtznm = where(objnm eq 'quartz', ntqz)
for ii=0, ntqz-1 do begin
startqtz = strmid(obnm[qtznm[ii]], 0,4)
endqtz = strmid(obnm[qtznm[ii]], 5,9)
print, 'startqtz: ', startqtz
print, 'endqtz: ', endqtz

for i=long(startqtz), long(endqtz) do $;
	spawn, 'cp -vf /raw/mir7/'+date+'/'+prefix[0]+'.'+strt(i, f='(I04)')+'.fits '+caldir
endfor

print, '****************************************************'
print, 'COPYING THAR EXPOSURES...'
print, '****************************************************'
;now to copy the Calibration ThAr exposures over to the same directory:
tharnums = where((objnm eq 'thar') AND (propid eq 'Calib'), numthars)
for ii=0, numthars-1 do begin
	spawn, 'cp -vf /raw/mir7/'+date+'/'+prefix[0]+'.'+strt(obnm[tharnums[ii]])+'.fits '+caldir
endfor

print, '****************************************************'
print, 'COPYING BIAS EXPOSURES...'
print, '****************************************************'
;now to copy the Calibration bias exposures over to the same directory:
biasnums = where((objnm eq 'bias') AND (propid eq 'Calib'), numbias)
for ii=0, numbias-1 do begin
	spawn, 'cp -vf /raw/mir7/'+date+'/'+prefix[0]+'.'+strt(obnm[biasnums[ii]])+'.fits '+caldir
endfor

print, '****************************************************'
print, 'COPYING IODINE EXPOSURES...'
print, '****************************************************'
;now to copy the Calibration iodine exposures over to the same directory:
iodnums = where((objnm eq 'iodine') AND (propid eq 'Calib'), numiod)
for ii=0, numiod-1 do begin
	spawn, 'cp -vf /raw/mir7/'+date+'/'+prefix[0]+'.'+strt(obnm[iodnums[ii]])+'.fits '+caldir
    spawn, 'cp -vf /tous/mir7/fitspec/'+date+'/a'+prefix[0]+'.'+strt(obnm[iodnums[ii]])+'.fits '+caldir
endfor

print, '****************************************************'
print, 'COPYING MASTER FLAT FILES...'
print, '****************************************************'
;now to copy the Calibration summed flat fields over to the same directory:
spawn, 'cp -vf /tous/mir7/flats/chi'+date+'.fiberflat.fits '+caldir
spawn, 'cp -vf /tous/mir7/flats/chi'+date+'.narrowflat.fits '+caldir
spawn, 'cp -vf /tous/mir7/flats/chi'+date+'.slicerflat.fits '+caldir
spawn, 'cp -vf /tous/mir7/flats/chi'+date+'.slitflat.fits '+caldir
endif;KW:nocalib
;stop

;****************************************************************************
;				NOW TO SORT THROUGH THE PLANS/PROPOSALS
;****************************************************************************

;proparr = ['N-11B-0602	      ', $  
; 'CN-10B-0043	      ', $  
; 'Vand-01           ', $  
; 'Roch-2            ', $  
; 'SUNY-11	          ', $  
; 'CN-11B-63	      ', $  
; 'CN-11B-21', $  
; 'CN-11B-15', $  
; 'N-11B-301', $  
; 'N-11B-0302', $  
; 'Van-08/N-144', $  
; 'SUNY-10', $  
; 'VAN-11B-0003', $  
; 'CN-21	', $  
; 'N-11B-0424', $  
; 'N-11B-0301']  
;   
; planarr = [15, $ ;N-11B-0602	        
; 0, $;CN-10B-0043	        
; 3, $;Vand-01             
; 0, $;Roch-2              
; 0, $;SUNY-11	            
; 0, $;CN-11B-63	        
; 13, $;CN-11B-21  
; 8, $;CN-11B-15  
; 0, $;N-11B-301  
; 16, $;N-11B-0302  
; 10, $;Van-08/N-144  
; 0, $;SUNY-10  
; 4, $;VAN-11B-0003  
; 12, $;CN-21	  
; 17, $;N-11B-0424  
; 0];N-11B-0301
propid2 = strt(propid)
propsrt = sort(propid2)
proparr1 = uniq(propid2[propsrt])
proparr = propid2[propsrt[proparr1]]
print, '*****************************************************'
print, ' OBJECTS FROM THE FOLLOWING PLANS WERE OBSERVED: '
print, '*****************************************************'
;print, transpose(proparr)
c1idx = where(proparr eq 'CPS', c1ct)
if c1ct then proparr=rm_elements(proparr, c1idx)
c2idx = where(proparr eq 'Calib', c2ct)
if c2ct then proparr = rm_elements(proparr, c2idx)
print, transpose(proparr)
c3idx = where(proparr eq '', c3ct)
if c3ct then proparr = rm_elements(proparr, c3idx)
c4idx = where(proparr eq 'Calib11', c4ct)
if c4ct then proparr = rm_elements(proparr, c4idx)
c5idx = where(proparr eq '35', c5ct)
if c5ct then proparr = rm_elements(proparr, c5idx)
c6idx = where(proparr eq '33', c6ct)
if c6ct then proparr = rm_elements(proparr, c6idx)

line=''

if n_elements(proparr) gt 0 then begin
;now to cycle through the proposals transferring people's files
;to their directories:
for pi=0, n_elements(proparr)-1 do begin
	;pilocs = where(strt(propid) eq strt(proparr[pi]) and strt(objnm) ne 'iodine', plct)
    ;not sure why I blocked out the i2s in the line above. Changed on 20111202 ~MJG:
	pilocs = where(strt(propid) eq strt(proparr[pi]), plct)
	print, 'pi is: ', strt(pi), '  |  plct is: ', strt(plct)
	;if they exist, copy them to the appropriate directory:
	if plct gt 0 then begin
		plandir = qodir+date+'_planid_'+strt(proparr[pi])+'/'
		spawn, 'mkdir '+plandir
	    if ~keyword_set(nocopy) then begin
		print, '********************************************'
		print, 'COPYING FILES TO PLAN '+proparr[pi]+' DIRECTORY...'
		print, '********************************************'
		for flnbr=0, plct-1 do begin
			spawn, 'cp -vf /raw/mir7/'+date+'/'+prefix+'.'+obnm[pilocs[flnbr]]+'.fits '+plandir
			spawn, 'cp -vf /tous/mir7/fitspec/'+date+'/r'+prefix+'.'+obnm[pilocs[flnbr]]+'.fits '+plandir
			spawn, 'cp -vf /tous/mir7/fitspec/'+date+'/a'+prefix+'.'+obnm[pilocs[flnbr]]+'.fits '+plandir
		endfor
		endif;KW(nocopy)
		if ~keyword_set(nolog) then begin
			print, '********************************************'
			print, 'CREATING PLAN-SPECIFIC LOGSHEET ...'
			print, '********************************************'
			;LogSheetFileName
			lsfnm = '/tous/mir7/logsheets/20'+strmid(date, 0, 2)+'/'+date+'.log'
			print, lsfnm
			;PlanSpecificLogSheetFileName
			pslsfn = plandir+date+'_planid_'+proparr[pi]+'.log'
			print, pslsfn
			openr, 1, lsfnm
			openw, 2, pslsfn
			for ln=0, 8 do begin
			readf, 1, line
			printf, 2, line
			endfor ;0-->8 (The logsheet header)
			while ~EOF(1) do begin
			readf, 1, line
			vars = strsplit(line, ' ', /extract)
			if n_elements(vars) ge 7 then begin
			ppid = strt(vars[7])
			if ppid eq 'Calib' or ppid eq proparr[pi] then printf, 2, line
			endif
			endwhile ;~EOF(1)
			free_lun, 1
			free_lun, 2
		endif;KW:nolog
	endif ;plct > 0
endfor ;pi=0 -> #proparr

print, '********************************************'
print, 'NOW BALLING THINGS UP...'
print, '********************************************'
if ~keyword_set(noballs) then begin
for pi=0, n_elements(proparr)-1 do begin
  pilocs = where(strt(propid) eq strt(proparr[pi]) and strt(objnm) ne 'iodine', plct)
  print, 'pi is: ', strt(pi), '  |  plct is: ', strt(plct)
  ;if they exist, copy them to the rightful directory:
  if plct gt 0 then begin
	 planflz = date+'_planid_'+strt(proparr[pi])+'/'
	 plantar = date+'_planid_'+strt(proparr[pi])+'.tgz'
	 spawn, 'cd '+qodir+ ' ; tar -cvzf '+plantar+' '+planflz
  endif
endfor ;pi
calflz = date+'_cals/'
caltar = date+'_cals.tgz'
spawn, 'cd '+qodir+' ; tar -cvzf '+caltar+' '+calflz
endif;KW:noballs

print, '********************************************'
print, 'NOW INSERTING THE TARFILE NAMES INTO '
print, 'THE MYSQL DATABASE FOR QUE OBSERVERS...'
print, '********************************************'
PRINT, 'the proparr is: ', proparr
mysqldate = '20'+strmid(date, 0, 2)+'-'+strmid(date, 2, 2)+'-'+strmid(date, 4, 2)
;Mysql Get Script ID File Name (mgsidfn):
mgsidfn = '/tous/mir7/queueobs/mysqlscripts/'+date+'_GetScriptID.txt'
openw, 2, mgsidfn
printf, 2, "SELECT  scripts.script_id  FROM scripts WHERE scripts.script_date_start = '"+mysqldate+"';"
close, 2
spawn, 'sshexo mysql -u automark -paut0marker chiron < '+mgsidfn, scriptidarr
scriptid = scriptidarr[1]
calball = '/'+date+'/'+date+'_cals.tgz'
for pa=0, n_elements(proparr)-1 do begin
  planid = strt(long(proparr[pa]))
  ;Mysql Get Plan File Name (mgplanfn):
  mgplanfn = '/tous/mir7/queueobs/mysqlscripts/'+date+'_GetObjPlan'+planid+'.txt'
  openw, 2, mgplanfn
  printf, 2, "SELECT nos_objects.object_id FROM nos_objects LEFT JOIN nos_results ON nos_objects.script_id = nos_results.script_id AND nos_objects.object_id = nos_results.object_id LEFT JOIN object_tarfiles ON nos_objects.object_id = object_tarfiles.object_id AND object_tarfiles.script_id = IFNULL(nos_objects.script_id,0) LEFT JOIN objects on nos_objects.object_id = objects.object_id WHERE nos_objects.script_id = "+scriptid+" AND nos_results.nosresults_status = 'Observed' AND object_tarfiles.object_id IS NULL AND objects.plan_id = "+strt(long(proparr[pa]))
  close, 2
  spawn, 'sshexo mysql -u automark -paut0marker chiron < '+mgplanfn, objects
  if strt(objects[0]) ne '' then begin
	objectarr = objects[1:*]
	planball = '/'+date+'/'+date+'_planid_'+planid+'.tgz'
	for obji=0, n_elements(objectarr)-1 do begin
	  objectid = strt(objectarr[obji])
	  mgobjfn = '/tous/mir7/queueobs/mysqlscripts/'+date+'_Plan'+planid+'_Obj'+objectid+'.txt'
	  openw, 2, mgobjfn
	  printf, 2, "INSERT IGNORE INTO object_tarfiles (object_id, script_id, otarfile_tarfile, otarfile_calibrations ) VALUES( "+objectid+", "+scriptid+", '"+planball+"', '"+calball+"')"
	  close, 2
	  spawn, 'sshexo mysql -u automark -paut0marker chiron < '+mgobjfn
	endfor;mysql object loop
  endif;objects ne ''
endfor;mysql plan loop
endif;proparr gt 0

end;chi_que_distrib.pro