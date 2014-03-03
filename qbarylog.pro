pro qbarylog,logfile, test=test, baryDir=baryDir, prefix=prefix 

; FUNCTION: Calculate Barycentric Correction for
; 	    Stellar spectra.

;  METHOD: 1.) Obtain Object Name, FLUX-WIEGHTED mean epoch, etc. from
;              online logsheets.
;          2.) Star  positions are stored in: /mir7/bary/ctio_st.dat,  hip.dat,
;              kother.ascii (ktranslation.dat rarely used to find HIP numbers)
;          3.) Drive qbary.pro, which makes calculation.
;          4.) Append qbcvel.ascii file


;Create: NOV-7-93	ECW
;Modified: JUL-95	CMc.   Modified to Drive New Barycentric Routines 
;Modified; JAN-96       ECW    Modified to get input info from online 
;			       logsheets.  Currently setup to run on 
;                              the machines: hodge,quark,coude.
;Modified; Jan-97       GWM    cleaned up, made identical to Keck version
;Modified; Nov-97       GWM    Auto-read log lines; Invoke hipparcos catalog
;Modified; Nov-01   CMc    To accomodate 1000+ image seismology runs
;                               & general cleanup
;Modified; May-02   CMc    To include remove acc. calc. + other improvements
;Modified; Feb-05   CMc:
; --    To execute bcvel.ascii backups into bcvel_backup/ directory
; --    Convert From kbarylog to barylog for CTIO.
; --    "HR" designated stars do not have BC's computed.
;
;Modified; Mar-09	JMB		
;	log files should be passed with complete paths, this script will determine
;	the directory the file lives in and the run prefix
;
;Modified; Jan-10	JMB
;	updated to examine multiple log-file date formats and notice wether the UT
;	date changed during the night.  If it did, later on we increment the day
;	for any obs which have times < 14 hours.
;
;Modified; Feb-10	JMB
;	updated strndgen to format integers between 0-9999 as 4 digit strings
;	with leading zeros.
;
;Modfied; May-10	JMB
;	updated log-file date format code to notice when date  didn't change but
;	log file date format _appeared_ to have two dates (ie 2010 May 11/11)
;
;Modified: 20111012 ~MJG
;	the formating of what is printed is cleaner and more robust. Prior to this 
;	modification the routine was failing to include the full filename now that
;	the image prefixes are longer.
;
;Modified: 20111111 ~MJG
;	added a few items to the ThAR and list of items to exclude.
;
; At least 1 log file is required
;
if ( n_elements(logfile) eq 0 ) then begin
	print, "qBaryLog:  No logfiles found to process.  Returning"
	stop
endif

; HANDLE PARAMETER DEFAULTS & FILE PATHS
if ( ~ keyword_set(baryDir) ) then baryDir = '/tous/mir7/bary/'  ; DF revision Mar2012
if ( keyword_set(test) ) then baryFile = 'qbcvel.test.ascii' else baryFile = 'qbcvel.ascii'

bcFile = baryDir + baryFile
structFile = 'ctio_st.dat'

; Verify that the neccessary files exist
if (n_elements(file_search(structFile)) eq 0) then message,structFile+ ' is missing.'
if (n_elements(file_search(bcFile)) eq 0) then message,bcFile+ ' is missing.'

;VARIABLE DECLARATIONS:
noerror=1	& chk=''		& du=''		& dum=''	& dummy=''	&  req=''  
log=''		& logline=''	& obtype='' & mid=''
year=0		& month=0		& day=0		& hr=0		& sc=0		& epoch=0

skiplist= ['WIDEFLAT','WIDE','WF','W','WIDE-FLAT','WIDE/FLAT','W/F', $
           'WIDEFLATS','WIDES','WFS','WIDE-FLATS','WIDE/FLATS', $
           'WIDE_FLAT','FLAT','JUNK', 'FOCUS','quartz','QUARTZ','Quartz',$
           'QTZ', 'qtz', 'Qtz','FLATFIELD','SATURATED','SKIP', 'BIAS', 'DARK']
iodlist = ['IODINE','I','I2','IOD','QTZIODINE', 'iodine']
thorlist = ['TH-AR','TH_AR','THAR','TH/AR','THORIUM','THORIUM-ARGON','hear', 'th-ar']
daylist = ['DAY','DS','DAYSKY','DAY-SKY','DAY_SKY']
moonlist = ['moon', 'Moon','MOON']
skylist = ['SKY', 'DARK','BIAS','narrowflat','narrflat','SUN','Sun','sun','SKYTEST']
months = { Jan: '01', Feb: '02', Mar: '03', Apr: '04', May: '05', Jun: '06', Jul: '07', Aug: '08', Sep: '09', Oct: '10', Nov: '11', Dec: '12' }
monthTags = tag_names(months)

;DEFINE STRUCTURES
maxsize = 5000                  ; 
log = {log, object: '', hour: '', min: '', sec: '', type: ''}
temp = {bcvel, filename:'', object:'', cz:0.d0, mjd:0.d0, ha:0.d0, type:''} ;Temporary structure for results
temp = replicate(temp[0],maxsize)

tempbcfile = strarr(maxsize)    ;temporary storage of ascii results: 200 lines

print,'     *************************************************************'
print,'     **      THE CTIO BARYCENTRIC CORRECTIONS PROGRAM           **'
print,'     **                                                         **'
print,'     ** You input:                                              **'
print,'     ** 		LOGSHEET: ,' + logFile +               '           **'
print,'     **                                                         **'
print,'     ** Output to: '  +  bcfile  +      '                       **'
print,'     **                                                         **'
print,'     *************************************************************'

logFileParts = stregex(logFile,'^(.*/)([0-9]+\.log)$',/EXTRACT,/SUBEXPR)
if ( n_elements(logFileParts) ne 3 or logFileParts[2] eq '' ) then message, "New logfile format or Unable to locate a valid log file: "+logFile

logDir = logFileParts[1]
logFile = logFileParts[2]

logfileorig = logfile

;READ LOGSHEET
print,'reading: ',logFileParts[0]

restore,structFile & cat = dum   ;RESTORE COORD. STRUCTURE

; spit out the file and grab the UT Date line
spawn, "cat '" + logdir + logfile + "'" , output
output = cull(output)           ; remove any blank lines
logline = (output(where(getwrds(output,0) eq 'UT')))[0] ; 1 element

if strmid(logfile,0,1) ne '1' then begin   ; DF revision Mar2012
if ( strlen(logline) ne 0 ) then stop,'no proper format for logline' 

; We have several possible log file date formats, below we have regexes for all of those
;	which describe log-files with multiple dates
;
   flip = 0                     ; check for a changing UT
   dateRegexes = [ $
		{	regex: 'UT Date: ([0-9]{4})[^a-zA-Z]+([a-zA-Z]{3})[^0-9]+([0-9]+)\/([0-9]{4})[^a-zA-Z]+([a-zA-Z]{3})[^0-9]+([0-9]+)', $ ; changing year
			flip: 1, namedMonth: 1, checkDay: 0 }, $
		{	regex: 'UT Date: ([0-9]{4})[^a-zA-Z]+([a-zA-Z]{3})[^0-9]+([0-9]+)\/([a-zA-Z]{3})[^0-9]+([0-9]+)', $ ; changing month
			flip: 1, namedMonth: 1, checkDay: 0 }, $
		{	regex: 'UT Date: ([0-9]{4})[^a-zA-Z]+([a-zA-Z]{3})[^0-9]+([0-9]+)\/([0-9]+)', $ ; changing day, same month
			flip: 1, namedMonth: 1, checkDay: 1 }, $
		{	regex: 'UT Date: ([0-9]{4})[^a-zA-Z]+([a-zA-Z]{3})[^0-9]+([0-9]+)[^\/]', $ ; single ut date
			flip: 0, namedMonth: 1, checkDay: 0 }, $
		{	regex: 'UT Date: ([0-9]{4})[ \t]+([0-9]{2})\/([0-9]{2})-([0-9]{2})\/([0-9]{2})', $ ; numeric instead of string dates
			flip: 1, namedMonth: 0, checkDay: 0 } $
	]
   i = 0                        ;
   dateMatch = 0                ;
   while ( i lt n_elements(dateRegexes) and dateMatch eq 0 ) do begin
	if ( stregex(logline,dateRegexes[i].regex,/BOOLEAN) ) then begin
		dateReg = dateRegexes[i].regex
		dateMatch = 1
		flip = dateRegexes[i].flip
		namedMonth = dateRegexes[i].namedMonth
		checkDay = dateRegexes[i].checkDay
	endif
	i = i + 1
     end

   if ( ~dateMatch ) then stop,'UT date was in an unknown format'

; Extract the YYYY MMM DD from the date
   logFileDate = stregex(logline,dateReg,/EXTRACT,/SUBEXPR)
   year = logFileDate[1]
   m = where(monthTags eq strupcase(logFileDate[2]))
   if ( namedMonth ) then month = months.(m) else month = m
   day = logFileDate[3]
   if ( checkDay && logFileDate[4] eq logFileDate[3] ) then flip = 0 ; somehow the format was two day even though there was only 1 day (ie May 11/11)
endif ;reading old format logfile for dates

if ( strmid(logfile,0,1) eq '1' ) then begin ; DF revision Mar2012
   flip = 1  ; add a day to get correct UT
   year = '20'+strmid(logfile,0,2)
   m = strmid(logfile,2,2) 
   month = m
   day = strmid(logfile,4,2) 
endif

ans='y'
;Save Backup copy of bcvel.ascii
if ( ~ keyword_set(test) ) then begin
    strdate =strtrim(strcompress(day),2)+'_'+strtrim(strcompress(month),2)+$
      '_'+strtrim(strcompress(year),2)
    command = 'cp "'+bcFile+'" "'+barydir+'bcvel_backup/'+baryFile+'_'+strdate+'"'
;    spawn,command
endif  else print,'NOT Backing up!  Since this version of bary is in test mode!'

; Open the log file and get the prefix from the header area,
; continue through the whole file getting all of the entries so that we can
; make sure this log sheet wasn't already run.
openr, LOGHEAD,logDir+logFile,/GET_LUN
;prefix = ''
;while ( eof(LOGHEAD) eq 0 and prefix eq '' ) do begin
;	readf,LOGHEAD,logLine
;	if ( strpos(logLine,'prefix:') gt -1 ) then begin
;		prefix = stregex(logLine,'prefix:.([a-zA-Z0-9]+)',/EXTRACT,/SUBEXPR)
;		prefix = prefix[1]
;		break
;	endif
;endwhile
;close,LOGHEAD
;if ( prefix eq '' ) then message, "Unable to find the run prefix (ie qa04)"
print,day


num = 0
skip = 0                        ;Reset to NOT skip
strindgen = strtrim(string(indgen(10000),FORMAT='(I4.4)'),2)

if ( keyword_set(test) ) then begin
	print,' '
	print,'  Filename      Object     I2 in?   Time   Barycentric Corr. (m/s)'
	
	print,'---------------------------------------------------------'
endif

openr,logune,logdir+logfile,/get_lun
;LOOP THROUGH EACH LINE ON THE LOGSHEET
WHILE eof(logune) eq 0 do begin ;check for end of file (eof = 1)
    readf,logune,logline        ;read one line in the logsheet

	;Read the first four entries on the line.
    recnum = strtrim(getwrd(logline[0],0),2) ;record number
		
    log.object = strtrim(strupcase(getwrd(logline,1)),2) ;object name
    if log.object eq 'LHS2627' then log.object = 'LHS2726'
    first2 = strupcase(strmid(log.object,0,2))
    celltest = strtrim(strupcase(getwrd(logline,2)),2) ;Was cell in?
    strtime = strtrim(getwrd(logline,3),2) ;time from logsheet
    linelen = (strlen(logline))[0] ; first element only
    temptest = strpos(strupcase(logline),'TEMP') ; test for word "Template"

	;Construct reduced filename
    filename = prefix + '.' + recnum

	;Guarantee that this is really an observation of something useful
    IF ((celltest eq 'Y' or celltest eq 'N') and $ ;Was cell specified?
        (select(skiplist,log.object) ne 1) and $ ;Not wide flat nor skiplist?
        select(strindgen,recnum)) and  (linelen gt 1)  THEN BEGIN 

        if (celltest eq 'Y') then log.type='o' else log.type='t'
        if select(iodlist,log.object) then begin
            log.type='i' & log.object='iodine'
        endif
        if select(thorlist,log.object) then begin
            log.type='u' & log.object='th-ar'
        endif
        if select(daylist,log.object) then begin
            log.type='s' & log.object='day_sky'
        endif
        if select(moonlist,log.object) then begin
            log.type='o' & log.object='moon'
        endif
        if select(skylist,log.object) then begin
            log.type = 'u'            
        endif

        if temptest[0] ne -1 and log.type ne 't' then begin ; Error Trap
            print,'****WARNING:  Possible Template Detected: '
            print,logline
            print,'But observation type is not "t".  Was I2 cell in?'
            help,log.type
        endif


		;ENTERING TIME RETRIEVAL SECTION
        colon1 = strpos(strtime,':') ;position of first colon
        colon2 = rstrpos(strtime,':') ;postiion of last colon

        strhour = strmid(strtime,0,colon1) ;string version of UT hour of obs.
        strminute = strmid(strtime,colon1+1,colon2-colon1-1) ;UT minute
        strsecond = strmid(strtime,colon2+1,10) ;UT second

        hour = float(strhour)
        minutes = float(strminute) + float(strsecond)/60.d0

;fischer added next 2 lines to try to adjust for two UT's 
        if flip eq 1 then begin 
           dd=day
           if hour lt 14. then dd=dd+1
        endif
        if flip eq 0 then dd=day
        jdUTC = jdate([year,month,dd,hour,minutes])
        mjd = jdUTC-2440000.d0  ; modified JD

		; Fix lengths of output strings
        len = (strlen(filename))[0]
        if len lt 9 then for jj=0,9-len-1 do filename = filename + ' '
        obj = log.object
        len = (strlen(obj))[0]
        if len lt 8 then for jj=0,8-len-1 do obj = ' '+obj 
        if strlen(strminute) eq 1 then strminute = '0'+strminute 
        if strlen(strsecond) eq 1 then strsecond = '0'+strsecond 
        strtime = strhour+':'+strminute+':'+strsecond
        len = (strlen(strtime))[0]
        if len lt 9 then for jj=0,9-len-1 do strtime = strtime + ' '
                                ;
        cz = 0.0d0 
        ha = 0.d0
        filename = filename

        IF select(['o','t'],log.type) then begin ;need barycentric correction

;       LOOKUP COORDINATES: lookup.pro takes starname (log.object) and finds
;         coords, equinox, proper motion, parallax 
;         coords=[ra,dec], pm=[ra_motion, dec_motion]

           print,filename, ' ', log.object
           if first2 ne 'MO' then begin
              if first2 ne 'HR' then begin ; SKIP B STARS, MOON (no B.C. for B*s)
                 lookup,log.object,coords,epoch,pm,parlax,radvel,hip=hip,$
                        barydir=barydir,cat=cat,tyc=tyc

                 if abs( coords(0)) eq 99.d0 then begin ;Logsheet star not found
                    coords = [0.d0,0.d0]                ;force ra and dec = 0. :no object found
                    pm     = [0.d0,0.d0]                ;dummy proper motion
                    epoch = 2000.d0                     ;dummy epoch
                 endif else begin 

                    qbary,jdUTC,coords,epoch,czi,obs='ctio',pm = pm,$
                          barydir=barydir, ha=ha
                    cz = rm_secacc(czi,pm,parlax,mjd)
                 endelse
              endif             ; else print,'Skipping Bstar'
           endif  ;skipping Moon
        ENDIF                   ; 

		if ( keyword_set(test) ) then begin
			;Print Status to Screen
			stcz = strtrim(string(cz),2)
			len = (strlen(stcz))[0]
			if len lt 7 then for jj=0,7-len-1 do stcz = ' '+stcz
			infoline = '|  '+filename+' |  '+obj+' |  '
			infoline = infoline + celltest+'  | '+strtime+' | '+stcz+' |'
			k = (strlen(infoline))[0]-1
			dashln = '-'
			for p = 1,k do dashln = dashln+'-'
			
			print,infoline &  print,dashln
		endif

		;Store results to Structure
        temp[num].filename = filename
        temp[num].object = log.object
        temp[num].cz = cz
        temp[num].mjd = mjd
        temp[num].ha = ha
        temp[num].type = log.type
        num=num+1
    ENDIF
END                             ;while

;STORE RESULTS IN BCVEL.ASCII ?
temp = temp[0:num-1]            ;trim temp structure array
if ( keyword_set(test) ) then begin
	print,' '
	ans = ' '
	read,'Did all the printed results (above) look OK? (y or n)?',ans
endif else ans = 'Y'

if strupcase(ans) eq 'Y' then begin
    get_lun,une                 ;get Logical Unit Number
    openu,une,bcfile,/append    ;open bcvel file for writing
;    form = '(A9,3X,A10,1X,D11.3,1X,D12.6,1X,F7.3,1X,A1)'; modified for
;            fn, object, cts, mjd, ha, type
    form = '(2X, A-16,A12, D12.3,D13.6,F8.3,A2)' ; long names
    
    print,'Printing results to '+bcfile+' ...'
    for j=0,num-1 do begin
        fn = temp[j].filename
        ob = temp[j].object
        cz = temp[j].cz
        mjd = temp[j].mjd
        ha = temp[j].ha
        type = temp[j].type
        printf,une,format=form,fn,ob,cz,mjd,ha,type
        print,format=form,fn,ob,cz,mjd,ha,type
        ;stop
        
    end
    free_lun,une
    print,'Done with '+logfileorig
    print,' '
    comm=' '
    if num gt 55 then comm = 'You can have your pancakes now.'
    if num lt 30 then comm = 'Warm up the maple syrup.'
    print,' You observed ',strcompress(num),' stars. '+comm
    print,' '
    print,'It''s been a long night. Go to sleep!  Sweet Dreams.'
endif else begin
    print,' '
    print,'Make necessary changes (to logsheet?) and start again.'
    print,'The file ',bcfile,' was not affected. Exiting ...'
endelse
close,/all
;free_lun,logune 

end


