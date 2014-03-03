;+
;	NAME: SORTING_HAT
;
;	PURPOSE: To sort files according to binning and slit pair with ThAr
;				to run reduction code for extraction
;
; Modes for the sorting hat: 
; 	narrow: narrow slit (3x1) fast r/o (templates and I2 
;				alpha Cen mode) R=137,000
;	slicer: slicer (3x1) fast r/o R=80,000
;	slit:  slit (3x1) fast r/o  R=96,000
;	fiber: fiber (4x4) slow r/o R=27,000
;
; KEYWORDS:
;
;	OPTIONAL KEYWORDS: 
;  
;	REDUCE: runs reduce_ctio (reduce and get thar soln before running iod2fits)
;          if running reduce, need to pass in array of flats
;
;  IOD2FITS: matches thar solutions to correct observations and writes 
;				in fits format skip is an array of obnm that don-t need 
;				fits files skip=['2268'] thar_soln is the wavelength array
;				soln (if a matching thar not taken this night)
;
;	DOPPLER: a keyword indicating to run the Doppler code on the observations
;				taken with the input slicer position, run number, and date. 
;
;	DOPTAG:	
;
;  END_CHECK: checks to see that all 3x1 binned observations with input slit have been
;             reduced (iodspec), fits (fitspec), thar_soln
;
;	EXAMPLES: 
;		sorting_hat, '110311', run='rqa31', /narrow, /reduce
;		sorting_hat, '110311', run='rqa31', /narrow, /doppler, doptag='dd'
;       sorting_hat, '120313', run='achi120313', mode='narrow',/reduce, obsnm=[1008,1090, 1091] 
;  
; MODIFICATION HISTORY:
;	20110313 - Modified from old code to work with new dual amp readout ~DAF
; 	20110321 - Narrow slit now ignores objects with the name 'junk' ~MJG
;	20110331 - Now uses quartz for order finding if neither alpha cen nor a bstar
;					are present. ~MJG 
;  20110414 - Fixed a bug when processing files at the start of a new run ~MJG
;  20110808 - structured  (AT)
;  20120419 - Modified to work with Torrent 4 amp readout. Added 1x1 binning as an option. ~MJG
;  20121023 - Added the ThAr filename to redpar to be written to the FITS header ~MJG
;-
;
pro sorting_hat, night, run=run, iod2fits=iod2fits, reduce=reduce, $
doppler=doppler, doptag=doptag, end_check=end_check, skip=skip, $
thar_soln=thar_soln, getthid=getthid, mode = mode, obsnm=obsnm, $
		bin11 = bin11, flatsonly=flatsonly, tharonly=tharonly


angstrom = '!6!sA!r!u!9 %!6!n'
ctparfn = -1
spawn, 'pwd', pwddir
case pwddir of
   '/home/matt/projects/CHIRON/QC': ctparfn = '/home/matt/projects/CHIRON/REDUCTION/ctio.par'
   '/home/matt/projects/CHIRON/REDUCTION': ctparfn = '/home/matt/projects/CHIRON/REDUCTION/ctio.par'
   '/tous/CHIRON/REDUCTION': ctparfn = '/tous/CHIRON/REDUCTION/ctio.par'
   '/tous/CHIRON/QC': ctparfn = '/tous/CHIRON/REDUCTION/ctio.par'
endcase
if ctparfn eq -1 then begin
  print, '******************************************************'
  print, 'You must be running things from a different directory.'
  print, 'Your current working directory is: '
  print, pwddir
  print, 'ctparfn has not set. '
  print, 'Either changed your working directory, or modify the case'
  print, 'statement above this line.'
  print, '******************************************************'
  stop
endif

redpar = readpar(ctparfn)
redpar.imdir = night+'/'  ; pass night into redpar
redpar.date = night
redpar.versiond=systime()

if ~keyword_set(run) then begin
  imdir = redpar.rootdir+redpar.rawdir+redpar.imdir
  l = strlen(imdir)
  tmp = file_search(imdir+'*.fits', count = count) ; look for the data files
  if count eq 0 then begin 
      print, 'The run name could not be determined!'
      stop
  endif
; check the QA prefix: no less than 5 files per night!
  sel = where(strmid(tmp,l,2) eq 'qa')
  if n_elements(sel) gt 5 then run = strmid(tmp[sel[0]],l,4) else run = 'chi'+night
endif ;run not specified 

if strpos(run,'.') lt 0 then run=run+'.' ; add the point
 redpar.prefix = run

print, 'SORTING_HAT: night '+night+' run: '+run

;   Modes keyword
if ~keyword_set(mode) then begin 
    print, 'MODE is not defined. Returning from sorting_hat'
    return
endif

modeidx = (where(mode eq redpar.modes))[0] ; which mode?
if keyword_set(bin11) then modeidx += 4
redpar.mode = modeidx  ; pass current mode to other programs
if modeidx lt 0 then begin
    print, 'Error: unrecognized mode. Returning from sorting_hat'
    return
 endif

;logpath = redpar.rootdir+redpar.logdir+'20'+strmid(night, 0, 2)+'/'
logpath = redpar.logdir+'20'+strmid(night, 0, 2)+'/'
redpar.logdir=logpath
logsheet = redpar.rootdir+logpath+night+'.log'

iodspec_path = redpar.rootdir+redpar.iodspecdir+redpar.imdir
fits_path = redpar.rootdir+redpar.fitsdir+redpar.imdir
if ~file_test(fits_path) then spawn, 'mkdir '+fits_path
thid_path = redpar.rootdir+redpar.thiddir
thid_path2 = redpar.rootdir+redpar.thidfiledir
pretag = redpar.prefix_tag

readcol,logsheet, skip=9, obnm, objnm, i2, mdpt, exptm, bin, slit, f='(a5, a13, a4, a14, a8, a3, a6)'

;now to expand the quartz items in the logsheet:
   print, 'obnm before is: ', obnm
qcombs = where(strlen(obnm) gt 4)
for qi=0, n_elements(qcombs)-1 do begin
	qinit = strmid(obnm[qcombs[qi]], 0,4)
	qfini = strmid(obnm[qcombs[qi]], 5,4)
	ncombs = long(qfini) - long(qinit) + 1
	nobnm = lindgen(ncombs) + long(qinit)
	nobjnm = strarr(ncombs)+objnm[qcombs[qi]]
	ni2 = strarr(ncombs)+i2[qcombs[qi]]
	nmdpt = strarr(ncombs)+mdpt[qcombs[qi]]
	nexptm = strarr(ncombs)+exptm[qcombs[qi]]
	nbin = strarr(ncombs)+bin[qcombs[qi]]
	nslit = strarr(ncombs)+slit[qcombs[qi]]

	if qcombs[qi] eq 0 then begin
	  obnm = [strt(nobnm, f='(I04)'), obnm[(qcombs[qi]+1):*]]
	  objnm = [nobjnm, objnm[(qcombs[qi]+1):*]]
	  i2 = [ni2, i2[(qcombs[qi]+1):*]]
	  mdpt = [nmdpt, mdpt[(qcombs[qi]+1):*]]
	  exptm = [nexptm, exptm[(qcombs[qi]+1):*]]
	  bin = [nbin, bin[(qcombs[qi]+1):*]]
	  slit = [nslit, slit[(qcombs[qi]+1):*]]
	endif
	if qcombs[qi] ne n_elements(obnm)-1 and qcombs[qi] ne 0 then begin
	  obnm = [obnm[0:(qcombs[qi]-1)], strt(nobnm, f='(I04)'), obnm[(qcombs[qi]+1):*]]
	  objnm = [objnm[0:(qcombs[qi]-1)], nobjnm, objnm[(qcombs[qi]+1):*]]
	  i2 = [i2[0:(qcombs[qi]-1)], ni2, i2[(qcombs[qi]+1):*]]
	  mdpt = [mdpt[0:(qcombs[qi]-1)], nmdpt, mdpt[(qcombs[qi]+1):*]]
	  exptm = [exptm[0:(qcombs[qi]-1)], nexptm, exptm[(qcombs[qi]+1):*]]
	  bin = [bin[0:(qcombs[qi]-1)], nbin, bin[(qcombs[qi]+1):*]]
	  slit = [slit[0:(qcombs[qi]-1)], nslit, slit[(qcombs[qi]+1):*]]
	endif
	if qcombs[qi] eq n_elements(obnm)-1 then begin
	  obnm = [obnm[0:(qcombs[qi]-1)], strt(nobnm, f='(I04)')]
	  objnm = [objnm[0:(qcombs[qi]-1)], nobjnm]
	  i2 = [i2[0:(qcombs[qi]-1)], ni2]
	  mdpt = [mdpt[0:(qcombs[qi]-1)], nmdpt]
	  exptm = [exptm[0:(qcombs[qi]-1)], nexptm]
	  bin = [bin[0:(qcombs[qi]-1)], nbin]
	  slit = [slit[0:(qcombs[qi]-1)], nslit]
	endif
	
	qcombs[qi:*] += (ncombs - 1)
endfor
   print, 'obnm after is: ', obnm
ut = gettime(mdpt) ; floating-point hours, >24h in the morning
;stop

;**************************************************************
; ******* REDUCE *******************
;REDUCING THE DATA:	
;**************************************************************
   if keyword_set(reduce) then begin
		xsl=where(bin eq redpar.binnings[modeidx] and slit eq redpar.modes[modeidx],n_modes)
;                if xsl[0] lt 0 then stop, 'Sorting_hat: no files found! Stopping.' 
                if n_modes gt 0 then begin

		obnm1=obnm[xsl]
		objnm1=objnm[xsl]
		
      ;only create the median bias from if set in ctio.par:
      if redpar.biasmode eq 0 then begin
		  ;restore the log structure:
		  restore, redpar.rootdir+redpar.logstdir+'20'+strmid(redpar.date, 0, 2)+'/'+redpar.date+'log.dat'
		  ;now create the median bias frames if need be:
		  if redpar.modes[modeidx] ne 'fiber' then begin
			 fname = redpar.rootdir+redpar.biasdir+$
			  redpar.date+'_bin31_normal_medbias.dat'
			  print, 'Now testing median bias frame filename: ', fname 
			 ;if ~file_test(fname) then begin
			    print, '3x1 normal...'
				chi_medianbias, redpar = redpar, log = log, /bin31, /normal
			    print, '4x4 normal...'
				chi_medianbias, redpar = redpar, log = log, /bin44, /normal
			    print, '1x1 normal...'
				chi_medianbias, redpar = redpar, log = log, /bin11, /normal
			    print, '3x1 fast...'
				chi_medianbias, redpar = redpar, log = log, /bin31, /fast
			    print, '1x1 fast...'
				chi_medianbias, redpar = redpar, log = log, /bin11, /fast
			 ;endif
		  endif;nonfiber median bias frame check/make
		  fnamef = redpar.rootdir+redpar.biasdir+$
		  redpar.date+'_bin44_normal_medbias.dat'
		  if ~file_test(fnamef) then begin
			 chi_medianbias, redpar = redpar, log = log, /bin44, /normal
		  endif;fiber median bias frame check/make
      endif
      
		flatindx=where(objnm1 eq 'quartz',num_flat)

 		if num_flat gt 0 then begin ; process dash in the logfile flat numbers
		  tmp=[0]
		  for ii=0,num_flat-1 do begin  ; fabricate flat fields
			dum=obnm1[flatindx[ii]]
			if strlen(dum) gt 5 then begin
				gf=fix(strmid(dum,0,4))  &  gl=fix(strmid(dum,5,4))
				diff=gl-gf+1
				tmp=[tmp,gf+indgen(diff)]
			endif else tmp=[tmp,dum]
		  endfor
		  flatset=tmp[1:*]  ; throw out the dummy
               endif else begin
                 print, 'Sorting-hat: no flat files found. Returning.'
                 return
              endelse

		thariodindx=where(objnm1 eq 'thar' or objnm1 eq 'iodine',num_thariod)
		print, '******************************'
		print, 'THORIUM ARGON AND IODINE OBSERVATIONS TO BE PROCESSED: '
		print, obnm1[thariodindx]
                thar = fix(obnm1[thariodindx])

		starindx=where(objnm1 ne 'iodine' and objnm1 ne 'thar' $
			and objnm1 ne 'focus' and objnm1 ne 'junk' and objnm1 ne 'dark' $
                               and objnm1 ne 'bias', num_star)
      if keyword_set(flatsonly) then starindx = where(objnm1 eq 'quartz')                               
      if keyword_set(flatsonly) then thar = 0
		star = fix(obnm1[starindx]) ; file numbers
		if keyword_set(obsnm) then star = obsnm
		if keyword_set(tharonly) then star = 0


		if redpar.debug ge 2 then print, 'Sorting-HAT: before calling reduce_ctio'
		reduce_ctio, redpar, mode, flatset=flatset, star=star, thar=thar, date=night
		endif ;n_modes > 0
    endif  ;reduce

;**************************************************************
; ******* ThAr processing *******************	
;**************************************************************
 if keyword_set(getthid) then begin
		xsl=where(bin eq redpar.binnings[modeidx] and slit eq redpar.modes[modeidx],n_modes)

   if n_modes gt 0 then begin
	   obnm1=obnm[xsl]  &   objnm1=objnm[xsl]  
	   tharindx=where(objnm1 eq 'thar',num_thar)
	   print, '******************************'
	   print, 'THORIUM ARGON TO BE PROCESSED: '
	   print, obnm1[tharindx]
	   thar = obnm1[tharindx] ; string array

	   if keyword_set(thar_soln) then thidfile =  thid_path2+thar_soln+'.thid' else begin 
		  if strmid(run,0,2) eq 'qa' then begin 
			 findthid, night, redpar, thidfile,run=run 
		  endif else  findthid, night, redpar, thidfile
		  if thidfile eq 'none' or thidfile eq '' then begin 
			 print, 'No previous THID files found, returning. Type ".c"'
			 stop
		  endif 
	   endelse ; thar_soln  

	   print, 'Initial THID file: '+thidfile
	   restore, thidfile
	   initwvc = thid.wvc 

	   print, 'Ready to go into AUTOTHID'   
	   !p.multi=[0,1,1]
	   for i=0,num_thar-1 do begin 
		  isfn = iodspec_path+pretag+run+thar[i]
		  print, 'ThAr file to fit: ', isfn
		  rdsk, t, isfn,1
		  ;NEW, AUTOMATED WAY OF DOING THINGS:
		  print, 'ThAr obs is: ', thar[i], ' ', strt((1d2*i)/(num_thar-1),f='(F8.2)'),'% complete.'
		  ;stop

		  rawfn = redpar.rootdir+redpar.rawdir+redpar.imdir+'/'+run+thar[i]+'.fits'
		  header = headfits(rawfn)
		  if strt(fxpar(header, 'COMPLAMP')) ne 'TH-AR' then begin
			 print, 'WARNING! NO TH-AR LAMP IN FOR: '
			 print, rawfn
			 print, 'TYPE THE IDL COMMAND: '
			 print, "chi_junk, date='"+redpar.date+"', seqnum='"+thar[i]+"', reason = 'No ThAr Lamp.', /chi_q, /log"
			 print, 'TO GET RID OF IT.'
			 stop
		  endif else begin
			 auto_thid, t, initwvc, 6., 6., .8, thid, awin=10, maxres=4, /orev
			 ;for fiber, narrow and regular slit modes:
			 ;thid, t, 64., 64.*[8797d,8898d], wvc, thid, init=initwvc, /orev 
			 ;for slicer mode:
			 ;thid, t, 65., 65.*[8662.4d,8761.9d], wvc, thid, init=initwvc, /orev 
			 
			 if thid.nlin lt 700d then begin
			   print, 'CRAPPY FIT TO THE THAR! INTERVENTION NEEDED!'
			   print, 'ONLY '+strt(thid.nlin)+' GOOD LINES FOUND!'
			   stop
			 endif
			 
			 fnm = thid_path2+pretag+run+thar[i]
			 fsuf = '.thid'
			 if file_test(fnm+fsuf) then spawn, $
			 'mv '+fnm+'.thid '+nextname(fnm,fsuf)
			 save, thid, file=fnm+fsuf

			 mkwave, w, thid.wvc
			 w = reverse(w,2) ; increase with increasing order numver
			 fnm = thid_path+'ctio_'+pretag+run+thar[i]
			 fsuf = '.dat'
			 if file_test(fnm+fsuf) then spawn, $
			 'mv '+fnm+'.dat '+nextname(fnm,fsuf)
			 save, w, file=fnm+fsuf
		  endelse
	   endfor
	endif;n_modes > 0
endif ; getthid

;**************************************************************
;******* Write FITS files for reduced data ************
; from the input logsheet, find all observations matching the selected mode
;**************************************************************
if keyword_set(iod2fits) then begin
     x1=where(bin eq redpar.binnings[modeidx] and slit eq redpar.modes[modeidx] $
	     and objnm ne 'junk' and objnm ne 'dark' $
            and objnm ne 'focus' and objnm ne 'bias',n_found)
 
     tharindx=where((objnm eq 'thar') and (bin eq redpar.binnings[modeidx]) and (slit eq redpar.modes[modeidx]),  num_thar)
;     if x1[0] lt 0 or num_thar eq 0 then stop, 'Sorting_hat: no matching observations or ThAr for iod2fits. Stop'
if ( (n_found gt 0) and (num_thar gt 0)) then begin

      print,'Number of '+mode+' observations: ',n_found
      print, 'ThAr files: ', obnm[tharindx]
 
;  thar file is defined on input?
    if keyword_set(thar_soln) then begin
           restore, thid_path2+thar_soln+'.thid' ; thid structure
            mkwave, w, thid.wvc
            w = reverse(w,2) ; increase with increasing order number    
       endif else begin             
; get all wavelength solutions for this night and this mode,e UT of ThAr exposures 
;       wavfiles = thid_path+'ctio_'+run+obnm[tharindx]+'.dat' ; string array of wavelength solution file names  
         thidfiles = thid_path2+pretag+run+obnm[tharindx]+'.thid' ; string array of wavelength solution file names  
         ;stop
         wavut = ut[tharindx] ; time of ThAr exposures
; check existence of wavelength solutions, stop if not found
        for k=0,num_thar-1 do begin
            res = file_search(thidfiles[k], count=count)
            if count eq 0 then stop, 'Missing THID file '+ thidfiles[k]
        endfor 

         restore, thidfiles[0] ; w, first solution of the night
         mkwave, w, thid.wvc
         w = reverse(w,2) ; increase with increasing order numver
         ww = dblarr(num_thar,n_elements(w[*,0]),n_elements(w[0,*]))
         ww[0,*,*] = w
         for k=1,num_thar-1 do begin ; all other solutions in ww array
           restore, thidfiles[k]
           mkwave, w, thid.wvc
           w = reverse(w,2) ; increase with increasing order numver
           ww[k,*,*] = w
         endfor
    endelse ;thar_soln 

     for i=0,n_found-1 do begin	
			obnm[i]=strtrim(obnm[x1[i]])
			nxck=0
			if keyword_set(skip) then xck=where(obnmx1[[i]] eq skip,nxck) 
			if nxck eq 0 then begin
				rdsk,sp,iodspec_path+pretag+run+obnm[x1[i]],1   
				rdsk,hd,iodspec_path+pretag+run+obnm[x1[i]],2   
				sz=size(sp)  &   ncol=sz[1]    &    nord=sz[2]
				spec=fltarr(2,ncol,nord)
            if ~keyword_set(thar_soln) then begin ; find closest ThAr
                      ut0 = ut[x1[i]]
                      timediff = abs(ut0 - wavut)
                      sel = (where(timediff eq min(timediff)))[0]  
                      w = ww[sel,*,*]
                      ;save the ThAr filename to write
                      ;to the FITS header a few lines later:
                      thidfile_name = thidfiles[sel]
             endif
                ;the zeroth dimension of spec is the wavelength solution:
				spec[0,*,*]=w
				;the first dimension of spec is the actual spectrum:
				spec[1,*,*]=sp
				outfile=pretag+run+obnm[i]+'.fits'
				;*******************************************************
				;now to add reduction code info to fits headers:
				;*******************************************************
				;the number of tags in the redpar structure:
				nt = n_tags(redpar)
				tnms = string(tag_names(redpar), format='(A-8)')
				endhd = hd[-1]
				hd = hd[0:n_elements(hd)-2]
				for ii=0, nt-1 do begin
				  ;print, 'line: ', i
				  remlen = 78 - strlen(tnms[ii]+' = ')
				  vals = redpar.(ii)
				  val = strt(vals[0])
				  for j=1, n_elements(vals)-1 do begin
					 val += ', '+strt(vals[j])
					 print, 'j is: ', j, 'val is now: ', val
				  endfor
					 hd = [hd, tnms[ii]+'= '+"'"+string(val+"'", format='(A-'+strt(remlen)+')')]
				endfor
				hd = [hd, string('THARFNAM', format='(A-8)')+'= '+"'"+string(thidfile_name+"'", format='(A-'+strt(remlen)+')')]
				hd = [hd,endhd]
				
				;now change the NAXIS and NAXISn values to reflect the reduced data:
				specsz = size(spec)
				fxaddpar, hd, 'NAXIS', specsz[0], 'Number of data axes'
				fxaddpar, hd, 'NAXIS1', specsz[1], 'Axis 1 length: 0=wavelength, 1=spectrum'
				fxaddpar, hd, 'NAXIS2', specsz[2], 'Axis 2 length: extracted pixels along each echelle order'
				fxaddpar, hd, 'NAXIS3', specsz[3], 'Axis 3 length: number of echelle orders extracted'
				fxaddpar, hd, 'RESOLUTN', thid.resol, 'Resolution determined from the ThAr.'
				fxaddpar, hd, 'THIDNLIN', thid.nlin, 'Number of ThAr lines used for wav soln.'
				
				print, 'now writing: ', outfile, ' ', strt(i/(n_found - 1d)*1d2),'% complete.'
				writefits,fits_path+outfile, spec,hd
				
			   if redpar.debug ge 1 and redpar.debug le 2 then begin
				 fdir = redpar.plotsdir + 'fits/'
				 spawn, 'mkdir '+fdir
				 fdir = redpar.plotsdir + 'fits/' + redpar.date
				 spawn, 'mkdir '+fdir
				 fdir = fdir + '/halpha'
				 spawn, 'mkdir '+fdir
				 fname = fdir+'/'+'halpha'+redpar.prefix+obnm[i]
				 if file_test(fname+'.eps') then spawn, 'mv '+fname+'.eps '+nextnameeps(fname+'_old')+'.eps'
				 ps_open, fname, /encaps, /color
				 !p.multi=[0,1,1]
			   endif;debug plot fname and dirs
			   
			   if redpar.debug ge 1 then begin
				 ;only plot if debug is greater than 0:
				 plot, spec[0,*,39], spec[1,*,39], /xsty, /ysty, $
				 xtitle='Wavelength['+angstrom+']', ytitle='Flux', $
				 title=outfile, yran=[0,1.1*max(spec[1,*,39])]
			   endif;debug plotting
			   if redpar.debug ge 1 and redpar.debug le 2 then begin
				 ps_close
				 spawn, 'convert -density 200 '+fname+'.eps '+fname+'.png'
			   endif ;ps_close & png
			endif
     endfor  ;  iod2fits
     endif; num_thar and n_found > 0
  endif ;iod2fits

;***** DOPPLER ********	
if keyword_set(doppler) then begin
		x1=where(bin eq redpar.binnings[modeidx] and $
		slit eq redpar.modes[modeidx] and $
		objnm ne 'quartz' and objnm ne 'junk',n_found)

		if n_found gt 0 then print,'Number of '+mode[modeidx]+' observations: ',n_found
; Oct 15 2011: plug the Doppler:
                print, 'Doppler is not yet operational in sorting_hat. Returning.'
                return
	endif ; doppler
	
;**************************************************************
; ******* END-CHECK *********************
;THE END CHECK TO MAKE SURE EVERYTHING HAS BEEN PROCESSED:
;**************************************************************
 if keyword_set(end_check) then  begin
	x1=where(bin eq redpar.binnings[modeidx] and slit eq redpar.modes[modeidx] and objnm ne 'quartz',n_check)
        if x1[0] lt 0 then begin
          print, 'Sorting_hat: no files found! returning'
          return
        endif

        	for k=0,n_check-1 do begin
			if objnm[x1[k]] eq 'thar' then begin
				fthar=file_search(thid_path+'*'+obnm[x1[k]]+'*',count=thar_count)
				if thar_count eq 0 then print, objnm[x1[k]]+' '+obnm[x1[k]]+' has no ThAr soln '
			endif else begin
				fiod=file_search(iodspec_path+'*'+obnm[x1[k]]+'*',count=iod_count)
				if iod_count eq 0 then print, objnm[x1[k]]+' '+obnm[x1[k]]+'  has no iodspec'
				
				ffits=file_search(fits_path+'*'+obnm[x1[k]]+'*',count=fits_count)
				if fits_count eq 0 then print, objnm[x1[k]]+' '+obnm[x1[k]]+'  has no fitspec'
			endelse 
		endfor       
      endif ; end_check



end ;sorting_hat.pro
