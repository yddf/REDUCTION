pro reduce_ctio4k,  $
   imdir=imdir,$
   imprefix=imprefix,$
   flatset=flatset,$
   thariod=thariod,$
   bstar_ind=bstar_ind, $
   star=star,$
   qbcvel=qbcvel, $                 ; location of qbcvel.ascii file
   xwid=xwid, $
   dpks=dpks

; Batch job to drive raw reduction for CTIO
; DAF 17-May-08
; 12-Mar-2011 DF revised for the e2v-4k chip 

; INPUT from sorting_hat
; imdir: 090102
; imprefix: qa01 
; flatset = [indices] 
; thar = [thar1, iod1, iod2]
; bstar_ind = 2053
; star = [indices]


; FOURTH: RUN REDUCTION CODE
; 8. create the logsheet with logmaker.pro
; 9. run barylog.pro to get the barycentric corrections
; 10. run the raw reduction code, reduce_ctio4k.pro
; 11. run the wavelength calibration, thid.pro
; 12. if flag = 1 (or "yes") then generate fits files, setup_fits.pro 

redpar = readpar('ctio.par') ; reduction parameters, directories

;1. CCD-Image Input Directory Path, where raw images are stored.
	indir=redpar.rootdir+redpar.rawdir+strcompress(string(imdir),/rem)+'/' ;e.g. /mir7/raw/090102/
;	indir='/mir7/raw/'+strcompress(string(imdir),/rem)+'/' ;e.g. /mir7/raw/090102/

;2. Prefix added to FITS headers:  qa01.
	prefix=imprefix   ; e.g. 'qa04'   
	if strpos(prefix,'.') lt 0 then prefix=prefix+'.'


;3. Record numbers for thar and iodine (don-t need sky subtraction)
    if keyword_set(thariod) then begin 
      threc=thariod
      xgd=where(threc ne 0)
      threc=threc[xgd]
      if n_elements(threc) eq 0 then threc = [-1]
    endif

;4. Flats. MUST be given on input!
    wrecint = flatset
    nrecf = n_elements(wrecint)
    recnums = strtrim(string(wrecint,format='(I4.4)'),2)  ;convert to strings
    flatfnums = prefix + recnums        
    flatfnames = indir + prefix + recnums + '.fits'  ;array of WF files 
;    if (numsets gt 0) then widerec = wrecint  
    
;5. Record number of Well-Exposed spectrum, -1 to use flat field

	ord_finder_rec = bstar_ind

;6. OUTPUT Directory Path and Output Prefix (ie, tapename for run)
    outdir= redpar.rootdir + redpar.iodspecdir

;7.  Record number of Stellar spectra here:
	recint = star

	;names of raw fits files 
    nrec = n_elements(recint)
    recnums = strtrim(string(recint, FORMAT='(I4.4)'),2)  ; convert to string with leading zeros    
    print,recnums
    spnums = prefix + recnums
    spfnames = indir + prefix + recnums +'.fits'

	; string array of spectrum file names
	outprefix = 'r' +  prefix
	outfnames= outdir + outprefix  + recnums 
    
; Make sure to include The IDL Directory Path which contains all reduction programs
;	pathSep = PATH_SEP(/SEARCH_PATH)
;	oldPath = !path
;    !path = '/mir7/pro/' + pathSep + !path

; Order-Finding Exposure: strong exposure, such as iodine or bright star(B star)
	recint = ord_finder_rec
	recnums = strtrim(string(recint,format='(I4.4)'),2) ;convert to strings with leading zeros
		  
	ordfnums = prefix + recnums   
	ordfname = indir + prefix + recnums + '.fits'   

;THORIUMS:  Insert record numbers to reduce  here:
      if keyword_set (thariod) then begin
	thnrec = n_elements(threc)
	threcnums = strtrim(string(threc,format='(I4.4)'),2) ;convert to strings 
	thspfnames = indir + prefix + threcnums + '.fits'
	thoutfnames = outdir + outprefix  + threcnums  
     endif else threcnums = 'none'
  

	print,''
	print,'    ****ECHOING PARAMETER VALUES FROM REDUCE_CTIO4k****'
	print,'If values are incorrect stop program and reset in REDUCE_CTIO4K.'
	print,' '
	print,'SPECTRA:'
	print,spnums
	print,' '
	print,'FLATS:'
	print,flatfnums
	print,' '
	print,'DEFAULT ORDER FILE:'
	print,ordfnums
	print,' '
	print, 'THORIUM/IODINE: '
	print, prefix + threcnums

        stop, 'Debug stop before flats, .c to continue'
;       if redpar.debug then stop, 'Debug stop before flats, .c to continue'

; CRUNCH WIDEFLATS
;	 ADDFLAT, flatfnames,prefix, redpar  ; dual amp r/o new chip 
;        if redpar.debug then stop, 'Debug stop afrer flats, .c to continue'

;FIND DEFAULT ORDER LOCATIONS. 
	if keyword_set(ordfname) then CTIO_DORD, prefix, ordfname, redpar, 0., iodfn=iodfn

; GET FLAT FIELD
        rdsk, flt, prefix+'.sum'                    ;get the previously determined flat
	rdsk,orc,prefix + '.ord'			;restore default order locs
       stop, 'REDUCE_CTIO4k: before getting flat' 
       flat = getflat(flt,orc,xwid)

;REDUCE ThAr,Iodines (RED)
        if keyword_set(thariod) then begin
	numthar=n_elements(threc)
	 FOR i=0,numthar-1 do begin
		CTIO_SPEC,prefix,thspfnames[i], thoutfnames[i],redpar, xwid=xwid, flat=flat, /thar,/nosky
	 ENDFOR
	 CATCH, /CANCEL ; clear the catch block in case we bailed on the last one
	endif

;STELLAR SPECTRA REDUCTION (RED)
	FOR i=0,nrec-1 do begin
     	CTIO_SPEC,prefix,spfnames(i),outfnames(i),redpar, xwid=xwid, flat=flat, /cosmics
     ENDFOR


end
