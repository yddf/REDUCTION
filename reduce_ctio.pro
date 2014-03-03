pro reduce_ctio,  redpar, mode, flatset=flatset, thar=thar, $                   
   order_ind=order_ind,  star=star, date=date

; Batch job to drive raw reduction for CTIO
; Based on 12-Mar-2011 DF revised for the e2v-4k chip 
; Oct 15, 2011 AT
; revised 20120405 ~MJG

; redpar = readpar('file.par') contains imdir, prefix, and other directories
; mode:  one of: narrow, slicer, slit, fiber
; flatset = [indices] 
; thar = [thar1, iod1, iod2] calibration files
; order_ind =  order definition file, set -1 to use summed flat instead
; star = [indices] ; stellar spectra to reduce
;date: yymmdd (used to create directories)


;2. Prefix added to FITS headers:  
	prefix=redpar.prefix   ; e.g. 'chi111003'   
	if strpos(prefix,'.') lt 0 then prefix=prefix+'.' ; add the point

; First check if all relevant input parameters are specified

   if ~keyword_set (order_ind) then begin
       print, 'REDUCE_CTIO: ORDER definition is not given at input, use flat instead'
       order_ind = -1
    endif 

; Identify the mode
   modeidx = redpar.mode

;1. CCD-Image Input Directory Path, where raw images are stored.
    indir=redpar.rootdir+redpar.rawdir+redpar.imdir ;e.g. /mir7/raw/090102/
;6. OUTPUT Directory Path and Output Prefix (ie, tapename for run)
    outdir= redpar.rootdir + redpar.iodspecdir + date + '/'
    if ~file_test(outdir) then spawn, 'mkdir '+outdir

; Try to read from the disk previously saved flats
if ~keyword_set (flatset) then begin
     name = redpar.rootdir+redpar.flatdir+prefix+mode+'.flat'
     tmp = file_search(name, count=flatcount)
     if flatcount eq 0 then begin
       print, 'REDUCE_CTIO: FLATS are not given at input, not found on disk, returning.'
       return
     endif else begin
       print, 'REDUCE_CTIO: reading previously saved flat from disk' 
       rdsk, sum, name, 1 ; restore saved flat
       flatfnums='SUM'        
    endelse 
 endif else begin ; flats are given
    nrecf = n_elements(flatset)
    recnums = strtrim(string(flatset,format='(I4.4)'),2)  ;convert to strings
    flatfnums = prefix + recnums        
    flatfnames = indir + prefix + recnums + '.fits'  ;array of flat-field files 
 endelse 

;7.  Record number of Stellar spectra here:
    nrec = n_elements(star)
    recnums = strt(star, f='(I4.4)')  ; convert to string with leading zeros    
    spnums = prefix + recnums
    spfnames = indir + prefix + recnums +'.fits'
   ; string array of spectrum file names
    outprefix = redpar.prefix_tag +  prefix
    outfnames= outdir + outprefix  + recnums 
    
; Order-Finding Exposure: strong exposure, such as iodine or bright star(B star)
      if order_ind ge 0 then begin
	recint = order_ind
	recnums = strtrim(string(recint,format='(I4.4)'),2) ;convert to strings with leading zeros   
	ordfname = indir + prefix + recnums + '.fits'   
     endif else ordframe='FLAT'

;THORIUMS:  Insert record numbers to reduce  here:
;3. Record numbers for thar and iodine (don-t need sky subtraction)
      if keyword_set(thar) then threc = thar else threc = -1 
      if threc[0] ge 0 then begin
	thnrec = n_elements(threc)
	threcnums = strtrim(string(threc,format='(I4.4)'),2) ;convert to strings 
	thspfnames = indir + prefix + threcnums + '.fits'
	thoutfnames = outdir + outprefix  + threcnums  
      endif else threcnums = 'none'
  
	print,''
	print,'    ****ECHOING PARAMETER VALUES FROM REDUCE_CTIO****'
	print,'If values are incorrect stop the program'
	print,' '
	print,'SPECTRA:'
	print,spnums
	print,' '
	print,'FLATS:'
	print,flatfnums
	print,' '
	print,'DEFAULT ORDER FILE:'
	print, order_ind
	print,' '
	print, 'THORIUM/IODINE: '
	print, prefix + threcnums
if redpar.debug ge 2 then print, 'REDUCE_CTIO: press .C to continue' 
if redpar.debug ge 2 then stop
; CRUNCH  FLATS
       name = redpar.rootdir+redpar.flatdir+prefix+mode+'.sum'          
       if keyword_set(flatset) then begin
;          if redpar.debug then stop, 'REDUCE_CTIO: debug stop before flats, .c to continue'
          ADDFLAT, flatfnames,sum, redpar, im_arr  ; crunch the flats (if redpar.flatnorm=0 then sum = wtd mean)
          if (size(sum))[0] lt 2 then stop ; no data!

          wdsk, sum, name, /new
          print, 'REDUCE_CTIO: summed flat is written to '+name  
;          if redpar.debug then stop, 'Debug stop after flats, .c to continue'
       endif else begin
         print, 'Using previously saved flat '+name 
         rdsk, sum, name, 1  ; get existing flat from disk
         bin = redpar.binnings[modeidx] ; set correct binning for order definition
         redpar.binning = [fix(strmid(bin,0,1)), fix(strmid(bin,2,1))]
         print, 'The binning is ', redpar.binning
       endelse

;FIND DEFAULT ORDER LOCATIONS.  
;SLICERFLAT=1 means use narrow slit to define slicer order locations
if redpar.slicerflat eq 0 or mode ne 'slicer' then begin
	if order_ind ge 0 then begin
	  ctio_dord, ordfname, redpar, orc, ome 
	endif else ctio_dord, ordfname, redpar, orc, ome, image=sum
	
	name = redpar.rootdir+redpar.orderdir+prefix+mode+'.orc'
	wdsk, orc, name, /new
	print, 'REDUCE_CTIO: order location is written to '+name  
	;         if redpar.debug then stop, 'Debug stop after order location, .c to continue'
endif else begin
	name = redpar.rootdir+redpar.orderdir+prefix+'narrow.orc'
	rdsk, orc, name, 1
	orc[0,*] += redpar.dpks[modeidx]
	;now subtract 2 outer orders since the slicer is larger than the slits:
	redpar.nords -= 2
	orc = orc[*,1:(redpar.nords - 1)]
	;stop
endelse

; GET FLAT FIELD
        xwid = redpar.xwids[modeidx]
;       if redpar.debug then stop, 'REDUCE_CTIO: debug stop before getting flat' 
        flat = getflat(sum, orc, xwid, redpar, im_arr=im_arr)
        name = redpar.rootdir+redpar.flatdir+prefix+mode+'.flat'
        fitsname = redpar.rootdir+redpar.flatdir+prefix+mode+'flat.fits'

        wdsk, flat, name, /new
        rdsk2fits, filename=fitsname, data = flat
        print, 'REDUCE_CTIO: extracted flat field is written to '+name  
        FF = flat[*,*,0] ; the actual flat
        if redpar.debug ge 2 then stop, 'Debug stop after flat field, .c to continue'

;REDUCE ThAr,Iodines (RED)
    if keyword_set(thar) then begin
		numthar=n_elements(threc)
	 	FOR i=0,numthar-1 do begin
	   		redpar.seqnum = strt(threcnums[i])
			CTIO_SPEC,prefix,thspfnames[i], thoutfnames[i],redpar, orc,xwid=xwid, flat=ff,/nosky
	 	ENDFOR
	 	CATCH, /CANCEL ; clear the catch block in case we bailed on the last one
	endif

PRINT, 'RIGHT BEFORE STELLAR CALL TO CTIO_SPEC'
if redpar.debug gt 1 then STOP
;STELLAR SPECTRA REDUCTION (RED)
	FOR i=0,nrec-1 do begin
	   redpar.seqnum = recnums[i]
     	CTIO_SPEC,prefix,spfnames[i],outfnames[i],redpar, orc, xwid=xwid, flat=ff, /cosmics
     ENDFOR
end
