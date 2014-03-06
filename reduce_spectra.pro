;+
;
;  NAME: 
;     reduce_spectra
;
;  PURPOSE: To reduce the YDDF data
;   
;
;  CATEGORY:
;      YDDF
;
;  CALLING SEQUENCE:
;      reduce_spectra, redpar, mode, flatset=flatset, star=star, thar=thar, date=date
;
;  KEYWORD PARAMETERS:
;		date: the date the file is from in yymmdd format
;		flatset: (optional) specify the sequence numbers of files to use
;			for flat=fielding
;		flatname: (optional) specify the name of the already-summed flat 
;			image to use for flat fielding
;		order_ind: (optional) the order indeces
;    
;  EXAMPLE:
;      reduce_spectra, redpar, mode, flatset=flatset, star=star, thar=thar, date=date
;
;  MODIFICATION HISTORY:
;        c. Matt Giguere 2014.03.05 11:26:24
;
;-
pro reduce_spectra, redpar, $
date=date, $
flatset=flatset, $
flatname=flatname, $
order_ind=order_ind, $
star=star, $
thar=thar

;restore the image prefix from the redpar structure:
prefix=redpar.prefix   ; e.g. 'yddf111003'   
if strpos(prefix,'.') lt 0 then prefix=prefix+'.' ; add the point

;check if all relevant input parameters are specified

if ~keyword_set (order_ind) then begin
   print, 'reduce_spectra: ORDER definition is not given at input, use flat instead'
   order_ind = -1
endif 

; Identify the mode
modeidx = redpar.mode

;restore the mode name:
mode = redpar.modes[modeidx]

;CCD-Image Input Directory Path, where raw images are stored.
indir=redpar.rootdir+redpar.rawdir+redpar.imdir 

;OUTPUT Directory Path and Output Prefix 
outdir= redpar.rootdir + redpar.fitsdir + date + '/'
if ~file_test(outdir) then spawn, 'mkdir '+outdir

;THREE OPTIONS FOR FLAT FIELDING
;1. restore flat from flatname provided:
if keyword_set(flatname) then begin
   flatset = 1
   sum  = mrdfits(flatname)
   flatfnums = 'SUM'
endif;kw(flatname)

;2. Try to read from the disk previously saved flats
if ~keyword_set (flatset) then begin
   flatname = redpar.rootdir+redpar.flatdir+prefix+mode+'.flat'
   if ~file_test(flatname) then begin
	  print, 'reduce_spectra: FLATS are not given at input, not found on disk.'
	  stop
   endif else begin
	  print, 'reduce_spectra: reading previously saved flat from disk' 
	  sum = mrdfits(flatname)
	  flatfnums='SUM'        
   endelse 
endif 

;3. create flat from provided set of flat seqnums
if ~keyword_set(flatname) and keyword_set(flatset) then begin
   nrecf = n_elements(flatset)
   recnums = strtrim(string(flatset,format='(I4.4)'),2)  ;convert to strings
   flatfnums = prefix + recnums        
   flatfnames = indir + prefix + recnums + '.fits'  ;array of flat-field files 

   flatname = redpar.rootdir+redpar.flatdir+prefix+mode+'_sum.fits'          
   ADDFLAT, flatfnames,sum, redpar, im_arr  
   if (size(sum))[0] lt 2 then stop ; no data!
   mwrfits, sum, flatname
   stop
   print, 'reduce_spectra: summed flat is written to '+flatname  
endif 

;Record number of Stellar spectra here:
nrec = n_elements(star)
recnums = strt(star, f='(I4.4)')  ; convert to string with leading zeros    
spnums = prefix + recnums
spfnames = indir + prefix + recnums +'.fits'
; string array of spectrum file names
outprefix = redpar.prefix_tag +  prefix
outfnames= outdir + outprefix  + recnums +'.fits'

;Order-Finding Exposure: strong exposure, such as iodine or bright star(B star)
if order_ind ge 0 then begin
   recint = order_ind
   recnums = strtrim(string(recint,format='(I4.4)'),2) ;convert to strings with leading zeros   
   ordfname = indir + prefix + recnums + '.fits'   
endif else ordframe='FLAT'

;THORIUMS:  Insert record numbers to reduce  here:
;Record numbers for thar and iodine (do not need sky subtraction)
if keyword_set(thar) then threc = thar else threc = -1 
if threc[0] ge 0 then begin
   thnrec = n_elements(threc)
   threcnums = strtrim(string(threc,format='(I4.4)'),2) ;convert to strings 
   thspfnames = indir + prefix + threcnums + '.fits'
   thoutfnames = outdir + outprefix  + threcnums  
endif else threcnums = 'none'

print,''
print,'    ****ECHOING PARAMETER VALUES FROM reduce_spectra****'
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
if keyword_set(thar) then begin
print, prefix + threcnums
endif else print, 'NONE'
if redpar.debug ge 2 then print, 'reduce_spectra: press .C to continue' 
if redpar.debug ge 2 then stop

;FIND DEFAULT ORDER LOCATIONS. If order finding image is passed in
;use that. Otherwise, used the sum flat field:
if order_ind ge 0 then begin
   ;restore the input order finding image:
   ordfndim = getimage(ordfname, redpar)
   fords,ordfndim, orc, ome, redpar
endif else fords, sum, orc, ome, redpar

;write the newfound order locations to file:
name = redpar.rootdir+redpar.orderdir+prefix+mode+'.fits'
mwrfits, orc, name
print, 'reduce_spectra: order location is written to '+name  
if redpar.debug ge 2 then stop, 'Debug stop after order location, .c to continue'

; GET FLAT FIELD
xwid = redpar.xwids[modeidx]
if redpar.debug ge 2 then stop, 'reduce_spectra: debug stop before getting flat' 
flat = getflat(sum, orc, xwid, redpar, im_arr=im_arr)
fitsname = redpar.rootdir+redpar.flatdir+prefix+mode+'flat.fits'
mwrfits, flat, fitsname

print, 'reduce_spectra: extracted flat field is written to '+fitsname
FF = flat[*,*,0] ; the actual flat
if redpar.debug ge 2 then stop, 'Debug stop after flat field, .c to continue'

;REDUCE ThAr,Iodines (RED)
if keyword_set(thar) then begin
   numthar=n_elements(threc)
   FOR i=0,numthar-1 do begin
	  redpar.seqnum = strt(threcnums[i])
	  EXTRACT_SPEC,thspfnames[i], thoutfnames[i],redpar, orc, flat=ff,/nosky
   ENDFOR
endif

PRINT, 'RIGHT BEFORE STELLAR CALL TO EXTRACT_SPEC'
if redpar.debug gt 1 then STOP
;STELLAR SPECTRA REDUCTION (RED)
FOR i=0,nrec-1 do begin
   redpar.seqnum = recnums[i]
   EXTRACT_SPEC,spfnames[i],outfnames[i],redpar, orc, flat=ff, /cosmics
ENDFOR
end;reduce_spectra.pro