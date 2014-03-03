pro getspec,im,orc,xwd,spec,sky = sky, cosmics = cosmics, diff = replace, $
optspec=optspec, gain=gain, ron=ron, redpar = redpar
;Subroutine extracts spectra from all orders described in orc. 
; im (input array (# columns , # rows)) image from which orc and were
;   determined and from which spectrum is to be extracted.
; orc (input array (# coeffs , # orders)) polynomial coefficients (from FORDS)
;   that describe the location of complete orders on the image.
; xwd (input scalar) fractional extraction width (from GETXWD)
; spec (output array (# columns , # orders)) final extracted spectrum from im.
; optspec  returns optimally extracted spectrum, spec is boxcar-extracted 
;Calls GETARC
;24-Oct-89 JAV	Create.
;01-Nov-89 GBB	Modified to allow no background subtraction.
;10-Nov-89 JAV  Cleaned up background subtraction logic.
;03-Dec-89 JAV	Added fractional extraction width to argument list.
;14-Dec-89 JAV	Fixed checks for swath off edge of spectrum.
;19-Jan-89 JAV	Fixed coefficient calculation in 'Arc Off Edge of Image' tests.
;23-Jan-89 JAV	Really fixed 'Arc Off Edge if Image' tests.
;06-Jun-90 JAV	Added argument to GETARC call so total counts are returned.
;04-Sep-90 JAV	Fixed background subtraction bug; backgd/pixel is subtracted
;		 from spectrum counts/pixel BEFORE conversion to total counts.
;13-Sep-90 JAV	Added user specified extraction width logic ($hamxwd stuff).
;18-Apr-92 JAV	Updated global variable list/interpretations. Changed xwd
;		 logic.
;22-Sep-92 ECW  Added test to determine how to extend orders on high and
;		low ends of image depending on value of xwd.
;12-Jun-01 JTW  Added cosmic ray removal machinery
;17-May-08 DAF  Adapted for CTIO
;
;
;
;BUG: IDL will halt with an error if image is indexed with an invalid index.
;     ANA used to return the last element in the array. New logic is needed
;     here to handle orders that are partially off the chip.
; 26-Oct-2011 AT No common block. gain and RON are passed by the
; caller, not checked on entry. 


if n_params() lt 4 then begin
  print,'syntax: getspec,im,orc,xwd,spec[,sky[,cosmics]'
  retall
endif

;  trace,25,'GETSPEC: Entering routine.'

;Define useful quantities.
im=double(im)
ncol = n_elements(im(*,0))				;# columns in image
nrow = n_elements(im(0,*))				;# rows in image
ncoef = n_elements(orc(*,0))				;# polyn. coeffs
nord =  n_elements(orc(0,*))				;# full orders in orc
ix = findgen(ncol)					;column indicies
spec = dblarr(ncol,nord)				;init spectrum
orcend = dblarr(ncoef,nord+2)				;init extended orcs

;GETARC needs order location coefficients (orc) on both sides of arc swath to 
;  be extracted. In order to extract the lowest and highest orders specified
;  in orc, we need to extend orc one extra order on each end. We shall do so
;  by linearly extrapolating the last two orc on each end.
;Extend orc on the low end. Check that requested swath lies on image.
orclo = 2*orc(*,0) - orc(*,1)				;extrapolate orc

coeff = orc(*,0)					;central coefficients
y=poly(ix,coeff) - xwd				;edge of arc

yoff = where(y lt 0,noff)				;pixels off low edge
if noff gt 0 then begin				;check if on image
  ;GETARC will reference im(j) where j<0. These array elements do not exist.
  if redpar.debug ge 2 then print,'GETSPEC: Top order off image in columns [' $
  + strtrim(string(yoff(0)),2) + ',' $
  + strtrim(string(yoff(noff-1)),2) + '].'
endif

;Extend orc on the high end. Check that requested swath lies on image.
orchi = 2*orc(*,nord-1) - orc(*,nord-2)		;extrapolate orc


coeff = orc(*,nord-1)					;central coefficients
y=poly(ix,coeff) + xwd				;edge of arc

yoff = where(y gt nrow-1,noff)			;pixels off high edge
if noff gt 0 then begin			
  ;   GETARC will reference im(j) where j > ncol*nrow-1. These array elements do
  ;     not exist.
  print,'GETSPEC: Bottom order off image in columns [' $
  + strtrim(string(yoff(0)),2) + ',' $
  + strtrim(string(yoff(noff-1)),2) + '].'
endif;noff gt 0

;Define an order set (orcend) extended one extra order on either end.
for n = 1,nord do orcend(*,n) = orc(*,n-1)
orcend(*,0) = orclo
orcend(*,nord+1) = orchi

;Now loop through orders extracting spectrum and maybe subtracting background.

if redpar.debug ge 2 then print,'GETSPEC: Extracting spectrum.'


stop
if n_elements(sky) eq 0 then sky=im*0.
if keyword_set(cosmics) then remove_cosmics, im, orc, xwd, sky, spec = optspec, cosmics = replace, mask = mask, fwhm = seeing, gain=gain, ron=ron

imsz = size(im)
maskim = dblarr(imsz[1], imsz[2])
orcsz = size(orcend)
ybarr = dblarr(imsz[1], orcsz[2])
ytarr = dblarr(imsz[1], orcsz[2])

; boxcar extraction
!p.multi=[1,1,2]
for onum=1,nord do begin				;loop thru orders
  ;extract counts/pixel
  if keyword_set(redpar) then begin
  	;add variable order width:
     if redpar.slcrxtrawid[0] gt 0 and onum gt redpar.slcrxtrawid[1] and redpar.mode eq 1 then $
     	xwd = redpar.xwids[redpar.mode] + redpar.slcrxtrawid[0] else xwd = redpar.xwids[redpar.mode]
	 getarc, im, orcend, onum, xwd, arc, pix, debug = redpar.debug, ybi, yti
	 for i=0, imsz[1]-1 do maskim[i,ybi[i]:yti[i]] += (255d - 255d / nord * onum)
  endif else getarc, im, orcend, onum, xwd, arc, pix, ybi, yti
  ybarr[*, onum-1] = ybi
  ytarr[*, onum-1] = yti
  ;store total counts
  spec[*,onum-1] = double(arc) * pix			
endfor

if redpar.debug ge 1 and redpar.debug le 2 then begin
  fdir = redpar.plotsdir + 'arcs/'
  spawn, 'mkdir '+fdir
  fdir = redpar.plotsdir + 'arcs/' + redpar.date
  spawn, 'mkdir '+fdir
  fname = fdir+'/'+'arcs'+redpar.prefix+redpar.seqnum
  if file_test(fname+'.eps') then spawn, 'mv '+fname+'.eps '+nextnameeps(fname+'_old')+'.eps'
  ps_open, fname, /encaps, /color
  !p.multi=[0,1,1]
endif;debug plots

if redpar.debug ge 1 then begin
  ;only plot if debug is greater than 0:
  display, im
  for i=0, nord-1 do oplot, ybarr[*,i], col=250
  for i=0, nord-1 do oplot, ytarr[*,i], col=120
endif
;for i=0, nord-1 do begin
; display, im
; oplot, ybarr[*,i], col=250
; oplot, ytarr[*,i], col=120
; plot, spec[*,i], /xsty, /ysty
 ;stop
;endfor
;plot, spec
if redpar.debug ge 1 and redpar.debug le 2 then begin
  ps_close
  spawn, 'convert -density 200 '+fname+'.eps '+fname+'.png'
endif

;if keyword_set(redpar) then begin
;  if redpar.debug ge 2 then print, 'JUST FINISHED GET ARC. DISPLAYED IS THE MASK'
;  if redpar.debug ge 1 then begin
;	 if debug ge 1 and debug le 2 then begin
;		fdir = redpar.plotsdir + 'getspec/'
;		spawn, 'mkdir '+fdir
;		fdir = redpar.plotsdir + 'getspec/' + redpar.date
;		spawn, 'mkdir '+fdir
;		ps_open, nextnameeps(fdir+'/'+'mask'), /encaps, /color
;	 endif;debug plots
;	 display, maskim
;	 if redpar.debug ge 1 and redpar.debug le 2 then ps_close
;	 endif;debug>1
;	 
;	 if debug ge 1 and debug le 2 then begin
;		ps_open, nextnameeps(fdir+'/'+'im'), /encaps, /color
;	 endif;debug plots
;	 display, maskim
;	 if redpar.debug ge 1 and redpar.debug le 2 then begin
;		ps_close
;		spawn, 'convert -density 200 '+qfn+'.eps '+qfn+'.png'
;	 endif;ps_close
;  endif;debug>1
;  if redpar.debug ge 2 then stop
;endif;redpar passed in
;  if keyword_set(cosmics) then begin
;    if keyword_set(optspec) then spec=optspec ; use optimally extracted spectrum
;    spec(0) = seeing
;  endif


if redpar.debug ge 1 then print,'GETSPEC: Spectrum extracted - about to return to caller.'
if redpar.debug ge 2 then stop

return
end
