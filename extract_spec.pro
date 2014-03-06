;+
;
;  NAME: 
;     extract_spec
;
;  PURPOSE: Extract the spectrum and divide by flat if specified
;
;  CATEGORY:
;      YDDF
;
;  CALLING SEQUENCE:
;
;      extract_spec
;
;  INPUTS:
;
;  KEYWORD PARAMETERS:
;    
;  EXAMPLE:
;      extract_spec
;
;  MODIFICATION HISTORY:
;        c. Matt Giguere 2014.03.05 12:27:23
;
;-
pro extract_spec, $
spfname, $
outfname, $
redpar, $
orc, $
cosmics = cosmics, $
flat = flat, $
nosky = nosky

angstrom = '!6!sA!r!u!9 %!6!n'
loadct, 39, /silent
usersymbol, 'circle', /fill, size_of_sym = 0.5

DEBUG=redpar.debug
xwid = redpar.xwids[redpar.mode]

print,'extract_spec: Entering routine.'
print,'spfname=',spfname

;Read the image file
im = getimage(spfname, redpar, header=header)  
if (size(im))[0] lt 2 then begin
   print, 'Image is not found. Returning from extract_spec.'
   stop
endif

sz = size(im)		
ncol = sz[1]				;# columns in image
nrow = sz[2]				;# rows in image
szf = size(flat)		
ncolf = szf[1]				;# columns in image
nrowf = szf[2]				;# rows in image

if ncol ne ncolf  then begin
   print, 'extract_spec: HALT! Your image is not the same size as your flat!'
   stop
endif

;OLD SCHOOL WAY OF FLAT FIELDING:
if keyword_set(flat) and redpar.flatnorm eq 2 then spec = im/flat 

;EXTRACT SPECTRUM
if not keyword_set(thar) then begin
   getspec, im, orc, xwid, spec, sky=sky, $
	  cosmics=cosmics, optspec=optspec, $
	  diff=replace, gain=redpar.gain, ron=redpar.ron, $
	  redpar = redpar
endif else begin
   ; ThAR - no cosmic removal
   getspec,im,orc,xwid,spec, gain=redpar.gain, ron=redpar.ron, $
	  redpar = redpar
endelse

;save the original spec
spec_o = spec

if redpar.debug ge 1 then begin
   print, '***********************************************'
   print, 'RIGHT BEFORE FLAT-FIELDING...'
   print, '***********************************************'
endif

if redpar.debug ge 2 then stop

; flat-field correction
if keyword_set(flat) and redpar.flatnorm le 1 then spec = double(spec)/flat else $
	print, 'extract_spec: WARNING: no flat-field correction!'
i=0
specsz = size(spec)
nords = specsz[2]

if redpar.debug ge 1 and redpar.debug le 2 then begin
!p.multi=[0, 1, 3]
fdir = redpar.plotsdir + 'extracts/'
spawn, 'mkdir '+fdir
fdir = redpar.plotsdir + 'extracts/' + redpar.imdir
spawn, 'mkdir '+fdir
fdir = redpar.plotsdir + 'extracts/' + redpar.imdir + redpar.seqnum
spawn, 'mkdir '+fdir
endif;debug plots

for i=0, nords-1 do begin
   if redpar.debug ge 1 and redpar.debug le 2 then begin
	  fname = fdir+'/'+redpar.prefix+redpar.seqnum+'_Ord'+strt(i)
	  if file_test(fname) then spawn, 'mv '+fname+' '+nextnameeps(fname+'_old')
	  ps_open, fname, /encaps, /color
   endif;debug plots

   if redpar.debug ge 1 then begin
	  plot, spec_o[*,i], title=redpar.prefix+redpar.seqnum+' Order '+strt(i)+' Extracted', /xsty, /ysty, ytitle='Flux'
	  plot, flat[*,i], title=redpar.date+' '+redpar.modes[redpar.mode]+' Mode Order '+strt(i)+' Flat', /xsty, /ysty, ytitle='Flux'
	  plot, spec[*,i], title=redpar.prefix+redpar.seqnum+' Order '+strt(i)+' Spec/Flat', /xsty, /ysty, $
	  xtitle='Dispersion Direction [pix]', ytitle='Flux'
   endif

   if redpar.debug ge 1 and redpar.debug le 2 then begin
   ps_close
   spawn, 'convert -density 200 '+fname+'.eps '+fname+'.png'
   endif
endfor


print,'extract_spec: Saving extracted spectrum to ' + outfname
spec=rotate(spec,7)
print, header
mwrfits, spec, outfname, header
end;extract_spec.pro