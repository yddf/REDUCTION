;+
; PURPOSE:Procedure to read a CCD image from the disk
; Subtract bias, correct for non-linearity, trim
;
; INPUT:  FILE: file name (in the data directory) 
;			 REDPAR: initial parameter file from (*.par)
;
; OUTPUT: image array, float, not rotated
;
; KEYWORDS: 
;         header: will contain the FITS header
;         geom: is the structure returned by chip_geometry.pro 
;
;CREATION:
; AT Oct 6, 2011
;
;MODIFICATION HISTORY:
;reworked to handle 4 amps and other problems ~MJG & DAF 201203-04
;--------------------------------------------------------
;-

function getimage, file, redpar, header=header, geom=geom

name = file ; for full name  

tmp = findfile(name, count=c)
if (c eq 0) then begin ; file not found
print, 'File '+name+' is not found, return ZERO'
return, 0
endif

im = double(readfits(name,header))

binsz = strt(fxpar(header, 'xbinning'))+strt(fxpar(header, 'ybinning'))

if redpar.biasmode eq 0 then begin
	fname = redpar.rootdir+redpar.biasdir+redpar.date+'_bin'+binsz+'_'+strt(rdspd)+'_medbias.dat'
	restore, fname

	;now subtract the median bias frame:
	im = im - bobsmed
	print, 'GETIMAGE: SUBTRACTED MEDIAN BIAS FRAME:'
	print, fname
endif

redpar.gain = 1d ; save gain for further use


; remember the binning in redpar
redpar.binning = [fxpar(header, 'xbinning'), fxpar(header, 'ybinning')]

im=rotate(im,3) 
return, im
end


