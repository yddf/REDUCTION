pro addflat, flatfiles, sum, redpar, im_arr
; median-co-adds flat-field files, result in SUM
compile_opt idl2

numwf=n_elements(flatfiles) & owidefiles = flatfiles
print, "Entering ADDFLAT routine Nflats= ",numwf

sum = 0
im = getimage(flatfiles[0], redpar, header=header) ; read the first image
if (size(im))[0] lt 2 then return ; file not found

geom = chip_geometry(flatfiles[0], hdr=header, redpar=redpar) ; returns CCD geometry in a structure 
;stop
sz = size(im)
nc=sz[1]  &  nr=sz[2]
im_arr1=dblarr(nc,nr,numwf)
im_arr1[*,*,0]=im
imnsz = size(im)
swidth = 50L
gdfltidx = dblarr(numwf)
normvalarr = 0d

ctwf=0   
fspot = 0 ;index the im array data cube
for j = 0, numwf-1 do begin
	im = getimage(flatfiles[j], redpar, header=header, geom=geom) ; read  image, correct for bias and non-linearity        
	if redpar.flatnorm eq 1 then begin
		imswath = im[(sz[2]/2d - swidth):(sz[2]/2d + swidth),*]
		imswmed = median(imswath, dimen=1, /double)
		normval = max(imswmed)
		print,'flat #: ', strt(j),' max ADU/pixel: ', round(normval), ' minimum flat ADU/pixel: ',strt(round(redpar.minflatval))
		if normval ge redpar.minflatval then ctwf++
		if normval ge redpar.minflatval then gdfltidx[j] = 1
		if normval ge redpar.minflatval then normvalarr = [normvalarr, normval]
	endif 
	if redpar.flatnorm eq 0 then normval = 1d
	print, 'fspot is: ', fspot
	if redpar.flatnorm le 1 then im_arr1[*,*,fspot] = im/normval
   if normval ge redpar.minflatval then fspot++
endfor

;now to remove the spots for images that had too few counts:
if ctwf lt numwf then begin
  print, 'WARNING! You had flats that had too few counts! Now excluding them!'
  print, strt(ctwf)+' out of '+strt(numwf)+' are being used.'
  print, 'The flats being used are: '
  printt, flatfiles[where(gdfltidx eq 1)]
  print, 'The flats NOT being used are: '
  printt, flatfiles[where(gdfltidx ne 1)]
  stop
  im_arr = im_arr1[*,*,0:(ctwf-1)]
endif else im_arr = im_arr1
if redpar.flatnorm eq 1 then im_arr *= mean(normvalarr[1:*])

print, 'ADDFLAT: calculating median flat...'
sum = dblarr(nc,nr)
for ncol=0,nc-1 do begin
  for nrow=0,nr-1 do begin
	 sum[ncol,nrow]=median(im_arr[ncol,nrow,*])
  endfor
endfor


badp = where(sum le 0, nbadp)   ;find pixels le 0, and counts them
if nbadp gt 0 then     sum[badp] = 1.0           ;and fix them ?? is this a fix?
print, 'ADDFLAT: Now leaving routine.'
end
