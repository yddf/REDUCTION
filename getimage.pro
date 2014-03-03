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

; test if the file exists
;   name = redpar.rootdir + redpar.rawdir + file ; omit rootdir and datdir
   name =  file ; for full name  

   tmp = findfile(name, count=c)
   if (c eq 0) then begin ; file not found
     print, 'File '+name+' is not found, return ZERO'
     return, 0
    endif

    im = float(readfits(name,header))
    ; returns CCD geometry in a structure 
    if ~keyword_set(geom) then geom = chip_geometry(name, hdr=header, redpar=redpar)
	namps = n_elements(strsplit(sxpar(header, 'amplist'), ' '))

     if (geom.status ne 'OK') then begin ; wrong data
        print, 'Wrong header in file '+name
        return, 0
     endif

if namps eq 2 then begin
; This code will work for two-amplifier mode only!
; no lower amps 
       imleft = im[geom.image_trim.upleft[0]:geom.image_trim.upleft[1],*]
       biasleft = im[geom.bias_trim.upleft[0]:geom.bias_trim.upleft[1],*]
       gainleft = geom.gain.upleft
       
       imright = im[geom.image_trim.upright[0]:geom.image_trim.upright[1],*]
       biasright = im[geom.bias_trim.upright[0]:geom.bias_trim.upright[1],*]
       gainright = geom.gain.upright
       
       ron2 = [variance(biasleft), variance(biasright)] ; estimate of readout noise
       redpar.ron = sqrt(mean(ron2))  ; redaout noise in ADU

 ; Subtract the bias
	sziml=size(imleft)  &   szbl=size(biasleft) 
	if sziml[2] ne szbl[2] then stop, 'BIAS and image have different number of lines!' 
	for k=0,sziml[2]-1 do imleft[*,k] = imleft[*,k]- median(biasleft[*,k])

	szimr=size(imright)  &   szbr=size(biasright) 
	if szimr[2] ne szbr[2] then stop,  'BIAS and image have different number of lines!' 
	for k=0,szimr[2]-1 do imright[*,k] = imright[*,k]- median(biasright[*,k])

; Non-linearity correction
	mode = where(strtrim(geom.readout_speed) eq redpar.readmodes, count) ; readout mode
	if count eq 0 then begin 
    	print, 'Unknown readout mode! Return ZERO'
    	return, 0
	endif
    nlc = redpar.nlc[*,mode]  
    gains = redpar.gains[*,mode]

	if strtrim(geom.readout_speed) eq 'fast' then nlc = redpar.nlc[*,0]
	if strtrim(geom.readout_speed) eq 'normal' then nlc = redpar.nlc[*,1]
	imleft = imleft * (1d - nlc[0] * imleft)
	imright = imright * (1d - nlc[1] * imright)*gains[1] 

	if redpar.debug ge 2 then begin 
   		print, 'READ OUT SPEED IS: ', geom.readout_speed
   		print, 'BINNING: ', geom.bin.row, geom.bin.col
   		print, 'RON noise [left,right]: ', sqrt(ron2)
	endif
  	gainleft=redpar.gains[0,mode]  & gainleft = gainleft[0]
  	gainright=redpar.gains[1,mode]  & gainright = gainright[0]
  	imleft=imleft*gainleft
  	imright=imright*gainright
  	im=[imleft, imright]  ; join the two parts

; Trim the image?
 	sz = size(im)
 	if (redpar.xtrim[0] eq 0) and (redpar.xtrim[1] eq 0) then xtrim=[0,sz[1]-1] else xtrim=redpar.xtrim/geom.bin.row 
 	if (redpar.ytrim[0] eq 0) and (redpar.ytrim[1] eq 0) then ytrim=[0,sz[2]-1] else ytrim=redpar.ytrim/geom.bin.col 
 	im = im[xtrim[0]:xtrim[1], ytrim[0]:ytrim[1]]

endif ;2 amps
      
if namps eq 4 then begin
	   if redpar.biasmode eq 0 then begin
		   ;If the median bias frame option is set in ctio.par use this method:
		   rdspd = geom.readout_speed
		   if strt(sxpar(header, 'CCDSUM')) eq '3 1' then binsz = '31'
		   if strt(sxpar(header, 'CCDSUM')) eq '1 1' then binsz = '11'
		   if strt(sxpar(header, 'CCDSUM')) eq '4 4' then binsz = '44'
		   fname = redpar.rootdir+redpar.biasdir+redpar.date+'_bin'+binsz+'_'+strt(rdspd)+'_medbias.dat'
		   restore, fname
		   ;First subtract the median overscan from each quadrant. This part floats around over the course of the night, 
		   ;but the substructure doesn't change:
		   ;1. subtract median value from upper left quadrant (both image and overscan region):
		   idx = [0L, geom.bias_full.upleft[1], geom.image_trim.upleft[2], geom.image_trim.upleft[3]]
		   im[idx[0]:idx[1], idx[2]:idx[3]] -= $
			 median(im[geom.bias_trim.upleft[0]:geom.bias_trim.upleft[1], geom.bias_trim.upleft[2]:geom.bias_trim.upleft[3]])
		   ;2. now do the same for the upper right quadrant:
		   idx = [geom.bias_full.upright[0], fxpar(header, 'NAXIS1')-1, geom.image_trim.upright[2], geom.image_trim.upright[3]]
		   im[idx[0]:idx[1], idx[2]:idx[3]] -= $
			 median(im[geom.bias_trim.upright[0]:geom.bias_trim.upright[1], geom.bias_trim.upright[2]:geom.bias_trim.upright[3]])
		   ;3. and the bottom left quadrant:
		   idx = [0L, geom.bias_full.botleft[1], geom.image_trim.botleft[2], geom.image_trim.botleft[3]]
		   im[idx[0]:idx[1], idx[2]:idx[3]] -= $
			 median(im[geom.bias_trim.botleft[0]:geom.bias_trim.botleft[1], geom.bias_trim.botleft[2]:geom.bias_trim.botleft[3]])
		   ;4. now the bottom right:
		   idx = [geom.bias_full.botright[0], fxpar(header, 'NAXIS1')-1, geom.image_trim.botright[2], geom.image_trim.botright[3]]
		   im[idx[0]:idx[1], idx[2]:idx[3]] -= $
			 median(im[geom.bias_trim.botright[0]:geom.bias_trim.botright[1], geom.bias_trim.botright[2]:geom.bias_trim.botright[3]])
		   ;now subtract the median bias frame:
		   im = im - bobsmed
		   print, 'GETIMAGE: SUBTRACTED MEDIAN BIAS FRAME:'
		   print, fname
	   endif

       imupleft = im[geom.image_trim.upleft[0]:geom.image_trim.upleft[1],geom.image_trim.upleft[2]:geom.image_trim.upleft[3]]
       biasupleft = im[geom.bias_trim.upleft[0]:geom.bias_trim.upleft[1],geom.bias_trim.upleft[2]:geom.bias_trim.upleft[3]]
       gainupleft = geom.gain.upleft
       
       imupright = im[geom.image_trim.upright[0]:geom.image_trim.upright[1],geom.image_trim.upright[2]:geom.image_trim.upright[3]]
       biasupright = im[geom.bias_trim.upright[0]:geom.bias_trim.upright[1],geom.bias_trim.upright[2]:geom.bias_trim.upright[3]]
       gainupright = geom.gain.upright

       imbotleft = im[geom.image_trim.botleft[0]:geom.image_trim.botleft[1],geom.image_trim.botleft[2]:geom.image_trim.botleft[3]]
       biasbotleft = im[geom.bias_trim.botleft[0]:geom.bias_trim.botleft[1],geom.bias_trim.botleft[2]:geom.bias_trim.botleft[3]]
       gainbotleft = geom.gain.botleft
       
       imbotright = im[geom.image_trim.botright[0]:geom.image_trim.botright[1],geom.image_trim.botright[2]:geom.image_trim.botright[3]]
       biasbotright = im[geom.bias_trim.botright[0]:geom.bias_trim.botright[1],geom.bias_trim.botright[2]:geom.bias_trim.botright[3]]
       gainbotright = geom.gain.botright

      ron2 = [variance(biasupleft), variance(biasupright), variance(biasbotleft), variance(biasbotright)] ; estimate of readout noise
      redpar.ron = sqrt(mean(ron2))  ; readout noise in ADU

if redpar.biasmode eq 1 then begin
; If the median overscan option is set in ctio.par, then subtract the bias using this method:
	szimupl=size(imupleft)  &   szbupl=size(biasupleft) 
	if szimupl[2] ne szbupl[2] then stop, 'BIAS and image have different number of lines!' 
	for k=0,szimupl[2]-1 do imupleft[*,k] = imupleft[*,k]- median(biasupleft[*,k])

	szimupr=size(imupright)  &   szbupr=size(biasupright) 
	if szimupr[2] ne szbupr[2] then stop, 'BIAS and image have different number of lines!' 
	for k=0,szimupr[2]-1 do imupright[*,k] = imupright[*,k]- median(biasupright[*,k])

	szimbotl=size(imbotleft)  &   szbbotl=size(biasbotleft) 
	if szimbotl[2] ne szbbotl[2] then stop, 'BIAS and image have different number of lines!' 
	for k=0,szimbotl[2]-1 do imbotleft[*,k] = imbotleft[*,k]- median(biasbotleft[*,k])

	szimbotr=size(imbotright)  &   szbbotr=size(biasbotright)
	if szimbotr[2] ne szbbotr[2] then stop, 'BIAS and image have different number of lines!' 
	for k=0,szimbotr[2]-1 do imbotright[*,k] = imbotright[*,k]- median(biasbotright[*,k])
endif

; Non-linearity correction
	mode = where(strtrim(geom.readout_speed) eq redpar.readmodes, count) ; readout mode
	if count eq 0 then begin 
   	print, 'Unknown readout mode! Return ZERO'
   	return, 0
	endif
    nlc = redpar.nlc[*,mode]  
    gains = redpar.gains[*,mode]
    redpar.gain = gains[0] ; save gain for further use
;if (redpar.debug) then print, 'Non-inearity coefs: ', nlc 

;gains = gains/gains[0] ; force left gain to one

	if strtrim(geom.readout_speed) eq 'fast' then  nlc = redpar.nlc[*,0]
	if strtrim(geom.readout_speed) eq 'normal' then nlc = redpar.nlc[*,1]

; nlc not defined for Torrent 
;	imleft = imleft * (1d - nlc[0] * imleft)
;	imright = imright * (1d - nlc[1] * imright)*gains[1] 

if redpar.debug ge 2 then begin 
   print, 'READ OUT SPEED IS: ', geom.readout_speed
   print, 'BINNING: ', geom.bin.row, geom.bin.col
   print, 'RON noise [upleft,upright, botleft, botright]: ', sqrt(ron2)
endif

  gainupleft=redpar.gains[1,mode]  & gainupleft = gainupleft[0]
  gainupright=redpar.gains[2,mode]  & gainupright = gainupright[0]
  gainbotleft=redpar.gains[0,mode]  & gainbotleft = gainbotleft[0]
  gainbotright=redpar.gains[3,mode]  & gainbotright = gainbotright[0]

  imupleft=imupleft*gainupleft
  imupright=imupright*gainupright
  imbotleft=imbotleft*gainbotleft
  imbotright=imbotright*gainbotright
;print, redpar.gains[*,mode]
;stop
  im=[[imbotleft, imbotright],[imupleft, imupright]]  ; join the four parts

; Trim the image?
 sz = size(im)
 if (redpar.xtrim[0] eq 0) and (redpar.xtrim[1] eq 0) then xtrim=[0,sz[1]-1] else xtrim=redpar.xtrim/geom.bin.row 
 if (redpar.ytrim[0] eq 0) and (redpar.ytrim[1] eq 0) then ytrim=[0,sz[2]-1] else ytrim=redpar.ytrim/geom.bin.col 
 im = im[xtrim[0]:xtrim[1], ytrim[0]:ytrim[1]]

endif ;quad amp r/o
      

; remember the binning in redpar
  redpar.binning = [geom.bin.row, geom.bin.col]

 im=rotate(im,1) 		

  return, im

end
  

