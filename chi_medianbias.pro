;+
;
;  NAME: 
;     chi_medianbias
;
;  PURPOSE: 
;    To create and store median bias frames for the various modes for 
;	 bias subtraction
;
;  CATEGORY:
;      CHIRON
;
;  CALLING SEQUENCE:
;
;      chi_medianbias
;
;  INPUTS:
;
;  OPTIONAL INPUTS:
;
;  OUTPUTS:
;
;  KEYWORD PARAMETERS:
;    
;  EXAMPLE:
;      chi_medianbias, redpar = redpar, log = log, /bin31, /normal
;
;  MODIFICATION HISTORY:
;        c. Matt Giguere 2012.04.17 04:09:19 PM
;
;-
pro chi_medianbias, $
help = help, $
postplot = postplot, $
redpar = redpar, $
bin11 = bin11, $
bin31 = bin31, $
bin44 = bin44, $
normal = normal, $
fast = fast, $
log = log, $
bobsmed = bobsmed

if keyword_set(bin31) then binsz='31'
if keyword_set(bin11) then binsz='11'
if keyword_set(bin44) then binsz='44'
if keyword_set(normal) then rdspd = 'normal'
if keyword_set(fast) then rdspd = 'fast'

;locate Bias OBServations (bobs):
bobs = where(strt(log.object) eq 'bias' and strt(log.ccdsum) eq '3 1' and strt(log.speedmod) eq 'normal', bobsct)

if keyword_set(bin11) and keyword_set(normal) then begin
  bobs = where(strt(log.object) eq 'bias' and strt(log.ccdsum) eq '1 1' and strt(log.speedmod) eq 'normal', bobsct)
endif

if keyword_set(bin11) and keyword_set(fast) then begin
  bobs = where(strt(log.object) eq 'bias' and strt(log.ccdsum) eq '1 1' and strt(log.speedmod) eq 'fast', bobsct)
endif

if keyword_set(bin31) and keyword_set(fast) then begin
  bobs = where(strt(log.object) eq 'bias' and strt(log.ccdsum) eq '3 1' and strt(log.speedmod) eq 'fast', bobsct)
endif

if keyword_set(bin44) and keyword_set(normal) then begin
  bobs = where(strt(log.object) eq 'bias' and strt(log.ccdsum) eq '4 4' and strt(log.speedmod) eq 'normal', bobsct)
endif
print, 'Number of bias frames in '+binsz+' '+rdspd+' mode: ', bobsct
if bobsct gt 2 then begin
bcube = dblarr(long(log[bobs[0]].naxis1), long(log[bobs[0]].naxis2), bobsct)
for i=0, bobsct-1 do begin
  biasim = double(readfits(log[bobs[i]].filename, hd))
  geom = chip_geometry(hdr=hd)
  ;1. subtract median value from upper left quadrant (both image and overscan region):
  idx = [0L, geom.bias_full.upleft[1], geom.image_trim.upleft[2], geom.image_trim.upleft[3]]
  biasim[idx[0]:idx[1], idx[2]:idx[3]] -= $
  	median(biasim[geom.bias_trim.upleft[0]:geom.bias_trim.upleft[1], geom.bias_trim.upleft[2]:geom.bias_trim.upleft[3]])
  ;2. now do the same for the upper right quadrant:
  idx = [geom.bias_full.upright[0], log[bobs[0]].naxis1-1, geom.image_trim.upright[2], geom.image_trim.upright[3]]
  biasim[idx[0]:idx[1], idx[2]:idx[3]] -= $
  	median(biasim[geom.bias_trim.upright[0]:geom.bias_trim.upright[1], geom.bias_trim.upright[2]:geom.bias_trim.upright[3]])
  ;3. and the bottom left quadrant:
  idx = [0L, geom.bias_full.botleft[1], geom.image_trim.botleft[2], geom.image_trim.botleft[3]]
  biasim[idx[0]:idx[1], idx[2]:idx[3]] -= $
  	median(biasim[geom.bias_trim.botleft[0]:geom.bias_trim.botleft[1], geom.bias_trim.botleft[2]:geom.bias_trim.botleft[3]])
  ;4. now the bottom right:
  idx = [geom.bias_full.botright[0], log[bobs[0]].naxis1-1, geom.image_trim.botright[2], geom.image_trim.botright[3]]
  biasim[idx[0]:idx[1], idx[2]:idx[3]] -= $
  	median(biasim[geom.bias_trim.botright[0]:geom.bias_trim.botright[1], geom.bias_trim.botright[2]:geom.bias_trim.botright[3]])

  ;now save the bias image to a cube:
  bcube[*,*,i] = biasim
endfor
bobsmed = median(bcube, /double, dimen=3)
fname = redpar.rootdir+redpar.biasdir+redpar.date+'_bin'+binsz+'_'+rdspd+'_medbias.dat'
save, bobsmed, filename=fname
print, 'Median bias frome filename saved as: ', fname
endif
end;chi_medianbias.pro