pro ctio_dord,  ordfname, redpar, orc,ome, image=image
;  Determines  the default order locations for a particular
;  spectrograph setting. 
;
;INPUTS:
;   ORDFNAME   (input string) Filename of FITS file to be used
;   REDPAR      parameter structure. Current mode is passed there.
;   IMAGE      optional image to use for order location, e.g. summed flat
;
;OUTPUTS:
;     ORC  (array (# coeffs , # orders))] coefficients from the
;          polynomial fits to the order peaks.
;     OME  (optional output vector (# orders))] each entry gives the mean of the
;           absolute value of the difference between order locations and the polynomial
;           fit to these locations.
; 6-Oct-11 AT  Added redpar as argument to use getimage.pro 


if n_params() lt 3 then begin
  print,'syntax: ctio_dord,ordfname, redpar,orc[,ome].'
  retall
end

  print,'CTIO_DORD: Entering routine.'

  if ~keyword_set(image) then begin ; read order-location image from the disk
     image = getimage(ordfname, redpar, header=head)  
     if (size(image))[0] lt 2 then begin
       print, 'CTIO_DORD: Image is not found. Returning from all'
       retall
     endif
  endif

  swid = 32 ; hard-coded!!

;stop
fords,image,swid,orc, ome, redpar	;find order location coeffs

;  print,'CTIO_DORD: Saving order locations to ' + prefix + '.ord'
;  comment = 'Order location fit coefficients.'  
;  wdsk,orc,prefix + '.ord',comment,/new		;  save ORCs to disk
  
  return
end
