pro mkwave,w,wvc,ord,col1,col2,vel=vel
;Program to build a wavelength scale from coefficients from solvarr.
; w (output vector) wavelength scale for order ord, columns col1-col2.
; wvc (input vector) wavelength information from fndwvc.pro
; ord (input scalar) absolute order number for which to find wavelengths.
; col1 (input scalar) starting column number of segment.
; col2 (input scalar) ending column number of segment.
; vel= (input scalar) radial velocity (in km/s) shift to apply to
;   wavelength scale (positive velocities increase wavelengths).
;04-Aug-92 JAV	Create.
;15-Aug-92 JAV	Now assume fit is vs. column/100, not column. Also now
;		 fitting cubic vs. column/100.
;19-Sep-92 JAV	Added huge new code section to handle version 2 wvc format.
;		 Can now build entire wavelength array, if requested.
;15-Jul-93 JAV,GB  Version 2.2; cross terms added; fits vs. order/100 now;
;		 bug in order polynomial fixed; increased nord to 2.
;03-Nov-93 GB	Added vel= keyword and logic.
;06-Dec-93 JAV	Documented last change.
;07-Jun-94 JAV	Negative ord, col1, and col2 are now equivalent to not
;		specifying the arguments at all.
;15-Jul-93 JAV,GB  Version 2.3; made cols/1000
;23-Jan-97 JAV	Return wavelength array in double precision.
;27-Nov-99 JAV	Allow 4 or 6 cross-terms (through version 2.5).

if n_params() lt 2 then begin
  print,'syntax: mkwave,w,coefs [,ord,col1,col2,vel=]'
  retall
endif

if n_elements(ord)  eq 0 then  ord = -1			;set default value
if n_elements(col1) eq 0 then col1 = -1			;set default value
if n_elements(col2) eq 0 then col2 = -1			;set default value

if col1 lt 0 then col1 = 0				;begin in first col

;Handle version 1.0 .wvc files.
  if n_elements(wvc) eq 6 then begin
    message,/info,'This is apparently an old wvc format.'
    message,/info,'You should convert to the new format.'
    if n_params() lt 3 then $
      message,'You must specify ord as third argument when using old wvc.'
    if col2 lt 0 then $
      col2 = long(wvc(5)) - 1			;build to last column
    len = col2 - col1 				;desired number of columns
    ic = (col1 + dindgen(len)) / 100d0		;column indicies / 100
    nl = (wvc(0) + wvc(4)*ord) $
      + ic*(wvc(1) + ic*(wvc(2) + ic*wvc(3))) ;calc. n*lambda efficiently
    w = nl / ord				;recover wavelengths
    return
  endif

;Assume we have new (version 2) format.
;Extract info from wvc information block.
  vers = wvc(0)					;file structure version number
  ncol = wvc(1)					;number of columns in spectrum
  nord = wvc(2)					;number of orders in spectrum
  obase = wvc(3)				;base order in spectrum
  fill = wvc(4:6)				;filler in reserved space
  ncross = fix(wvc(7))				;number of cross-terms
  coldeg = fix(wvc(8))				;degree of column fit poly
  orddeg = fix(wvc(9))				;degree of order fit poly
  coeff = wvc(10:*)				;fit coefficients

;Build column indicies/100.0 for polynomial construction.
  if col2 lt 0 then $
    col2 = long(ncol) - 1			;build to last column
  len = col2 - col1 + 1 			;desired number of columns
  ic = (col1 + dindgen(len)) / 100d0		;column indicies / 100
  if vers ge 2.3 then ic = ic / 10.d0		;column indicies / 1000 (2.3 up)

;Set the desired order range.
  if ord lt 0 then begin			;true: build entire array
    obeg = obase				;begin with zeroeth order
    oend = obase + nord - 1			;end with last order
  endif else begin				;else: only one order
    obeg = ord					;begin with the order we want
    oend = ord					;end with the order we want
  endelse
  no = oend - obeg + 1				;number orders to build

;Loop thru orders, building wavelength array.
  w = dblarr(len,no)				;init wavelength array
  for i=0,no-1 do begin				;loop thru build orders
    order = (obeg + i) / 100d0			;current order
    nlc = dblarr(len)				;init n*lambda column piece
    for j=coldeg,1,-1 do begin			;loop back thru col coeffs
      nlc = ic * (nlc + coeff(j))		;build polynomial
    endfor
    nlo = dblarr(len)				;init n*lambda order piece
    for j=coldeg+orddeg,coldeg+1,-1 do begin	;loop back thru order coeffs
      nlo = order * (nlo + coeff(j))		;build polynomial
    endfor
    if vers lt 2.2 then begin			;no cross-terms
      nlx = 0.0
    endif else begin				;cross-terms added
      inx = coldeg + orddeg +1			;start of cross-terms
      case ncross of
        4: nlx = order*ic *( coeff(inx) + ic*coeff(inx+1) $
	       + order*coeff(inx+2) + ic*order*coeff(inx+3) )
        6: nlx = order*ic *( coeff(inx) + ic*coeff(inx+1) $
	       + order*coeff(inx+2) + ic*order*coeff(inx+3) $
               + ic^2*coeff(inx+4) + order^2*coeff(inx+5))
        else: begin
                print, 'mkwave: unexpected number of cross-terms: ' $
                     + strtrim(ncross, 2)
                return
              end
      endcase
    endelse
    nl = coeff(0) + nlc + nlo + nlx		;build n*lambda
    w(*,i) = nl / double(obeg + i)		;recover wavelengths
  endfor

  if keyword_set(vel) then begin
    w = w * (1.0 + vel/2.9979246e5)		;apply velocity shift
  endif

  return
end
