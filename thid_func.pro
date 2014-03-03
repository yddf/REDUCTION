function thid_func, dummy, par
;Echelle wavelength fitting function.
;Inputs:
; dummy (ignored) order and pixel number passed by common block
; par (vector(15)) coefficients of order*lambda fit:
;   [const, X^[1,2,...,maxordx], Y^[1,2,...,maxordy] $
;         , X*Y,X^2*Y,X*Y^2,X^2*Y^2,X^3*Y,X*Y^3]
;Returns:
; predicted order*wavelength at order and pixel values passed by common block.
;History:
; 28-Jun-99 Valenti  Wrote.
; 24-Nov-99 Valenti  Generalized to allow arbitrary maximum order.

if n_params() lt 2 then begin
  print, 'syntax: wair = thid_func(dummy, par)'
  return, 0
endif

common thid_common, thid_order, thid_pixel, thid_maxordx, thid_maxordy

;Constant term.
  mlam = par(0)

;X terms.
  base = 1
  xterms = 0d0
  for i=base+thid_maxordx-1, base, -1 do begin
    xterms = thid_pixel * (xterms + par(i))
  endfor
  mlam = mlam + xterms

;Y terms.
  base = thid_maxordx + 1
  yterms = 0
  for i=base+thid_maxordy-1, base, -1 do begin
    yterms = thid_order * (yterms + par(i))
  endfor
  mlam = mlam + yterms

;Cross terms.
  base = thid_maxordx + thid_maxordy + 1
  cross = thid_pixel     * thid_order     * par(base)    $
        + thid_pixel^2.0 * thid_order     * par(base+1)  $
        + thid_pixel     * thid_order^2.0 * par(base+2)  $
        + thid_pixel^2.0 * thid_order^2.0 * par(base+3)  $
        + thid_pixel^3.0 * thid_order     * par(base+4)  $
        + thid_pixel     * thid_order^3.0 * par(base+5)
  mlam = mlam + cross

;Return result.
  return, mlam

end
