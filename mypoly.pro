;+
;
;  NAME: 
;     mypoly
;
;  PURPOSE: 
;   
;
;  CATEGORY:
;      UTILITIES
;
;  CALLING SEQUENCE:
;
;      mypoly
;
;  INPUTS:
;
;  OPTIONAL INPUTS:
;
;  OUTPUTS:
;
;  OPTIONAL OUTPUTS:
;
;  KEYWORD PARAMETERS:
;    
;  EXAMPLE:
;      mypoly
;
;  MODIFICATION HISTORY:
;        c. Matt Giguere 2011.02.25 06:36:31 PM
;
;-
function mypoly, X, P
  P = double(P)
  OUT = dblarr(n_elements(x)) + P[0]
  np = n_elements(P)
	 ;plot, x, out, ps=8, /ysty
    ;stop
    ;P[0] = Constant/ offset
    ;P[1] = linear term
    ;P[2] = qudratic term
    ;.
    ;.
    ;.
    ;P[N-1] = Constant /x-offset (must be included!!)
    
  for i=1, np-2 do begin
	 OUT += P[i]*(X-P[np-1])^i
	 ;plot, x, out, ps=8, /ysty
	 ;print, 'i is: ', i
	 ;print, 'P[i] is: ', P[i]
	 ;stop
  endfor
  return, OUT
end;mypoly.pro