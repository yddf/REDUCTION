; Convert UT midpoint time to float hours
; input; string array in 'hh:mm:ss' format
; output: ut, >24h after midnight
function gettime, mdpt
  n = n_elements(mdpt)
  ut = fltarr(n)
  for i=0,n-1 do begin
    s = mdpt[i]
    hh = fix(strmid(s,0,2)) & mm = fix(strmid(s,3,2)) & sec = fix(strmid(s,6,2))
    ut[i] = hh + mm/60. + sec/3600.
    if ut[i] lt 12. then ut[i] += 24. ; past midnight
  endfor
  return, ut
end


