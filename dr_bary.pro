pro dr_bary

ff=file_search('/tous/mir7/logsheets/2012/12*.log',count=num)
for i=0,num-1 do print, i, ' ',ff[i]
stop
for i=0,num-1 do begin 
   pfx='chi'+strmid(ff[i],26, 6)
   print,ff[i]
   qbarylog,ff[i],baryDir='/tous/mir7/bary/', prefix=pfx
endfor 

barystruct_dbl,observ='ctio'

end  ; pro 
