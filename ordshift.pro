; Check vertical order shift in the spectrum in the specified order
pro ordshift, im, orc, order, xwid

poff = xwid/2+1  ; half-width
n = 2*poff + 1   ; pixels across order

offset = findgen(n)- poff ; vertical offset 

  sz = size(im)				;variable info block
  ncol = sz(1)				;number of cols in image
  nrow = sz(2)				;number of rows in image

x = findgen(ncol)
y = poly(x,orc(*,order))
prof = fltarr(n)
for k=0,n-1 do for i=0,ncol-1 do prof[k]+=im[i,y[i]+offset[k]+0.5]

center = total(offset*prof)/total(prof)
print, 'Order center [pixels]: ', center
plot, offset, prof

;stop

end
