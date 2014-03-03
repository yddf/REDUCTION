; Top-level script for CHIRON data reduction
; input: night like '111003'
; Oct 18, 2011 AT
pro allreduce, night

;modes=['slit','slicer', 'narrow', 'fiber']
;modes=['slit', 'narrow', 'fiber']
;modes=['slicer','fiber']
modes=['slicer']

for i=0, n_elements(modes)-1 do sorting_hat, night, mode=modes[i], /reduce, /getthid, /iod2fits
;for i=0, n_elements(modes)-1 do sorting_hat, night, mode=modes[i], /reduce,  /iod2fits
end

