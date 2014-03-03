;+
;
;  NAME: 
;     chi_barystruct
;
;  PURPOSE: To convert qbcvel.ascii into qbcvel.dat
;
;  CATEGORY:
;      CHIRON
;
;  INPUTS:
;		ASCIIFN: The filename of the ascii file.
;		STRUCTFN: The filename of the idl data structure
;
;  EXAMPLE:
;      chi_barystruct, asciifn='/tous/mir7/bary/qbcvel.ascii', $
;		structfn='/tous/mir7/bary/qbcvel.dat'
;
;  MODIFICATION HISTORY:
;        c. Matt Giguere 2013.11.13 14:36:53
;
;-
pro chi_barystruct, $
asciifn = asciifn, $
structfn = structfn

if ~keyword_set(asciifn) then asciifn = '/tous/mir7/bary/qbcvel.ascii'
if ~keyword_set(structfn) then structfn = '/tous/mir7/bary/qbcvel.dat'

nobs = file_lines(asciifn)
bcat_i = create_struct('obsnm', '', 'objnm', '', $
	'bc', 0d, 'jd', 0d, 'ha', 0d, 'obtype', '')
bcat = replicate(bcat_i, nobs)

openr, flun, asciifn, /get_lun
line=''
for i=0L, nobs-1 do begin
  readf, flun, line
  vars = strsplit(line, ' ', /extract)
  if n_elements(vars) gt 6 then stop
  bcat[i].obsnm = vars[0]
  bcat[i].objnm = vars[1]
  bcat[i].bc = vars[2]
  bcat[i].jd = vars[3]
  bcat[i].ha = vars[4]
  bcat[i].obtype = vars[5]
endfor;loop through lines
free_lun, flun
save, bcat, filename=structfn
end;chi_barystruct.pro