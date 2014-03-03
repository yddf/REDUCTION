pro ck_snr

qtzs=readfits('achi120313.1019.fits') ; slit obs of qtz
qtzn=readfits('achi120313.1008.fits') ; ; narrow slit obs of qtz 

snr_arr=fltarr(16,4)  & slit_rat=fltarr(16)  & narr_rat=fltarr(16) 

print,'        Narrow       Narrow    Slit       Slit' 
print,'Order   sqrt(cts)    stddev    sqrt(cts)  stddev ' 
form1='(a4, a8, a12, a10, a10)'
for i=0,15 do begin  ; checking snr in 16 orders
    pixoffset=100*i
    pixst=500+pixoffset
    pixend=800+pixoffset
    xarr=indgen(301)
    ord=i+16
	snipn=qtzn[pixst:pixend,ord]
	coefn=poly_fit(xarr,snipn,1)
	snipn=median(snipn)*snipn/poly(xarr,coefn)   ;trend removed
	qn=median(snipn)
	stdn=stddev(snipn)
	strsqrtn=strcompress(string(fix(sqrt(qn))),/rem)
	strsdsn=strcompress(string(fix(qn/stddev(snipn))),/rem)
	
	snips=qtzs[pixst:pixend,ord]
	coefs=poly_fit(xarr,snips,1)
	snips=median(snips)*snips/poly(xarr,coefs)    ;trend removed
	qs=median(snips)
	stds=stddev(snips)
	strsqrts=strcompress(string(fix(sqrt(qs))),/rem)
	strsdss=strcompress(string(fix(qs/stddev(snips))),/rem)

	print,strcompress(string(ord),/rem), strsqrtn, strsdsn, strsqrts, strsdss, format=form1

snr_arr[i,0] = strsqrtn
snr_arr[i,1] = strsdsn
snr_arr[i,2] = strsqrts
snr_arr[i,3] = strsdss
narr_rat[i] = (qn/stddev(snipn))/sqrt(qn)
slit_rat[i] = (qs/stddev(snips))/sqrt(qs)

endfor

print,'Ratio of stddev to sqrt for narrow: ', median(narr_rat)
print,'Ratio of stddev to sqrt for slit: ', median(slit_rat)

end 


