; Rough estimate of rms noise in extracted FF spectrum
; input: 1D order
pro noise, order, ron=ron, xwid=xwid, gain=gain

;gain = 2.59 ; el/ADU hard-wired
if ~keyword_set(gain) then gain = 2.3 ; el/ADU hard-wired

n = n_elements(order)
x = findgen(n)

nmin = 2*(n/6) ; even number
nmax = 2*(n/3)-1 ; odd number  

c = poly_fit(x,order, 5)
y = poly(x,c) ; smoothed order
d = (order - y)/y
d1 = d[nmin:nmax] ; central part

;plot, d1, xs=1, ys=1

d0 = mean(order[nmin:nmax])
sig2 = variance(d1) 

sig2phot = 1./(d0*gain) 

if keyword_set(ron) and keyword_set(xwid) then begin ; add readout noise
  rnoise2 = (ron)^2*xwid/d0^2
  print, 'RON contribution: ', sqrt(rnoise2)
  sig2phot += rnoise2
endif

chi2 = sig2/sig2phot


;stop

n1 = n_elements(d1)
power = abs( (fft(d1))[0:n1/2])^2
freq = findgen(n1/2)/n1

plot, freq, power, xtitle='Frequency, pix^-1', ytitle='Power'

; high-frequency noise
sig2hf = 4.*total(power[n1/4:n1/2])

fmt = '(A,10F8.4)'
print, 'Mean signal: ', d0
print, 'RMS relative [TOT, HF, photon]: ', sqrt([sig2, sig2hf, sig2phot]), f=fmt
print, 'CHI2 [all, HF]: ', chi2, sig2hf/sig2phot, f=fmt

;print, 'RMS white-noise: ', sqrt(sig2hf)

;stop

end
