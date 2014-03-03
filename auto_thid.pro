pro auto_thid, spec, init, ordx, ordy, othresh, thid $
             , awin=awin, maxres=maxres, minamp=minamp $
             , nclip=nclip, plot=plot, orev=orev, $
             kwnord = kwnord
;Automatically identify thorium lines in an extracted echelle spectrum
;of a thorium lamp, using a previous solution as an initial guess.
;
;Input:
; spec (array(npix,nord)) extracted echelle spectrum of a thorium lamp
; init (vector) wavelength coefficients (e.g. thid.wvc) that will be used
;  to construct an initial guess for the wavelength scale.
; maxres (scalar) maximum allowed residual when 
; ordx (scalar) degree polynomial to fit vs. pixel number
; ordy (scalar) degree polynomial to fit vs. order number
; othresh (scalar) maximum residual (in pixels) used during outlier rejection
; awin= (scalar) half-width (in pixels) of window when fitting line centers
; maxres= (scalar) maximum allowed residual (in pixels) when finding lines
; minamp= (scalar) minimum counts (in ADU) in peaks of lines to find and fit
; nclip= (scalar) number of points to clip each iteration
; plot= (scalar) diagnostic plot level (0=none, 1=summary, 2=more, etc.)
; /orev (switch) indicates that orders are in descending order.
; kwnord: a keyword that sets the numbers of orders. This was needed for
;	running the CHIRON data. 
;
;Output:
; thid (structure) information about lines used in fit
;
;Notes:
; The directory containing the program thid.pro must also contain an IDL
; save file named "thid_data.xdr", which in turn must contain the variables:
;        th_wair     air wavelengths (in Angstroms)
;        th_wnum     wavenumbers (in 1/cm)
;        th_inten    line relative intensities
; On non-Unix systems, the line data file may be placed in the current
; working directory.
;
;History:
; 05-Jul-99 Valenti  Coding of initial version completed.
; 24-Nov-99 Valenti  Replaced fspline with interpol, but meaintained speed.
;                    Output version 2.4 -> 2.5 allows maxord up to 6.
; 01-Sep-00 Valenti  Allow staring guess to be specified as polynomial fit
;                    vs. column number to one or more orders (e.g. from disp).
; 24-May-02 Valenti  Created auto_ version from original thid.pro.
; 20110909: MJG ~ kwnord added

if n_params() lt 6 then begin
  print, 'syntax: auto_thid, spec, init, ordx, ordy, othresh, thid'
  print, '                [, awin=, maxres= ,minamp=, nclip= ,plot= ,/orev ]'
  return
endif

common thid_common, thid_order, thid_pixel, thid_maxordx, thid_maxordy

;Defualt values for optional parameters.
  if n_elements(awin) eq 0 then awin = 5	;half-width (pixels)
  if n_elements(maxres) eq 0 then maxres = 1.0	;max residual (pixels)
  if n_elements(minamp) eq 0 then minamp = 5.0	;min amplitude (ADU)
  if n_elements(nclip) eq 0 then nclip = 10	;points to clip each iter
  if n_elements(plot) eq 0 then plot = 0	;diagnostic plot level

;Set default location for line data file.
  help, calls=calls				;get stack of IDL calls
  call = calls(0)				;extract call to this routine
  ipos1 = strpos(call, '/')			;beginning of unix path
  ipos2 = rstrpos(call, '/')			;end of unix path
  if ipos1 ge 0 and ipos2 ge 0 then begin
    dir = strmid(call, ipos1, ipos2-ipos1+1)	;look in dir with thar.pro
  endif else begin
    dir = ''					;look in current directory
  endelse

;Internal parameters.
  maxlist = 10000				;maximum lines to identify
  estsig = 1.0					;estimated line sigma (pixels)
  dmax = 5.0					;inital plot maximum
  maxordx = 6					;maximum polynomial in x
  maxordy = 6					;maximum polynomial in y
  maxcross = 6

;Check that requested polynomial orders are allowed.
  if ordx lt 0 or ordx gt maxordx then begin
    print, 'ordx (' + strtrim(ordx, 2) + ') must be between 0 and ' $
         + strtrim(maxordx, 2)
  endif
  if ordy lt 0 or ordy gt maxordy then begin
    print, 'ordy (' + strtrim(ordy, 2) + ') must be between 0 and ' $
         + strtrim(maxordy, 2)
  endif

;Load parameters shared with thid_func.pro
  thid_maxordx = maxordx
  thid_maxordy = maxordy

;Useful sizes.
  sz = size(spec)				;variable info block
  if sz(0) ne 2 then begin			;true: not an array
    message, /info, 'spectrum must be a two-dimensional array'
    return
  endif
  npix = sz(1)					;number of pixels
  nord = sz(2)					;number of orders

;Extract information from initial guess.
  obase = init(3)				;base order number

;Useful variables.
  plist = indgen(npix)				;list of pixel indices
  olist = obase + indgen(nord)			;list of orders
  if keyword_set(orev) then begin
    olist = rotate(olist, 5)
  endif
  thmax = fltarr(nord)				;line intensity scaling
  thmin = fltarr(nord)				;line intensity scaling
  m_wave = dblarr(maxlist)			;marked wavelength
  m_ord = fltarr(maxlist)			;marked order
  m_pix = fltarr(maxlist)			;marked pixel
  m_wid = fltarr(maxlist)			;marked gaussian sigma
  m_flag = intarr(maxlist)			;marked line flags
  par = dblarr(1+maxordx+maxordy+maxcross)	;fit coefficients
  nm = 0					;number of marked lines

;Read line data from IDL save file.
  restore, dir + 'thid_data.xdr'
  i_argon=where(th_inten lt 0, n_argon)
  if(n_argon gt 0) then th_inten(i_argon)=-th_inten(i_argon)

specsz = size(spec)
nord = specsz[2]
;Preprocess observed spectrum to enhance plot visibility.
  s = fltarr(npix, nord)		;init normalized spectrum
  for i=0, nord-1 do begin		;loop thru orders
    tmp = spec(*,i)
    contf, tmp, cfit, /bot, frac=0.5
    tmp = tmp - cfit
    tmp = tmp / median(abs(tmp))
    s(*,i) = tmp
  endfor

;Initial guess for wavelength scale.
  mkwave, wave, init
  if keyword_set(orev) then wave = rotate(wave, 7)

; AT: match input wavelength to the actual binning (number of pixels)
    if ((size(wave))[1] ne npix) then begin ; match the dimensions of previous solution
       print, 'Re-formatting the initial wavelength solution!'
       tmp = dblarr(npix,nord)
       for i=0,nord-1 do tmp[*,i] = congrid(wave[*,i], npix)
       wave = tmp
    endif

;Save current plot state.
  if plot ge 1 then begin
    xsav = !x
    ysav = !y
    psav = !p

;Set desired plot characteristics.
    !p.charsize=1.5
    !x.margin=[8,2]
    !y.margin=[4,1]

;Set up the color table.
    loadct, 0, /silent
    colblu = !d.table_size-2
    colgrn = !d.table_size-3
    colorg = !d.table_size-4
    colred = !d.table_size-5
    colyel = !d.table_size-6
    tvlct, ct_r, ct_g, ct_b, /get
    ct_r(colblu) = 128
    ct_g(colblu) = 128
    ct_b(colblu) = 255
    ct_r(colgrn) = 0
    ct_g(colgrn) = 200
    ct_b(colgrn) = 0
    ct_r(colorg) = 255
    ct_g(colorg) = 128
    ct_b(colorg) = 0
    ct_r(colred) = 255
    ct_g(colred) = 0
    ct_b(colred) = 0
    ct_r(colyel) = 255
    ct_g(colyel) = 255
    ct_b(colyel) = 0
    tvlct, ct_r, ct_g, ct_b
  endif

;Automatically associate lines.
  nm = 0
  for i=0, nord-1 do begin
    onum = olist(i)
    wbeg = min(wave(*,i), max=wend)
    iwhr = where(th_wair ge wbeg and th_wair le wend, nwhr)
    for j=0, nwhr-1 do begin
      th_wav = th_wair(iwhr(j))
      th_pix = (interpol(findgen(npix), wave(*,i), th_wav))(0)
      ipix = round(th_pix)
      if ipix gt awin and ipix lt (npix-awin-1) then begin
        xfit = ipix + dindgen(2*awin+1) - awin
        yfit = s(ipix-awin:ipix+awin,i)
        bg = min(yfit)
        amp = max(yfit, imax) - bg
        if amp gt minamp and imax ne 0 and imax ne 2*awin then begin
          gpar = double([amp, ipix, estsig, bg])
;          gfit = gauss_fit(xfit, yfit, gpar) ; uses home-made  gauss_fit and curfit, slow
          gfit = gaussfit(xfit, yfit, gpar, nterm=4)
          if abs(gpar(1) - th_pix) le maxres and $
             gpar(2) gt 0 then begin
            m_wave(nm) = th_wav
            m_ord(nm) = olist(i)
            m_pix(nm) = gpar(1)
            m_wid(nm) = gpar(2)
            m_flag(nm) = 1
            nm = nm + 1
          endif
        endif
      endif
    endfor
  endfor
  print, strtrim(nm, 2) + ' lines found initially'

  if plot ge 1 then begin
    print, 'Half-width of window to use when fitting lines: ' $
         + strtrim(string(awin, form='(i9)'),2) + ' pixels'
    print, 'Maximum allowed residual when associating lines: ' $
         + strtrim(string(maxres, form='(f9.1)'),2) + ' pixels'
    print, 'Minimum amplitude of lines to fits: ' $
         + strtrim(string(minamp, form='(f9.1)'),2) + ' ADU'
  endif

;Decide on cross-terms (use all that are OK).
  xok = intarr(maxcross)
  xlist = ''
  nxok = 0
  if ordx ge 1 and ordy ge 1 and nm ge ordx+ordy+nxok+2 then begin
    xok(0) = 1
    xlist = xlist + 'XY '
    nxok = nxok + 1
  endif
  if ordx ge 2 and ordy ge 1 and nm ge ordx+ordy+nxok+2 then begin
    xok(1) = 1
    xlist = xlist + 'XXY '
    nxok = nxok + 1
  endif
  if ordx ge 1 and ordy ge 2 and nm ge ordx+ordy+nxok+2 then begin
    xok(2) = 1
    xlist = xlist + 'XYY '
    nxok = nxok + 1
  endif
  if ordx ge 2 and ordy ge 2 and nm ge ordx+ordy+nxok+2 then begin
    xok(3) = 1
    xlist = xlist + 'XXYY '
    nxok = nxok + 1
  endif
  if ordx ge 3 and ordy ge 1 and nm ge ordx+ordy+nxok+2 then begin
    xok(4) = 1
    xlist = xlist + 'XXXY '
    nxok = nxok + 1
  endif
  if ordx ge 1 and ordy ge 3 and nm ge ordx+ordy+nxok+2 then begin
    xok(5) = 1
    xlist = xlist + 'XYYY '
    nxok = nxok + 1
  endif
  nxok = total(xok)
  if nxok gt 0 and plot ge 0 then print, 'Cross terms: ' + xlist

; ----- Begin outlier rejection loop.
  repeat begin

;Identify associated lines.
    nfree = ordx + ordy + nxok + 1
    m_w = double(m_wave(0:nm-1))
    m_p = double(m_pix(0:nm-1))
    m_o = double(m_ord(0:nm-1))
    m_f = double(m_flag(0:nm-1))
    igd = where(m_f ne 0, ngd)
    if ngd lt nfree then begin
      message, 'too few lines (' + strtrim(ngd, 2) + ') to fit'
    endif
    m_w = m_w(igd)
    m_p = m_p(igd)
    m_o = m_o(igd)
    m_f = m_f(igd)
    m_ml = m_o * m_w

;Load common block used by fitting routine.
    thid_order = m_o / 1d2
    thid_pixel = m_p / 1d3

;Use Marquardt to fit marked lines.
    dummy = fltarr(ngd)
    sig = replicate(1.0, ngd)
    par = dblarr(maxordx + 1 + maxordy + maxcross)
    par(0:1) = poly_fit(m_p/1d3, m_ml, 1, /double)
    dpar = [ 1d1, replicate(1, maxordx) $
                , replicate(1, maxordy) $
                , replicate(1, maxcross) ]
    if ordx lt maxordx then dpar(ordx+1:maxordx) = 0.0
    if ordy lt maxordy then dpar(ordy+maxordx+1:maxordx+maxordy) = 0.0
    iwhr = where(xok eq 0, nwhr)
    if nwhr gt 0 then dpar(maxordx+maxordy+1+iwhr) = 0.0
    mlfit = marq('thid_func', dummy, m_ml, sig, par, dpar, trace=(plot ge 2))

;Calculate wavelengths for every pixel.
    dummy = fltarr(npix)
    thid_pixel = dindgen(npix) / 1d3
    for i=0, nord-1 do begin
      thid_order = replicate(olist(i), npix) / 1d2
      wave(*,i) = thid_func(dummy, par) / olist(i)
    endfor

;Report residuals.
    resid = fltarr(ngd)
    xpix = findgen(npix)
    for iord=0, nord-1 do begin
      ord = olist(iord)
      iwhr = where(m_o eq ord, nwhr)
      if nwhr gt 0 then begin
        resid(iwhr) = m_p(iwhr) - interpol(xpix, wave(*,iord), m_w(iwhr))
      endif
    endfor
    if plot gt 0 then begin
      print, strtrim(ngd, 2) + ' good lines' $
           + ';  Residuals' $
           + ': Min=' + strtrim(string(min(resid), form='(f10.2)'), 2) $
           + ', Max=' + strtrim(string(max(resid), form='(f10.2)'), 2) $
           + ', RMS=' + strtrim(string(stddev(resid), form='(f10.2)'), 2) $
           + ' pixels'
    endif

;Calculate residuals.
    if not keyword_set(resid) then begin
      message, 'attempt to reject outliers without fitting line locations'
    endif
    m_w = double(m_wave(0:nm-1))
    m_p = double(m_pix(0:nm-1))
    m_o = double(m_ord(0:nm-1))
    m_f = double(m_flag(0:nm-1))
    resid = fltarr(nm)
    xpix = findgen(npix)
    for iord=0, nord-1 do begin
      ord = olist(iord)
      iwhr = where(m_o eq ord, nwhr)
      if nwhr gt 0 then begin
        resid(iwhr) = m_p(iwhr) - interpol(xpix, wave(*,iord), m_w(iwhr))
      endif
    endfor

;Reject worst remaining outliers exceeding threshold.
    good = (m_f gt 0)
    igd = where(good eq 1)
    iout = where(abs(resid(igd)) gt othresh, nout)
    if nout gt 0 then begin
      if nout gt nclip then nc = nclip else nc = 1
      iout = igd(iout)
      isort = sort(abs(resid(iout)))		;sort residuals
      ibad = iout(isort(nout-nc:nout-1))	;worst outliers
      good(ibad) = 0				;mark as bad
    endif
    m_flag(igd) = 2
    ibad = where(good eq 0, nbad)
    if nbad gt 0 then m_flag(ibad) = 0
    if plot ge 2 then print, strtrim(ngd, 2) + ' good lines'

;Diagnostic plot.
    if plot ge 2 then begin
      !p.multi = [0, 1, 2]
      plot, m_p(igd), resid(igd), /ps, xsty=3, ysty=3, chars=1.5 $
          , xtit='Pixel Number', ytit='Residual (pixels)'
      oplot, !x.crange, [0,0], co=colred
      plot, m_o(igd), resid(igd), /ps, xsty=3, ysty=3, chars=1.5 $
          , xtit='Order Number', ytit='Residual (pixels)'
      oplot, !x.crange, [0,0], co=colred
      medres = replicate(1e10, nord)
      for i=0,nord-1 do begin
        iwhr = where(m_o(igd) eq olist(i), nwhr)
        if nwhr gt 0 then medres(i) = median(resid(igd(iwhr)))
      endfor
      oplot, olist, medres, max=1e9, co=colblu
      !p.multi = 0
      print, 'AUTO-THID: press any key to continue, q to quit'
      if ngd lt 7d2 then junk = get_kbrd(1) else junk=' '
      if junk eq 'q' then retall
    endif

; ----End of outlier rejection loop.
  endrep until nout eq 0

; final report
     print, strtrim(ngd, 2) + ' good lines;  Residuals' $
           + ': Min=' + strtrim(string(min(resid[igd]), form='(f10.2)'), 2) $
           + ', Max=' + strtrim(string(max(resid[igd]), form='(f10.2)'), 2) $
           + ', RMS=' + strtrim(string(stddev(resid[igd]), form='(f10.2)'), 2) $
           + ' pixels'

;Diagnostic plot.
  if plot ge 1 then begin
    !p.multi = [0, 1, 2]
    plot, m_p(igd), resid(igd), /ps, xsty=3, ysty=3, chars=1.5 $
        , xtit='Pixel Number', ytit='Residual (pixels)'
    oplot, !x.crange, [0,0], co=colred
    plot, m_o(igd), resid(igd), /ps, xsty=3, ysty=3, chars=1.5 $
        , xtit='Order Number', ytit='Residual (pixels)'
    oplot, !x.crange, [0,0], co=colred
    medres = replicate(1e10, nord)
    for i=0,nord-1 do begin
      iwhr = where(m_o(igd) eq olist(i), nwhr)
      if nwhr gt 0 then medres(i) = median(resid(igd(iwhr)))
    endfor
    oplot, olist, medres, max=1e9, co=colblu
    !p.multi = 0
     print, 'AUTO-THID: press any key to continue, q to quit'
    if ngd lt 2d2 then junk = get_kbrd(1) else junk=' '
    if junk eq 'q' then retall
  endif

;2D diagnostic plot of spectrum and fit.
  if plot ge 1 then begin
    image = imgexp(s, plist, olist, xs, ys, xran, yran, pos=dev_pos)
    image = imgscl(image, min=0, max=dmax, top=!d.table_size-7)
    erase
    tv, image, /dev, dev_pos(0), dev_pos(1) $
             , xsize=dev_pos(2)-dev_pos(0) $
             , ysize=dev_pos(3)-dev_pos(1)
    plot, [0,1], /noerase, /xsty, /ysty, /dev, /nodata $
        , pos=dev_pos, xr=xran, yr=yran, yticks=yticks $
        , xtit='Pixel Number', ytit='Order Number'
    for i=0, nord-1 do begin
      onum = olist(i)
      wbeg = min(wave(*,i), max=wend)
      iwhr = where(th_wair ge wbeg and th_wair le wend, nwhr)
      if nwhr gt 0 then begin
        th_wav = th_wair(iwhr)
        th_pix = interpol(findgen(npix), wave(*,i), th_wav)
        th_int = alog10(th_inten(iwhr))
        thmax(i) = max(th_int)
        thmin(i) = min(th_int)
        th_int = (th_int - thmin(i)) / (thmax(i) - thmin(i))
        th_int = 0.8*th_int + 0.1			;map into 0.1-0.9
        for j=0, nwhr-1 do begin
          color = colblu
          ii = where(iwhr(j) eq i_argon, n_argon)
          if(n_argon gt 0) then  color=colyel
          if nm gt 0 then begin
            if min(abs(m_wave - th_wav(j)), imin) lt 0.0001 then begin
              case m_flag(imin) of
                0: color = colred
                1: color = colorg
                2: color = colgrn
                else: message, 'Unknown flag value'
              endcase
            endif
          endif
          oplot, th_pix(j)+[0,0], onum+[-0.45,0.45]*th_int(j) $
               , co=color,th=1+(th_int(j) gt 0.5)
        endfor
      endif
    endfor

;Restore plot state.
    !p = psav
    !y = ysav
    !x = xsav
    loadct, 0
  endif

;Construct wavelength coefficients.
  vers = 2.5
  fill = 0.0
  wvc = [ vers, npix, nord, obase, fill, fill, fill $
        , maxcross, maxordx, maxordy, par ]

;Construct line statistics structure.
  if nm gt 0 then begin
    m_w = double(m_wave(0:nm-1))
    m_p = float(m_pix(0:nm-1))
    m_o = fix(m_ord(0:nm-1))
    m_f = m_flag(0:nm-1)
    fwhm = m_wid(0:nm-1) * sqrt(alog(256))
    resid = fltarr(nm)
    xpix = findgen(npix)
    for iord=0, nord-1 do begin
      ord = olist(iord)
      iwhr = where(m_o eq ord, nwhr)
      if nwhr gt 0 then begin
        resid(iwhr) = m_p(iwhr) - interpol(xpix, wave(*,iord), m_w(iwhr))
      endif
    endfor
    good = (m_f gt 0)
    igd = where(good eq 1, ngd)

;Report focus information.
    print, 'Median line FWHM = ' $
         + strtrim(string(median(fwhm), form='(f9.2)'), 2) $
         + ' pixels'
    resol = median( 0.5*(wave(1:npix-1,*) + wave(0:npix-2,*))  $
                       /(wave(1:npix-1,*) - wave(0:npix-2,*))) $
          / median(fwhm)
    resol = abs(resol)
    print, 'Median resolution = ' + strtrim(round(resol), 2)

    if ngd gt 1 then rmspix = stddev(resid(igd)) else rmspix = -1.0
    rmsang = rmspix * median(wave(1:npix-1,*) - wave(0:npix-2,*))
    thid = { created: systime() $
           , vers: vers         $
           , wvc: wvc           $
           , nlin: ngd          $
           , wair: m_w(igd)     $
           , order: m_o(igd)    $
           , pixel: m_p(igd)    $
           , resid: resid(igd)  $
           , rms_pix: rmspix    $
           , rms_ang: rmsang    $
           , fwhm: fwhm(igd)    $
           , resol: resol       }

  endif else begin
    lines = 0
  endelse

end
