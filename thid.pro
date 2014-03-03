pro thid, spec, obase, mlam, wvc, thid, init=init, orev=orev $
        , inpord=inpord, inpcoef=inpcoef, fitdeg=fitdeg
;Identify thorium lines in an extracted echelle spectrum of a thorium lamp.
;
;Input:
; spec (array(npix,nord)) extracted echelle spectrum of a thorium lamp
; obase (scalar) lowest order number in spectrum
; mlam (vector(2)) approximate range of order number times wavelength in
;   Angstroms, covered by the spectrum
; init= (vector(21)) wavelength coefficients that will be used to construct
;   an initial guess for the wavelength scale.
; /orev (switch) indicates that orders are reversed (i.e. in descending
;   order. Obase is still interpreted as the lowest number order.
; inpord= (scalar or vector(niord) list of orders with polynomial wavelength
;   fit coefficients in inpcoef.
; inpcoef= (vector(nicoef) or array(nicoef,niord)) polynomial coefficients
;   giving wavelength versus pixel number (starting at 0) for each order
;   listed in inpord.
; fitdeg= (scalar) degree of polynomial that should be used to interpolate
;   and/or extrapolate m*lambda vs. order number (m), using the input
;   wavelength fits in inpord= and inpcoef=.
;Output:
; wvc (vector(21)) coefficients of wavelength fit in standard form
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
; 01-Sep-00 Valenti  Allow starting guess to be specified as polynomial fit
;                    vs. column number to one or more orders (e.g. from disp).

if n_params() lt 4 then begin
  print, 'syntax: thid, spec, obase, mlam, wvc'
  print, '      [,thid ,/orev ,init= ,inpord= ,inpcoef=, fitdeg=]'
  return
endif

common thid_common, thid_order, thid_pixel, thid_maxordx, thid_maxordy

;Defualt values for optional parameters.
  if n_elements(fitdeg) eq 0 then fitdeg = 1	;degree to fit input vs. order

;Set default location for line data file.
  help, calls=calls				;get stack of IDL calls
  call = calls(0)				;extract call to this routine
  ipos1 = strpos(call, '/')			;beginning of unix path
  ipos2 = strpos(call, '/', /reverse_search)	;end of unix path
  if ipos1 ge 0 and ipos2 ge 0 then begin
    dir = strmid(call, ipos1, ipos2-ipos1+1)	;look in dir with thar.pro
  endif else begin
    dir = ''					;look in current directory
  endelse

;Internal parameters.
  maxlist = 10000				;maximum lines to identify
  dmax = 5.0					;inital plot maximum
  mwin = 15.0					;mark line window (pixels)
  maxordx = 6					;maximum polynomial in x
  maxordy = 6					;maximum polynomial in y
  maxcross = 6

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
  if keyword_set(init) then begin	;initial guess in output wvc format

;Initial solution from existing 2D wavelength fit coefficients.
    mkwave, wave, init
    if keyword_set(orev) then wave = rotate(wave, 7)

;stop
; AT: match input wavelength to the actual binning (number of pixels)
    if ((size(wave))[1] ne npix) then begin ; match the dimensions of previous solution
       print, 'Re-formatting the initial wavelength solution!'
       tmp = dblarr(npix,nord)
       for i=0,nord-1 do tmp[*,i] = congrid(wave[*,i], npix)
       wave = tmp
    endif

  endif else begin			;initial guess from single order fits
    wave = dblarr(npix, nord)
    if keyword_set(inpord) or keyword_set(inpcoef) then begin

;Initial solution from inpord= and inpcoef= keywords.
      if not keyword_set(inpord) or not keyword_set(inpcoef) then begin
        message, /info, 'must supply both inpord= and inpcoef='
        return
      endif
      ipix = dindgen(npix)
      niord = n_elements(inpord)
;Check that reduced spectrum contains all orders with input fits.
      for i=0, niord-1 do begin
        relord = inpord(i) - obase
        if relord lt 0 or relord ge nord then begin
          message, /info, 'orders in inpord= must be in range ' $
                 + strtrim(obase,2) + '-' + strtrim(obase+nord-1, 2)
          return
        endif
      endfor
;If fit was provided for only one order, scale by ratio of order number.
      if niord eq 1 then begin
        wiord = poly(ipix, inpcoef)
        for i=0, nord-1 do begin
          wave(*,i) = wiord * double(inpord(0)) / (obase+i)
        endfor
;If fit provided for multiple orders, linearly interpolate vs. order.
      endif else begin
        nterms = (size(inpcoef))(1)
        coco = dblarr(fitdeg+1,nterms)
        for j=0, nterms-1 do begin		;fit inpcoef vs. order
          coco(*,j) = poly_fit(inpord, inpord*reform(inpcoef(j,*)) $
                              ,fitdeg, /double)
        endfor
        for i=0, nord-1 do begin
          ordcoef = dblarr(nterms)
          for j=0, nterms-1 do ordcoef(j) = poly([obase+i], coco(*,j))
          wave(*,i) = poly(ipix, ordcoef) / (obase+i)
        endfor
      endelse

;Construct constant dispersion solution from input m*lambda range.
    endif else begin			;initial guess from m*lambda range
      mlam_vec = mlam(0) + (mlam(1) - mlam(0)) * dindgen(npix) / (npix - 1)
      for i=0, nord-1 do wave(*,i) = mlam_vec / olist(i)
    endelse
  endelse

;List of valid keyboard commands.
  valid_cmd = ['a','h','f','m','o','p','q','r','t','w','x']

;Save current plot state.
  xsav = !x
  ysav = !y
  psav = !p

;Set desired plot characteristics.
  !p.charsize=1.5
  !x.margin=[7,2]
  !y.margin=[4,1]

;Set up the color table.
  erase
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

;Begin command loop (exit when user types "x").
  cmd = 'r'
  while cmd ne 'x' and cmd ne 'q' do begin

;Handle pending plot refresh commands.
    if cmd eq 'r' then begin
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
;          th_pix = interpol(findgen(npix), congrid(wave(*,i),npix), th_wav)
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
    endif

;Get next command.
    repeat begin
      print, form='(a,$)', 'cmd: '
      cmd = get_kbrd(1)
      icmd = where(valid_cmd eq cmd, valid)
      if valid then begin
        print, cmd
      endif else begin
        print, cmd + " (invalid, type 'h' for list of commands)"
      endelse
    endrep until valid

;Help handler.
    if cmd eq 'h' then begin
      print, '  a - automatically mark lines'
      print, '  h - print command summary'
      print, '  f - fit marked lines to get wavelength solution'
      print, '  m - mark/unmark lines'
      print, '  o - remove outliers'
      print, '  p - print list of marked lines'
      print, '  r - refresh plot'
      print, '  q - quit program'
      print, '  t - set display threshold'
      print, '  w - set line location window'
      print, '  x - exit program'
    endif

;Display threshold handler.
    if cmd eq 't' then begin
      read, 'Display threshold is currently ' + string(dmax, '(f10.2)') $
          + ', enter new value: ', dmax
      dmax = dmax > 1				;should be at least 1
      cmd = 'r'					;redisplay spectrum
    endif

;Line location window handler.
    if cmd eq 'w' then begin
      read, ' line window is currently +/-' + string(mwin, '(f10.2)') $
          + ' pixels, enter new value: ', mwin
      mwin = mwin > 2				;must be at least 2 pixels
    endif

;Line marking handler.
    if cmd eq 'm' then begin
      done = 0
      repeat begin

;Get atlas line.
        print, form='(a,$)', 'Mark an isolated atlas line: '
        cursor, cx1, cy1, /DOWN
        if cx1 lt min(!x.crange) or cx1 gt max(!x.crange) then done = 1
        if cy1 lt min(!y.crange) or cy1 gt max(!y.crange) then done = 1
        if done eq 0 then begin
          cord = round(cy1)
          junk = min(abs(olist - cord), iord)
          cwav1 = wave(round(cx1), iord)
          junk = min(abs(th_wair - cwav1), iwav1)
          xpix = (interpol(findgen(npix), wave(*,iord), th_wair(iwav1)))(0)
          lstr = (alog10(th_inten(iwav1)) - thmin(iord)) $
               / (thmax(iord) - thmin(iord))
          lstr = 0.8*lstr + 0.1				;map into 0.1-0.9
          oplot, xpix+[0,0], cord+[-0.5,0.5]*lstr $
               , co=colorg, th=1+(lstr gt 0.5)
          print, 'lcen=' + string(cwav1, '(f10.4)') $
               + ', order=' + strtrim(cord, 2)

;Get observed line.
;          wait, 0.2
          print, form='(a,$)', 'Mark corresponding observed line: '
          cursor, cx2, cy2, /DOWN
          bad = 0
          if cx2 lt min(!x.crange) or cx2 gt max(!x.crange) then bad = 1
          if cy2 lt min(!y.crange) or cy2 gt max(!y.crange) then bad = 1
          if round(cy2) ne cord then bad = 1
          if bad eq 0 then begin
            junk = max(s(cx2-mwin:cx2+mwin,iord), imin)
            ipk = round(cx2) - mwin + imin
            xfit = ipk + dindgen(2*mwin+1) - mwin
            yfit = s(ipk-mwin:ipk+mwin,iord)
            gfit = gaussfit(xfit, yfit, gpar, nterm=4)
            ipix2 = gpar(1)
            oplot, ipix2+[0,0], cord+[-0.5,0.5], co=colorg
            print, 'pixel=' + string(ipix2, '(f10.2)')
          endif
          if bad eq 1 then begin
            print, '(line not marked)'
            oplot, xpix+[0,0], cord+[-0.45,0.45]*lstr $
                 , co=colblu, th=1+(lstr gt 0.5)
          endif else begin
            m_wave(nm) = th_wair(iwav1)
            m_ord(nm) = cord
            m_pix(nm) = ipix2
            m_wid(nm) = gpar(2)
            m_flag(nm) = 1
            nm = nm + 1
          endelse
          wait, 0.2
        endif
      endrep until done eq 1
      print, '(done marking lines)'
    endif

;Print marked lines handler.
    if cmd eq 'p' then begin
      if nm gt 0 then begin
        buff = strarr(nm)
        isort = sort(m_wave(0:nm-1))
        m_wave(0:nm-1) = m_wave(isort)
        m_ord(0:nm-1) = m_ord(isort)
        m_pix(0:nm-1) = m_pix(isort)
        m_wid(0:nm-1) = m_wid(isort)
        m_flag(0:nm-1) = m_flag(isort)
        for i=0, nm-1 do begin
          buff(i) = string(i, form='(i4)') + ': ' $
                  + strtrim(string(m_wave(i), form='(f10.4)'), 2) + ',' $
                  + string(m_ord(i), form='(i5)') + ',' $
                  + string(m_pix(i), form='(f9.2)')
        endfor
;        hprint, buff
        print, buff
      endif else begin
        print, 'No marked lines'
      endelse
    endif

;Automatic line association handler.
    if cmd eq 'a' then begin
      if nm eq 0 then estsig = 2.5/2.35 else estsig = total(m_wid) / nm
      h = histogram(s, min=-20, max=20)
      xh = -20 + findgen(n_elements(h))
      junk = gaussfit(xh, h, gpar, nterm=3)
      minamp = 5.0 * gpar(2)
      buff = ''
      awin = round(5.0 * estsig) > 2
      read, 'Enter half width of fit window in pixels (5*sig=' $
          + strtrim(awin, 2) + '):  ', buff
      if strtrim(buff, 2) ne '' then awin = round(float(buff))
      awin = awin > 2
      maxres = 0.5 * gpar(2)
      read, 'Enter maximum allowed residual in pixels (sig/2=' $
          + strtrim(string(maxres, form='(f9.2)'),2) + '):  ', buff
      if strtrim(buff, 2) ne '' then maxres = float(buff)
      read, 'Enter minimum line height to fit in ADU (5*sig=' $
          + strtrim(string(minamp, form='(f9.1)'),2) + '):  ', buff
      if strtrim(buff, 2) ne '' then minamp = round(float(buff))
      nm = 0
      for i=0, nord-1 do begin
        onum = olist(i)
        wbeg = min(wave(*,i), max=wend)
        iwhr = where(th_wair ge wbeg and th_wair le wend, nwhr)
        for j=0, nwhr-1 do begin
          th_wav = th_wair(iwhr(j))
          th_pix = (interpol(findgen(npix), wave(*,i), th_wav))(0)
;          th_pix = (interpol(findgen(npix), congrid(wave(*,i),npix), th_wav))(0)
          ipix = round(th_pix)
          if ipix gt awin and ipix lt (npix-awin-1) then begin
            xfit = ipix + dindgen(2*awin+1) - awin
            yfit = s(ipix-awin:ipix+awin,i)
            bg = min(yfit)
            amp = max(yfit) - bg
            if amp gt minamp then begin
              gpar = double([amp, ipix, estsig, bg])
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
      cmd = 'r'
    endif

;Fit marked lines handler.
    if cmd eq 'f' then begin
      if nm eq 0 then begin
        print, 'Mark lines before fitting wavelengths'
      endif else if nm eq 1 then begin
        junk = min(abs(olist - m_ord(0)), iord)
        wcur = (interpol(wave(*,iord), findgen(npix), m_pix(0)))(0)
        dmlam = m_ord(0) * (m_wave(0) - wcur)
        for i=0, nord-1 do wave(*,i) = wave(*,i) + dmlam / olist(i)
        cmd = 'r'
      endif else begin
        m_w = double(m_wave(0:nm-1))
        m_p = double(m_pix(0:nm-1))
        m_o = double(m_ord(0:nm-1))
        m_f = double(m_flag(0:nm-1))
        igd = where(m_f ne 0, ngd)
        m_w = m_w(igd)
        m_p = m_p(igd)
        m_o = m_o(igd)
        m_f = m_f(igd)
        m_ml = m_o * m_w

;Determine what order polynomials can be constrained.
        ordxok = intarr(maxordx+1)
        for i=0, n_elements(ordxok)-1 do begin
          bsiz = npix / float(i+1)
          h = histogram(m_p, min=0, max=npix-1, bin=bsiz)
          if min(h) gt 0 then ordxok(i) = 1
        endfor
        nordxok = total(ordxok)
        ordyok = intarr(maxordy+1)
        for i=0, n_elements(ordyok)-1 do begin
          bsiz = nord / float(i+1)
          h = histogram(m_o, min=obase, max=obase+nord-1, bin=bsiz)
          if min(h) gt 0 then ordyok(i) = 1
        endfor
        nordyok = total(ordyok)

;Choose polynomial order in X.
        if nordxok eq 1 then begin
          iwhr = where(ordxok eq 1)
          ordx = iwhr(0)
          print, 'Polynomial order in X: ', strtrim(ordx, 2)
        endif else begin
          list = ''
          for i=0, nordxok-1 do begin
            if ordxok(i) eq 1 then begin
              list = list + ',' + strtrim(i, 2)
            endif
          endfor
          print, form='(a,$)', 'Enter polynomial order in X {' $
                             + strmid(list, 1, 99) + '}:'
          read, '  ', ordx
        endelse

;Choose polynomial order in Y.
        if nordyok eq 1 then begin
          iwhr = where(ordyok eq 1)
          ordy = iwhr(0)
          print, 'Polynomial order in Y: ', strtrim(ordy, 2)
        endif else begin
          list = ''
          for i=0, nordyok-1 do begin
            if ordyok(i) eq 1 then begin
              list = list + ',' + strtrim(i, 2)
            endif
          endfor
          print, form='(a,$)', 'Enter polynomial order in Y {' $
                             + strmid(list, 1, 99) + '}:'
          read, '  ', ordy
        endelse

;Decide on cross-terms (use all that are OK).
        xok = intarr(maxcross)
        xlist = ''
        nxok = 0
        if ordx ge 1 and ordy ge 1 and ngd ge ordx+ordy+nxok+2 then begin
          xok(0) = 1
          xlist = xlist + 'XY '
          nxok = nxok + 1
        endif
        if ordx ge 2 and ordy ge 1 and ngd ge ordx+ordy+nxok+2 then begin
          xok(1) = 1
          xlist = xlist + 'XXY '
          nxok = nxok + 1
        endif
        if ordx ge 1 and ordy ge 2 and ngd ge ordx+ordy+nxok+2 then begin
          xok(2) = 1
          xlist = xlist + 'XYY '
          nxok = nxok + 1
        endif
        if ordx ge 2 and ordy ge 2 and ngd ge ordx+ordy+nxok+2 then begin
          xok(3) = 1
          xlist = xlist + 'XXYY '
          nxok = nxok + 1
        endif
        if ordx ge 3 and ordy ge 1 and ngd ge ordx+ordy+nxok+2 then begin
          xok(4) = 1
          xlist = xlist + 'XXXY '
          nxok = nxok + 1
        endif
        if ordx ge 1 and ordy ge 3 and ngd ge ordx+ordy+nxok+2 then begin
          xok(5) = 1
          xlist = xlist + 'XYYY '
          nxok = nxok + 1
        endif
        nxok = total(xok)
        if nxok gt 0 then print, 'Cross terms: ' + xlist

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
        mlfit = marq('thid_func', dummy, m_ml, sig, par, dpar, /trace)

;Calculate wavelengths for every pixel.
        dummy = fltarr(npix)
        thid_pixel = dindgen(npix) / 1d3

; AT 5-Oct-2011
;       wave = fltarr(npix,nord) ; re-define wave
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
        print, 'Residuals:  Min=' $
             + strtrim(string(min(resid), form='(f10.2)'), 2) $
             + ',  Max=' + strtrim(string(max(resid), form='(f10.2)'), 2) $
             + ',  RMS=' + strtrim(string(stddev(resid), form='(f10.2)'), 2) $
             + ' pixels'
        cmd = 'r'
      endelse
    endif

;Outlier handler.
    if cmd eq 'o' then begin
      if not keyword_set(resid) then begin
        print, 'Must fit marked lines before removing outliers'
      endif else begin
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

        good = (m_f gt 0)
        print, 'Draw box containing outliers to reject'
        done = 0
        repeat begin
          igd = where(good eq 1)
          plot, m_p(igd), resid(igd), /ps, xsty=3, ysty=3, chars=1.5 $
              , xtit='Pixel Number', ytit='Residual (pixels)'
          cursor, cx1, cy1, wait=3
          if cx1 lt min(!x.crange) or cx1 gt max(!x.crange) then done = 1
          if cy1 lt min(!y.crange) or cy1 gt max(!y.crange) then done = 1
          if done eq 0 then begin
            cx2 = cx1
            cy2 = cy1
            repeat begin
              cx2o = cx2
              cy2o = cy2
              cursor, cx2, cy2, wait=2
              oplot, [cx1, cx1, cx2o, cx2o, cx1] $
                   , [cy1, cy2o, cy2o, cy1, cy1], co=!p.background
              oplot, m_p(igd), resid(igd), /ps
              oplot, [cx1, cx1, cx2, cx2, cx1] $
                   , [cy1, cy2, cy2, cy1, cy1], co=!p.color
              cursor, junk, junk, /nowait
            endrep until !err eq 0
            ibad = where(m_p gt (cx1<cx2) and m_p lt (cx1>cx2) $
                     and resid gt (cy1<cy2) and resid lt (cy1>cy2), nbad)
            if nbad gt 0 then good(ibad) = 0
          endif
        endrep until done eq 1

        print, 'Draw box containing outliers to reject'
        done = 0
        wait, 0.2
        repeat begin
          igd = where(good eq 1, ngd)
          plot, m_o(igd), resid(igd), /ps, xsty=3, ysty=3, chars=1.5 $
              , xtit='Order Number', ytit='Residual (pixels)'
          medres = replicate(1e10, nord)
          for i=0,nord-1 do begin
            iwhr = where(m_o(igd) eq olist(i), nwhr)
            if nwhr gt 0 then medres(i) = median(resid(igd(iwhr)))
          endfor
          oplot, olist, medres, max=1e9, co=colblu
          cursor, cx1, cy1, wait=3
          if cx1 lt min(!x.crange) or cx1 gt max(!x.crange) then done = 1
          if cy1 lt min(!y.crange) or cy1 gt max(!y.crange) then done = 1
          if done eq 0 then begin
            cx2 = cx1
            cy2 = cy1
            repeat begin
              cx2o = cx2
              cy2o = cy2
              cursor, cx2, cy2, wait=2
              oplot, [cx1, cx1, cx2o, cx2o, cx1] $
                   , [cy1, cy2o, cy2o, cy1, cy1], co=!p.background
              oplot, m_o(igd), resid(igd), /ps
              oplot, [cx1, cx1, cx2, cx2, cx1] $
                   , [cy1, cy2, cy2, cy1, cy1], co=!p.color
              cursor, junk, junk, /nowait
            endrep until !err eq 0
            ibad = where(m_o gt (cx1<cx2) and m_o lt (cx1>cx2) $
                     and resid gt (cy1<cy2) and resid lt (cy1>cy2), nbad)
            if nbad gt 0 then good(ibad) = 0
          endif
        endrep until done eq 1
        m_flag(igd) = 2
        ibad = where(good eq 0, nbad)
        if nbad gt 0 then m_flag(ibad) = 0
        print, strtrim(ngd, 2) + ' good lines'
        cmd = 'r'
      endelse
    endif

;End of command loop.
  endwhile

;Restore plot state.
  !p = psav
  !y = ysav
  !x = xsav

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
