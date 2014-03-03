function biasTrim, inBias
	return, [ $
		long(inBias[0]) + 3, $ ; col start
		long(inBias[1]) - 3, $ ; col end
		long(inBias[2]), $ ; row start
		long(inBias[3]) $ ; row end
	]
END

function chip_geometry, fileName, hdr=hdr, redpar=redpar
	; given a CTIO fits file, return a structure which contains the following:
	;	bias region (full according to the controller)
	;	bias region (trimmed by a few column pixels in beginning and end)
	;	image region (full according to the controller)
	;	image region (trimmed by 10 pixels...in the header already)
	; 20110330 - Modified to output the readout_speed to correct for
	;		non-linearity ~MJG
	; 20120218 - Added torrent
	; 20130828 - added ccd_full ~MJG
	;
	
	; If only a fileName was passed, open it and grab the header.
	;	If the header was passed, there is no need to get the info ourselves.
	;
	if ( ~keyword_set(hdr) ) then begin
		f = readfits(fileName,hdr,/SILENT)
	endif
	if keyword_set(redpar) then begin
	  date = redpar.date
	endif else begin
	  dkw = strt(sxpar(hdr,'date'))
	  date = strmid(dkw,2,2)+strmid(dkw,5,2)+strmid(dkw,8,2)
	endelse
	
	namps = n_elements(strsplit(sxpar(hdr, 'amplist'), ' '))
	; Set up the results structure
	;	We could have more amp readouts than just left and right,
	;	for the moment that would constitute an error though so
	;	we will just leave upleft and upright for now.
	
	if namps eq 2 then begin
	results = { 				$
		status: 'error', 		$
		controller: 'unknown',	$
		ccd: 'unknown',    $
		amps: ['upleft','upright'], $
		bias_trim: {upleft: [0,0,0,0], upright: [0,0,0,0]},		$
		bias_full: {upleft: [0,0,0,0], upright: [0,0,0,0]},		$
		image_trim: {upleft: [0,0,0,0], upright: [0,0,0,0]},	$
		image_full: {upleft: [0,0,0,0], upright: [0,0,0,0]},	$
		gain:	{upleft: 0.0, upright: 0.0},					$
		read_noise: {upleft: 0.0, upright: 0.0}, 				$
		bin: { row: 1,	col: 1 },								$
		readout_speed: 'unknown'    $
	}
	endif
	
	if namps eq 4 then begin
	results = { 				$
		status: 'error', 		$
		controller: 'unknown',	$
		ccd: 'unknown',    $
		amps: ['upleft','upright', 'botleft', 'botright'], $
		bias_trim: {upleft: [0,0,0,0], upright: [0,0,0,0], botleft: [0,0,0,0], botright: [0,0,0,0]},		$
		bias_full: {upleft: [0,0,0,0], upright: [0,0,0,0], botleft: [0,0,0,0], botright: [0,0,0,0]},		$
		image_trim: {upleft: [0,0,0,0], upright: [0,0,0,0], botleft: [0,0,0,0], botright: [0,0,0,0]},	$
		image_full: {upleft: [0,0,0,0], upright: [0,0,0,0], botleft: [0,0,0,0], botright: [0,0,0,0]},	$
		ccd_full: {upleft: [0,0,0,0], upright: [0,0,0,0], botleft: [0,0,0,0], botright: [0,0,0,0]},	$
		gain:	{upleft: 0.0, upright: 0.0, botleft: 0.0, botright: 0.0},					$
		read_noise: {upleft: 0.0, upright: 0.0, botleft: 0.0, botright: 0.0}, 				$
		bin: { row: 1,	col: 1 },								$
		readout_speed: 'unknown'    $
	}
	endif
	

	;
	amps = { upleft: '21', upright: '22' } ;, lowleft: '11', lowright: '12' }
	validAmpList = '21 22' ; we could be off if we use a different arrangement of amps
	
	; first determine if we are using the controller (1/1/2010) or the old
	;	we will just pick keywords unique to old and new
	;	(old 'GTINDEX', new 'SPEEDMOD')
	XSTART = sxpar(hdr,'XSTART',count=early2008)
	GTINDEX = sxpar(hdr,'GTINDEX',count=oldCount)
	SPEEDMOD = sxpar(hdr,'SPEEDMOD',count=newCount)
	DETECTOR = sxpar(hdr,'DETECTOR',count=det_tag)
	DHEINF = sxpar(hdr,'DHEINF',count=control_tag)

	if ((strpos(DHEINF, 'MNSN') ge 0) and (strpos(DHEINF, 'torrent') lt 0)) then controller='mnsn' else controller='unknown' 
	if strpos(DHEINF, 'torrent') ge 0 then controller = 'torrent'
	if strpos(DETECTOR,'4K') ge 0 then ccd = 'E2V-4K'
	
	if ( early2008 eq 1 ) then begin
		; early 2008 format of headers, chip layout same as pre-monsoon controller
		results.controller = 'arcon'
		results.ccd = '200 (Blanco 2k, 24micron)'
		keys = { $
			amplist: 'AMPLIST',		$
			bias: 'BSEC',			$
			image_full: 'DSEC',		$
			image_trim: 'TSEC',		$
			gain: 'GTGAIN',			$
			read_noise: 'GTRON',	$
			bin:	 'CCDSUM'		$
		}
	endif else if ( oldCount eq 1 AND newCount eq 0 and ccd ne 'E2V-4K') then begin
		; using the old controller (pre-monsoon)
		results.controller = 'arcon'
		results.ccd = '200 (Blanco 2k, 24micron)'
		keys = { $
			amplist: 'AMPLIST',		$
			bias: 'BSEC',			$
			image_full: 'DSEC',		$
			image_trim: 'TSEC',		$
			gain: 'GTGAIN',			$
			read_noise: 'GTRON', 	$
			bin: 'CCDSUM'			$
		}
	endif else if ( oldCount eq 0 AND newCount eq 1 and ccd ne 'E2V-4K') then begin
		; using the 'new' monsoon controller
		results.controller = 'new'
		results.ccd = '200 (Blanco 2k, 24micron, mnsn)'
		keys = { $
			amplist: 'AMPLIST',		$
			bias: 'BSEC',			$
			image_full: 'DSEC',		$
			image_trim:	'TSEC',		$
			gain: 'GAIN',			$
			read_noise: 'RON',		$
			bin: 'CCDSUM'			$
		}	
	endif else if ( ccd eq 'E2V-4K' AND controller eq 'mnsn') then begin
		; using the MNSN controller with the 4K e2v detector
		results.controller = 'mnsn'
		results.ccd = '201 (e2v 4k, 15micron)'
		keys = { $
			amplist: 'AMPLIST',		$
			bias: 'BSEC',			$
			image_full: 'DSEC',		$
			image_trim:	'TSEC',		$
			gain: 'GAIN',			$
			read_noise: 'RON',		$
			bin: 'CCDSUM'			$
		}
		results.readout_speed = SPEEDMOD
	endif else if (ccd eq 'E2V-4K' AND controller eq 'torrent') then begin
		; using the torrent controller with the 4K e2v detector
		results.controller = 'torrent'
		results.ccd = '201 (e2v 4k, 15micron)'
		keys = { $
			amplist: 'AMPLIST',		$
			bias: 'BSEC',			$
			ccd_full: 'SCSEC',      $
			image_full: 'DSEC',		$
			image_trim:	'TSEC',		$
			gain: 'GAIN',			$
			read_noise: 'RON',		$
			bin: 'CCDSUM'			$
		}
		results.readout_speed = SPEEDMOD
		amps = { upleft: '21', upright: '22' , botleft: '11', botright: '12' }
		validAmpList = '11 12 21 22' ; we could be off if we use a different arrangement of amps
	endif else begin
		; don't know what we are using, bail
		return, results
	endelse

	;
	; Next, check to see if we have the correct amps.  We should have
	; 21 (upper left) and 22 (upper right).  Anything else and we should
	; exit out.
	;
	usedAmps = sxpar(hdr,keys.amplist,count=hasAmpList)
	if ( hasAmpList ne 1 OR STRTRIM(usedAmps,2) ne validAmpList ) then return, results
	
	;
	; Now get the info we are looking for
	;	pixel values will be in the format [1045:1098,1:1200]
	;	by splitting on : and , we will get col_beg,col_end,row_beg,row_end
	;	we also include [ and ] because sometimes the parameters have them
	;	as if they were an IDL array and it messed with the extraction.
	;
	;	All pixel values are 1-based but the IDL arrays are all 0 based,
	;	so we will subtract 1 from every pixel value.
	;
	
	; bias
	results.bias_full.upleft = long(strsplit(sxpar(hdr,keys.bias + amps.upleft),'[],:',/EXTRACT)) -1.
	results.bias_full.upright = long(strsplit(sxpar(hdr,keys.bias + amps.upright),'[],:',/EXTRACT)) -1.
	if namps eq 4 then begin
	results.bias_full.botleft = long(strsplit(sxpar(hdr,keys.bias + amps.botleft),'[],:',/EXTRACT)) -1.
	results.bias_full.botright = long(strsplit(sxpar(hdr,keys.bias + amps.botright),'[],:',/EXTRACT)) -1.
	endif
	
	; bias trim: trim off 3 pixels on both sides of actual bias
	results.bias_trim.upleft = biasTrim(results.bias_full.upleft)
	results.bias_trim.upright = biasTrim(results.bias_full.upright)
	if namps eq 4 then begin
	results.bias_trim.botleft = biasTrim(results.bias_full.botleft)
	results.bias_trim.botright = biasTrim(results.bias_full.botright)
	endif
	
	; full ccd area
	if namps eq 4 then begin
	results.ccd_full.upleft = long(strsplit(sxpar(hdr,keys.ccd_full+amps.upleft),'[],:',/EXTRACT)) -1.
	results.ccd_full.upright = long(strsplit(sxpar(hdr,keys.ccd_full+amps.upright),'[],:',/EXTRACT)) -1.
	results.ccd_full.botleft = long(strsplit(sxpar(hdr,keys.ccd_full+amps.botleft),'[],:',/EXTRACT)) -1.
	results.ccd_full.botright = long(strsplit(sxpar(hdr,keys.ccd_full+amps.botright),'[],:',/EXTRACT)) -1.
	endif

	; full image area
	results.image_full.upleft = long(strsplit(sxpar(hdr,keys.image_full+amps.upleft),'[],:',/EXTRACT)) -1.
	results.image_full.upright = long(strsplit(sxpar(hdr,keys.image_full+amps.upright),'[],:',/EXTRACT)) -1.
	if namps eq 4 then begin
	results.image_full.botleft = long(strsplit(sxpar(hdr,keys.image_full+amps.botleft),'[],:',/EXTRACT)) -1.
	results.image_full.botright = long(strsplit(sxpar(hdr,keys.image_full+amps.botright),'[],:',/EXTRACT)) -1.
	endif
	
	; trimmed image area (useable image)
	results.image_trim.upleft = long(strsplit(sxpar(hdr,keys.image_trim+amps.upleft),'[],:',/EXTRACT)) -1.
	results.image_trim.upright = long(strsplit(sxpar(hdr,keys.image_trim+amps.upright),'[],:',/EXTRACT)) -1.
	if namps eq 4 then begin
	results.image_trim.botleft = long(strsplit(sxpar(hdr,keys.image_trim+amps.botleft),'[],:',/EXTRACT)) -1.
	results.image_trim.botright = long(strsplit(sxpar(hdr,keys.image_trim+amps.botright),'[],:',/EXTRACT)) -1.
	endif
   if long(date) lt 120502L then begin
	  ;fix binning problem that was fixed on 120502. See emails to/from
	  ;Marco on that date or CHIRON Blog
	  results.image_trim.upleft[1] += 1
	  results.image_trim.botleft[1] += 1
	  results.image_trim.upright[0] -= 1
	  results.image_trim.botright[0] -= 1
	endif
	
	; gain
	results.gain.upleft = sxpar(hdr,keys.gain+amps.upleft)
	results.gain.upright = sxpar(hdr,keys.gain+amps.upright)
	if namps eq 4 then begin
	results.gain.botleft = sxpar(hdr,keys.gain+amps.botleft)
	results.gain.botright = sxpar(hdr,keys.gain+amps.botright)
	endif
	;
	; GAIN OVERRIDE
	;	Gain test have derived the following values which currently do not correspond to the
	;	values supplied in the headers.
	;
	if ( results.controller eq 'old' ) then begin
		results.gain.upleft = 1.38
		results.gain.upright = 1.23
	endif else if ( results.controller eq 'new' ) then begin
		results.gain.upleft = 2.12
		results.gain.upright = 1.96
	endif else if ( results.controller eq 'mnsn' ) then begin
		results.gain.upleft = 2.5
		results.gain.upright = 2.6
	endif
	if ( results.controller eq 'torrent' ) then begin
		results.gain.upleft = 1.30
		results.gain.upright = 1.250
		results.gain.botleft = 1.264
		results.gain.botright = 1.257
	endif
	if ( results.controller eq 'new' ) then begin
		; At the current time (March 1, 2010) the image trim region for the upper left amp is off by 1 pixel
		;  as if the offsets are given by a zero-based system rather than one based.  We add the pixel
		;  back here.  Revisit this if the header changes.
		
		; Now, another change...this has been corrected, so as long as we aren't starting
		;	on column 20, we will skip adding the offset.  We should really correct all
		;	of the headers and remove this kludge
		;
		if ( results.image_trim.upleft[0] eq 20 ) then results.image_trim.upleft[0:1] += 1
		;results.image_trim.upright[0:1] += 1
	endif
	
	; binning (should be 1x1)
	bin = sxpar(hdr,keys.bin)
	bin = stregex(bin,'([0-9]+)[ ]+([0-9]+)',/EXTRACT,/SUBEXPR)
	results.bin.row = bin[1]
	results.bin.col = bin[2]
	
	; read noise
	results.read_noise.upleft = sxpar(hdr,keys.read_noise+amps.upleft)
	results.read_noise.upright = sxpar(hdr,keys.read_noise+amps.upleft)
	if namps eq 4 then begin
	results.read_noise.botleft = sxpar(hdr,keys.read_noise+amps.botleft)
	results.read_noise.botright = sxpar(hdr,keys.read_noise+amps.botleft)
	endif
	
	results.status = 'OK'
	return, results
END

