;+
;
;  NAME: 
;     yddf_logmaker
;
;  PURPOSE: Given the date, this routine will generate a logsheet
;
;  CATEGORY:
;      CHIRON
;
;  CALLING SEQUENCE:
;
;      yddf_logmaker, date
;
;  INPUTS:
;
;  OPTIONAL INPUTS:
;
;  OUTPUTS:
;
;  OPTIONAL OUTPUTS:
;
;  KEYWORD PARAMETERS:
;    
;  EXAMPLE:
;      yddf_logmaker
;
;  MODIFICATION HISTORY:
;        c. Matt Giguere 2014.03.14 17:06:41
;
;-
pro yddf_logmaker, date

;make string and/or chop excess whitespace for input date:
date = strt(date)

;make a list of all the files in the raw directory:
rawdir = '/raw/yddf/'+date+'/'
spawn, 'ls -1 '+rawdir, filelist

printt, filelist
;the logsheet filename and directory:
lgdir = '/tous/yddf/logsheets/20'+strmid(date, 0, 2)+'/'

;make the directory if it doesn't exist for the current year:
if ~file_test(lgdir) then spawn, 'mkdir '+lgdir
lgfn = lgdir+'yddf'+date+'.log'

;open the logsheet for writing:
openw, wlun, lgfn, /get_lun

;add the header information
printf, wlun, ' YALE DOPPLER DIAGNOSTIC FACILITY (YDDF) EXPOSURE LOG '
printf, wlun, ''
printf, wlun, '-------------------------------------------------------'+$
	'-----------------------------------'
printf, wlun, 'DATE: ', date
printf, wlun, 'Number of files: ', strt(n_elements(filelist))
printf, wlun, '-------------------------------------------------------'+$
	'-----------------------------------'
printf, wlun, '  OBS #   |  OBJECT NAME       |    DATE & TIME    | EXPTIME | BIN 

;now loop through all the files, writing the contents to a logsheet:
for i=0, n_elements(filelist)-1 do begin
hd = headfits(rawdir+filelist[i])
fnarr = strsplit(filelist[i], '.', /extract)
printf, wlun, $
	string(fnarr[1], format='(I-10)'), '|', $;the sequence number
	string(fxpar(hd, 'object'), format='(A20)'), '|', $ ;object name
	string(fxpar(hd, 'DATE-OBS'), format='(A19)'), '|', $ ;observation date
	string(fxpar(hd, 'EXPTIME'), format='(A9)'), '|', $ ;exposure time
	' ',string(strt(fxpar(hd, 'XBINNING')), f='(A1)'), 'x', $
	string(strt(fxpar(hd, 'YBINNING')), f='(A1)') ;binning

endfor
close, wlun
print, 'Your logsheet has been created. It can be found at: '
print, lgfn
end;yddf_logmaker.pro