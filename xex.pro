; CHIRON data browsing GUI
; AT Oct 20, 2011
pro xex_event, ev

common xexcomm, xbase,redpar,mode,logs,stars,thars,nightidx,staridx,tharidx,ordidx,log,sp,w,flags

  type = tag_names(ev, /structure) ; name of the {ev} structure

  case type of
;--------------------
'WIDGET_BUTTON': begin
 widget_control,ev.id, get_value=value 
 case value of
   'Read': read_sp ; print, 'Read'
   'Plot': plot_sp ; print, 'Plot'
   '2Dspec': plot_2d ; print, 'Plot'
   'Exit': widget_control, ev.top, /destroy
   '>':   nextorder, 1   
   '<':   nextorder, -1
 else:
 endcase ; buttons
 end      ; button

'WIDGET_COMBOBOX': begin
 case ev.id of
  widget_info(xbase,find_by_uname='night_sel'):  begin
      nightidx = ev.index
      print, 'Selected night: '+logs[nightidx]
      getlog
      getstars
   end  ; night_sel  
  widget_info(xbase,find_by_uname='mode_sel'):   begin
     mode = redpar.modes[ev.index]
     redpar.mode = ev.index
     print, 'current mode: ', mode
     getstars
  end ; mode_sel
  widget_info(xbase,find_by_uname='star_sel'): begin
    staridx = ev.index
    showlog
    read_sp
;    print, 'star=',staridx
;    stop
 end ; star_sel
 widget_info(xbase,find_by_uname='order_sel'):  begin
    ordidx = ev.index
    plot_sp
 end
 else:
  endcase
;  showdat
end       ; combobox

;--------------------------- Edit events
else: begin
 widget_control,ev.id, get_uvalue=uvalue 
 if (uvalue eq 'Paredit') then begin
  widget_control,ev.id, get_value=value 
 case ev.id of

  widget_info(xbase,find_by_uname='order_label'):  ordidx = fix(value)
;  widget_info(xbase,find_by_uname='type_sys'):  sysel.type = trim(value)
;  widget_info(xbase,find_by_uname='sep_sys'):  sysel.sep = float(value)
  else:
endcase
 endif
end
 
endcase 
;----------------------
 
end  ; xex_event
;-----------------------------------------------------------
pro xex

common xexcomm

print, 'Entering XEX GUI. Type >common xex to access the variabled from command line' 
print, 'common xexcomm, xbase,redpar,mode,logs,stars,thars,nightidx,staridx,tharidx,ordidx,log,sp,w,flags'

redpar = readpar('ctio.par') ; main parameters
logdir = redpar.rootdir+redpar.logdir
;logs = file_search(logdir+'*.log', count=nlogs)
logs = file_search(logdir+'[0-9]*.log', count=nlogs)

if nlogs eq 0 then return ; No logfiles
print, nlogs,' logsheets found in '+logdir
logs = strmid(logs, strlen(logdir)) ; remove the path
logs = strmid(logs,0,6) ; remove .log

ord = sort(logs)
logs = logs[reverse(ord)] ; reverse time order, last night first
print, 'Last night: '+logs[0]
nightidx = 0 ; current night
staridx = 0  ; current star 
ordidx = 28L  ; current order
redpar.mode= 1 ; initial mode: slit
redpar.nords = 40 ; set the order number

order_string = string(indgen(redpar.nords),'(I2)')

getlog  ; read the log
;getstars ; read the stars


  xbase = widget_base(title='CHIRON DATA GUI',/column)
  widget_control, xbase, /realize
  
  row1 = widget_base(xbase, /row) 
  row2 = widget_base(xbase, /row) 
  row3 = widget_base(xbase, /row) 
  row4 = widget_base(xbase, /row) 

; ------ top row ------
;  night_label =  widget_label(row1,uname='night_label',xsize=100,value=logs[nightidx])
  night_sel = widget_combobox(row1,uname='night_sel',xsize=80,value=logs)

  mode_sel = widget_combobox(row1,uname='mode_sel',xsize=80,value=redpar.modes)
  widget_control, mode_sel, set_combobox_select=redpar.mode ; pre-select the mode

  star_label =  widget_label(row1,uname='star_label',xsize=50,value='Star')
  star_sel = widget_combobox(row1,uname='star_sel',xsize=200,value=['none'])

  thar_label =  widget_label(row1,uname='thar_label',xsize=50,value='ThAr')
  thar_sel = widget_combobox(row1,uname='thar_sel',xsize=80,value=['none'])

; ------ 2nd row: log-file data
  logline_label = widget_label(row2,uname='logline_label',xsize=400,value='---')

; -------  3rd row: status data
  but0 = widget_button(row3, xsize=60, value='Read',sensitive=1)
  data_label = widget_label(row3,uname='data_label',xsize=200,value='---')

; --------- bottom row: controls

  but1 = widget_button(row4, xsize=60, value='<',sensitive=1)
  order_sel = widget_combobox(row4,uname='order_sel',xsize=80,value=order_string)
  widget_control, order_sel, set_combobox_select=ordidx
  but2 = widget_button(row4, xsize=60, value='>',sensitive=1)

  but3 = widget_button(row4, xsize=60, value='Plot',sensitive=1)
  but4 = widget_button(row4, xsize=60, value='2Dspec',sensitive=1)
  but5 = widget_button(row4, xsize=60, value='Exit',sensitive=1)

 xmanager, 'xex', xbase, /no_block
 getstars

end
;----------------------------------------------------------------------
pro getlog  ; read log-file for current night
common xexcomm
 logname = redpar.rootdir+redpar.logdir + logs[nightidx]+'.log'
 readcol, logname, skip=9, obnm, objnm, i2, mdpt, exptm, bin, slit,  f='(a5, a13, a4, a14, a8, a3, a6)'
 log = {obnm:obnm,objnm:objnm,i2:i2,mdpt:mdpt,exptm:exptm,bin:bin,slit:slit} 

; Find image prefix and remember it in redpar.prefix
  imdir = redpar.rootdir+redpar.rawdir+logs[nightidx]+'/'
  l = strlen(imdir)
  tmp = file_search(imdir+'*.fits', count = count) ; look for the data files
  if count eq 0 then begin
     print, 'No data files found for this night!'
     return
   endif
; check the QA prefix
;stop
  sel = where(strmid(tmp,l,2) eq 'qa')
  if n_elements(sel) gt 5 then run = strmid(tmp[sel[0]],l,4) else run = 'chi'+logs[nightidx]
   if strpos(run,'.') lt 0 then run=run+'.' ; add the point
   redpar.prefix = run
   print, 'Prefix = '+redpar.prefix
end
;----------------------------------------------------------------------
pro getstars   ; find stars for current night and mode
common xexcomm

  modeidx = redpar.mode ; current mode
;  xsl=where(log.bin eq redpar.binnings[modeidx] and log.slit eq redpar.modes[modeidx],n_modes) 
   xsl=where(log.slit eq redpar.modes[modeidx],n_modes) ; ignore binning

   objnm1=log.objnm[xsl]
   sel=where(objnm1 ne 'quartz' and objnm1 ne 'iodine' and objnm1 ne 'thar' $
			and objnm1 ne 'focus' and objnm1 ne 'junk' and objnm1 ne 'dark' $
                               and objnm1 ne 'bias', num_star)
   print, 'Stars found :',num_star

   if num_star gt 0 then stars = xsl[sel] else stars=[0]
   staridx = 0 ; first star for the night

   sel=where(objnm1 eq 'thar',num_thar)
   print, 'ThAr found :',num_thar
   if num_thar gt 0 then thars = xsl[sel] else thars=[0]
   tharidx = 0

  widget_control, widget_info(xbase,find_by_uname='star_label'), set_value='N**'+string(num_star,'(I3)')
  if num_star gt 0 then names = log.obnm[stars]+':'+log.objnm[stars] else names = ['none']
  widget_control, widget_info(xbase,find_by_uname='star_sel'), set_value=names

  widget_control, widget_info(xbase,find_by_uname='thar_label'), set_value='Th '+string(num_thar,'(I2)')
  if num_thar gt 0 then names = log.obnm[thars] else names = ['none']
  widget_control, widget_info(xbase,find_by_uname='thar_sel'), set_value=names

  showlog

end
;----------------------------------------------------------------------
pro nextorder, step
common xexcomm
 ordidx = (ordidx + step) 
 if ordidx ge redpar.nords then ordidx = redpar.nords -1
 ordidx = ordidx > 0
; widget_control, widget_info(xbase,find_by_uname='order_label'), set_value=string(ordidx,'(I2)')
 widget_control, widget_info(xbase,find_by_uname='order_sel'), set_combobox_select=ordidx

 plot_sp
; print, ordidx
end
;----------------------------------------------------------------------
pro showlog
common xexcomm

 i = stars[staridx] ; pointer to the star record 
; print, 'i= ',i
;obnm, objnm, i2, mdpt, exptm, bin, slit,  f='(a5, a13, a4, a14, a8, a3, a6)'
 logstring = log.obnm[i]+'   '+log.objnm[i]+'  '+log.i2[i]+'    '+log.mdpt[i] +'   '+ $
   log.exptm[i]+'    '+log.bin[i]+' '+log.slit[i] 
 widget_control, widget_info(xbase,find_by_uname='logline_label'), set_value=logstring
end
;----------------------------------------------------------------------
pro read_sp  ; read selected spectrum
common xexcomm

   filename = redpar.rootdir+redpar.fitsdir+'r'+redpar.prefix+log.obnm[stars[staridx]]+'.fits'
   tmp = file_search(filename, count=count)
   if count eq 0 then begin
      print, 'File '+filename+'  is not found!'
      return
   endif
   sp = readfits(filename)
  widget_control, widget_info(xbase,find_by_uname='data_label'), set_value=filename
;   sp = tmp[1,*,*]
;   w = tmp[0,*,*]
;stop
end
;----------------------------------------------------------------------
pro plot_sp  ; plot selected order
common xexcomm
; remove stupid last column!
ncol = n_elements(sp[0,*,0])
sp[1,ncol-1,*] = sp[1,ncol-2,*]

plot, sp[0,*,ordidx], sp[1,*,ordidx], xs=1
end
;----------------------------------------------------------------------
pro plot_2d  ; plot 2D spectrum
common xexcomm

sz = size(sp[0,*,*])
xs = sz[2] ; number of columns
while xs gt 1000 do xs=xs/2 ;reduce pixels
ys = sz[3] ; number of orders
;window, 0, xs=xs, ys=10*ys ; 8 pixels per order
tvscl, congrid(reform(sp[1,*,*]), xs, 10*ys)
;stop
end
;----------------------------------------------------------------------
;----------------------------------------------------------------------
;----------------------------------------------------------------------
;----------------------------------------------------------------------

