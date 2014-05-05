;+
;
;  NAME: 
;     printt
;
;  PURPOSE: Print the transpose of an array
;   
;
;  CATEGORY:
;      UTILITIES
;
;  CALLING SEQUENCE:
;
;      printt
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
;      printt
;
;  MODIFICATION HISTORY:
;        c. Matt Giguere 2012.03.10 12:20:43 PM
;
;-
pro printt, $
help = help, $
arg

if keyword_set(help) then begin
	print, '*************************************************'
	print, '*************************************************'
	print, '        HELP INFORMATION FOR printt'
	print, 'KEYWORDS: '
	print, ''
	print, 'HELP: Use this keyword to print all available arguments'
	print, ''
	print, ''
	print, ''
	print, '*************************************************'
	print, '                     EXAMPLE                     '
	print, "IDL>"
	print, 'IDL> '
	print, '*************************************************'
	stop
endif

print, transpose(arg)
end;printt.pro