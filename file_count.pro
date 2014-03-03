pro file_count, dir=dir, prefix, missing
	; Get all of the fits files in the specified directory,
	; look for any gaps between the sequentially numbered files.
	; return a list of the files which are missing.
	;
	print,'Your logsheet has been made with files from '+dir
	if ( strmid(dir,0,/REVERSE_OFFSET) ne '/' ) then dir = dir + '/'
	
	; Get all of the fits files in the directory
	;
	allFitsFiles=File_Search(dir+prefix+'*.fits',count=nFiles)
	if nFiles eq 0 then stop, 'no fits files found for ' + dir+prefix + '!'
	
	; make sure that the files we found are formatted like: qa04.nnnn.fits or qa04_nnnn.fits
    ;
    obsFiles = where(stregex(allFitsFiles,'/([0-9a-zA-Z]+).([0-9]+)\.fits$',/BOOLEAN))
    if ( n_elements(obsFiles) eq 0 ) then begin
    	print, 'No files found with name format "/([0-9a-zA-Z]+).([0-9]+)\.fits$"'
    	stop
    endif
    obsFiles = allFitsFiles[obsFiles]
    nObs = n_elements(obsFiles)

	obsNum=intarr(nObs)
	
	; Collect all of the sequence numbers
	;
	for i=0,nObs-1 do begin
	   f = stregex(obsFiles[i],"([0-9]+)\.fits$",/EXTRACT,/SUBEXPR)
	   obsNum[i] = fix(f[1])
	endfor
	
	; Get the range of sequence numbers and create an array without any holes
	; starting at the minimum sequence number that was found.
	;
	totalSeqNums=max(obsNum)-min(obsNum)+1
	if ( totalSeqNums eq nObs ) then begin
		; nothing to do, we have the proper number of files in the range min->max
		missing = -1
		return
	endif
	
	seqNums=indgen(totalSeqNums)+min(obsNum)
	
	; locate the missing sequence numbers
	;
	nMissing = totalSeqNums - nObs
	missing=intarr(nMissing)
	for i=0,totalSeqNums-1 do begin
	   x=where(obsNum eq seqNums[i],nx)
	   if nx eq 0 then missing=[missing,seqNums[i]]
	endfor
	
	print,''
	print,'You have '+nMissing+' missing files: '
	print,''
	print,'The following files are missing from the subdirectory: ',dir
	print, ''
	print, prefix + '.' + missing + '.fits'
	print, ''
	print, 'Please make a note of the missing files in the logsheet'
	print, ''
end
