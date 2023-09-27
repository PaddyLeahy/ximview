; -----------------------------------------------------------------------------
;
;  Copyright (C) 2007-2013   J. P. Leahy
;
;
;  This file is part of Ximview and of HEALPix
;
;  Ximview and HEALPix are free software; you can redistribute them and/or modify
;  them under the terms of the GNU General Public License as published by
;  the Free Software Foundation; either version 2 of the License, or
;  (at your option) any later version.
;
;  Ximview and HEALPix are distributed in the hope that they will be useful,
;  but WITHOUT ANY WARRANTY; without even the implied warranty of
;  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
;  GNU General Public License for more details.
;
;  You should have received a copy of the GNU General Public License
;  along with HEALPix; if not, write to the Free Software
;  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
;
;
; -----------------------------------------------------------------------------
; Module PARSE_INPUT: decodes input parameter and converts to image
;
; J. P. Leahy 2008
;
; Contents:
;:
; parse_input_event  Exists solely to delete the optional plot.
; bb2nan             traps apparent bad bits
; bad2nan            Converts HP bad values to NaN
;                      also squashes random bad values
; hp_check           Sees if data is HEALPix
; get_image          Wrapper that reads both primary HDU and extensions
; get_fits           Opens, decodes and reads a FITS file
; form_label         Extracts info from header to form data labels
; parse_input        Main routine
;
; Standard IDL documentation block follows declaration of parse_input
;
PRO parse_input_event, event
COMPILE_OPT IDL2, HIDDEN
WIDGET_CONTROL, event.TOP, /DESTROY
END
;
PRO bb2nan, array
;
; Kludge to trap bad bits
;
SWITCH SIZE(array, /TYPE) OF
    4:                          ; Single precision
    6: BEGIN                    ; complex single
        nan  = !values.F_NAN
        BREAK
    END
    5:                          ; Double precision
    9: BEGIN                    ; complex double
        nan  = !values.D_NAN
        BREAK
    END
    ELSE: RETURN
ENDSWITCH

bad_bits = array GT 1e30 OR array LT -1e30
bad_bits = WHERE(bad_bits, nbb)
IF nbb GT 0 THEN BEGIN
   PRINT, nbb, $
          FORMAT="('Found',I3,' likely bad bits: setting to NaN')"
   array[bad_bits] = nan
ENDIF
       
END
;
PRO bad2nan, array, bad, RESTORE = restore
;
; Sets values of array that have HEALPix bad values to NaN
; On restore, undoes this.
;
COMPILE_OPT IDL2, HIDDEN
ON_ERROR, 2
hpx_sbadval = -1.6375E30
hpx_dbadval = -1.6375D30

SWITCH SIZE(array, /TYPE) OF
    4:                          ; Single precision
    6: BEGIN                    ; complex single
        flag = hpx_sbadval
        nan  = !values.F_NAN
        BREAK
    END
    5:                          ; Double precision
    9: BEGIN                    ; complex double
        flag = hpx_dbadval
        nan  = !values.D_NAN
        BREAK
    END
    ELSE: RETURN
ENDSWITCH

IF KEYWORD_SET(restore) THEN BEGIN
   IF bad[0] NE -1 THEN array[bad] = flag 
ENDIF ELSE BEGIN
    bad = WHERE(array EQ flag) 
    IF bad[0] NE -1 THEN array[bad] = nan
ENDELSE

END

PRO hp_check, order, header, dsize, healpix, cut4, nside
;
; Checks header and data size to see if image could be one or several
; Healpix datasets (including CUT4 format)
;
; Inputs:
;    Order:  Parse_Input input parameter: if set implies HEALPix.
;    header: Image header
;    dsize:  SIZE structure for dataset
;
; Outputs:
;    Healpix: True if dataset is healpix (including CUT4)
;    Cut4:    True if cut4 format
;    Nside:   HEALPix size parameter
;
COMPILE_OPT IDL2, HIDDEN

healpix = KEYWORD_SET(order)
cut4 = 0B

pixtype = SXPAR(header, 'PIXTYPE', COUNT = count)
hphead  = STRCMP(pixtype, 'HEALPIX', 7, /FOLD_CASE)
healpix = healpix || hphead
IF (count GT 0) && healpix && ~hphead THEN MESSAGE, /INFORMATIONAL, $
          'Assuming HEALPix per input parameters, despite header'

IF ~healpix THEN BEGIN
    test = SXPAR(header,'CTYPE1', COUNT = count)
    IF count GT 0 THEN RETURN ; it's an image
ENDIF

IF healpix || count EQ 0 THEN BEGIN
    nside = SXPAR(header, 'NSIDE')
    nsh = NPIX2NSIDE(dsize.DIMENSIONS[0])
    IF nside EQ 0 THEN nside = nsh
ENDIF
IF nside NE nsh THEN BEGIN
    IF healpix THEN BEGIN       ; Check for cut4 file
                                ; Check if we have required columns
        test = SXPAR(header, 'TTYPE1')
        cut4 = STRCMP(test, 'PIXEL', 5, /FOLD_CASE)
        test = SXPAR(header, 'TTYPE2')
        cut4 = cut4 && STRCMP(test,'SIGNAL', 5, /FOLD_CASE)
                                ; Check for required values of
                                ; optional keywords (if present)
        test = SXPAR(header, 'GRAIN', COUNT = count)
        IF count GT 0 THEN cut4 = cut4 && test EQ 1
        test = SXPAR(header, 'INDXSCHM', COUNT = count) ;
        IF count GT 0 THEN cut4 = $
          cut4 && STRCMP(test, 'EXPLICIT', 8, /FOLD_CASE)
        test = SXPAR(header, 'OBS_NPIX')
        IF test NE 0 THEN cut4 = cut4 && dsize.DIMENSIONS[0] EQ test
        IF ~cut4 THEN nside = -1
    ENDIF ELSE nside = -1       ; could be gridded HEALPix
ENDIF
IF healpix && nside EQ -1 THEN MESSAGE, $
  'Image was claimed to be HEALPix but number of pixels is invalid'

healpix = healpix || (count EQ 0 && nside GT 0)

END

PRO get_image, file, header, naxis, col, verbose, start, exten, data
;
; Wrapper to read image: covers both primary HDU and extensions.
;
; We could swap to using LUN not file name but no obvious benefit.
;
; Modified 27/06/2023: /noscale in READFITS, apply scaling externally,
; in order to correctly find blank pixels and convert to NaNs
;
COMPILE_OPT IDL2, HIDDEN

nax = SXPAR(header,'NAXIS*')
ntot = 1
FOR i = 2,naxis-1 DO ntot *= nax[i]

ncol = N_ELEMENTS(col)
IF ncol EQ 0 THEN BEGIN
    ncol = ntot
    col = INDGEN(ntot) + 1
ENDIF
                                ; Check that requested plane(s) exists
IF ntot LT MAX(col) THEN MESSAGE, 'Requested image slice(s) not found'

doall = ntot EQ 1

IF doall THEN BEGIN
   data = READFITS(file, header, EXTEN_NO = exten, /NOSCALE) 
   fits_scale, data, header
ENDIF ELSE BEGIN
    ptrs = PTRARR(ncol, /ALLOCATE_HEAP)
    
; Find wanted slices = coord (zero-relative) on last dimension
    lastdim = naxis-3
    subdim = ntot / nax[naxis-1] ; number of image planes per slice
    IF lastdim EQ 0 THEN BEGIN
        slices = col-1
        sl2get = slices
        nslice = ncol
    ENDIF ELSE BEGIN
        ipix = ARRAY_INDICES(nax[2:*], col-1, /DIMENSION)
        PRINT, 'IPIX:', ipix
        slices = ipix[lastdim,*]
        sl2get = slices[UNIQ(slices, SORT(slices))]
        nslice = N_ELEMENTS(sl2get)
    ENDELSE
    PRINT, 'Slices:', slices 
    PRINT, 'Image slices to read:', sl2get
    FOR ii = 0, nslice-1 DO BEGIN
                                ; This is a poor way to access the data,
                                ; but avoids unnecessary virtual memory usage
                                ; as long as nplanes = 1. 
                                ;Really, lower-level I/O would be better.
        slin = sl2get[ii]
        col_offset = subdim*slin + 1
        planes = WHERE(slices EQ slin, nplanes)  
        data = READFITS(file, header, NSLICE = slin, EXTEN_NO = exten, $
                       /NOSCALE)
        fits_scale, data, header
        IF subdim GT 1 THEN BEGIN
            data = REFORM(data, nax[0], nax[1], subdim, /OVERWRITE)
            jpix = col[planes] - col_offset
            FOR jj = 0, nplanes -1 DO *ptrs[planes[jj]] = data[*,*,jpix[jj]]
        ENDIF ELSE *ptrs[planes[0]] = TEMPORARY(data) 
    ENDFOR
    FOR ii=0,ncol-1 DO *ptrs[ii] = REFORM(*ptrs[ii],/OVERWRITE)
    IF ncol EQ 1 THEN BEGIN ; retrieve data from heap, at least temporarily
        data = TEMPORARY(*ptrs[0])
        PTR_FREE, ptrs
    ENDIF ELSE data = ptrs  
ENDELSE

text = exten GT 0 ? 'Image extension' : 'file'
IF verbose THEN MESSAGE, /INFORMATIONAL, text + ' read at' + $
  STRING(SYSTIME(1) - start)
IF verbose THEN HELP, /MEMORY

END

PRO get_fits, name, col, exten, file, data, header, howto, xtype, $
              verbose, start
;
; Reads a file named "name", "name.fits" or "name.FITS" and extracts
; either a set of slices from an image or a set of columns from a
; table, according to which it encounters first. It only looks in the
; main HDU and the first extension unless explicitly requested.
;
; Inputs:
;    name:   file name or root
; Input/Output
;    col:    list of required image planes or data columns
;            set to list in file if unset on input
;    exten:  FITS extension. (0 = Primary HDU). If unset look in 0, then 1
; Outputs
;    file:    completed and trimmed file name
;    data:    binary data, or pointer(s) to it
;    header:  FITS header constructed from primary header + extension
;             header (if appropriate).
;    howto:   data format... see header for parse_input
;    xtype:   string containing either "IMAGE" or "BINTABLE"
;    verbose: print diagnostics
;    start:   time (for diagnostics)
;
COMPILE_OPT IDL2, HIDDEN
ON_ERROR, 2
ncol = N_ELEMENTS(col)
col_set = ncol GT 0

; see if there is a file of this name:
file = STRTRIM(name,2)
info = FILE_INFO(file)
IF ~info.EXISTS THEN BEGIN
    file = file+'.fits'
    info = FILE_INFO(file)
    IF ~info.EXISTS THEN BEGIN
        file = STRTRIM(name,2)+'.FITS'
        info = FILE_INFO(file)
        IF ~info.EXISTS THEN MESSAGE, "Can't find file named " + name
    ENDIF
ENDIF

;error_status = 0                ; For debugging
CATCH, error_status              ; In case it's not a fits file
IF error_status THEN BEGIN
   PRINT, !error_state.msg
    CATCH, /CANCEL
    HELP, /LAST_MESSAGE
    IF N_ELEMENTS(unit) GT 0 THEN FREE_LUN, unit
    MESSAGE, 'Cannot read '+file+'. May not be FITS'
ENDIF

OPENR,  unit, file, /GET_LUN, ERROR = err
IF err NE 0 THEN MESSAGE, 'Cannot open file ' + file

MESSAGE, 'Reading '+file, /INFORMATIONAL
howto = 1

FXHREAD, unit, primhdr, status
IF status NE 0 THEN MESSAGE, 'Problem reading primary FITS header'

                                ; Is there more than just a header?
naxis = SXPAR(primhdr,'NAXIS')

; At this point exten may not be set

; widgets set exten = -1 if not set:
exten_set = N_ELEMENTS(exten) GT 0 && exten GE 0 
IF ~exten_set || exten EQ 0 THEN BEGIN
    IF  naxis GT 0 THEN BEGIN   ; Data is in primary HDU
        header = primhdr
        get_image, file, header, naxis, col, verbose, start, 0, data
;        header = primhdr
        IF SIZE(data,/TYPE) EQ 10 THEN howto = 3
        exten = 0
        FREE_LUN, unit
        RETURN                  ; all done!
    ENDIF ELSE IF exten_set THEN MESSAGE, 'No data in primary HDU'
ENDIF

; If we reach here the primary HDU was empty or unwanted. Read the
; extension

IF ~exten_set THEN exten = 1

IF naxis GT 0 THEN BEGIN
                    ; no routine available to advance to end of 
                    ; data so close & reopen file
    FREE_LUN, unit
    OPENR, unit, file, /GET_LUN
    hdu_read = 0
ENDIF ELSE hdu_read = 1

; Move forward to wanted HDU:
n2move = exten - hdu_read
IF n2move GT 0 THEN status = FXMOVE(unit, n2move)
IF status NE 0 THEN MESSAGE, $
    STRING(exten, FORMAT = "('Can''t find extension number',I3)")

FXHREAD, unit, xtnhdr, status
FREE_LUN, unit
IF status NE 0 THEN MESSAGE, 'Problem reading extension header'
xtype = SXPAR(xtnhdr,'XTENSION')

IF xtype EQ 'BINTABLE' OR xtype EQ 'A3DTABLE' THEN BEGIN
                                ; binary table
    xtype = 'BINTABLE'          ; replace synonym
    FXBOPEN, unit, file, exten, xtnhdr
                                ; Check that requested column exists
    msg = "Requested column(s) not found"
    nfields = SXPAR(xtnhdr, 'TFIELDS')

    ttype = SXPAR(xtnhdr,'TTYPE*', COUNT = count)
                                ; Check for first column = pixels
                                ; (doesn't count)
    IF count GE 1 THEN cut4 = STRCMP(ttype[0], 'PIXEL', 5, /FOLD_CASE)

    IF ~col_set THEN col = INDGEN(nfields - cut4) + 1 + cut4 ELSE $
      IF SIZE(col,/TYPE) EQ 7 THEN BEGIN ; We have a list of names...
        ttype = STRTRIM(STRUPCASE(ttype), 2)
        col   = STRTRIM(STRUPCASE(col), 2)
        match, ttype, col, colnum, COUNT= count
        IF count LT ncol THEN MESSAGE, msg
        col = colnum + 1  ; From now on use numbers not names
    ENDIF ELSE BEGIN
        IF nfields - cut4 LT MAX(col) THEN MESSAGE, msg
        IF cut4 then col += 1
    ENDELSE
    IF cut4 AND col[0] NE 1 THEN col = [1, col]
    ncol = N_ELEMENTS(col)

    IF ncol EQ 1 THEN FXBREADM, unit, col, data  ELSE BEGIN
        data = PTRARR(ncol)
        FXBREADM, unit, col, POINTERS=data, PASS_METHOD='POINTER'
                                ; Convert to 1-D arrays since the only
                                ; option here (so far) is HEALPix arrays.
        npix = N_ELEMENTS(*data[0])
        FOR i=0,ncol-1 DO *data[i] = REFORM(*data[i], npix, /OVERWRITE)

        howto = 3
    ENDELSE
    IF verbose THEN MESSAGE, /INFORMATIONAL, $
      STRING('Binary extension read at', SYSTIME(1) - start)
    IF verbose THEN HELP, /MEMORY
    FXBCLOSE, unit

    dsize = SIZE(data)
    IF dsize[0] GT 1 && ncol EQ 1 THEN $
      data = REFORM(data,dsize[dsize[0]+2],/OVERWRITE)
                                ; We only read one column
    nmand = 8
ENDIF ELSE IF xtype EQ 'IMAGE   ' THEN BEGIN ; Image extension
    naxis = SXPAR(xtnhdr,'NAXIS')
    get_image, file, header, naxis, col, verbose, start, 1, data
    xtnhdr = header
    nmand = 5+naxis
    IF SIZE(data,/TYPE) EQ 10 THEN howto = 3
ENDIF ELSE MESSAGE, 'Cannot process extension type '+xtype

; Remove SIMPLE, BITPIX, NAXIS, EXTEND and END cards from primary header and
; add rest to extension header immediately after the mandatory
; keywords.
SXDELPAR, primhdr, 'EXTEND'
lp = N_ELEMENTS(primhdr) - 2
IF lp GE 3 THEN BEGIN
    primhdr = primhdr[3:lp+1] ; Keep END card for the moment
      ; Strip keywords from primary that recur in extension:
    keys = STRMID(primhdr, 0, 8)
    pactive = WHERE (keys NE 'COMMENT ' AND keys NE 'HISTORY ' AND $
                     keys NE 'HIERARCH', nprim)
ENDIF ELSE nprim = 0

IF nprim EQ 0 THEN header = TEMPORARY(xtnhdr) ELSE BEGIN
    pkeys = keys[pactive]
    lext = N_ELEMENTS(xtnhdr)
    keys = STRMID(xtnhdr[nmand:lext-2], 0, 8)
    eactive = WHERE (keys NE 'COMMENT ' AND keys NE 'HISTORY ' AND $
                     keys NE 'HIERARCH', next)
    ekeys = keys[eactive]
    idx = -1L                   ; dummy starter value
    FOR i=0,nprim-1 DO BEGIN
        null = WHERE(pkeys[i] EQ ekeys, hit)
        IF hit GT 0 THEN SXDELPAR, primhdr, pkeys[i]
    ENDFOR
    lp = N_ELEMENTS(primhdr) - 1
    IF lp GT 0 THEN BEGIN
        lastactive = nmand+eactive[next-1]
        header = [xtnhdr[0:lastactive], '', 'COMMENT ' + $
    '===== The following lines are from the original primary header =====', $
        primhdr[0:lp-1], 'COMMENT ' + $
    '================ End of lines from primary header ==================', $
        '', xtnhdr[lastactive+1:*]]
    ENDIF ELSE header = TEMPORARY(xtnhdr)
ENDELSE
IF verbose THEN MESSAGE, /INFORMATIONAL, $
  STRING('Headers merged at', SYSTIME(1) - start)
  
END

PRO form_label, path, file, header, checksum, types, ncol, col, $
                xtype, exten, naxis, name, colab, namestr, polcode, freqcode
;
; Extracts as much information as possible from the fits header, and
; records it in the namestr structure, the name string (concatenation of 
; namestr), labels for each column (colab), and polarization codes.
; 
COMPILE_OPT IDL2, HIDDEN
ON_ERROR, 2

scodes = ['YX','XY','YY','XX', 'LR', 'RL', 'LL','RR','','I','Q','U','V', $
          'PPOL', 'FPOL', 'PANG', 'SPINDEX', 'TAU', 'RM']
          ; Values after 'V' are non-standard AIPS convention.
polcode  = REPLICATE(0, ncol)
freqcode = STRARR(ncol)
colab = STRARR(ncol)

telescope = SXPAR(header, 'TELESCOP', COUNT = count, /SILENT)
IF count EQ 0 THEN telescope = ''
instrument = SXPAR(header, 'INSTRUME', COUNT = count, /SILENT)
IF count EQ 0 THEN instrument = ''
IF instrument EQ telescope THEN instrument = ''

object = SXPAR(header, 'OBJECT', COUNT = count, /SILENT)
IF count EQ 0 THEN object = ''
creator = SXPAR(header, 'CREATOR', COUNT = count, /SILENT)
IF count EQ 0 THEN creator = ''

IF xtype EQ 'BINTABLE' THEN BEGIN
    ttype = (SXPAR(header, 'TTYPE*', COUNT = count))[col-1]
    
                                ; Remove underscores & shorten
    FOR i=0,ncol-1 DO BEGIN
        ttype[i] = STRJOIN(STRSPLIT(ttype[i],'_',/EXTRACT),' ')
        ipol  = STREGEX(ttype[i], 'polari(s|z)ation|stokes', /FOLD_CASE, $
                        LENGTH = len)
        IF ipol GE 0 THEN BEGIN
            first = STRMID(ttype[i],0,ipol)
            last  = STRMID(ttype[i],ipol+len)
            icode = WHERE(scodes EQ STRTRIM(first,2))
            IF icode EQ -1 THEN  icode = WHERE(scodes EQ STRTRIM(last,2))
            IF icode NE -1 THEN polcode[i] = icode - 8
            ttype[i] = first + 'Pol' + last
        ENDIF ELSE BEGIN
            itemp = STREGEX(ttype[i], 'temperature', /FOLD_CASE, LENGTH = len)
            IF itemp GE 0 THEN ttype[i] = STRMID(ttype[i],0,itemp) + 'Temp' $
              + STRMID(ttype[i],itemp+len)
        ENDELSE
    ENDFOR
    colab = STRTRIM(ttype,2)
ENDIF

fkind = ''
funit = ''
fval = 0
freq_set = 0B
freq_active = 0B
freq_axis = -1
stokes = ''
pol_set = ~ARRAY_EQUAL(polcode,INTARR(ncol))

IF xtype EQ 'IMAGE' && naxis GT 2 THEN BEGIN  ; Process extra axes in header
    nax = SXPAR(header,'NAXIS*', COUNT = count)
    IF count LT naxis THEN nax = [nax, REPLICATE(1,naxis-count)]
    
; Find planes requested
; The axes of the image plane are defined as being 1 & 2. Note that these
; may be, e.g. frequency & RA, not necessarily RA & Dec.
;
; "active" array has one value per non-image dimension. Value is true if the 
; selected image slices ("columns") have different values on that axis, 
; otherwise false: that is, inactive axes give a global description of the 
; data; active axes distinguish between different columns.
;
; If a non-image dimension has a size > 1 but is not active, it might become
; active if another channel from the same axis is loaded later. So we use
; such axes as the tab label if there are no active axes.
;
; colcoord(i,j) gives the (zero-relative) coordinate of column i on axis j
; (again, excluding image axes).

    colcoord = INTARR(ncol,naxis-2)
    active = BYTARR(naxis)      
    block = REPLICATE(1,naxis-2)
    FOR i = 2,naxis-2 DO block[i-1] *= nax[i]*block[i-2]
    list = col-1
    FOR i=naxis-1,2,-1 DO BEGIN
        colcoord[0,i-2] = list / block[i-2]
        list = list - colcoord[*,i-2]*block[i-2]
        active[i] = MIN(colcoord[*,i-2]) LT MAX(colcoord[*,i-2])
        IF (~active[i]) && (nax[i-2] GT 1) THEN active[i] = 2B
    ENDFOR
 
    rvals = SXPAR(header,'CRVAL*',COUNT=count)
    IF count LT naxis THEN rvals = count GT 0 $
      ? [rvals,REPLICATE(1.,naxis-count)] : REPLICATE(1.,naxis)
    rpixs = SXPAR(header,'CRPIX*',COUNT=count)
    IF count LT naxis THEN rpixs = count GT 0 $
      ? [rpixs,REPLICATE(1.,naxis-count)] : REPLICATE(1,naxis)
    delts = SXPAR(header,'CDELT*',COUNT=count)
    IF count LT naxis THEN delts = count GT 0 $
      ? [delts,REPLICATE(1.,naxis-count)] : REPLICATE(1.,naxis)
    units = SXPAR(header,'CUNIT*',COUNT=count)
    IF count LT naxis THEN units = count GT 0 $
      ? [units,REPLICATE('',naxis-count)] : REPLICATE('',naxis)

                                ; Search for frequency info
                                ; NB: non-linearities in frequency axis can
                                ; be specified via WCS info but this is not
                                ; yet handled.
    fcodes = ['FREQ', 'WAVE', 'AWAV', 'WAVN', 'ENER', $
              'VRAD', 'VOPT', 'ZOPT', 'VELO', 'BETA']
    defun  = ['Hz', 'm', 'm', 'm-1', 'J', 'm s-1', 'm s-1', '', 'm s-1', '']
    FOR i=2,naxis-1 DO BEGIN
        match = WHERE(STRCMP(types[i],fcodes,4),hit)
        IF hit GT 0 THEN BEGIN
            freq_set = 1B
            freq_active = active[i] EQ 1B
            freq_axis = i
            fkind = types[i]
            funit = units[i] NE '' ? units[i] : defun[match]
            funit = STRTRIM(funit,2)
            IF ~freq_active THEN BEGIN
                fval = delts[i]*(colcoord[0,i-2]-rpixs[i]+1) + rvals[i]
                prec = delts[i] NE 0 ? ALOG10(fval /ABS(delts[i])) > 2 : 4
                freq = numunit(fval, funit, PRECISION = prec)
                IF funit EQ '' THEN freq = fkind+': '+freq
            ENDIF
            BREAK
        ENDIF
    ENDFOR
 
                                  ; Search for stokes info
    match = WHERE(STRCMP(types,'STOKES',6),hit)
    IF hit GT 0 THEN BEGIN
        pol_set = 1B
        pol_active = active[match] EQ 1B
        IF ~pol_active THEN BEGIN
            match = match[hit-1]
            istokes = ROUND( delts[match]*(colcoord[0,match-2] $
                              - rpixs[match]+1) + rvals[match] )
            stokes = scodes[istokes+8]
            polcode[*] = istokes
        ENDIF
    ENDIF

; Now redefine active as a list of active dimensions
    many = WHERE(active EQ 2B, nmany) ; many pixels but only one selected
    active = WHERE(active EQ 1B, nactive)
    IF ncol EQ 1 THEN BEGIN
        active = many ; no genuinely active dimensions
        nactive = nmany
    ENDIF

;  make column labels.
    IF ncol GT 1 || nmany GT 0 THEN BEGIN ; At least one dimension > 2 is active
        colco = TRANSPOSE(colcoord[*, active-2])
        da = REBIN(delts[active],nactive,ncol)
        ra = REBIN(rpixs[active],nactive,ncol)
        va = REBIN(rvals[active],nactive,ncol)
        vals = da*(colco - (ra - 1)) + va
        FOR i = 0,ncol-1 DO BEGIN
            FOR j = 0,nactive-1 DO BEGIN
                IF STRCMP(types[active[j]],'STOKES',6) THEN BEGIN
                    istokes = ROUND(vals[j,i])
                    polcode[i] = istokes
; Don't include stokes in the label at this point
                ENDIF ELSE IF active[j] EQ freq_axis THEN BEGIN
                    prec = da[j] NE 0 ? ALOG10(vals[j,i]/ABS(da[j])) + 1 : 4
                    freqcode[i] = numunit(vals[j,i], funit, PRECISION = prec >2)
                    IF funit EQ '' THEN freqcode[i] = fkind+': '+freqcode[i]
; Don't include frequency in the label at this point
                ENDIF ELSE BEGIN
                    colab[i] += ' ' + types[active[j]] + ' '
                                       ; precision required
                    prec = da[j] GT 0 ? ALOG10(vals[j,i]/da[j]) > 2 : 4
                    colab[i] += numunit(vals[j,i], units[j], PRECISION = prec)
                ENDELSE
            ENDFOR
            colab[i] = STRTRIM( STRCOMPRESS(colab[i]), 2)
        ENDFOR
    ENDIF
ENDIF

IF ~freq_set THEN BEGIN                      ; Try non-standard keywords:
    freq = SXPAR(header,'FREQ', COUNT = hit)
    IF hit GT 0 THEN BEGIN
        fkind = 'FREQ'
        ftype = SIZE(freq,/TYPE)
        IF ftype NE 7 THEN BEGIN ; not a string: try to guess units
            IF freq LT 1E6 THEN freq = numunit(freq,' GHz')+' (?)' $
            ELSE freq = numunit(freq,' Hz')+' (?)'
        ENDIF
    ENDIF ELSE freq = '' 
ENDIF
IF N_ELEMENTS(freq) EQ 0 THEN freq = ''

IF ~freq_active THEN freqcode[*] = freq

IF ~pol_set THEN BEGIN                      ; Try non-standard keyword:
    stokes = SXPAR(header,'STOKES', COUNT = hit)
    IF hit GT 0 THEN BEGIN
        istokes = WHERE(scodes EQ STRTRIM(stokes,2),hit)
        IF hit GT 0 THEN polcode[*] = istokes - 8
     ENDIF ELSE stokes = ''
ENDIF
IF N_ELEMENTS(stokes) EQ 0 THEN stokes = ''

namestr = {path: path, file: file,  $
           telescope: telescope, instrument: instrument, $
           creator: creator, object: object, freq: freq, stokes: stokes, $
           extension: exten, extension_type: xtype, spectral_type: fkind, $
           header: PTR_NEW(), checksum: checksum, name: ''}


namestr.FILE       = file       ? file : 'Online data'
namestr.TELESCOPE  = telescope  ?  ' ' + STRTRIM(telescope, 2)  : ''
namestr.INSTRUMENT = instrument ?  ' ' + STRTRIM(instrument, 2) : ''
namestr.CREATOR    = creator    ?  ' Created by ' + STRTRIM(creator, 2) : ''
namestr.OBJECT     = object     ?  ' ' + STRTRIM(object, 2)  : ''
namestr.FREQ       = freq       ? '  ' + STRTRIM(freq, 2)   : ''
namestr.STOKES     = stokes     ? ' Stokes ' + STRTRIM(stokes, 2) : ''

name = namestr.FILE + ': '
FOR i=2,7 DO name += namestr.(i)

; If we have definite polarizations for some columns, assume the rest
; are Stokes I (gets round CMB habit of labelling I as "temperature").
test = WHERE(polcode NE 0, npol, COMPLEMENT=not_set)
IF npol GT 0 && npol LT ncol THEN polcode[not_set] = 1

test = WHERE(colab EQ '', nblank)
IF nblank NE 0 THEN colab[test] = 'Map '+STRTRIM(STRING(test+1),2)

END

PRO parse_input, thing, ORDER = order, PROJ = proj, TEMPORARY = temporary, $
                 data, header, name, howto, colab, scale_par, namestr, $
                 polcode, freqcode, GET_SCALE = get_scale, COLUMN = col, $
                 EXTENSION = exten, PLOT = do_plot, VERBOSE = verbose
;+
; NAME:
;       PARSE_INPUT
;
; PURPOSE:
;       Decodes input parameter (thing): is it
;
;       (a) A stack of 2-D images
;       (b) A stack of HEALPix (HP) arrays
;       (c) A file name for a FITS file containing images or HP arrays
;       (d) A file name minus the ".fits" or ".FITS" ending
;       (e) A structure containing header + set of 2-D images (one per tag)
;       (f) A structure containing header + a set of HP arrays (one per tag)
;       (g) A pointer to any of the above
;       (h) An array of pointers to 2-D images or HP arrays.
;
;       Headers are assumed to be more or less FITS.
;       Returns extracted data, header, and information from the
;       header. Optionally gets robust statistics of the data for use
;       in scaling. Optionally, plots histograms to check statistics.
;
; CATEGORY:
;       Input/Output, Widget
;
; CALLING SEQUENCE:
;
;       PARSE_INPUT, Thing, Data, Header, Name, Howto, Colab, $
;                    Scale_par, Namestr, Polcode
;
; INPUTS:
;       Thing:     An IDL variable that somehow specifies the wanted data
;
; KEYWORD PARAMETERS:
;       ORDER:     'RING' or 'NESTED' (default: header value or 'RING')
;
;       PROJ:      'GRID', 'NPOLE', 'SPOLE' if this is a healpix array
;                  (default: 'GRID'), or if it is already gridded but
;                  header is missing.
;
;       COLUMN:    Requested table column(s) or image plane(s). (Below
;                  referred to simply as maps). Default: all
;                  First col is #1.
;                  Pixel column in cut4 files does not count.
;
;       EXTENSION: FITS extension containing the data
;                  Default: look in 0, then 1.
;
;       GET_SCALE: If set, get statistics of the data
;
;       PLOT:      Set to launch a widget containing plots of the
;                  histogram of each map around the noise level (only
;                  if /GET_SCALE is also set).
;
;       TEMPORARY: Input array should be overwritten to save space
;
;       VERBOSE:   Print diagnostics
;
; OUTPUTS:
;       Data:      2+D image or array of pointers to 2+D image
;                  (ie 2D or 3D if more than one map requested).
;
;       Header:    FITS-like header describing data (may not fully
;                  conform to FITS standard)
;
;       Name:      Text string to use as suggested title
;
;       Howto:     = 0: use the original variable (= 2D image)
;                  = 1: data is the 2D image
;                  = 2: data is a pointer to a 2D image
;                  = 3: data is an array of pointers to 2D images.
;
; OPTIONAL OUTPUTS:
;      Colab:      Array of strings usable as labels for individual maps
;
;      Scale_par:  [4, N_map] array containing min, max, mode, standard
;                  deviation in each map. (Set if /GET_SCALE is set).
;
;      Namestr:    Structure containing header elements suitable for
;                  constructing the "Name" and "Colab" strings.
;
;      Polcode:    Array of code numbers describing polarization each
;                  map: standard FITS "STOKES" codes plus:
;                     5 = linearly polarized intensity
;                     6 = Fractional linear polarization
;                     7 = Polarization angle
;                     8 = Spectral Index
;                     9 = Optical depth
;                    10 = Rotation Measure
;
; SIDE EFFECTS:
;      IF /GET_SCALE and /PLOT are set, a widget is launched showing
;      the histograms of each map around the mode, overplotted with
;      the parabolic fit used to find the mode.
;
;      Thing may be overwritten if /TEMPORARY is set, and if not, HEALPix
;      bad values will be converted to NaNs.
;
; EXAMPLE:
;               dummy = FINDGEN(100,100)
;               PARSE_INPUT, dummy, data, header, name, howto
;
;      Returns howto = 0, data is undefined, header is a minimal FITS
;      header, and name = "Online data:"
;
;
; MODIFICATION HISTORY:
;       Written by:     J. P. Leahy, January 2008
;       March 2008:     Added CUT4 support
;       April 2008:     Bug fixes
;       August 2013:    Bug fixes, form_label extracted from main program &
;                       updated.
;       Sept  2020:     Added trapping of bad bits
;-
COMPILE_OPT IDL2, HIDDEN
ON_ERROR, 2

start = SYSTIME(1)
ispointer = 0
path = ''  &  file = ''
xtype = 'IMAGE'
do_plot = KEYWORD_SET(do_plot)  ; Flag to plot histogram of data.
IF N_ELEMENTS(proj) EQ 0 THEN proj = ''

temporary = KEYWORD_SET(temporary)
; Set temporary: automatically true if there is nothing to return it to...
intype = SIZE(thing, /TYPE)      ; 10 = pointer
temporary = temporary || (intype NE 10 && ~ARG_PRESENT(thing))

verbose = KEYWORD_SET(verbose)
tsize = SIZE(thing,/STRUCTURE)
ncol = N_ELEMENTS(col)
col_set = ncol GT 0
IF ncol GT 1 THEN BEGIN
    col = col[UNIQ(col, BSORT(col))] ; sort & eliminate duplicates
    ncol = N_ELEMENTS(col)
ENDIF
IF col_set && MIN(col) LT 1 THEN MESSAGE, 'Column numbers start at 1'

REDO:

CASE tsize.TYPE OF
    0: MESSAGE, 'Undefined image parameter'

    7: BEGIN ; String... should be a file name
        name = ispointer ? (*thing)[0] : thing[0]

        get_fits, name, col, exten, file, data, header, howto, xtype, $
          verbose, start
                        ; Strip directory info from filename (if any)
        sep = PATH_SEP()
        IF sep EQ '\' THEN sep = '\\'
        firstchar = STRSPLIT(file, '.*'+sep, /REGEX)
        path = STRMID(file, 0, firstchar-1)
        file = STRMID(file, firstchar)
        temporary = 1
        ncol = N_ELEMENTS(col)
        col_set = ncol GT 0
        IF ~col_set THEN MESSAGE, 'Internal error: no columns read from file'
    END

    8: BEGIN                    ; Structure or pointer to structure
        IF tsize.N_DIMENSIONS GT 1 THEN $
          MESSAGE, 'Expected image but found array of structures'

                                ; Check it has header and (potential) data:
        header = ispointer ? *thing.(0) : thing.(0)
        hdtype = SIZE(header, /TYPE)
        badhd = hdtype NE 7     ; Should be a FITS header, ie strings.
        IF ~badhd THEN BEGIN    ; check it is a FITS-ish header
            simple = SXPAR(header,'SIMPLE')
            IF ~simple THEN BEGIN ; Not primary header
                xtype = SXPAR(header,'XTENSION',COUNT=count)
                exten = 1 ; arbitrarily define extension number
                badhd = count EQ 0
            ENDIF ELSE exten = 0
        ENDIF
        ntags = N_TAGS(ispointer ? *thing : thing)
        IF badhd || ntags LT 2 THEN $
          MESSAGE, 'Unknown structure type supplied as input'
        nfields = ntags - 1     ; Number of tags with data
        tagnames = TAG_NAMES(ispointer ? *thing : thing)
                                ; Check for cut4 file
        IF nfields GE 2 THEN $
          cut4 = STRCMP(tagnames[1],  'PIXEL', 5, /FOLD_CASE) && $
                 STRCMP(tagnames[2], 'SIGNAL', 6, /FOLD_CASE) ELSE cut4 = 0

                                ; Have a look at first item to see
                                ; what type it is:
        dsize = SIZE( ispointer ? *thing.(1): thing.(1) ,/STRUCTURE) 
        datatype = dsize.TYPE
        IF datatype EQ 10 THEN BEGIN
            howto = 2
            dsize = SIZE(*(ispointer ? *thing.(1) : thing.(1)), /STRUCTURE)
        ENDIF ELSE howto = 1
        wasone =  howto EQ 1

        is_image = 0B
        IF ntags EQ 2 && dsize.N_DIMENSIONS GT 1 THEN BEGIN
            is_image = 1B
                                ; See if this is a healpix dataset:
            hp_check, order, header, dsize, healpix, cut4, nside
            npix = dsize.N_ELEMENTS
            nax  = dsize.DIMENSIONS
            mappix = healpix ? nax[0] : nax[0] * nax[1]
            nfields =  npix / mappix
        ENDIF

; Find which columns we are supposed to read by interpreting the
; "col" input parameter
                                ; Potential message if things go wrong:
        msg = 'Input structure does not contain all requested columns'
        IF ~col_set THEN BEGIN
            col = INDGEN(nfields) + 1
            col_set = 1B
        ENDIF ELSE IF ~is_image && SIZE(col, /TYPE) EQ 7 THEN BEGIN
                                ; list of names...
            tagnames = STRTRIM(STRUPCASE(tagnames), 2)
            col      = STRTRIM(STRUPCASE(col), 2)
            match, tagnames[1:*], col, colnum, COUNT = count
            IF count LT ncol THEN MESSAGE, msg
            col = colnum + 1    ; From now on use numbers not names
        ENDIF ELSE IF cut4 then col += 1

        IF cut4 AND col[0] NE 1 THEN col = [1, col]
        ncol = N_ELEMENTS(col)
        IF nfields-cut4 LT MAX(col) THEN MESSAGE, msg

; Read in first column:
        data = ispointer ? *thing.(col[0]): thing.(col[0])
        
        IF ncol GT 1 THEN BEGIN ; Store data via a pointer array
            tmp = TEMPORARY(data)
            IF is_image THEN BEGIN ; one 2+D image
                                ; Delete original now, as we need yet
                                ; another copy...
                IF temporary THEN BEGIN
                    IF ispointer THEN PTR_FREE, thing
                    ispointer = 0B
                ENDIF
                IF howto EQ 2 THEN BEGIN ; swap to non-heap variable
                    tmp2 = temporary ? TEMPORARY(*tmp) : *tmp
                    IF temporary THEN PTR_FREE, tmp
                    tmp = TEMPORARY(tmp2)
                ENDIF

                tmp = REFORM(tmp, mappix, nfields, /OVERWRITE)
                data = PTRARR(ncol, /ALLOCATE_HEAP)
                FOR i=0,ncol-1 DO BEGIN
                    tmp2 = tmp[*,col[i]-1]
                    IF ~healpix THEN tmp2 = REFORM(tmp2, nax[0], nax[1], $
                                                   /OVERWRITE)
                    *data[i] = TEMPORARY(tmp2)
                ENDFOR
                tmp = 0
            ENDIF ELSE BEGIN ; One col per tab
                IF howto EQ 2 THEN BEGIN
                    data = PTRARR(ncol)
                    data[0] = tmp
                    FOR i=1,ncol-1 DO $
                      data[i]  = ispointer ? *thing.(col[i]) : thing.(col[i])
                ENDIF ELSE BEGIN
                    data = PTRARR(ncol, /ALLOCATE_HEAP)
                    *data[0] = TEMPORARY(tmp)
                    FOR i=1,ncol-1 DO $
                      *data[i] = ispointer ? *thing.(col[i]) : thing.(col[i])
                 ENDELSE
            ENDELSE
            howto = 3
        ENDIF
        IF temporary THEN BEGIN ; Delete thing
; NB: if structure contained pointers, their heap variables have been
; copied to data, unless is_image in which case any pointer was freed above.
            IF ispointer THEN PTR_FREE, thing
            thing = 0
        ENDIF
        IF wasone THEN temporary = 1B ; we can delete data heap variable(s)
    END

    10: BEGIN ; Pointer
        IF ispointer THEN MESSAGE, 'Cannot parse chains of pointers'
        ndims = tsize.N_DIMENSIONS
        CASE ndims OF
            0: BEGIN ; Pointer to something else: find out what
                tsize = SIZE(*thing,/STRUCTURE)
                ispointer = 1
                GOTO, REDO
            END
            1: BEGIN ; Should be array of pointers to 1- or 2-D datasets
                ncol = tsize.DIMENSIONS[0]
                testcol = INDGEN(ncol) + 1
                IF col_set && ARRAY_EQUAL(col,testcol) EQ 0B THEN MESSAGE, $
                  'Please trim data array to required size'
                col = testcol
                col_set = 1B
                T0 = SIZE(*thing[0])
                SWITCH T0[0] OF
                    0: MESSAGE, $
                      'Pointer arrays should point directly to bulk data!'
                    1: ; Same as 2:
                    2: BEGIN
                        FOR i=1,tsize.N_DIMENSIONS DO BEGIN
                            T = SIZE(*thing[i])
                            IF ARRAY_EQUAL(T0[0:T0[0]], T[0:T[0]]) EQ 0B $
                              THEN MESSAGE, 'Mismatched array sizes in input'
                        ENDFOR
                        tsize = SIZE(*thing[0], /STRUCTURE)
                        MKHDR, header, tsize.TYPE, $
                          [tsize.DIMENSIONS[0:T0[0]-1], ndims]
                        data = thing
                        howto = 3
                        exten = 0
                        BREAK
                    END
                    ELSE: MESSAGE, 'Each pointer should point to 1 or 2-D data'
                ENDSWITCH
            END
            ELSE: MESSAGE, $
              'Expected image but found multi-dimensional array of pointers'
        ENDCASE
    END
    ELSE: BEGIN                 ; Array (hopefully) of numbers
                ; Make sure we have only been passed the planes needed

        MKHDR, header, tsize.TYPE, tsize.DIMENSIONS[0:tsize.N_DIMENSIONS-1]
                                ; See if this is a healpix dataset:
        exten = 0
        hp_check, order, header, tsize, healpix, cut4, nside
        npix  = tsize.N_ELEMENTS
        nax   = tsize.DIMENSIONS
        mappix = healpix ? nax[0] : nax[0] * nax[1]
        IF mappix EQ 0 THEN MESSAGE, $
          'Cannot process input: 1-D and does not seem to be HEALPix'

        nblock =  npix / mappix
        testcol = LINDGEN(nblock) + 1
        IF col_set && ~ARRAY_EQUAL(col, testcol) THEN MESSAGE, $
                'Please trim data array to required size'
        col = testcol           ; Assume that all these data passed are wanted.
        ncol = nblock
        col_set = 1B
        IF ncol GT 1 THEN BEGIN
                                ; swap to non-heap variable, split into cols,
                                ; & restore to heap.
            IF temporary THEN BEGIN
                tmp = TEMPORARY(ispointer ? *thing : thing)
                IF ispointer THEN PTR_FREE, thing
            ENDIF ELSE tmp = ispointer ? *thing : thing
            tmp = REFORM(tmp, mappix, nblock, /OVERWRITE)
            data = PTRARR(ncol, /ALLOCATE_HEAP)
            FOR i=0,ncol-1 DO BEGIN
                tmp2 = tmp[*,col[i]-1]
                IF ~healpix THEN tmp2 = REFORM(tmp2, nax[0], nax[1], $
                                               /OVERWRITE)
                *data[i] = TEMPORARY(tmp2)
            ENDFOR
            tmp = 0
            howto = 3
        ENDIF ELSE IF ispointer THEN BEGIN
            data = thing
            howto = 2
        ENDIF ELSE IF temporary THEN BEGIN
            data = TEMPORARY(thing)
            howto = 1
        ENDIF ELSE howto = 0
    END
ENDCASE

IF verbose THEN BEGIN
    MESSAGE, /INFORMATIONAL, $
    STRING('Got '+xtype+' data at', SYSTIME(1) - start)
    HELP, /MEMORY
ENDIF
naxis  = SXPAR(header,'NAXIS')

CASE howto OF
    0: dsize = SIZE(thing, /STRUCTURE)
    1: dsize = SIZE(data, /STRUCTURE)
    2: dsize = SIZE(*data, /STRUCTURE)
    3: dsize = SIZE(*data[0], /STRUCTURE)
ENDCASE
IF dsize.N_DIMENSIONS EQ 0 THEN $
  MESSAGE, 'Found scalar when image expected'

; What kind of data do we have?
IF howto EQ 3 THEN BEGIN
    type = INTARR(ncol)
    FOR i=0,ncol-1 DO type[i] = SIZE(*data[i], /TYPE)
ENDIF ELSE type = dsize.TYPE

utype = type[UNIQ(type, SORT(type))]
match, utype, [7, 8, 10, 11], suba, subb
badtypes = ['STRING', 'STRUCTURE', 'POINTER', 'OBJREF']
IF subb[0] NE -1 THEN MESSAGE, 'Cannot process data type ' + badtypes[subb]

; Find safe data type for means etc:
match, utype, [4, 5, 6, 9], suba, subb
IF subb[0] NE -1 THEN stype = MAX(utype[suba]) ELSE stype = 5

; See if this is a healpix dataset, if not already done:
IF N_ELEMENTS(healpix) EQ 0 THEN $
  hp_check, order, header, dsize, healpix, cut4, nside
  
checksum32, BYTE(header), checksum

IF healpix THEN BEGIN           ; convert HEALPix bad values to NaN
                                ; BEFORE trying to scale image!
    IF howto GE 2 THEN BEGIN
        IF ~temporary THEN  badarr = PTRARR(ncol,/ALLOCATE_HEAP)
        FOR icol = 0L, ncol-1 DO BEGIN
            dptr = data[icol]
            IF ~temporary THEN BEGIN
                badptr = badarr[icol]
                bad2nan, *dptr, *badptr
            ENDIF ELSE bad2nan, *dptr
        ENDFOR 
    ENDIF ELSE IF howto EQ 1 THEN bad2nan, data ELSE bad2nan, thing, bad_idx
ENDIF

; Trap bad bits
IF howto GE 2 THEN BEGIN
   FOR icol = 0L, ncol-1 DO BEGIN
      dptr = data[icol]
      bb2nan, *dptr
   ENDFOR 
ENDIF ELSE IF howto EQ 1 THEN bb2nan, data ELSE bb2nan, thing

IF KEYWORD_SET(get_scale) THEN BEGIN
                                ; Analyse each channel for scaling purposes:
    scale_par = MAKE_ARRAY(4, ncol, TYPE = stype)
    nblock = ncol - cut4

    IF do_plot THEN BEGIN ; Launch widget for plot
        swap_lut, dummy1, dummy, old_graph
        ttit = 'Data Histogram'
        IF nblock GT 1 THEN ttit += 's'
        tlb = WIDGET_BASE(/COLUMN, TITLE = ttit)
        row  = WIDGET_BASE(tlb, /ROW)
        draw   = LONARR(nblock)
        FOR imap = 0, nblock-1 DO BEGIN
            icol = col[imap + cut4]
            base = WIDGET_BASE(row, /COLUMN)
            void = WIDGET_LABEL(base, VALUE = $
                                'Column '+STRTRIM(STRING(icol),2) )
            draw[imap] = WIDGET_DRAW(base, XSIZE=300, YSIZE=300, RETAIN = 2)
        ENDFOR
        done = WIDGET_BUTTON(tlb, VALUE = 'Done')
        WIDGET_CONTROL, tlb, /REALIZE
    ENDIF

    FOR imap = 0, nblock-1 DO BEGIN
        IF do_plot THEN BEGIN
            WIDGET_CONTROL, draw[imap], GET_VALUE = window
            WSET, window
        ENDIF
        scaling_params, thing, data, imap+cut4, healpix, howto, col, $
          ar, mode, sdev, PLOT = do_plot
        PRINT, imap, ar, mode, sdev, FORMAT = $
       "('Column',I2,': Min & MAX:',2E11.3,', Mode:',E11.3,', est RMS:',E10.3)"
        scale_par[0,imap] = [ar, mode, sdev]
    ENDFOR

    IF verbose THEN MESSAGE, /INFORMATIONAL, $
      STRING('Scaling parameters derived', SYSTIME(1) - start)
    IF verbose THEN HELP, /MEMORY
    IF do_plot THEN BEGIN
                                ; Allow user to get rid of plot!
        XMANAGER, 'parse_input', tlb, /NO_BLOCK
        restore_lut, dummy, old_graph
    ENDIF
ENDIF

IF healpix THEN BEGIN           ; Now convert to grid
    IF cut4 THEN BEGIN          ; must be howto = 3
        tmp = cut4grid(*data[0], data[1:*], header, order, proj, $
                       VERBOSE = verbose)
        IF temporary THEN PTR_FREE, data
        data = TEMPORARY(tmp)
        col = col[1:*] - 1      ; We have eliminated col 0 = PIXEL
        ncol = ncol - 1
    ENDIF ELSE CASE howto OF
        0: data = hpgrid(thing, header, order, proj, VERBOSE = verbose)
        1: data = hpgrid(TEMPORARY(data), header, order, proj, VERBOSE=verbose)
        2: data = temporary ? hpgrid(TEMPORARY(*data), header, order, proj, $
                                     VERBOSE = verbose) $
                            : hpgrid(*data, header, order, proj, $
                                     VERBOSE = verbose)
        3: BEGIN
            tmp = hpgrid(data, header, order, proj, VERBOSE = verbose)
            IF temporary THEN PTR_FREE, data
            data = TEMPORARY(tmp)
        END
    ENDCASE
                ; Restore healpix bad values in original array if required
    IF ~temporary THEN BEGIN
        IF howto EQ 0 THEN bad2nan, thing, bad_idx, /RESTORE
        IF howto GE 2 THEN BEGIN
            FOR icol = 0L, ncol-1 DO BEGIN
                dptr = thing[icol]
                badptr = badarr[icol]  
                bad2nan, *dptr, *badptr, /RESTORE
            ENDFOR 
            PTR_FREE, badarr
        ENDIF
    ENDIF
    
    IF dsize.N_DIMENSIONS GT 1 THEN howto = 3 ; we get a pointer array

    IF howto LT 3 THEN BEGIN
        dsize = SIZE(data, /STRUCTURE)
        howto = 1
    ENDIF ELSE dsize = SIZE(*data[0], /STRUCTURE)
    
    temporary = 1B

    IF verbose THEN MESSAGE, /INFORMATIONAL, $
      STRING('Converted to HEALPix at', SYSTIME(1) - start)

ENDIF ELSE BEGIN                ; definitely not HEALPix array
    IF dsize.N_DIMENSIONS EQ 1 THEN $
      MESSAGE, 'Found 1-D array when image expected'
ENDELSE

IF verbose THEN HELP, /MEMORY

ndims = dsize.N_DIMENSIONS
extras = ndims GT 2

; If we are allowed to, strip data down to 2D (but this should never happen)
IF ~col_set THEN BEGIN
    MESSAGE, /INFORMATIONAL, 'Internal problem: col unset'
    CASE howto OF
        0:                      ; Do nothing: never temporary
        1: IF extras THEN data = data[*,*]
        2: IF temporary && extras THEN *data = *data[*,*]
    ENDCASE
    IF temporary OR howto EQ 1 THEN extras = 0
ENDIF

; Interpret & maybe modify header:

types = SXPAR(header,'CTYPE*',COUNT = count)
IF count LT naxis THEN BEGIN
    extra = STRING(INDGEN(naxis-count)+count+1, FORMAT = "('AXIS',I2)")
    types = count GT 0 ? [types,extra] : extra
ENDIF ELSE IF count GT naxis THEN naxis = count
types = STRTRIM(types,2)

IF ~healpix && proj NE '' THEN BEGIN
                                ; May be already gridded HEALPix
    ngrid = dsize.DIMENSIONS[0]
    nside = STRCMP(proj,'GRID', 4, /FOLD_CASE) ? ngrid/5 : ngrid/4
    npix = NSIDE2NPIX(nside)
    IF npix EQ -1 OR ngrid NE dsize.DIMENSIONS[1] THEN BEGIN
        MESSAGE, /INFORMATIONAL, 'Projection set to "'+proj+'" but image wrong size'
    ENDIF ELSE BEGIN            ; Update header
        header = extras ? grid_header(header, nside, proj, ndims-2, $
                                      dsize.DIMENSIONS[2:ndims-1]) $
                        : grid_header(header, nside, proj, 0, 0)
    ENDELSE
ENDIF

IF verbose THEN MESSAGE, /INFORMATIONAL, $
  STRING('Data ready', SYSTIME(1) - start)
IF verbose THEN HELP,/MEMORY

;
; Process header to get name, colab, namestr, polcode, freqcode:
;
form_label, path, file, header, checksum, types, ncol, col, $
            xtype, exten, naxis, name, colab, namestr, polcode, freqcode

IF verbose THEN MESSAGE, /INFORMATIONAL, $
  STRING('Labels set', SYSTIME(1) - start)
END
