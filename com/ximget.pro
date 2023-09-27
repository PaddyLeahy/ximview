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
FUNCTION ximget, tab, header, NO_COPY = no_copy, STRUCTURE = structure, $
                 POINTER = pointer
;+
; NAME: 
;       XIMGET
;
; PURPOSE: 
;       Retrieves an image displayed on a Ximview tab, making it
;       available on the IDL command line.
;
; CATEGORY:
;       Widget helpers
;
; CALLING SEQUENCE: 
;
;        Result = XIMGET(Tab, [header], [/STRUCTURE], [/NO_COPY], [/POINTER])
;
; INPUTS:
;       tab:     The tab containing the required image. Specify as
;                either the index number (starting with 0 on the left)
;                or the tab title. 
;
; KEYWORD PARAMETERS:
;       STRUCTURE: return IDL structure rather than simple array, with fields
;                  - HDR (FITS header)
;                  - IMG (2-D image array).
;
;       POINTER: returns pointer to image instead of the image itself.
;
;       NO_COPY: Set to save memory by transfering the image rather
;                than copying it. Does not help if STRUCTURE or POINTER
;                is set, so is ignored in that case.
;
; OUTPUTS:
;       Result is by default a 2-D array, of whatever size and type is stored 
;       in XIMVIEW. If the requested tab shows an RGB image, a 3-D
;       byte array is returned, with third dimension for RGB. 
;       If /STRUCTURE is set, result is a structure with two
;       fields: .HDR contains the FITS header, .IMG contains the image
;       If POINTER is set, returned value or .IMG is a pointer
;       (pointer triplet for RGB). 
;
; OPTIONAL OUTPUTS:
;       header:  FITS header describing the image (redundant if
;       /STRUCTURE is set). 
;
; SIDE EFFECTS:
;       If NO_COPY is set, the specified tab on the XIMVIEW widget is
;       deleted.  
;
; RESTRICTIONS:
;       Does not work for tabs containing RGB or HSV images.
;
;       Specifying /NO_COPY when there is only one tab is ignored,
;       since the last tab cannot be deleted. 
;
; EXAMPLE:
;       Display some WMAP map data with XIMVIEW and copy the total
;       intensity image to the command line environment:
;
;            XIMVIEW, 'wmap_band_iqumap_r9_3yr_K_v2', '*', COL=[1,2,3]
;            stokes_I = XIMGET(0)
;
; MODIFICATION HISTORY:
;       Written by:      J. P. Leahy, Feb 2008
;                        header, structure & pointer options added
;                        2013
;-
COMPILE_OPT IDL2
ON_ERROR, 2

make_struct = KEYWORD_SET(structure)
pointer = KEYWORD_SET(pointer)
get_head = ARG_PRESENT(header) || make_struct

; Find Ximview

test = WIDGET_INFO(/MANAGED)
IF test[0] EQ 0 THEN MESSAGE, 'No widgets currently being managed'

FOR i = 0,N_ELEMENTS(test)-1 DO BEGIN
    uname = WIDGET_INFO(test[i], /UNAME)
    IF uname EQ 'XIMVIEW' THEN GOTO, GOTIT
ENDFOR

MESSAGE, 'Ximview is not currently being managed'

GOTIT:

top = test[i]
WIDGET_CONTROL, top, GET_UVALUE = state
labels = STRTRIM(get_tab_uvals(state.TABS), 2)

; Find the appropriate tab
IF SIZE(tab,/TYPE) EQ 7 THEN BEGIN ; we have a tab name not number
    itab = WHERE(STRCMP(labels, STRTRIM(tab, 2), /FOLD_CASE), ntab)
    IF itab[0] EQ -1 THEN MESSAGE, 'Tab ' + tab + ' not found'
    IF ntab GT 1 THEN $
        MESSAGE, /INFORMATIONAL, 'Ambiguous tab name, returning first'
    itab = itab[0]
ENDIF ELSE itab = tab

screens = (*state.TABARR).SCREEN
iscreen = WHERE(itab[0] EQ screens)
IF iscreen EQ -1 THEN MESSAGE, 'Request to get non-existent tab or blink tab: quitting'

str = (*state.TABARR)[iscreen]

mono = PTR_VALID(str.BYTE_PTR)
IF ~mono THEN BEGIN
    MESSAGE, /INFORMATIONAL, 'Tab '+labels[itab]+' is an RGB image'
    IF pointer THEN MESSAGE, /INFO, 'Returning RGB triplet of pointers' $
    ELSE MESSAGE, /INFO, 'Returning byte image cube, third dimension is RGB'
    dims = SIZE(*str.RGB[0], /DIMENSIONS)
ENDIF

no_copy = KEYWORD_SET(no_copy) 
IF no_copy THEN BEGIN
    IF N_ELEMENTS(screens) LE 1 THEN BEGIN
        no_copy = 0B
        MESSAGE, /INFORMATIONAL, 'Overriding /NO_COPY to preserve last tab'
    ENDIF ELSE IF make_struct THEN BEGIN
        no_copy = 0B
        MESSAGE, /INFORMATIONAL, $
          'Overriding /NO_COPY as it does not save space for structure output'
    ENDIF ELSE IF pointer THEN BEGIN
        no_copy = 0B
        MESSAGE, /INFORMATIONAL, $
          'Overriding /NO_COPY as there is no advantage for pointer output'
    ENDIF ELSE IF ~mono THEN BEGIN
        no_copy = 0B
        MESSAGE, /INFORMATIONAL, $
          'Overriding /NO_COPY for RGB image as may have unexpected consequences'
    ENDIF
ENDIF

; get header
IF get_head THEN BEGIN
    map_seq = str.MAP_SEQ
    IF map_seq EQ -1 THEN MESSAGE, 'Header requested but this RGB image has no header'

    filestr = (*state.FILES)[str.MAP_SEQ-1]
    header = *filestr.HEADER
    col = str.COLUMN

    test = SXPAR(header, 'SIMPLE', COUNT = hit)
    IF hit EQ 0 THEN SXADDPAR, header, 'SIMPLE', 'T', ' FITS format', 'BITPIX'

    ttypes = SXPAR(header, 'TTYPE*', COUNT = ntype)
    IF ntype GT 0 THEN BEGIN ; Converted from binary table extension
        units = SXPAR(header,'TUNIT*')
        SXADDPAR, header, 'BUNIT', units[col-1]
        SXADDHIST, 'Original HEALPix column label: '+ ttypes[col - 1], header
        istr = STRTRIM( STRING(INDGEN(ntype)+1), 2)
        SXDELPAR, header, 'TTYPE'+istr
        SXDELPAR, header, 'TUNIT'+istr
    ENDIF ELSE BEGIN  ; data was an image.
        naxis = SXPAR(header,'NAXIS')
        nax   = SXPAR(header,'NAXIS*')
        IF naxis GT 2 THEN BEGIN
                                     ; Find FITS pixel number of our column on axes 3:naxis
            ipix = ARRAY_INDICES(nax[2:*], col-1, /DIMENSION) + 1 
            rpixs = SXPAR(header,'CRPIX*',COUNT=count)
            IF count LT naxis THEN rpixs = count GT 0 $
              ? [rpixs,REPLICATE(1.,naxis-count)] : REPLICATE(1,naxis)
            FOR i = 2, naxis-1 DO BEGIN
                istr = STRTRIM(STRING(i+1), 2)
                SXADDPAR, header, 'NAXIS'+istr, 1
                SXADDPAR, header, 'CRPIX'+istr, rpixs[i] - ipix[i-2] + 1
            ENDFOR
        ENDIF
    ENDELSE
ENDIF

IF pointer THEN BEGIN
    ptr = mono ? str.IM_PTR : str.RGB
    IF make_struct THEN data = {HDR: header, IMG: ptr} ELSE data = ptr
ENDIF ELSE IF make_struct THEN BEGIN
    IF ~mono THEN BEGIN
        data = {HDR: header, IMG: BYTARR(dims[0], dims[1], 3, /NOZERO)}
        FOR i=0,2 DO data.IMG[0,0,i] = *str.RGB[i]
    ENDIF ELSE data = {HDR: header, IMG: *str.IM_PTR} 
ENDIF ELSE BEGIN
    IF ~mono THEN BEGIN
        data = BYTARR(dims[0], dims[1], 3, /NOZERO)
        FOR i=0,2 DO data[0,0,i] = *str.RGB[i]
    ENDIF ELSE IF no_copy THEN BEGIN
        ptr = str.IM_PTR
        data = TEMPORARY(*ptr)
    ENDIF ELSE data = *str.IM_PTR
ENDELSE

; Delete tab if NO_COPY set
IF no_copy THEN BEGIN
    newevent = {ID: 0L, TOP: 0L, HANDLER: 0L, TAB: itab}
    WIDGET_CONTROL, state.PAD2, SEND_EVENT = newevent
ENDIF

RETURN, data

END

