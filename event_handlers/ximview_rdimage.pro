; -----------------------------------------------------------------------------
;
;  Copyright (C) 2013   J. P. Leahy
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
; Module XIMVIEW_RDIMAGE
;
; J. P. Leahy August 2013
;
; Contains the event handlers for loading image files on the fly in Ximview.
;
PRO ximview_rdimage_event, event
; Handles events from _rdimage dialog
;
COMPILE_OPT IDL2, HIDDEN
ON_ERROR, 1
WIDGET_CONTROL, event.TOP, GET_UVALUE = info

tag = TAG_NAMES(event, /STRUCTURE_NAME)
CASE STRMID(tag,0,12) OF
    'FILESEL_EVEN': BEGIN
        CASE event.DONE OF
            0: RETURN
            1: BEGIN
                WIDGET_CONTROL, info.INDID, GET_VALUE = index
                WIDGET_CONTROL, info.TITID, GET_VALUE = title
                result = {cancel: 0B, index: index, title: title, $
                          file: event.VALUE}
            END
            2: result = {cancel: 1B}
        ENDCASE
        WIDGET_CONTROL, info.RETID, SET_UVALUE = result
        WIDGET_CONTROL, event.TOP, /DESTROY
    END
    'WIDGET_TEXT':             ; do nothing
    ELSE: MESSAGE, /INFORMATIONAL, 'Unexpected event type '+tag
ENDCASE
END

PRO ximview_rdimage, event
;  Reads an IDL-compatible image file and loads it in a new tab
;
COMPILE_OPT IDL2, HIDDEN
COMMON gr_global
ON_ERROR, 1

start = SYSTIME(1)
WIDGET_CONTROL, event.TOP, GET_UVALUE = state
verbose = state.VERBOSE
first = state.FIRST
WIDGET_CONTROL, event.ID, GET_UVALUE = save
query = WIDGET_BASE(GROUP_LEADER = event.TOP, /MODAL, /COLUMN, $
                    TITLE = 'Enter image file name')
base  = WIDGET_BASE(query, /ROW)
void  = WIDGET_LABEL(base, VALUE = 'Index No. (multi-image files):')
indID = WIDGET_TEXT(base, VALUE = '0', XSIZE = 3, /EDITABLE)
base2 = WIDGET_BASE(query, /ROW)
void  = WIDGET_LABEL(base2, VALUE = 'Title:')
titID = WIDGET_TEXT(base2, VALUE = '', XSIZE = 40, /EDITABLE)
fileID = CW_FILESEL(query, /IMAGE_FILTER)
info = {retID: event.ID, indID: indID, titID: titID, fileID: fileID}
WIDGET_CONTROL, query, SET_UVALUE = info

IF verbose THEN MESSAGE,/INFORMATIONAL, $
  STRING(SYSTIME(1)-start, "Launching chooser", $
         FORMAT = "(F7.3,' seconds: ', A)")
  
; Fire off widget and collect result
WIDGET_CONTROL, query, /REALIZE
XMANAGER, 'ximview_rdimage', query

WIDGET_CONTROL, event.ID, GET_UVALUE = result
WIDGET_CONTROL, event.ID, SET_UVALUE = save

                            ; Interpret results here
IF SIZE(result,/TYPE) NE 8 THEN RETURN                
IF result.CANCEL THEN RETURN

name = result.TITLE[0]
index = FIX(result.INDEX)
image = READ_IMAGE(result.FILE, red, green, blue, IMAGE_INDEX = index)

imsize = SIZE(image)
trucol = imsize[0] EQ 3
dims = imsize[1:imsize[0]]
IF trucol THEN dims = dims[WHERE(dims NE 3)]

; Make nominal FITS header,
mkhdr, header, imsize[imsize[0]+1], imsize[1:imsize[0]]

; Parameters usually set by parse_input:

; Strip directory info from filename (if any)
sep = PATH_SEP()
IF sep EQ '\' THEN sep = '\\'
firstchar = STRSPLIT(result.FILE, '.*'+sep, /REGEX)
path = STRMID(result.FILE, 0, firstchar-1)
file = STRMID(result.FILE, firstchar)

instrument = ''
tiff = QUERY_TIFF(result.FILE, info, GEOTIFF = geotiff)
IF tiff THEN BEGIN
    instrument = info.DESCRIPTION
    rots = [0,7,2,5,0,1,6,3,4]
    rotdir = rots[info.ORIENTATION]
    tn = TAG_NAMES(info)
    len = MAX(strlen(tn))
    fmt = "(A-"+STRTRIM(STRING(len),2)+",'   ')"
    tn = STRING(tn, FORMAT=fmt)
    SXADDHIST, 'TIFF info records:', header, location='END'
    FOR itag = 0,N_ELEMENTS(tn)-1 DO BEGIN
        item = STRING(tn[itag], info.(itag))
        SXADDHIST, item, header, location='END'
    ENDFOR
    IF SIZE(geotiff,/TYPE) EQ 8 THEN BEGIN
        tn = TAG_NAMES(geotiff)
        len = MAX(strlen(tn))
        fmt = "(A-"+STRTRIM(STRING(len),2)+",'   ')"
        tn = STRING(tn, FORMAT=fmt)
        SXADDHIST, 'GEOTIFF records:', header, location='END'
        FOR itag = 0,N_ELEMENTS(tn)-1 DO BEGIN
            item = STRING(tn[itag], geotiff.(itag))
            SXADDHIST, item, header, location='END'
        ENDFOR
    ENDIF 
ENDIF ELSE rotdir = 0

howto = 1
tablab = 'Map 1'
checksum32, BYTE(header), checksum
namestr = {path: path, file: file, telescope: '', instrument: instrument, $
           creator: '', object: '', freq: '', stokes: '', $
           extension: 0, extension_type: 'not FITS', spectral_type: '', $
           header: PTR_NEW(), checksum: checksum, name: ''}
newname = namestr.FILE + ': '
FOR i=2,7 DO newname += namestr.(i)

WIDGET_CONTROL, state.TABS, GET_UVALUE = mode

                                ; Turn off blinking, make tab 0
                                ; current, swap graphics state:
prep_screen, state, mode, tab_arr, old_graph
WIDGET_CONTROL, /HOURGLASS

ntab = N_ELEMENTS(tab_arr)
column = 1
str = tab_arr[0]
tabptr = state.TABARR ; pointer, as opposed to tab_arr = actual array.

lutstr  = *str.LUT
scaled = imsize[imsize[0]+1] EQ 1 ; Bytes already
lut_set = N_ELEMENTS(red) NE 0
IF lut_set THEN BEGIN
    lutstr.R = red  & lutstr.G = green  &  lutstr.B = blue
ENDIF

howto = 2      ; we will move data to the heap below.
temporary = 1B ; image data can be deleted
polcodes  = 0  ; no polarization info
freqcodes = 0  ; no frequency info

IF trucol THEN BEGIN
    rgbptr = PTRARR(3)
    inter = WHERE(imsize[1:3] EQ 3)
    CASE inter[0] OF
        0: FOR i=0,2 DO rgbptr[i] = PTR_NEW( REFORM(image[i,*,*]))
        1: FOR i=0,2 DO rgbptr[i] = PTR_NEW( REFORM(image[*,i,*]))
        2: FOR i=0,2 DO rgbptr[i] = PTR_NEW( REFORM(image[*,*,i]))
        ELSE: MESSAGE, 'Malformed array for true-color image'
    ENDCASE
    image = 0
                             ; TIFF files may have non-standard orientation
    IF rotdir NE 0 THEN $
      FOR i=0,2 DO *rgbptr[i] = ROTATE(TEMPORARY(*rgbptr[i]),rotdir)
    iptr = rgbptr[0]
ENDIF ELSE BEGIN
    IF rotdir NE 0 THEN image = ROTATE(TEMPORARY(image),rotdir)
    iptr = PTR_NEW(TEMPORARY(image))
ENDELSE

make_tabs, dummy, '', column, 0, 0, name, temporary, state, $
           ntab, iptr, header, newname, howto, tablab, namestr, polcodes, $
           freqcodes, ncol, T, line, extradims, mismatch

IF mismatch THEN RETURN

str = (*tabptr)[ntab]
IF trucol THEN BEGIN
    str.DECOMPOSED = 1
    str.ABSRANGE = [0,255]
    str.RGB = rgbptr
    str.UNIT = '  R   G   B'
    str.COLLAB = ['', '', '']
    str.IM_PTR   = PTR_NEW() ; assigned to rgbptr by make_tabs
    str.BYTE_PTR = PTR_NEW() ;
    ramp = LINDGEN(255)
    lutstr = {r: ramp, g: ramp, b: ramp, line: !P.color, absent: !P.background}
ENDIF ELSE BEGIN
    IF scaled THEN str.ABSRANGE = [0B,255B] $ ; BYTE_PTR assigned in make_tabs 
    ELSE BEGIN
        range = 'Full'
        *str.BYTE_PTR = scale_image(*iptr, range, ABS_RANGE = ar)
        str.ABSRANGE = ar
    ENDELSE
    str.UNIT  = ''
ENDELSE
str.RANGE = str.ABSRANGE
*str.LUT = lutstr

(*tabptr)[ntab] = str

title = state.TITLE
title_string = title.HEAD + form_unit(str.UNIT) + title.TAIL
WIDGET_CONTROL, state.READLAB, SET_VALUE = title_string
    
WIDGET_CONTROL, state.LABEL, SET_VALUE = topname
   
start_gscroll, ncol, state, T, ntab, mode, start

                                ; Draw initial screens:
fill_screens, ncol, ntab, tabptr, mode, state, start

; Enable blinking if possible
IF ncol + ntab GE 2 THEN BEGIN
    WIDGET_CONTROL, state.BLINK,  /SENSITIVE
    WIDGET_CONTROL, state.FRAMES, /SENSITIVE
    *state.BLINK_SEQ = INDGEN(ncol+ntab)
ENDIF

state.FIRST = 0B  

IF first THEN BEGIN
   start_msg = ['Click cursor to define centre of zoom/scroll, then:',$
                '      Mouse button 1 to drag image', $
                '                   2 to mark point and record value', $
                '      zoom with scroll wheel or buttons at screen left', '', $
                'Maxfit and Imstats buttons work around last marked point.']
   IF ~!ximview_expert THEN $
      ok = DIALOG_MESSAGE(start_msg,/INFORMATION,TITLE='Ximview', /CENTER)

    line = [line, title_string]
    PRINTF, state.LOGLUN, line, FORMAT="(A)"
    PRINT, line, FORMAT="(A)"
                                ; Turn on widget events
    str = (*state.TABARR)[0]
    WIDGET_CONTROL, str.DRAW, /SENSITIVE
    WIDGET_CONTROL, state.DISP, /SENSITIVE
    WIDGET_CONTROL, state.OV_BUTTON, /SENSITIVE
    WIDGET_CONTROL, state.IMSTATS, /SENSITIVE
    WIDGET_CONTROL, state.PEAKFIT, /SENSITIVE
ENDIF

restore_lut, *state.XIM_GRAPH, old_graph

WIDGET_CONTROL, event.TOP, SET_UVALUE = state
WIDGET_CONTROL, state.TABS,  SET_UVALUE = mode
WIDGET_CONTROL, state.TABS, /SENSITIVE

END
