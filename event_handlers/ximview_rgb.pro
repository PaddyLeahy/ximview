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
; Module XIMVIEW_RGB
;
; J. P. Leahy February 2008 Revised 2013
;
; Contains the event handlers & helper routines for 3-channel colour
; in Ximview.
;
PRO rgb_setup_event, event
; Handles modal widget setting rgb properties
;
COMPILE_OPT IDL2, HIDDEN
ON_ERROR, 1

WIDGET_CONTROL, event.TOP, GET_UVALUE = info

tag = TAG_NAMES(event, /STRUCTURE_NAME)
CASE STRMID(tag,0,12) OF
    'WIDGET_BUTTO': BEGIN
        WIDGET_CONTROL, event.ID, GET_VALUE = label
        ntab = N_ELEMENTS(info.RGBROW)
        CASE label OF
            'Clear selection': BEGIN
                FOR i=0,2 DO BEGIN
                    WIDGET_CONTROL, info.RGBROW[i], SET_UVALUE = -1
                    WIDGET_CONTROL, info.PBASE[i], SET_BUTTON = 0
                ENDFOR

                RETURN
            END
            'Done': BEGIN       ; Extract results from widget
                rgb = INTARR(3)
                FOR i=0,2 DO BEGIN
                    WIDGET_CONTROL, info.RGBROW[i], GET_UVALUE = tab
                    rgb[i] = tab
                ENDFOR
                WIDGET_CONTROL, info.CONVID, GET_VALUE = conv
                result = {cancel: 0B, rgb: rgb, convention: conv}
            END
            'Cancel': result = {cancel: 1B}
        ENDCASE
        WIDGET_CONTROL, info.RETID, SET_UVALUE = result
        WIDGET_CONTROL, event.TOP, /DESTROY
    END
    '': WIDGET_CONTROL, event.ID, SET_UVALUE = event.VALUE
ENDCASE

END

PRO rgb_setup, event, ntab, tablab, rgb, convention, cancel
;
; Launches modal widget to get parameters and returns results.
;
; Inputs:
;   event:  Event structure containing Widget IDs used for passing
;           results.
;   ntab:   Number of tabs (excluding blink tab if any)
;   tablab: Labels for each tab.
;
; Output:
;   rgb:        Tabs to be used for red, green, blue channels
;   convention: Convention for colour scaling
;   cancel:     True if operation aborted.
;
COMPILE_OPT IDL2, HIDDEN

WIDGET_CONTROL, event.ID, GET_UVALUE = save, /NO_COPY

; Set up widget
rgbrow = LONARR(3)
rgbstrs = ['Red  ', 'Green', 'Blue ']
vals = INDGEN(3)
pbase = LONARR(3)

query = WIDGET_BASE(GROUP_LEADER = event.TOP, /MODAL, $
                    /COLUMN, TITLE = 'Choose RGB channels')
void = WIDGET_LABEL(query, VALUE = 'Tab', /ALIGN_CENTER)

vals = INDGEN(3)
IF ntab LT 3 THEN vals[ntab:*] = -1
FOR i=0,2 DO BEGIN
    val = vals[i]
    rgbrow[i] = val NE -1 $
      ? CW_BGROUP(query, tablab[0:ntab-1], UVALUE = val, /ROW, /EXCLUSIVE, $
                  IDS = ids, SET_VALUE = val, LABEL_LEFT = rgbstrs[i]) $
      : CW_BGROUP(query, tablab[0:ntab-1], UVALUE = val, /ROW, /EXCLUSIVE, $
                  IDS = ids, LABEL_LEFT = rgbstrs[i])
    pbase[i] = WIDGET_INFO(ids[0], /PARENT)
ENDFOR

conventions = ['Saturate to white','Preserve colour for saturated pixels']
convID = CW_BGROUP(query, conventions, SET_VALUE = 0, /EXCLUSIVE, /COLUMN, $
                   LABEL_LEFT = 'Colour scaling:')
bbase = WIDGET_BASE(query, /ROW)
void = WIDGET_BUTTON(bbase, VALUE = 'Clear selection')
void = WIDGET_BUTTON(bbase, VALUE = 'Done')
void = WIDGET_BUTTON(bbase, VALUE = 'Cancel')

                                ; Structure with return ID, internal
                                ; IDs, and defaults
info = {retid: event.ID, rgbrow: rgbrow, pbase: pbase, convID: convID, $
        rgb: vals}
WIDGET_CONTROL, query, SET_UVALUE = info, /NO_COPY

WIDGET_CONTROL, query, /REALIZE
XMANAGER, 'rgb_setup', query

WIDGET_CONTROL, event.ID, GET_UVALUE = result, /NO_COPY
WIDGET_CONTROL, event.ID, SET_UVALUE = save, /NO_COPY

IF size(result,/TYPE) NE 8 THEN cancel = 1B ELSE cancel = result.CANCEL
IF cancel THEN RETURN

rgb = result.RGB
convention = result.CONVENTION
cancel = MAX(rgb) EQ -1

RETURN

END

PRO rgb_scale_event, event
;
; Interprets events from dialog launched by rgb_scale for colour-preserving
; case
;
COMPILE_OPT IDL2, HIDDEN
ON_ERROR, 1

WIDGET_CONTROL, event.TOP, GET_UVALUE = info

tag = TAG_NAMES(event, /STRUCTURE_NAME)
CASE tag OF
    'WIDGET_BUTTON': BEGIN
        WIDGET_CONTROL, event.ID, GET_VALUE = label
        CASE label OF
            'Accept': BEGIN       ; Extract results from widget
                trfunc = WIDGET_INFO(info.FUNCID, /DROPLIST_SELECT)
                WIDGET_CONTROL, info.WRAPID, GET_VALUE = wrap
                wrap = ([0, -1, 1])[wrap]
                WIDGET_CONTROL, info.ADJUSTID,  GET_VALUE = adjust
                adjust = FLOAT(STRSPLIT(adjust,', ',/EXTRACT, COUNT = count))
                IF count NE info.NCHAN THEN BEGIN
                    ncstr = STRTRIM(STRING(info.NCHAN),2)
                    ok = DIALOG_MESSAGE( /INFORMATION, $
                                         DIALOG_PARENT = event.TOP, $
               'Please enter exactly ' + ncstr + ' scaling adjustment numbers')
                    RETURN
                ENDIF
                WIDGET_CONTROL, info.RTEXT,  GET_VALUE = range
                WIDGET_CONTROL, info.BTEXT,  GET_VALUE = beta
                range = FLOAT(STRSPLIT( range,', ',/EXTRACT, COUNT = count))
                CASE count OF
                    1:  range = [0.0, range]
                    2:          ; should be OK
                    ELSE: BEGIN
                        ok = DIALOG_MESSAGE( /INFORMATION, $
                                             DIALOG_PARENT = event.TOP, $
                             'Enter 1 or 2 numbers for intensity range')
                        RETURN
                    END
                ENDCASE
                result = {cancel: 0B, adjust: adjust, trfunc: trfunc, $
                          range: range, wrap: wrap, beta: FLOAT(beta[0])}
            END
            'Cancel': result = {cancel: 1B}
        ENDCASE
        WIDGET_CONTROL, info.RETID, SET_UVALUE = result
        WIDGET_CONTROL, event.TOP, /DESTROY
    END
    '': WIDGET_CONTROL, event.ID, SET_UVALUE = event.VALUE
    ELSE: ; Droplist or text events need no action.
ENDCASE

END
PRO rgb_scale, event, tabarr, rgbtab, collab, result, cancel
;
; Launches a dialog widget to get info about scaling of the "averaged"
; channel for the colour-preserving case of RGB imaging
;
; Inputs:
;     event:  The event structure from ximview_rgb
;     tabarr: The array of structures describing each tab
;     rgbtab: index to tabarr for the RGB channels
;     collab: Label for each RGB channel
; Outputs:
;     result: Structure containing scaling parameters set by dialog
;

query = WIDGET_BASE(GROUP_LEADER = event.TOP, /MODAL, $
                    /COLUMN, TITLE = 'Set RGB Scaling')
void = WIDGET_LABEL(query, VALUE = $
                    'Data ranges after subtracting currently-set zero levels:')
rounding = FLTARR(2,3)
unit = ''
FOR i=0,2 DO IF rgbtab[i] GE 0 THEN BEGIN
    tab = tabarr[rgbtab[i]]
    IF unit EQ '' THEN BEGIN
        unit = tab.UNIT ; unit for user input below
        mult = tab.MULT
    ENDIF
    range = (tab.ABSRANGE - tab.ZERO)*tab.MULT
    range_str = STRING(range, FORMAT = "(2F8.2)")
    test = FLOAT(STRSPLIT( range_str,' ',/EXTRACT))
    rounding[0,i] = range - test/tab.MULT
    label = collab[i] + STRING(range, tab.UNIT, FORMAT = $
                               "(F8.2, ' to', F8.2, 1X, A)")
    void = WIDGET_LABEL(query, VALUE = label, /ALIGN_LEFT)
ENDIF

base = WIDGET_BASE(query, /ROW)
void = WIDGET_LABEL(base, VALUE='Rescale factors for each channel:')
void = WHERE (rgbtab NE -1, nchan)
adjstr = STRJOIN(STRTRIM(STRING(REPLICATE(1,nchan)), 2),',')
adjustID = WIDGET_TEXT(base, /EDITABLE, VALUE = adjstr, XSIZE = 16)

void = WIDGET_LABEL(query, VALUE = $
                    'Choose scaling to apply to colour-average image:')
funcID = WIDGET_DROPLIST(query, VALUE = scale_funs(-1), $
                         TITLE =  'Transfer function:')
wrapID = CW_BGROUP(query, ['No','Bright end only','Full'], /EXCLUSIVE, $
                   SET_VALUE = 0, LABEL_LEFT = 'Wrap colour table?', $
                   /ROW, /RETURN_NAME)

base = WIDGET_BASE(query, /ROW)
void   = WIDGET_LABEL(base,  VALUE ='Data range for scaling:')
rtext  = WIDGET_TEXT(base, /EDITABLE, XSIZE = 20, VALUE = '')
void   = WIDGET_LABEL(base, VALUE = unit)
base = WIDGET_BASE(query, /ROW)
void   = WIDGET_LABEL(base,  VALUE = 'Beta for Asinh scaling:')
btext  = WIDGET_TEXT(base, /EDITABLE, XSIZE = 10, VALUE = '')

bbase  = WIDGET_BASE(query, /ROW)
void   = WIDGET_BUTTON(bbase, VALUE = 'Accept')
void   = WIDGET_BUTTON(bbase, VALUE = 'Cancel')
                                ; Structure with return ID, & internal
                                ; IDs
info = {retid: event.ID, funcID: funcID, wrapID: wrapID, nchan: nchan, $
        adjustID: adjustID, rtext: rtext, btext: btext}

WIDGET_CONTROL, query, SET_UVALUE = info, /NO_COPY
WIDGET_CONTROL, funcID, SET_DROPLIST_SELECT = 0
WIDGET_CONTROL, event.ID, GET_UVALUE = save, /NO_COPY

; Fire off widget and collect result
WIDGET_CONTROL, query, /REALIZE
XMANAGER, 'rgb_scale', query

WIDGET_CONTROL, event.ID, GET_UVALUE = result, /NO_COPY
WIDGET_CONTROL, event.ID, SET_UVALUE = save, /NO_COPY

IF SIZE(result,/TYPE) NE 8 THEN cancel = 1B ELSE BEGIN
    cancel = result.CANCEL
    result.RANGE /= mult ; since we requested it in the matching units.
ENDELSE

END

PRO ximview_rgb, event
; Sets rgb options and creates a tab for a new RGB image
;
COMPILE_OPT IDL2, HIDDEN
COMMON gr_global
ON_ERROR, 1

start = SYSTIME(1)
WIDGET_CONTROL, event.TOP,  GET_UVALUE = state
WIDGET_CONTROL, state.TABS, GET_UVALUE = mode

                                ; Find tabs suitable as RGB channels,
                                ; i.e. not blink or RGB:
tabarr = *state.TABARR
screens = tabarr.SCREEN
good    = WHERE( PTR_VALID( tabarr.BYTE_PTR ), ntab)
iscreen = BSORT(screens[good])
screens = screens[good[iscreen]]
tablab = get_tab_uvals(state.TABS)
tablab = tablab[screens]
fmt = STRING(MAX(STRLEN(tablab)), FORMAT = "('(A',I2,')')")
tablab = STRING(tablab, FORMAT = fmt)
colstr = ['Red:   ','Green: ','Blue:  ']

                                ; Interpret requested action
button = 0B
tag = TAG_NAMES(event, /STRUCTURE_NAME)
CASE tag OF
    'WIDGET_BUTTON': BEGIN
        WIDGET_CONTROL, event.ID, GET_VALUE = label
        CASE label OF
            'Red-Green-Blue': BEGIN ; Launch modal widget to get info:

                is_mswin = STRCMP(!version.OS_FAMILY,'Windows', 4,/ FOLD_CASE)
                IF is_mswin THEN BEGIN ; Set a fixed-width font:
                    old_font = WIDGET_INFO(state.LABEL, /FONTNAME)
                    font = 'lucida console*10'
                    WIDGET_CONTROL, DEFAULT_FONT = font
                ENDIF

                rgb_setup, event, ntab, tablab, rgb, convention, cancel

                                ; Restore system default (?)
                IF is_mswin THEN WIDGET_CONTROL, DEFAULT_FONT = old_font

                IF cancel || ARRAY_EQUAL(rgb,[-1,-1,-1]) THEN RETURN
            END
            ELSE: MESSAGE, /INFORMATIONAL, 'Unknown button: ' + label
        ENDCASE
    END
    ELSE: MESSAGE, /INFORMATIONAL, 'Unknown event type received: ' + name
ENDCASE

prep_screen, state, mode, tabarr, oldgraph

rgbptr = PTRARR(3)
rgbtab = REPLICATE(-1, 3)
collab = STRARR(3)
FOR i=0,2 DO IF rgb[i] EQ -1 THEN rgbptr[i] = PTR_NEW() ELSE BEGIN
    rgbtab[i] = WHERE(screens[rgb[i]] EQ tabarr.SCREEN)
    rgbptr[i] = tabarr[rgbtab[i]].BYTE_PTR
    collab[i] = colstr[i] + tablab[rgb[i]]
ENDELSE

; Add entry to tab array
str = tabarr[rgbtab[0]]  ; ensures that it's not a funny tab.

lut = str.LUT
;
;          Following commented out because bad and absent will come out
;          at the right grey level naturally.
;
;bad_col = [(*lut).R[badcol], (*lut).G[badcol], (*lut).B[badcol]]
;
;abscol = [(*lut).R[!P.background], (*lut).G[!P.background], $
;          (*lut).B[!P.background]]
;                                ; Set colours for bad and absent pixels
;absent = WHERE(*valid EQ !P.background)
;IF absent[0] NE -1 THEN FOR i=0,2 DO IF rgbtab[i] NE -1 THEN $
;  (*rgbptr[i])[absent] = abscol[i]
;absent = 0

dims = state.IMSIZE[1:2]
nx = dims[0]  &  ny = dims[1]
npix = nx*ny
void = WHERE (rgbtab NE -1, nchan)

IF convention EQ 1 THEN BEGIN
; Lupton et al (2004, PASP) idea.
; We have to re-write byte images for the individual channels.

; Algorithm here does some scaling twice; could be sped up by 30% or
; so at the cost of larger memory overheads.
;
                                ; Launch scaling dialog
    IF is_mswin THEN WIDGET_CONTROL, DEFAULT_FONT = font
    rgb_scale, event, tabarr, rgbtab, collab, result, cancel
    IF is_mswin THEN WIDGET_CONTROL, DEFAULT_FONT = old_font
    IF cancel THEN BEGIN
        update_screen, state.TABARR, 0, mode, done
        GOTO, QUIT
    ENDIF

    WIDGET_CONTROL, /HOURGLASS
    IF state.VERBOSE THEN BEGIN
        MESSAGE, /INFORMATIONAL, 'Begin scaling byte channels'
        HELP,/MEMORY
    ENDIF

    adjust = [1,1,1]
    adjust[WHERE(rgbtab NE -1)] = result.ADJUST
    str.TRFUNC = result.TRFUNC
    str.WRAP   = result.WRAP
    str.ZERO   = 0.0
    str.BETA   = result.BETA
    str.RANGE  = result.RANGE

    r1 = str.RANGE[0]
 
    lchunk = 1024L^2 ; Loop overheads should be negligible for chunks this big
    IF npix LT lchunk THEN lchunk = npix
    nchunk = DIVUP(npix, lchunk)
    idc = LINDGEN(lchunk)
    last = npix - (nchunk-1L)*lchunk

    channel = PTRARR(3,/ALLOCATE_HEAP)

    rescale  = FLTARR(npix)
    maxval = 0.0
    FOR ichunk = 0, nchunk-1 DO BEGIN
        nidc = N_ELEMENTS(idc)
        sum = FLTARR(nidc)
        FOR i=0,2 DO IF rgbtab[i] NE -1 THEN BEGIN
            *channel[i] = FLTARR(nidc)
            ptr = tabarr[rgbtab[i]].IM_PTR
            IF PTR_VALID(ptr) THEN temp = (*ptr)[idc] ELSE BEGIN
                WIDGET_CONTROL, state.LABEL, GET_UVALUE=image, /NO_COPY
                temp = image[idc]
                WIDGET_CONTROL, state.LABEL, SET_UVALUE=image, /NO_COPY
            ENDELSE
            goodpix = WHERE(FINITE(temp), COMPLEMENT = badpix)
            IF goodpix[0] GE 0 THEN BEGIN
                temp = temp[goodpix]
                temp -= tabarr[rgbtab[i]].ZERO
                temp *= adjust[i]
                (*channel[i])[goodpix] = temp
                temp /= nchan
                sum[goodpix] += temp
            ENDIF
            IF badpix[0] GE 0 THEN sum[badpix] = !values.F_NAN
        ENDIF
        temp = 0  &  goodpix = 0  &  badpix = 0
        fsum = scale_image(sum, str.RANGE, str.WRAP, str.TRFUNC, $
                           str.ZERO, str.BETA)

        bad = WHERE(~FINITE(sum), nbad)
        IF nbad NE 0 THEN sum[bad]  = 255.0 ; Avoid div by bad

        sum -= r1
        zeros = WHERE(sum LE 0.0, nzero)
        IF nzero + nbad LT nidc THEN BEGIN
            IF nzero  NE 0 THEN sum[zeros] = 255.0 ; Avoid div by zero
            fsum = FLOAT(fsum) / sum
            FOR i=0,2 DO IF rgbtab[i] NE -1 THEN BEGIN
                *channel[i] -= r1
                *channel[i] *= fsum
                IF zeros[0] NE -1 THEN (*channel[i])[zeros] = 0.0
                maxval = MAX(*channel[i]) >  maxval
            ENDIF
        ENDIF
        IF nbad  NE 0 THEN fsum[bad]   = !values.F_NAN
        IF nzero NE 0 THEN fsum[zeros] = 0.
        rescale[idc] = fsum

        fsum = 0  &  zeros = 0  &  bad = 0

        IF ichunk LT nchunk - 2 THEN idc += lchunk ELSE $
          idc = LINDGEN(last) + (nchunk-1L)*lchunk
    ENDFOR
    PTR_FREE, channel
    idc = 0  &  sum = 0

    goodpix = WHERE(FINITE(rescale))
    IF goodpix[0] EQ -1 THEN BEGIN
        MESSAGE, /INFORMATIONAL, $
          'All pixels blanked on one or another colour channel: Abandoning RGB'
        WIDGET_CONTROL, state.TABS, /SENSITIVE
        RETURN
    ENDIF
    rescale[goodpix] *= (255. / maxval)
    FOR i = 0,2 DO IF rgbtab[i] NE -1 THEN BEGIN
        ptr = tabarr[rgbtab[i]].IM_PTR
        IF PTR_VALID(ptr) THEN temp = (*ptr) ELSE BEGIN
            WIDGET_CONTROL, state.LABEL, GET_UVALUE=image, /NO_COPY
            temp = image
            WIDGET_CONTROL, state.LABEL, SET_UVALUE=image, /NO_COPY
        ENDELSE
        goodpix = WHERE(FINITE(temp), ngood)
        IF goodpix[0] GE 0 THEN BEGIN
            temp = temp[goodpix]
            temp -= tabarr[rgbtab[i]].ZERO
            temp *= adjust[i]
            temp -= r1
            temp *= rescale[goodpix]
            bad = WHERE(~FINITE(temp))
            IF bad[0] NE -1 THEN temp[bad] = bad_col[i]
            temp >= 0.0
            (*rgbptr[i])[goodpix] = temp
        ENDIF
    ENDIF
    rescale = 0  &  temp = 0  &  goodpix = 0  &  bad = 0

    IF state.VERBOSE THEN BEGIN
        MESSAGE, /INFORMATIONAL, STRING(SYSTIME(1) - start, FORMAT = $
                                        "('Saved byte channels at',F7.2)")
        HELP, /MEMORY
    ENDIF
ENDIF ; ELSE BEGIN                ; End of no-saturate mode
;                                ; Set colours for bad pixels
;    borg = BYTARR(npix)
;    FOR i=0,2 DO IF rgbtab[i] NE -1 THEN $
;      borg AND= *rgbptr[i] EQ badcol
;    bad = WHERE(borg)
;    borg = 0
;    IF bad[0] NE -1 THEN FOR i=0,2 DO IF rgbtab[i] NE -1 THEN $
;      (*rgbptr[i])[bad] = bad_col[i]
;    bad = 0
;ENDELSE

str.TEMPORARY = 0B ; Don't delete byte arrays when tab is deleted.
make_rgb_tab, str, ntab, rgbptr, collab, 'RGB', state

ncol = 1 ; number of new tabs to make
start_gscroll, ncol, state, T, ntab, mode, start

tabarr = state.TABARR

                                ; Draw initial screens:
fill_screens, ncol, ntab, tabarr, mode, state, start

; Enable blinking with this tab
*state.BLINK_SEQ = INDGEN(ncol+ntab)

                                ;   Update units on readout label
mid = form_unit( (*tabarr)[0].UNIT)
title_string = state.TITLE.HEAD + mid + state.TITLE.TAIL
WIDGET_CONTROL, state.READLAB, SET_VALUE=title_string

WIDGET_CONTROL, state.TABS,  SET_UVALUE = mode

QUIT:

restore_lut, *state.XIM_GRAPH, oldgraph

WIDGET_CONTROL, state.TABS, /SENSITIVE


END
