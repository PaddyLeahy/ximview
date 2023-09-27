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
; Module XIMVIEW_POL
;
; J. P. Leahy March 2008
;
; Contains the event handlers for polarization display + helper
; routines:
;
; hsv_label:          Draws a colour label for the plot.
; set_pol_chan_event: Handles event from the set-up widget
; set_pol_chan:       Creates set-up widget and interprets results
; ximview_pol:        Creates polarization display
;
PRO hsv_label, str
;
; Draws a scale for a polarization hue-value plot (or HSV, but
; saturation is not represented).
;
; Inputs:
;     str:    Structure describing tab. Range parameter should
;             describe value scale. Hue scale assumed 0 to 180 deg.
;     badcol: colour index for bad values
;
COMPILE_OPT IDL2, HIDDEN
COMMON gr_global

WSET, str.SCALE_INDEX

absent = !P.background

DEVICE, /DECOMPOSED
ERASE, COLOR = absent*256^2+absent*256+absent

nxpix = !D.x_vsize - 60  &  nypix = !D.y_vsize - 18

chi   = REBIN(360.0 * FINDGEN(nxpix) / (nxpix-1), nxpix, nypix)
value = REBIN(REFORM(FINDGEN(nypix)/(nypix-1), 1, nypix, /OVERWRITE), $
              nxpix, nypix)

HSV2RGB, chi, 1.0, value, badcol, scaleptr, /HVONLY

x0 = 45  &  x1 = x0 + nxpix - 1
y0 = 15  &  y1 = y0 + nypix - 1
TV, [[[*scaleptr[0]]],[[*scaleptr[1]]],[[*scaleptr[2]]]], x0, y0, $
  TRUE = 3
PTR_FREE, scaleptr
                                ; label the scale
position = CONVERT_COORD([x0, x1], [y0, y1], /DEVICE, /TO_NORMAL)

xold = !X  &  yold = !Y
!Y.window = [position[1,0], position[1,1]]
!X.window = [position[0,0], position[0,1]]

r1 = str.RANGE[1]
ypow = 10.0^FLOOR( ALOG10(r1) )
yint = FLOOR(r1 / ypow) * ypow
xtl = 5.0 / nypix  &  ytl = 5.0 / nxpix
col = 256L^3 - 1
AXIS, XRANGE = [0,180], XSTYLE = 1, XTICKLEN = xtl, XTICKS = 4, $
  XTICKFORMAT = "(I0,'!9%!X')", COLOR = col
AXIS, YAXIS = 0, YRANGE = str.RANGE, YSTYLE = 1, YTICKLEN = ytl, $
  YMINOR = 1, YTICKINTERVAL = yint, COLOR = col
AXIS, XAXIS = 1, XTICKS = 1, COLOR = col

!X = xold  &  !Y = yold
DEVICE, DECOMPOSED = 0

END

PRO set_pol_chan_event, event
; Processes polarization channel dialog events
;
COMPILE_OPT IDL2, HIDDEN
ON_ERROR, 1

WIDGET_CONTROL, event.TOP, GET_UVALUE = info

tag = TAG_NAMES(event, /STRUCTURE_NAME)
CASE tag OF
    'WIDGET_BUTTON': BEGIN
        WIDGET_CONTROL, event.ID, GET_VALUE = label
        CASE label OF
            'Clear selection': BEGIN
                FOR i=0,2 DO BEGIN
                    WIDGET_CONTROL, info.POLROW[i], SET_UVALUE = -1
                    WIDGET_CONTROL, info.PBASE[i], SET_BUTTON = 0
                ENDFOR
;                WIDGET_CONTROL, info.IMINID, SET_VALUE = ''
                WIDGET_CONTROL, info.IBASE,  SENSITIVE = 0B
                WIDGET_CONTROL, info.PIBASE,  SENSITIVE = 1B
                RETURN
            END
            'Accept': BEGIN       ; Extract results from widget
                polchan = INTARR(3)
                FOR i=0,2 DO BEGIN
                    WIDGET_CONTROL, info.POLROW[i], GET_UVALUE = tab
                    polchan[i] = tab
                ENDFOR
                poltype = WIDGET_INFO(info.TYPEID, /DROPLIST_SELECT)
                WIDGET_CONTROL, info.PIMID, GET_VALUE = pimax
                pimax = FLOAT(pimax[0])
                WIDGET_CONTROL, info.FPMID, GET_VALUE = fpmax
                fpmax = FLOAT(fpmax[0])
                imin = 0
                result = {cancel: 0B, polchan: polchan, poltype: poltype, $
                          pimax: pimax, fpmax: fpmax, imin: imin}
            END
            'Cancel': result = {cancel: 1B}
        ENDCASE
        WIDGET_CONTROL, info.RETID, SET_UVALUE = result
        WIDGET_CONTROL, event.TOP, /DESTROY
    END
    'WIDGET_DROPLIST': BEGIN
        poltype = event.INDEX
        pol1lab = ['Stokes Q:     ', 'Pol intensity:', 'Pol fraction: ']
        pol2lab = ['Stokes U:     ', 'Pol angle:    ', 'Pol angle:    ']
        WIDGET_CONTROL, info.POLLAB[0], SET_VALUE = pol1lab[poltype]
        WIDGET_CONTROL, info.POLLAB[1], SET_VALUE = pol2lab[poltype]
        WIDGET_CONTROL, info.PIMID,  SET_VALUE = ''
        WIDGET_CONTROL, info.PIBASE, SENSITIVE = poltype EQ 0
    END
    '': BEGIN
        WIDGET_CONTROL, event.ID, SET_UVALUE = event.VALUE
        IF event.ID EQ info.POLROW[2] THEN BEGIN
            WIDGET_CONTROL, info.PIBASE, SENSITIVE = 0
            WIDGET_CONTROL, info.PIMID, SET_VALUE = ''
            WIDGET_CONTROL, info.IBASE,  /SENSITIVE
            sig3 = STRTRIM(STRING(3.0*info.SDEV),2)
;            WIDGET_CONTROL, info.IMINID, SET_VALUE = sig3
        ENDIF
    END
    'WIDGET_TEXT_CH': ; no action needed
    ELSE: MESSAGE, /INFORMATIONAL, 'Unexpected tag '+tag
ENDCASE
END

PRO set_pol_chan, event, tabarr, tablab, cancel, poltype, polchan, pimax, fpmax
; Launches modal dialog to get the image channels to use for polarization
;
COMPILE_OPT IDL2, HIDDEN
ON_ERROR, 1

ntab = N_ELEMENTS(tablab)
polcode = tabarr.POLCODE
sdev    = tabarr.SDEV
                                ; Try to identify Stokes channels
vals = [(WHERE(polcode EQ 2))[0], (WHERE(polcode EQ 3))[0], $
       (WHERE(polcode EQ 1))[0]]
extra = WHERE(vals EQ -1, nextra, COMPLEMENT = good)
IF nextra GT 0 THEN BEGIN
    nloose = ntab - 3 + nextra
    n2set = nloose < nextra
    IF n2set GT 0 THEN BEGIN
        loose = LINDGEN(ntab)
        FOR i=0,2-nextra DO loose = loose[WHERE(loose NE vals[good[i]])]
        vals[extra[0:n2set-1]] = loose[0:n2set-1]
    ENDIF
ENDIF
                                ; Set up widget
query = WIDGET_BASE(GROUP_LEADER = event.TOP, /MODAL, $
                    /COLUMN, TITLE = 'Choose Polarization channels')

polopt = ['Stokes Q and U','Polarized intensity and angle', $
              'Polarized fraction and angle']
typeID = WIDGET_DROPLIST(query, VALUE = polopt, $
                         TITLE =  'Polarization stored as:')

void = WIDGET_LABEL(query, VALUE = 'Tab', /ALIGN_CENTER)
polstr = ['Stokes Q:     ','Stokes U:     ','Stokes I:     ']
pollab = LONARR(3)  &  polrow = LONARR(3)  &  pbase  = LONARR(3)
FOR i=0,2 DO BEGIN
    val = vals[i]
    base = WIDGET_BASE(query, /ROW)
    pollab[i] = WIDGET_LABEL(base, VALUE = polstr[i], /ALIGN_LEFT)
    polrow[i] = val GE 0 $
      ? CW_BGROUP(base, STRTRIM(tablab,2), UVALUE = val, $
                  /ROW, /EXCLUSIVE, IDS = ids, SET_VALUE = val) $
      : CW_BGROUP(base, STRTRIM(tablab,2), UVALUE = val, $
                  /ROW, /EXCLUSIVE, IDS = ids)
    pbase[i] = WIDGET_INFO(ids[0], /PARENT)
ENDFOR

pibase = WIDGET_BASE(query, /ROW)
void = WIDGET_LABEL(pibase, VALUE= 'Max for pol intensity range:', /ALIGN_LEFT)
pimID = WIDGET_TEXT(pibase, VALUE = '', XSIZE = 20, /EDITABLE)
ibase = WIDGET_BASE(query, /ROW)
void = WIDGET_LABEL(ibase, VALUE = 'Max for frac. pol range:    ', /ALIGN_LEFT)
sig3 = vals[2] NE -1 ? STRTRIM(STRING(3.0*sdev[vals[2]]),2) : ''
fpmID = WIDGET_TEXT(ibase, VALUE = '1.0', XSIZE = 20, /EDITABLE)
IF vals[2] EQ -1 THEN WIDGET_CONTROL, ibase, SENSITIVE = 0 $
ELSE WIDGET_CONTROL, pibase, SENSITIVE = 0

iminID = 0 ; for later development

bbase  = WIDGET_BASE(query, /ROW)
void   = WIDGET_BUTTON(bbase, VALUE = 'Clear selection')
void   = WIDGET_BUTTON(bbase, VALUE = 'Accept')
void   = WIDGET_BUTTON(bbase, VALUE = 'Cancel')

                                ; Structure with return ID, & internal
                                ; IDs
info = {retID: event.ID, typeID: typeID, pollab: pollab, polrow: polrow, $
        pbase: pbase, pimID: pimID, fpmID: fpmID, pibase: pibase, $
        ibase: ibase, iminID: iminID, sdev: sdev}

WIDGET_CONTROL, query, SET_UVALUE = info, /NO_COPY
WIDGET_CONTROL, typeID, SET_DROPLIST_SELECT = 0
WIDGET_CONTROL, event.ID, GET_UVALUE = save, /NO_COPY

; Fire off widget and collect result
WIDGET_CONTROL, query, /REALIZE
XMANAGER, 'set_pol_chan', query

WIDGET_CONTROL, event.ID, GET_UVALUE = result, /NO_COPY
WIDGET_CONTROL, event.ID, SET_UVALUE = save, /NO_COPY

IF SIZE(result,/TYPE) NE 8 THEN cancel = 1B ELSE cancel = result.CANCEL
IF cancel THEN RETURN

polchan = result.POLCHAN
IF polchan[0] EQ -1 OR polchan[1] EQ -1 THEN BEGIN
    MESSAGE, /INFORMATIONAL, $
      'First two polarization channels required to make a display...returning'
    cancel = 1B
    RETURN
ENDIF

poltype = result.POLTYPE
pimax   = result.PIMAX
fpmax   = result.FPMAX
imin    = result.IMIN
                                ; Check that selected channels are
                                ; meaningful
freq0 = tabarr[polchan[0]].FREQCODE
scodes = ['YX','XY','YY','XX', 'LR', 'RL', 'LL','RR','','I','Q','U','V', $
          'PPOL', 'FPOL', 'PANG']
chanstr = ['1','2','3']
typecode = [[2,3,1],[5,7,1],[6,7,1]]
bbad = 0B
FOR i=0,2 DO IF polchan[i] NE -1 THEN BEGIN
    testcode = polcode[polchan[i]]
    bad = testcode NE typecode[i,poltype] && testcode NE 0
    IF bad THEN MESSAGE, /INFORMATIONAL, 'Channel' + chanstr[i] + $
      ': polarization: ' + scodes[typecode[i,poltype]+8] + ' expected but ' $
      + scodes[testcode+8] + ' found.'
    bbad = bbad || bad
    freq = tabarr[polchan[i]].FREQCODE
    bad = freq NE freq0
    IF bad THEN MESSAGE, /INFORMATIONAL, 'Frequency mismatch: channel ' $
      + chanstr[i] + ': ' + freq + '; channel 1: ' + freq0
    bbad = bbad || bad
ENDIF
IF bbad THEN BEGIN
    MESSAGE, /INFORMATIONAL, $
      'Revise channels, or set image polarization/frequency via Analysis menu'
    cancel = 1B
ENDIF

END

PRO ximview_pol, event
; Sets options for polarization display
;
COMPILE_OPT IDL2, HIDDEN
COMMON gr_global
ON_ERROR, 1

start = SYSTIME(1)
WIDGET_CONTROL, event.ID, GET_VALUE = label
WIDGET_CONTROL, event.TOP,  GET_UVALUE = state

                                ; Find tabs suitable as pol channels,
                                ; i.e. not blink or RGB:
tabarr = *state.TABARR
screens = tabarr.SCREEN
good    = WHERE( PTR_VALID( tabarr.BYTE_PTR ))
iscreen = BSORT(screens[good])
screens = screens[good[iscreen]]
tablab = get_tab_uvals(state.TABS)
tablab = tablab[screens]
fmt = STRING(MAX(STRLEN(tablab)), FORMAT = "('(A',I2,')')")
tablab = STRING(tablab, FORMAT = fmt)

is_mswin = STRCMP(!version.OS_FAMILY,'Windows', 4,/ FOLD_CASE)
IF is_mswin THEN BEGIN ; Set a fixed-width font:
    old_font = WIDGET_INFO(state.LABEL, /FONTNAME)
    font = 'lucida console*10'
    WIDGET_CONTROL, DEFAULT_FONT = font
ENDIF

set_pol_chan, event, tabarr[good[iscreen]], tablab, $
  cancel, poltype, polchan, pimax, fpmax

; Restore system default (?)
IF is_mswin THEN WIDGET_CONTROL, DEFAULT_FONT = old_font

IF cancel  || ARRAY_EQUAL(polchan,[-1,-1,-1]) THEN RETURN

WIDGET_CONTROL, state.TABS, GET_UVALUE = mode
prep_screen, state, mode, tabarr, oldgraph

good = WHERE( PTR_VALID( tabarr.BYTE_PTR))  ; after shift

WIDGET_CONTROL, /HOURGLASS
ntab  = N_ELEMENTS(tabarr)
collab = STRARR(3)
do_inten = polchan[2] NE -1
                                ; Look out for gores
sky = WHERE(*tabarr[good[0]].BYTE_PTR NE !P.background)

FOR i=0,1 DO BEGIN
    ptr = tabarr[good[polchan[i]]].IM_PTR
    IF PTR_VALID(ptr) THEN temp = (*ptr)[sky] ELSE BEGIN
        WIDGET_CONTROL, state.LABEL, GET_UVALUE=image, /NO_COPY
        temp = image[sky]
        WIDGET_CONTROL, state.LABEL, SET_UVALUE=image, /NO_COPY
    ENDELSE
    IF i EQ 0 THEN pol1 = TEMPORARY(temp) ELSE pol2 = TEMPORARY(temp)
ENDFOR

bptr_istemp = 0B

CASE label OF
    'Colour': BEGIN
        CASE poltype OF
            0: BEGIN ; Stokes Q & U
                collab[0] = 'Stokes Q: ' + tablab[polchan[0]]
                collab[1] = 'Stokes U: ' + tablab[polchan[1]]
                collab[2] = ''
                chi = ATAN( pol2, pol1)
                chi *= !radeg
                                ; Deal with negative values, if any
                neg = WHERE(chi LT 0., Nneg)
                IF Nneg GT 0 THEN chi[neg] += 360.0
                pol1 ^= 2
                pol2 ^= 2
                ppol  = TEMPORARY(pol1)
                ppol += pol2
                pol2 = 0
                ppol  = SQRT(TEMPORARY(ppol))
                maxpp = MAX(ppol,/NAN)
                minpp = 0.0
                pimax = pimax EQ 0.0 ? maxpp : MIN([maxpp,pimax])
                ppol <= pimax
            END
            1: BEGIN ; Pol intensity and angle
                collab[0] = 'Pol Int:  ' + tablab[polchan[0]]
                collab[1] = 'Pol Angle:' + tablab[polchan[1]]
                collab[2] = ''
                ppol = TEMPORARY(pol1)
                bptr = tabarr[good[polchan[0]]].BYTE_PTR
                chi  = TEMPORARY(pol2)
            END
            2: BEGIN ; Pol fraction and angle
                collab[0] = 'Frac Pol: ' + tablab[polchan[0]]
                collab[1] = 'Pol Angle:' + tablab[polchan[1]]
                collab[2] = ''
                ppol = TEMPORARY(pol1)
                bptr = tabarr[good[polchan[0]]].BYTE_PTR
                chi  = TEMPORARY(pol2)
            END
        ENDCASE
        IF do_inten THEN BEGIN
            new_label = 'HSV Pol'
            collab[2] = 'Stokes I: ' + tablab[polchan[2]]
            bptr = tabarr[good[polchan[2]]].BYTE_PTR
            IF poltype NE 2 THEN BEGIN
                iptr = tabarr[good[polchan[2]]].IM_PTR
                IF PTR_VALID(iptr) THEN temp = (*iptr)[sky] ELSE BEGIN
                    WIDGET_CONTROL, state.LABEL, GET_UVALUE=image, /NO_COPY
                    temp = image[sky]
                    WIDGET_CONTROL, state.LABEL, SET_UVALUE=image, /NO_COPY
                ENDELSE
                ppol /= temp
                bad = WHERE(~FINITE(ppol) OR ppol LT 0.0)
                IF bad[0] NE -1 THEN ppol[bad] = fpmax
                ppol <= fpmax
            ENDIF
            HSV2RGB, chi, ppol, (*bptr)[sky], badcol, rgbptr, satmin, satmax

            MESSAGE, /INFORMATIONAL, 'Saturation scale:'
            MESSAGE, /INFORMATIONAL, 'from fractional polarization of ' + $
              STRTRIM(STRING(satmin),2) + ' (grey)'
            MESSAGE, /INFORMATIONAL, 'to   fractional polarization of ' + $
              STRTRIM(STRING(satmax),2) + ' (saturated colours)'

            range = tabarr[good[polchan[2]]].RANGE
        ENDIF ELSE BEGIN
            new_label = 'HV Pol'
            IF poltype EQ 0 THEN BEGIN
                range = [0.0, pimax]
                scale = 254/pimax
                ppol *= scale
                ppol = BYTE(ppol)
                lift = WHERE(ppol GE badcol)
                ppol[lift] += 1B
                lift = 0
                bptr = PTR_NEW(*tabarr[good[0]].BYTE_PTR)
                (*bptr)[sky] = ppol
                bptr_istemp = 1B
            ENDIF ELSE range = tabarr[good[polchan[0]]].RANGE
            HSV2RGB, chi, 1.0, (*bptr)[sky], badcol, rgbptr, /HVONLY
        ENDELSE
                                ; Replace gores
        temp = *tabarr[good[0]].BYTE_PTR
        FOR i=0,2 DO BEGIN
            temp[sky] = *rgbptr[i]
            *rgbptr[i] = temp
        ENDFOR
        temp = 0

        str = tabarr[0]
        str.TRFUNC = 0 ; result.TRFUNC
        str.WRAP   = 0 ; result.WRAP
        str.ZERO   = 0.0
        str.BETA   = 0.0; result.BETA
        str.RANGE  = range
        str.TEMPORARY = 1B

        collab[0] = 'HSV'

        make_rgb_tab, str, ntab, rgbptr, collab, new_label, state

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
        
    END
;    'Vector:'
;    'LIC':
    ELSE: BEGIN
        MESSAGE, /INFORMATIONAL, 'Polarization option ' + label + $
          ' is not yet available'
                                ; Restore current tab to the one now showing
        current = WIDGET_INFO(state.TABS,/TAB_CURRENT)
        itab = WHERE(tabarr.SCREEN EQ current)
        WSET, tabarr[itab].WINDOW
        gscroll_newscreen, itab, tabarr, mode.ZOOM_FACTOR, mode.X_CENTRE, $
          mode.Y_CENTRE, mode.XHALF, mode.YHALF, done, mode.OVERVIEW

        *state.TABARR = tabarr
    END
ENDCASE

IF bptr_istemp THEN PTR_FREE, bptr

restore_lut, *state.XIM_GRAPH, oldgraph
 
WIDGET_CONTROL, state.TABS, /SENSITIVE

END
