; -----------------------------------------------------------------------------
;
;  Copyright (C) 2007-2020   J. P. Leahy
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
; Module XIMVIEW_SCALE
;
; Contains the event handlers and associated routines for changing the
; intensity scaling in Ximview.
;
PRO reset_scale_dialog, info, str, AUTO = auto
; Resets the scale dialog widget with data from a new screen, or set
; to defaults:
;
COMPILE_OPT IDL2, HIDDEN

full = ~KEYWORD_SET(auto)
IF SIZE(str,/TYPE) EQ 8 THEN BEGIN ; it's a structure
    IF full THEN BEGIN ; Don't set these if we've only autoscaled
        IF info.COLID NE -1 THEN $
          WIDGET_CONTROL, info.COLID, SET_DROPLIST_SELECT = str.COLTAB
        WIDGET_CONTROL, info.FUNCID, SET_DROPLIST_SELECT = str.TRFUNC
        WIDGET_CONTROL, info.WRAPID, SET_VALUE = ([1,0,2])[str.WRAP+1]
    ENDIF

    unset = ARRAY_EQUAL(str.ABSRANGE, [0.0, 0.0])
    abs_str = unset ? 'Values not yet returned' : $
      STRING(str.ABSRANGE*str.MULT, str.UNIT, FORMAT = $
             "(G10.4,' to ',G10.4,' ',A)")
    WIDGET_CONTROL, info.ABSLAB, SET_VALUE = abs_str

    unset = str.SDEV EQ 0.0
    mode_str = unset ? 'Autoscale to set' $
      : STRING(str.MODE*str.MULT, FORMAT="(G11.4)")
    WIDGET_CONTROL, info.MTEXT, SET_VALUE = mode_str
    sdev_str = unset ? 'Autoscale to set' $
      : STRING(str.SDEV*str.MULT, FORMAT="(G10.4)")
    WIDGET_CONTROL, info.SDTEXT, SET_VALUE = sdev_str
    range_str = STRING(str.RANGE*str.MULT, FORMAT="(G11.4,',',G11.4)")
    WIDGET_CONTROL, info.RTEXT, SET_VALUE = range_str
    beta_str = STRING(str.BETA*str.MULT, FORMAT="(G10.4)")
    WIDGET_CONTROL, info.BTEXT, SET_VALUE = beta_str
    zero_str = STRING(str.ZERO*str.MULT, FORMAT="(G11.4)")
    WIDGET_CONTROL, info.ZTEXT, SET_VALUE = zero_str
ENDIF ELSE BEGIN
    WIDGET_CONTROL, info.ABSLAB, SET_VALUE = 'no image set'
    WIDGET_CONTROL, info.MTEXT,  SET_VALUE = ''
    WIDGET_CONTROL, info.SDTEXT, SET_VALUE = ''
    WIDGET_CONTROL, info.RTEXT,  SET_VALUE = ''
    WIDGET_CONTROL, info.BTEXT,  SET_VALUE = ''
    WIDGET_CONTROL, info.ZTEXT,  SET_VALUE = ''
ENDELSE

END

PRO ximview_scale_event, event
; Processes events from the "set scaling" dialog
;
COMPILE_OPT IDL2, HIDDEN
ON_ERROR, 1

WIDGET_CONTROL, event.TOP, GET_UVALUE = info

tag = TAG_NAMES(event, /STRUCTURE_NAME)
CASE STRMID(tag,0,12) OF
    'WIDGET_BUTTO': BEGIN
        ntab = N_ELEMENTS(info.ISCREEN)
        WIDGET_CONTROL, event.ID, GET_VALUE = label
        CASE label OF
            'Select all': BEGIN
                vals = REPLICATE(1B,ntab)
                WIDGET_CONTROL, info.SEQROW, SET_VALUE = vals
                RETURN 
            END
            'Clear all': BEGIN
                vals = BYTARR(ntab)
                WIDGET_CONTROL, info.SEQROW, SET_VALUE = vals
                WIDGET_CONTROL, info.PBASE, SET_BUTTON = 0
                reset_scale_dialog, info, -1
                info.CURRENT = -1
                WIDGET_CONTROL, event.TOP, SET_UVALUE = info
                RETURN
            END
            'Autoscale': BEGIN
                                ; For the current selected tab
                                ; autoscale algorithm to calculate
                                ; zero and range
                IF info.CURRENT NE -1 THEN BEGIN
                    frames = info.CURRENT
                                ; NB we get "current frame" even if
                                ; not currently selected!

                ENDIF ELSE BEGIN
                    WIDGET_CONTROL, info.SEQROW, GET_VALUE = frames
                    frames = WHERE(frames EQ 1)
                    IF frames[0] EQ -1 THEN BEGIN
                        MESSAGE, /INFORMATIONAL, $
                          'Please select a tab for autoscaling'
                        RETURN
                    ENDIF
                    frames = info.ISCREEN[frames]
                ENDELSE
                                ; Unselect other frames to make it
                                ; clear we are only using one
                new = BYTARR(ntab)
                iframe = (WHERE(info.ISCREEN EQ frames[0]))[0]
                new[iframe] = 1B
                WIDGET_CONTROL, info.SEQROW, SET_VALUE = new
                WIDGET_CONTROL, info.PBASE, SET_BUTTON = 0
                but = WIDGET_INFO(info.PBASE,/CHILD)
                FOR i=1,iframe DO but = WIDGET_INFO(but,/SIBLING)
                WIDGET_CONTROL, but, SET_BUTTON = 1

                itab = WHERE((*info.TABARR).SCREEN EQ frames[0])
                str = (*info.TABARR)[itab]
                oldcol   = str.COLTAB
                oldfun   = str.TRFUNC
                oldrange = str.RANGE
                oldbeta  = str.BETA
                oldzero  = str.ZERO
                unset = str.SDEV EQ 0.0
                IF unset THEN BEGIN ; process image to get mode, rms
                    ; This might take some time
                    WIDGET_CONTROL, event.TOP, SENSITIVE = 0
                    WIDGET_CONTROL, /HOURGLASS

                    IF ~str.IM_PTR THEN BEGIN
                        howto = 1
                        imap = frames[0]
                        WIDGET_CONTROL, info.LABEL, GET_UVALUE = data, /NO_COPY
                    ENDIF ELSE BEGIN
                        howto = 3
                        data = str.IM_PTR
                        imap = 0
                    ENDELSE

                    IF info.DOPLOTID NE -1 THEN BEGIN
                        WIDGET_CONTROL, info.DOPLOTID, GET_VALUE = do_plot
                        WIDGET_CONTROL, info.DOPLOTID, /DESTROY
                        info.DOPLOTID = -1
                    ENDIF ELSE do_plot = 0B

                    IF do_plot[0] THEN BEGIN
                        swap_lut, dummy1, dummy, old_graph
                        draw = WIDGET_DRAW(info.RIGHT, XSIZE = 300, $
                                           YSIZE = 300, RETAIN = 2)
                        WIDGET_CONTROL, draw, GET_VALUE = plotwin
                        WSET, plotwin
                        WIDGET_CONTROL, event.TOP, SET_UVALUE = info
                    ENDIF

                    scaling_params, dummy, data, imap, 0, howto, $
                      info.ISCREEN, ar, zero, sdev, PLOT = do_plot

                    absmax = 0.1d0*MAX(ABS(ar))
                    IF str.MULT NE 1.0 THEN BEGIN ; Recover original unit:
                       test = NUMUNIT(str.MULT, str.UNIT, MULTIPLIER = mult, $
                                      OUT_UNIT = unit, /FORCE)
                    ENDIF ELSE unit = str.UNIT

                    test = NUMUNIT(absmax, unit, MULTIPLIER = mult, $
                                   OUT_UNIT = ounit, /FORCE)
                    str.UNIT = ounit
                    str.MULT = mult
                    
                    IF howto EQ 1 THEN $
                      WIDGET_CONTROL, info.LABEL, SET_UVALUE = data, /NO_COPY
                    IF do_plot THEN restore_lut, dummy, old_graph

                    WIDGET_CONTROL, event.TOP, /SENSITIVE
                ENDIF ELSE BEGIN
                    ar   = str.ABSRANGE
                    zero = str.MODE
                    sdev = str.SDEV
                ENDELSE

                str.COLTAB = WIDGET_INFO(info.COLID, /DROPLIST_SELECT)
                str.TRFUNC = WIDGET_INFO(info.FUNCID, /DROPLIST_SELECT)
                set_scale, [ar, zero, sdev], str
                reset_scale_dialog, info, str, /AUTO

; Avoid rounding errors on range:
                WIDGET_CONTROL, info.RTEXT,  GET_VALUE = range
                test = FLOAT(STRSPLIT( range,', ',/EXTRACT))
                info.ROUNDING = str.RANGE*str.MULT - test
                WIDGET_CONTROL, event.TOP, SET_UVALUE = info

                str.RANGE  = oldrange ; Don't update these until accepted.
                str.BETA   = oldbeta  ;
                str.ZERO   = oldzero  ;
                str.TRFUNC = oldfun   ;
                str.COLTAB = oldcol
                (*info.TABARR)[itab] = str


                RETURN
            END
            'Accept': BEGIN ; Extract results from widget
                WIDGET_CONTROL, info.SEQROW, GET_VALUE = frames
                frames = WHERE(frames EQ 1)
                screens = frames[0] EQ -1 ? -1 : info.ISCREEN[frames]
                coltab = info.COLID NE -1 ?  $
                  WIDGET_INFO(info.COLID,  /DROPLIST_SELECT) : 0
                trfunc = WIDGET_INFO(info.FUNCID, /DROPLIST_SELECT)
                WIDGET_CONTROL, info.WRAPID, GET_VALUE = wrap
                wrap = ([0, -1, 1])[wrap]
                WIDGET_CONTROL, info.RTEXT,  GET_VALUE = range
                WIDGET_CONTROL, info.BTEXT,  GET_VALUE = beta
                WIDGET_CONTROL, info.ZTEXT,  GET_VALUE = zero
                range = FLOAT(STRSPLIT( range,', ',/EXTRACT, COUNT = count))
                CASE count OF
                    0: range = 'Not set'
                    1:          ; should be OK
                    2: range = range + info.ROUNDING
                    ELSE: BEGIN
                        ok = DIALOG_MESSAGE( /INFORMATION, $
                                             DIALOG_PARENT = event.TOP, $
                             'Enter at most 2 numbers for intensity range')
                        RETURN
                    END
                ENDCASE
                zero =  VALID_NUM(zero[0]) ? FLOAT(zero[0]) : 0.0
                beta =  VALID_NUM(beta[0]) ? FLOAT(beta[0]) : 0.0
                IF trfunc EQ 1 && (beta EQ 0.0 || ~FINITE(beta)) THEN BEGIN
                        ok = DIALOG_MESSAGE( /INFORMATION, $
                                             DIALOG_PARENT = event.TOP, $
                             'Set Beta for Asinh scaling')
                        RETURN
                ENDIF
                result = {frames: screens, coltab: coltab, wrap: wrap, $
                          trfunc: trfunc, range: range, beta: beta, $
                          zero: zero, current: info.CURRENT}
            END
            'Cancel': result = {frames: -1}
        ENDCASE
        WIDGET_CONTROL, info.RETID, SET_UVALUE = result
        WIDGET_CONTROL, event.TOP, /DESTROY
    END
    'WIDGET_TEXT_':             ; no action needed
    '': IF event.ID EQ info.DEFROW THEN BEGIN
        screen = info.ISCREEN[event.VALUE]
        CASE event.SELECT OF
            0:                  ; Do nothing, exclusive button so every
                                ; deselection is followed by a selection

; Code for non-exclusive case:
;          IF event.VALUE EQ info.CURRENT THEN BEGIN
;                                ; Replace settings with details for
;                                ; next tab or nothing
;                WIDGET_CONTROL, info.DEFROW, GET_VALUE = frames
;                frames = WHERE(frames EQ 1)
;                IF frames[0] EQ -1 THEN BEGIN
;                    reset_scale_dialog, info, -1
;                ENDIF ELSE BEGIN
;                    iscreen = WHERE((*info.TABARR).SCREEN EQ frames[0])
;                    reset_scale_dialog, info, (*info.TABARR)[iscreen]
;                ENDELSE
;                info.CURRENT = frames[0]
;                WIDGET_CONTROL, event.TOP, SET_UVALUE = info
;           ENDIF
            1: BEGIN
                                ; Replace settings with current details for
                                ; this tab
                isc = WHERE((*info.TABARR).SCREEN EQ screen)
                reset_scale_dialog, info, (*info.TABARR)[isc]

                info.CURRENT = screen
                WIDGET_CONTROL, event.TOP, SET_UVALUE = info
            END
        ENDCASE
    ENDIF
    'WIDGET_DROPL': ; no action needed
    ELSE: MESSAGE, /INFORMATIONAL, 'Received unknown event type '+tag
ENDCASE
END

PRO ximview_scale, event
; Launches dialog to re-set the scale of one or more tabs
;
COMPILE_OPT IDL2, HIDDEN
COMMON gr_global
ON_ERROR, 1

WIDGET_CONTROL, event.TOP,  GET_UVALUE = state
WIDGET_CONTROL, state.TABS, GET_UVALUE = mode
tabarr = state.TABARR
tablab = get_tab_uvals(state.TABS)
first = state.FIRST

                                ; Select screens to scale. Exclude RGB
null = PTR_NEW()
good     = WHERE((*tabarr).BYTE_PTR NE null)
IF good[0] EQ -1 THEN RETURN ; No image to scale.

screens  = (*tabarr).SCREEN
current  = screens[0]
itab     = BSORT(screens)
iscreen  = screens[good]
iscreen  = iscreen[BSORT(iscreen)]
ntab     = N_ELEMENTS(iscreen)
icurrent = WHERE(iscreen EQ current, count)
IF count EQ 0 THEN BEGIN
    current = iscreen[0]
    icurrent = 0
ENDIF

WIDGET_CONTROL, state.TABS, SENSITIVE = 0 ; Avoids confusion

; see if this event carries extra info:
tags = TAG_NAMES(event)
vals = WHERE(tags EQ 'VALUE')
are_new = vals NE -1
IF are_new THEN BEGIN
    newtabs = event.VALUE ; list of newly-created tabs
    setup = BYTARR(N_ELEMENTS(newtabs))
    str = (*tabarr)[newtabs[0]]
    icurrent = WHERE(iscreen EQ newtabs[0])
    current = str.SCREEN
ENDIF ELSE str = (*tabarr)[0]

swap_lut, *state.XIM_GRAPH, str, oldgraph

; Set up widget
seqid   = INTARR(ntab)
seqrow  = iscreen
seqstrs = STRING(seqrow, FORMAT = "(I2)")
pbase   = LONARR(ntab)

is_mswin = STRCMP(!version.OS_FAMILY,'Windows', 4,/ FOLD_CASE)
IF is_mswin THEN BEGIN ; Set a fixed-width font:
    old_font = WIDGET_INFO(state.LABEL, /FONTNAME)
    font = 'lucida console*10'
    WIDGET_CONTROL, DEFAULT_FONT = font
ENDIF

init = BYTARR(ntab)
IF are_new THEN init[newtabs[0]] = 1 ELSE init[icurrent] = 1

done = 1                        ; Default if no tabs to be updated
newscale = BYTARR(ntab)

REDO:

query  = WIDGET_BASE(GROUP_LEADER = event.TOP, /MODAL, $
                    /COLUMN, TITLE = 'Set Intensity Scaling Parameters')
block  = WIDGET_BASE(query,/COLUMN, /ALIGN_TOP)
seqrow = CW_BGROUP(block, tablab[iscreen], UVALUE = -1, /ROW, /NONEXCLUSIVE, $
                   IDS = ids, LABEL_LEFT = 'Choose Tabs: ', $
                   SET_VALUE = init)
geom = WIDGET_INFO(query,/GEOMETRY)
DEVICE, GET_SCREEN_SIZE = screen
doscroll = geom.SCR_XSIZE+100 GT screen[0] 
IF doscroll THEN BEGIN
   WIDGET_CONTROL, seqrow, /DESTROY
   scroll_length = 600 < (geom.SCR_XSIZE - 100)
   seqrow = CW_BGROUP(block, tablab[iscreen], UVALUE = -1, /ROW, $
                      /NONEXCLUSIVE, IDS = ids, LABEL_LEFT = 'Choose Tabs: ', $
                      SET_VALUE = init, $
                      /SCROLL, X_SCROLL_SIZE=500, Y_SCROLL_SIZE=30)
ENDIF
pbase  = WIDGET_INFO(ids[0], /PARENT)
bbase  = WIDGET_BASE(block, /ROW, /ALIGN_RIGHT)
void   = WIDGET_BUTTON(bbase, VALUE = 'Select all')
void   = WIDGET_BUTTON(bbase, VALUE = 'Clear all')

IF doscroll THEN BEGIN
   IF icurrent EQ -1 THEN defrow = $
      CW_BGROUP(query, tablab[iscreen], /EXCLUSIVE, /ROW, $
      LABEL_LEFT = 'Make Current:', $
      /SCROLL, X_SCROLL_SIZE=500, Y_SCROLL_SIZE=30) $
   ELSE defrow = $
      CW_BGROUP(query, tablab[iscreen], /EXCLUSIVE, /ROW, $
                LABEL_LEFT = 'Make Current:', $
                SET_VALUE = icurrent, $
                /SCROLL, X_SCROLL_SIZE=500, Y_SCROLL_SIZE=30)
ENDIF ELSE BEGIN
   IF icurrent EQ -1 THEN defrow = $
      CW_BGROUP(query, tablab[iscreen], /EXCLUSIVE, /ROW, $
      LABEL_LEFT = 'Make Current:') $
   ELSE defrow = $
      CW_BGROUP(query, tablab[iscreen], /EXCLUSIVE, /ROW, $
                LABEL_LEFT = 'Make Current:', $
                SET_VALUE = icurrent)
ENDELSE
base1  = WIDGET_BASE(query,/ROW)
left   = WIDGET_BASE(base1,/COLUMN)
right  = WIDGET_BASE(base1,/COLUMN) ; reserved for plot

IF colmap THEN $
  colID  = WIDGET_DROPLIST(left, VALUE = colour_schemes(-1), $
                           TITLE = 'Colour table:     ') ELSE colID = -1

funcID = WIDGET_DROPLIST(left, VALUE = scale_funs(-1), $
                         TITLE =   'Transfer function:')
IF ~are_new THEN BEGIN
    winit = ([1,0,2])[str.WRAP+1]
    wrapID = CW_BGROUP(left, ['No','Bright end only','Full'], /EXCLUSIVE, $
                       SET_VALUE = winit, LABEL_LEFT = 'Wrap colour table?', $
                       /ROW, /RETURN_NAME)
ENDIF ELSE wrapID = CW_BGROUP(left, ['No','Bright end only','Full'], $
                              /EXCLUSIVE, LABEL_LEFT = 'Wrap colour table?', $
                              /ROW, /RETURN_NAME)

unset = ARRAY_EQUAL(str.ABSRANGE, [0.0, 0.0])
abs_str = unset ? 'Values not yet returned' : $
       STRING(str.ABSRANGE*str.MULT, str.UNIT, FORMAT = $
              "(G11.4,' to ',G11.4,' ',A)")
base   = WIDGET_BASE(left,/ROW)
void   = WIDGET_LABEL(base, VALUE = 'Full data range:     ', /ALIGN_LEFT)
abslab = WIDGET_LABEL(base, VALUE = abs_str,  /DYNAMIC_RESIZE, /ALIGN_LEFT)

unset = str.SDEV EQ 0.0
base   = WIDGET_BASE(left, /ROW)
void   = WIDGET_LABEL(base,  VALUE = 'Mode: ')
mode_str = unset ? 'Autoscale to set' $
                 : STRING(str.MODE*str.MULT, FORMAT="(G11.4)")
mtext = WIDGET_LABEL(base, /DYNAMIC_RESIZE, /ALIGN_LEFT, VALUE = mode_str)
void   = WIDGET_LABEL(base,  VALUE = 'Estimated rms: ')
sdev_str = unset ? 'Autoscale to set' $
                 : STRING(str.SDEV*str.MULT, FORMAT="(G10.4)")
sdtext = WIDGET_LABEL(base, /DYNAMIC_RESIZE, /ALIGN_LEFT, VALUE = sdev_str)

doplotID = unset ? CW_BGROUP(left, ['Plot histogram on autoscale?'], $
                            /NONEXCLUSIVE) : -1

unset = ARRAY_EQUAL(str.RANGE, [0.0, 0.0])
base   = WIDGET_BASE(left, /ROW)
void   = WIDGET_LABEL(base,  VALUE = 'Data range for scaling:')
range_str = unset ? '' : STRING(str.RANGE*str.MULT, FORMAT="(G11.4,',',G11.4)")
rtext  = WIDGET_TEXT(base, /EDITABLE, XSIZE = 23, VALUE = range_str)

; Take action to avoid rounding errors on range:
IF ~unset THEN BEGIN
    test = FLOAT(STRSPLIT( range_str,', ',/EXTRACT, COUNT = count))
    rounding = str.RANGE - test/str.MULT
ENDIF ELSE rounding = [0., 0.]

unset = str.BETA EQ 1.0
base   = WIDGET_BASE(left, /ROW)
void   = WIDGET_LABEL(base,  VALUE = 'Beta for Asinh scaling:')
beta_str = unset ? '' : STRING(str.BETA*str.MULT, FORMAT="(G10.4)")
btext  = WIDGET_TEXT(base, /EDITABLE, XSIZE = 12, VALUE = beta_str)

base   = WIDGET_BASE(left, /ROW)
void   = WIDGET_LABEL(base,  VALUE = 'Effective zero level:  ')
zero_str = STRING(str.ZERO*str.MULT, FORMAT="(G11.4)")
ztext  = WIDGET_TEXT(base, /EDITABLE, XSIZE = 12, VALUE = zero_str)

bbase  = WIDGET_BASE(query, /ROW, /ALIGN_CENTER)
;void   = WIDGET_BUTTON(bbase, VALUE = 'Clear tab selection')
autoID = WIDGET_BUTTON(bbase, VALUE = 'Autoscale')
void   = WIDGET_BUTTON(bbase, VALUE = 'Accept')
void   = WIDGET_BUTTON(bbase, VALUE = 'Cancel')

                                ; Structure with return ID, internal
                                ; IDs, and defaults
info = {top: event.TOP, retid: event.ID, seqrow: seqrow, pbase: pbase, $
        defrow: defrow, colID: colID, funcID: funcID, wrapID: wrapID, $
        abslab: abslab, mtext: mtext, sdtext: sdtext, doplotID: doplotID, $
        rtext: rtext, btext: btext, ztext: ztext, autoID: autoID, $
        tabarr: tabarr, label: state.LABEL, iscreen: iscreen, $
        current: current, rounding: [0., 0.], plot: -1L, right: right}

WIDGET_CONTROL, query,    SET_UVALUE = info, /NO_COPY
IF colmap THEN WIDGET_CONTROL, colID,    SET_DROPLIST_SELECT = str.COLTAB
WIDGET_CONTROL, funcID,   SET_DROPLIST_SELECT = str.TRFUNC
WIDGET_CONTROL, event.ID, GET_UVALUE = save, /NO_COPY

; Fire off widget and collect result
WIDGET_CONTROL, query, /REALIZE
XMANAGER, 'ximview_scale', query

WIDGET_CONTROL, event.ID, GET_UVALUE = result, /NO_COPY
WIDGET_CONTROL, event.ID, SET_UVALUE = save, /NO_COPY

                                ; Interpret result here
IF SIZE(result,/TYPE) NE 8 THEN frame = -1 ELSE frame = result.FRAMES

nframe = N_ELEMENTS(frame)

schemes = colour_schemes(-1)
scales  = scale_funs(-1)

IF nframe GT 0 THEN IF frame[0] EQ -1 THEN nframe = 0
IF nframe GT 0 THEN BEGIN
    WIDGET_CONTROL, /HOURGLASS
    current_str = (*tabarr)[itab[result.CURRENT]]
    mult0 = current_str.MULT  
ENDIF

IF first THEN BEGIN
    old_coltab = -1
    old_izero  = -1
    lut = 0
ENDIF ELSE BEGIN
    oldstr = (*tabarr)[itab[0]]
    old_coltab = oldstr.COLTAB
    old_izero = oldstr.IZERO
    lut = *oldstr.LUT
ENDELSE

FOR i=0,nframe-1 DO BEGIN
                                ; Indices in this loop:
                                ; i indexes tabs returned by widget
                                ; idx indexes tabs sent to widget
                                ; frame[i] indexes the frames=tabs
                                ; itab indexes tabarr, i.e.
                                ; tabarr[itab[i]] has .SCREEN = i
    resca = 0B  &  relut = 0B
    str = (*tabarr)[itab[frame[i]]]
    idx = WHERE(frame[i] EQ iscreen)
    
    IF str.BETA NE 1.0 THEN BEGIN
        range_str = STRING(str.RANGE*str.MULT, FORMAT="(G11.4,',',G11.4)")
        test = FLOAT(STRSPLIT( range_str,', ',/EXTRACT, COUNT = count))
        rounding = str.RANGE - test/str.MULT
    ENDIF ELSE rounding = [0., 0.]   
    
    zero = result.ZERO / mult0
    IF str.ZERO NE zero THEN BEGIN
        str.ZERO = zero
        resca = result.TRFUNC EQ WHERE(scales EQ 'Asinh')
    ENDIF
    beta = result.BETA / mult0
    IF str.BETA NE beta THEN BEGIN
        str.BETA = beta
        IF result.TRFUNC EQ WHERE(scales EQ 'Asinh') THEN resca = 1B
    ENDIF
    IF (str.COLTAB NE result.COLTAB)  || first THEN BEGIN
        relut = 1B
        str.COLTAB = result.COLTAB
    ENDIF
    IF str.WRAP NE result.WRAP THEN BEGIN
        resca = 1B
        str.WRAP = result.WRAP
    ENDIF
    IF str.TRFUNC NE result.TRFUNC THEN BEGIN
        resca = 1B
        str.TRFUNC = result.TRFUNC
    ENDIF
    range = result.RANGE

; Take action to avoid rounding errors on range:
    redo = SIZE(range, /TYPE) EQ 7 ; 'Not set' returned if range is not set
    IF redo THEN resca = 1B ELSE BEGIN
        range /= mult0
        IF ARRAY_EQUAL(str.RANGE, range+rounding) THEN $
          range = str.RANGE ELSE resca = 1B
    ENDELSE

    IF are_new THEN BEGIN
        thistab = WHERE(frame[i] EQ newtabs)
        IF thistab NE -1 THEN BEGIN
            resca = 1B
            setup[thistab] = 1B
        ENDIF
    ENDIF

    IF resca THEN BEGIN
        IF ~str.IM_PTR THEN BEGIN
            WIDGET_CONTROL, state.LABEL, GET_UVALUE = data, /NO_COPY
            *str.BYTE_PTR = scale_image(data, range, str.WRAP, str.TRFUNC, $
                                        zero, beta, ABS_RANGE = ar)
            WIDGET_CONTROL, state.LABEL, SET_UVALUE = data, /NO_COPY
        ENDIF ELSE *str.BYTE_PTR = scale_image(*str.IM_PTR, range, str.WRAP, $
                                       str.TRFUNC, zero, beta, ABS_RANGE = ar)

        str.ABSRANGE = ar
        str.RANGE    = range
        absmax = 0.1d0*MAX(ABS(ar))
        IF str.MULT NE 1.0 THEN BEGIN ; Recover original unit:
            test = NUMUNIT(str.MULT, str.UNIT, MULTIPLIER = mult, $
                           OUT_UNIT = unit, /FORCE)
        ENDIF ELSE unit = str.UNIT

        test = NUMUNIT(absmax, unit, MULTIPLIER = mult, $
                       OUT_UNIT = ounit, /FORCE)
        str.UNIT = ounit
        str.MULT = mult
        set_print_fmt, str
        
        newscale[idx] = 1B
    ENDIF

    izero = scale_image(zero, range, str.WRAP, str.TRFUNC, zero, beta)
                            
    IF izero NE str.IZERO && str.COLTAB EQ 4 THEN relut = 1B
    str.IZERO = izero
 
    IF N_ELEMENTS(*str.LUT) EQ 0 THEN relut = 1B
    
    IF relut THEN BEGIN
        IF str.COLTAB EQ old_coltab && $
          (str.COLTAB NE 4 || str.IZERO EQ old_izero) $
          THEN *str.LUT = lut ELSE BEGIN
            ximview_lut, str.COLTAB, str.IZERO, decomp
            str.DECOMPOSED = decomp
            TVLCT, r, g, b, /GET
            lut = {r: r, g: g, b: b, line: !P.color, absent: !P.background}
            *str.LUT = lut
            gscroll_setpar, /BLANK
        ENDELSE
         
        IF redraw_req THEN newscale[idx] = 1B
        IF state.GLOBAL_COLOUR THEN BEGIN
            (*tabarr).COLTAB = result.COLTAB
            FOR jtab=0,N_ELEMENTS(*tabarr) - 1 DO *(*tabarr)[jtab].LUT = lut
            IF redraw_req THEN newscale[*] = 1B
        ENDIF
    ENDIF
 
    IF resca THEN BEGIN
        nside = state.NSIDE
        IF nside NE 0 && state.IS_ASTROM THEN $
          fill_gores, nside, state.IMSIZE, state.ASTROM, str.BYTE_PTR
    ENDIF
    
    (*tabarr)[itab[frame[i]]] = str
    old_coltab = str.COLTAB
    old_izero  = str.IZERO
ENDFOR

IF are_new THEN BEGIN ; Make sure all new tabs have byte images:
    unset = WHERE(~setup)
    IF unset[0] NE -1 THEN BEGIN
       init[*] = 0B
       init[newtabs[unset[0]]] = 1B
       icurrent = itab[newtabs[unset[0]]]
       str = (*tabarr)[icurrent]
       current = str.SCREEN
       GOTO, REDO
    ENDIF
ENDIF

; Restore system default (?)
IF is_mswin THEN WIDGET_CONTROL, DEFAULT_FONT = old_font

FOR i = 0, ntab-1 DO IF newscale[i] THEN BEGIN
    j = itab[iscreen[i]]
    IF redraw_req THEN BEGIN
        lutptr = (*tabarr)[j].LUT
        TVLCT, (*lutptr).R, (*lutptr).G, (*lutptr).B
        !P.background = (*lutptr).ABSENT
        !P.color      = (*lutptr).LINE
    ENDIF
    update_screen, tabarr, j, mode, done
ENDIF

IF redraw_req THEN DEVICE, DECOMPOSED = 0B $
              ELSE DEVICE, DECOMPOSED = (*tabarr)[0].DECOMPOSED

mode.DONE = done
WIDGET_CONTROL, state.TABS,  SET_UVALUE = mode

IF are_new THEN BEGIN ; Make one of the new tabs current
    newevent = {ID: 0L, TOP: 0L, HANDLER: 0L, TAB: newtabs[0]}
    WIDGET_CONTROL, state.TABS, SET_TAB_CURRENT = newtabs[0]
    WIDGET_CONTROL, state.TABS, SEND_EVENT = newevent
ENDIF

                                ; Request another go if loading not finished
IF done EQ 0 THEN WIDGET_CONTROL, (*tabarr)[0].DRAW, TIMER = 0.5

restore_lut, *state.XIM_GRAPH, oldgraph

WIDGET_CONTROL, state.TABS, /SENSITIVE

END

