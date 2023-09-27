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
; Module XIMVIEW_BLINK
;
; J. P. Leahy February 2008
;             Updated August 2013: now updates unit as well as label.
;
; Contains the event handlers for blinking in Ximview.
;
PRO blink_setup_event, event
; Handles modal widget setting blink properties
;
COMPILE_OPT IDL2, HIDDEN
ON_ERROR, 1

WIDGET_CONTROL, event.TOP, GET_UVALUE = info

tag = TAG_NAMES(event, /STRUCTURE_NAME)
CASE STRMID(tag,0,12) OF
    'WIDGET_BUTTO': BEGIN
        WIDGET_CONTROL, event.ID, GET_VALUE = label
        ntab = N_ELEMENTS(info.seqrow)
        CASE label OF
            'Clear sequence': BEGIN
                FOR i=0,ntab-1 DO BEGIN
                    WIDGET_CONTROL, info.seqrow[i], SET_UVALUE = -1
                    WIDGET_CONTROL, info.pbase[i], SET_BUTTON = 0
                ENDFOR

                RETURN
            END
            'Done': BEGIN ; Extract results from widget
                WIDGET_CONTROL, info.TEXT, GET_VALUE = period
                WIDGET_CONTROL, info.scalID, GET_VALUE = do_scale
                WIDGET_CONTROL, info.seqrow[0], GET_UVALUE = seq
                FOR i=1,ntab-1 DO BEGIN
                    WIDGET_CONTROL, info.seqrow[i], GET_UVALUE=tab
                    seq = [seq, tab]
                ENDFOR
                result = {cancel: 0B, seq: seq, period: period, $
                          do_scale: do_scale}
            END
            'Cancel': result = {cancel: 1B, seq: info.SEQ, period: info.PERIOD}
        ENDCASE
        WIDGET_CONTROL, info.RETID, SET_UVALUE = result
        WIDGET_CONTROL, event.TOP, /DESTROY
    END
    'WIDGET_TEXT_': ; no action needed
    '': BEGIN
        WIDGET_CONTROL, event.ID, SET_UVALUE = event.VALUE
    END
ENDCASE

END

PRO blink_setup, event, ntab, tablab, seq, period, do_scale, cancel
;
; Launches modal widget to get blink parameters and returns results.
;
; Inputs:
;   event:  Event structure containing Widget IDs used for passing
;           results.
;   ntab:   Number of tabs (excluding blink tab if any)
;   tablab: Labels for each tab.
;
; Input/output:
;   seq:    Blink sequence array
;   period: time per tab while blinking (seconds)
;   do_scale: blink scale bar as well as image? 
;   cancel: true if operation aborted.
;
COMPILE_OPT IDL2, HIDDEN

WIDGET_CONTROL, event.ID, GET_UVALUE = save, /NO_COPY

nseq = N_ELEMENTS(seq)

; Set up widget
seqid = INTARR(ntab)
seqrow = INDGEN(ntab)
seqstrs = STRING(seqrow, FORMAT = "(I2)")
perstr  = STRING(period, FORMAT = "(F5.2)")
vals = nseq EQ ntab ? seq : [seq, REPLICATE(-1, ntab-nseq)]
pbase = LONARR(ntab)
query = WIDGET_BASE(GROUP_LEADER = event.TOP, /MODAL, $
                    /COLUMN, TITLE = 'Set Blink Parameters')
void = WIDGET_LABEL(query, VALUE = 'Tab',/ALIGN_CENTER)

FOR i=0,ntab-1 DO BEGIN
    seqrow[i] =  vals[i] EQ -1 ?           $
      CW_BGROUP(query, tablab[0:ntab-1], UVALUE = -1, /ROW, /EXCLUSIVE, $
                IDS = ids, LABEL_LEFT = 'Seq. '+seqstrs[i]) :           $
      CW_BGROUP(query, tablab[0:ntab-1], UVALUE = vals[i], /ROW, /EXCLUSIVE, $
                IDS = ids, SET_VALUE = vals[i], $
                LABEL_LEFT = 'Seq. ' + seqstrs[i])
    pbase[i] = WIDGET_INFO(ids[0], /PARENT)
ENDFOR

void = WIDGET_LABEL(query,VALUE = 'Blink period (s):')
text = WIDGET_TEXT(query, /EDITABLE, XSIZE = 6, VALUE = perstr)
scalID = CW_BGROUP(query, ['Yes', 'No'], LABEL_LEFT = 'Omit scale bar?', $
                   SET_VALUE = do_scale, UVALUE = do_scale, /EXCLUSIVE, /ROW)
bbase = WIDGET_BASE(query, /ROW)
void = WIDGET_BUTTON(bbase, VALUE = 'Clear sequence')
void = WIDGET_BUTTON(bbase, VALUE = 'Done')
void = WIDGET_BUTTON(bbase, VALUE = 'Cancel')

                                ; Structure with return ID, internal
                                ; IDs, and defaults
info = {retid: event.ID, seqrow: seqrow, pbase: pbase, $
        text: text, seq: seq, period: period, scalID: scalID}
WIDGET_CONTROL, query, SET_UVALUE = info, /NO_COPY

WIDGET_CONTROL, query, /REALIZE
XMANAGER, 'blink_setup', query

WIDGET_CONTROL, event.ID, GET_UVALUE = result, /NO_COPY
WIDGET_CONTROL, event.ID, SET_UVALUE = save, /NO_COPY

IF SIZE(result,/TYPE) NE 8 THEN cancel = 1B ELSE cancel = result.CANCEL
IF cancel THEN RETURN

seq    = result.SEQ
good   = WHERE(seq NE -1)
IF good[0] NE -1 THEN seq = seq[good] ELSE seq = -1

period = result.PERIOD
do_scale = result.DO_SCALE
                                ; Finally, clear any events scheduled
                                ; (i.e. pending timer events), because
                                ; we are starting over.
WIDGET_CONTROL, event.TOP, /CLEAR_EVENTS

END

PRO ximview_blink, event
; Sets blink options and switches blinking on or off, including
; creating or destroying the blink tab.
;
COMPILE_OPT IDL2, HIDDEN
COMMON gr_global
ON_ERROR, 1

start = SYSTIME(1)

WIDGET_CONTROL, event.TOP,  GET_UVALUE = state
WIDGET_CONTROL, state.TABS, GET_UVALUE = mode
tabarr = state.TABARR

ntab = WIDGET_INFO(state.TABS, /TAB_NUMBER)
IF mode.BLINK THEN ntab = ntab - 1

seq = *state.BLINK_SEQ
; First winnow sequence of any tabs that have been deleted:
mseq = MAX(seq)
IF mseq GE ntab THEN BEGIN
    iseq = WHERE(seq LT ntab)
    seq = seq[iseq]
    *state.BLINK_SEQ = seq
ENDIF

nseq     = N_ELEMENTS(seq)
IF nseq GT 0 THEN IF seq[0] EQ -1 THEN nseq = 0
oldblink = mode.BLINK
period   = mode.PERIOD
do_scale = mode.DO_SCALE

tablab = get_tab_uvals(state.TABS)

                                ; Interpret requested action
button = 0B
tag = TAG_NAMES(event, /STRUCTURE_NAME)
CASE tag OF
    'WIDGET_BUTTON': BEGIN
        WIDGET_CONTROL, event.ID, GET_VALUE = label
        CASE label OF
            'Blink on/off': BEGIN
                mode.BLINK = ~mode.BLINK
                button = 1B
            END
            'Blink setup': BEGIN ; Launch modal widget to get info:
                blink_setup, event, ntab, tablab, seq, period, do_scale, cancel
                IF cancel THEN RETURN
                *state.BLINK_SEQ = seq
                mode.PERIOD      = period
                mode.DO_SCALE    = do_scale
                nseq = N_ELEMENTS(seq)
                IF nseq GT 0 THEN IF seq[0] EQ -1 THEN nseq = 0
                mode.BLINK = nseq GT 1
            END
            ELSE: MESSAGE, /INFORMATIONAL, 'Unknown button: ' + label
        ENDCASE
    END
    'WIDGET_TIMER':  ; Blinking in progress, no action needed
    '': mode.BLINK = 0B ; anonymous event is a request to stop
    ELSE: MESSAGE, /INFORMATIONAL, 'Unknown event type received: ' + name
ENDCASE

IF mode.BLINK THEN BEGIN ; check that we can...
    IF nseq LE 1 THEN BEGIN
        IF ntab GT 1 AND button THEN BEGIN
                                ; launch setup widget to get info:
            blink_setup, event, ntab, tablab, seq, period, do_scale, cancel
            IF cancel THEN RETURN
            *state.BLINK_SEQ = seq
            mode.PERIOD      = period
            mode.DO_SCALE    = do_scale
            nseq = N_ELEMENTS(seq)
            IF nseq GT 0 THEN IF seq[0] EQ -1 THEN nseq = 0
            mode.BLINK = nseq GT 1
        ENDIF ELSE BEGIN
            mode.BLINK = 0B
            IF button THEN MESSAGE, /INFORMATIONAL, $
              'Blinking with one screen requested (should not be possible).'
        ENDELSE
    ENDIF
ENDIF

; Final check:
mseq = MAX(seq)
IF mseq GE ntab THEN BEGIN
    iseq = WHERE(seq LT ntab)
    seq = seq[iseq]
    *state.BLINK_SEQ = seq
    nseq = N_ELEMENTS(seq)
    MESSAGE, /INFORMATIONAL, STRING(nseq, FORMAT="(I1)") + $
      ' requested tabs available'
    IF mode.BLINK AND nseq LE 1 THEN BEGIN
        MESSAGE, /INFORMATIONAL, 'Blinking disabled'
        mode.BLINK = 0B
    ENDIF
ENDIF

;
swap_lut, *state.XIM_GRAPH, (*tabarr)[0], old_graph

; If blink state has changed, create or destroy blink tab
IF mode.BLINK NE oldblink THEN BEGIN
    IF mode.BLINK THEN BEGIN
                                ; Make sure all screens are up to date
        WSET, (*tabarr)[0].WINDOW
        gscroll_setpar, HIDDEN = 0
        gscroll_newscreen, 0, *tabarr, mode.ZOOM_FACTOR, mode.X_CENTRE, $
          mode.Y_CENTRE, mode.XHALF, mode.YHALF, done, mode.OVERVIEW
        mode.DONE = done
        IF done EQ 0B THEN WIDGET_CONTROL, (*tabarr)[0].DRAW, TIMER = 0.5

                                ; Now initialize blink tab:
        base  = WIDGET_BASE(state.TABS, TITLE='Blink', UVALUE = 'Blink', $
                            /COLUMN, XPAD = 0, YPAD = 0, SPACE = 0)
        mode.BBASE = base
        WSET, (*tabarr)[0].WINDOW
        draw  = WIDGET_DRAW(base, XSIZE = !D.x_vsize, YSIZE = !D.y_vsize, $
                            RETAIN = 0, EVENT_FUNC = 'ximview_scroll', $
                            /BUTTON_EVENTS, /MOTION_EVENTS)
        scale = WIDGET_DRAW(base, XSIZE = !D.x_vsize, YSIZE = 45, RETAIN = 0)
        WIDGET_CONTROL, state.TABS, SET_TAB_CURRENT = ntab
        WIDGET_CONTROL, draw, GET_VALUE = index
        mode.BWIN = index
        WIDGET_CONTROL, scale, GET_VALUE = index
        mode.BSWIN = index
    ENDIF ELSE BEGIN
        WIDGET_CONTROL, mode.BBASE, /DESTROY
        mode.BWIN = -1  &  mode.BSWIN = -1   & mode.BBASE = -1
                                ; Synchronize
        gscroll_setpar, /HIDDEN
        current = nseq GT 0 ? seq[0] : 0
        WIDGET_CONTROL, state.TABS, SET_TAB_CURRENT = current
        newevent = {ID: 0L, TOP: 0L, HANDLER: 0L, TAB: current}
        WIDGET_CONTROL, state.TABS, SEND_EVENT = newevent
    ENDELSE
ENDIF

IF mode.BLINK THEN BEGIN
                                ; Is the blink tab on top or not?
    current = WIDGET_INFO(state.TABS, /TAB_CURRENT)
    ontop = current EQ ntab

                                ; Do the actual blinking:
    fmt = STRING(MAX(STRLEN(tablab[seq])),FORMAT = "('(A',I2,')')")
    temp_title = STRING(tablab[seq[1]], FORMAT = fmt)
    s1 = WHERE((*tabarr).SCREEN EQ seq[1])
    from = (*tabarr)[s1].WINDOW
    WSET,  mode.BWIN
                                ; Update BLINK_SEQ just before loading
                                ; to give best chance that other
                                ; events (e.g. pix_print) will find
                                ; the right screen
    *state.BLINK_SEQ = SHIFT(seq, -1)
    IF ~redraw_req && ontop THEN BEGIN
        DEVICE, DECOMPOSED = (*tabarr)[s1].DECOMPOSED
        lutptr = (*tabarr)[s1].LUT
        TVLCT, (*lutptr).R, (*lutptr).G, (*lutptr).B
    ENDIF

    DEVICE, COPY = [0, 0, !D.x_vsize, !D.y_vsize, 0, 0, from]
    IF do_scale THEN BEGIN
        from = (*tabarr)[s1].SCALE_INDEX
        WSET, mode.BSWIN
        DEVICE, COPY = [0, 0, !D.x_vsize, !D.y_vsize, 0, 0, from]
    ENDIF
    WIDGET_CONTROL, mode.BBASE, BASE_SET_TITLE = temp_title

                                ; Print coord if we are on blink tab:
    IF ontop THEN BEGIN
                                    ;   Update units on readout label
        mid = form_unit( (*tabarr)[s1].UNIT)

        title_string = state.TITLE.HEAD + mid + state.TITLE.TAIL
        WIDGET_CONTROL, state.READLAB, SET_VALUE=title_string
    
        pix_print, state, 0, start
    ENDIF
    WIDGET_CONTROL, event.ID, TIMER = period
ENDIF

WIDGET_CONTROL, state.TABS, SET_UVALUE = mode

restore_lut, dummy, old_graph

END
