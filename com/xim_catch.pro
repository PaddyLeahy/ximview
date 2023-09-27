; -----------------------------------------------------------------------------
;
;  Copyright (C) 2016   J. P. Leahy
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
PRO xim_catch, top, state, mode
;+
; NAME: 
;       XIM_catch
;
; PURPOSE: 
;       routine to allow programs to interact with Ximview
;       widget, by returning the state and mode structures
;
; CATEGORY:
;       Widget helpers
;
; CALLING SEQUENCE: 
;
;        XIM_CATCH, top, state, mode
;
; INPUTS:
;       none
;
; OUTPUTS:
;       top         Widget ID of Ximview top-level base.
;       state       Top-level control structure for ximview
;       mode        structure specifying current state of widget
;       
;
; EXAMPLE:
;
;            XIMVIEW, 'wmap_band_iqumap_r9_3yr_K_v2', '*', COL=[1,2,3]
;            XIM_CATCH, top, state, mode
;                              ; If you are going to interact with the
;                              ; widget, freeze it right now:
;            WIDGET_CONTROL, state.TABS, SENSITIVE = 0
;
;  now you can really mess it up! Don't forget to set
;  sensitive=1 before you quit your code.
;
;            
; MODIFICATION HISTORY:
;       Written by:      J. P. Leahy, September 2020
;                        (extracted from ximgo)
;-
COMPILE_OPT IDL2
ON_ERROR, 2
start = SYSTIME(1)

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
WIDGET_CONTROL, state.TABS, GET_UVALUE = mode

END

