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
PRO swap_lut, new_graph, tabstr, old_graph
;
; Saves current graphics state and swaps it for the one set up in the
; tab control structure tabstr.
;
; Inputs:
;   new_graph:   Structure containing graphics state for ximview
;   tabstr:      Structure describing current tab

; Inputs in common gr_global
;   windev:      Device used for widget graphics ('X' or 'WIN')
;   colmap:      TRUE if device supports color maps (it should!)
;   redraw_req:  TRUE if a redraw is required to change colour scheme.
;                (will be false for visual classes "PseudoColor" or
;                "DirectColor"
; Output:
;   old_graph:   Structure containing details of former graphics
;   state.
;
COMPILE_OPT IDL2, HIDDEN
COMMON gr_global, windev, redraw_req, colmap, badcol, syscol, ibot, itop

ON_ERROR, 0
; Set common if needed
IF N_ELEMENTS(windev) EQ 0 THEN BEGIN
    system = STRUPCASE(!VERSION.os_family)
    windev = STRCMP(system, 'WINDOW', 6) ? 'WIN' : 'X'
ENDIF

; Remember previous colour setting
TVLCT, old_red, old_green, old_blue, /GET

old_device = !D.name

old_x = !X & old_y = !Y & old_Z = !Z
IF old_device NE windev THEN SET_PLOT, windev

; These parameters only apply for windows so get them *after* SET_PLOT
old_window = !D.window
DEVICE, GET_DECOMPOSED = old_decomposed

old_graph = {red: old_red, green: old_green, blue: old_blue, $
             x: old_x, y: old_y, z: old_z, P: !P, order: !ORDER, $
             window: old_window, decomposed: old_decomposed, device: old_device}

IF N_ELEMENTS(new_graph) EQ 0 THEN BEGIN
                     ; Set up graphics structure for ximview:
                     ; default system variables
    !P.CHARSIZE = 1.0
    !P.CHARTHICK = 1
    !P.clip = [0, 0, 639, 511, 0, 0]
    !P.FONT = -1  
    !P.LINESTYLE = 0
    !P.MULTI = 0
    !P.NOCLIP = 0
    !P.NOERASE = 0
    !P.NSUM = 0
    !P.POSITION = 0.0
    !P.PSYM = 0
    !P.REGION = 0.0
    !P.SUBTITLE = ''
    !P.SYMSIZE = 0.0
    !P.T = IDENTITY(4,/DOUBLE)
    !P.T3D = 0
    !P.THICK = 1.0
    !P.TITLE = ''
    !P.TICKLEN = 0.02
    !P.CHANNEL = 0
    
    !ORDER = 0
    
    !X.title = ''
    !X.type = 0
    !X.style = 0
    !X.ticks = 0
    !X.ticklen = 0.
    !X.thick = 0.
    !X.range = 0.
    !X.crange = 0.
    !X.s = [0d0, 1d0]
    !X.margin = [10., 3.]
    !X.omargin = 0.
    !X.window = 0.
    !X.region = 0.
    !X.charsize = 0.
    !X.minor = 0
    !X.tickv = 0d0
    !X.tickname = ''
    !X.gridstyle = 0
    !X.tickformat = ''
    !X.tickinterval = 0
    !X.tickunits = ''

    !Y = !X
    !Y.margin = [4., 2.]

    new_graph = {x: !X, y: !Y, P: !P}
ENDIF ELSE BEGIN
    !X = new_graph.X
    !Y = new_graph.Y
    !P = new_graph.P
    !ORDER = 0
ENDELSE

; Install requested colour table if there is one:
IF SIZE(tabstr, /TYPE) EQ 8 THEN BEGIN
    IF ~colmap THEN BEGIN
        MESSAGE, /INFORMATIONAL, $
          'Apparently this device does not support colour tables'
        RETURN
    ENDIF
    IF ~PTR_VALID(tabstr.LUT) THEN MESSAGE, $
      'Internal error: missing LUT structure on screen ' + $
      STRTRIM( STRING(tabstr.SCREEN), 2)
    IF N_ELEMENTS(*tabstr.LUT) EQ 0 THEN RETURN ; lut not set yet
    lut = *tabstr.LUT
    IF redraw_req THEN DEVICE, DECOMPOSED = 0B $
                  ELSE DEVICE, DECOMPOSED = tabstr.DECOMPOSED

    TVLCT, lut.R, lut.G, lut.B

    !P.background = lut.ABSENT
    !P.color      = lut.LINE
    new_graph.P = !P 
    WSET, tabstr.WINDOW
ENDIF

END
