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
PRO restore_lut, new_graph, old_graph
;
; Restore old graphics state
;
; Inputs:
;    new_graph: structure containing current graphics state
;    old_graph: structure containing old details created by swap_lut
;    
;
COMPILE_OPT IDL2, HIDDEN
COMMON gr_global

IF N_ELEMENTS(new_graph) NE 0 THEN BEGIN
    new_graph.X = !X
    new_graph.Y = !Y
    new_graph.P = !P
ENDIF

IF old_graph.WINDOW EQ -1 THEN WSET, -1 ELSE BEGIN
    DEVICE, WINDOW_STATE = ws
    IF ws[old_graph.WINDOW] EQ 1 THEN WSET, old_graph.WINDOW ELSE BEGIN
        WSET, -1
        IF old_graph.WINDOW NE 0 THEN MESSAGE, /INFORMATIONAL, $
          'Cannot find old window #'+STRTRIM(STRING(old_graph.WINDOW),2)
    ENDELSE
ENDELSE

DEVICE, DECOMPOSED = old_graph.DECOMPOSED
IF windev NE old_graph.DEVICE THEN SET_PLOT, old_graph.DEVICE

!X = old_graph.X  &  !Y = old_graph.Y  &  !Z = old_graph.Z
!ORDER = old_graph.ORDER
!P = old_graph.P

TVLCT, old_graph.RED, old_graph.GREEN, old_graph.BLUE

END
