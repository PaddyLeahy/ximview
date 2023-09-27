; -----------------------------------------------------------------------------
;
;  Copyright (C) 2007-2008   J. P. Leahy
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
FUNCTION rotcomb, a, b
; Looks up the combined effects of two calls to ROTATE
COMPILE_OPT IDL2, HIDDEN

lut = [[0, 1, 2, 3, 4, 5, 6, 7], $
       [1, 2, 3, 0, 7, 4, 5, 6], $
       [2, 3, 0, 1, 6, 7, 4, 5], $
       [3, 0, 1, 2, 5, 6, 7, 4], $
       [4, 5, 6, 7, 0, 1, 2, 3], $
       [5, 6, 7, 4, 3, 0, 1, 2], $
       [6, 7, 4, 5, 2, 3, 0, 1], $
       [7, 4, 5, 6, 1, 2, 3, 0]]
RETURN, lut[a,b]
END

FUNCTION gscroll_find_panel, iip, pixmap_panel, rot
;
; Sees if image panel ip is loaded with correct rotation
; Inputs:
;      iip: 1-D index to an image panel
;      pixmap_panel: array for the screen of interest
;      rot: required rotation.
; Returns:
;      1-D index to pixmap_panel array (ordered by distance from centre)
;
COMPILE_OPT IDL2, HIDDEN

COMMON GRID_GSCROLL

i = WHERE(pixmap_panel.idx EQ iip AND pixmap_panel.rot EQ rot)
IF i[0] EQ -1 && do_wrap && image_panel[iip].wrap GE 0 THEN $
; Check to see if the wrapped version of the panel is loaded:
    i = WHERE(pixmap_panel.idx EQ image_panel[iip].wrap AND $
              rotcomb(pixmap_panel.rot,image_panel[iip].wrap_rot) EQ rot)
RETURN, i[0]

END

