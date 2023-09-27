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
FUNCTION cellindex, x0, x1, y0, y1, nx, xoff, yoff, SORT=sort
; Returns as a 1D array the index values for the array elements in the
; array section (x0:x1)+xoff, (y0:y1)+yoff, optionally sorted with
; distance from (xoff,yoff)
;
COMPILE_OPT IDL2, HIDDEN

IF ~N_ELEMENTS(xoff) THEN xoff = 0
IF ~N_ELEMENTS(yoff) THEN yoff = 0

dx = x1 - x0 + 1 & dy = y1 - y0 + 1
IF dx GT nx THEN MESSAGE, 'Window too big'

xrange = LINDGEN(dx) + x0 + xoff
yrange = LINDGEN(dy) + y0 + yoff
xrange = REBIN(xrange, dx, dy, /SAMPLE)
yrange = REBIN(REFORM(yrange, 1, dy), dx, dy, /SAMPLE)

IF KEYWORD_SET(sort) THEN BEGIN
; Sort via distance from offset point:
    dist = SQRT(FLOAT((xrange-xoff)^2 + (yrange-yoff)^2))
    sort = BSORT(dist)
    xrange += yrange*nx
    RETURN, xrange[sort]
ENDIF ELSE RETURN, xrange + yrange*nx

END
