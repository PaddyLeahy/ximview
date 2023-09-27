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
PRO gscroll_refresh, xhalf, yhalf
;
; Re-draws current screen
;
COMPILE_OPT IDL2, HIDDEN

COMMON GRID_GSCROLL

; Decode zoom factor to integers:
zoom = old_zoomf GE 1.0
zfac = zoom ? FIX(old_zoomf) : FIX(1./old_zoomf)

WSET, screens[0].viewwin
gscroll_pix2view, xhalf, yhalf, zoom, zfac, screens[0].pixwin

END
