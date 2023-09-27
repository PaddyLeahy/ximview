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
FUNCTION set_ring_offset, nside
;
; Finds pixel offset for start of each ring in a RING-ordered HEALPix array.
;
iside = LINDGEN(nside)+1
; number of pixels in each ring:
nring = 4*[iside,REPLICATE(nside,2*nside-1),REVERSE(iside)]
; Get index for start pixel of each ring:
tot = TOTAL(nring,/CUMULATIVE)
; TOTAL should be /INTEGER but not available until version 6.1
tring = LONARR(4*nside)
tring[1] = tot                  ; Fills in whole array after [1].

RETURN, tring
END
