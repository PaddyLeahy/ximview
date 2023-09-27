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
FUNCTION set_xy_nest, nside, iface
;
; Sets up look-up table for one facet. Bit manipulation code is from
; HEALPix routine nest2ring.pro 
;
; Arguments: nside: usual HEALPix parameter
;            iface: index array for pixels in facet
;
; Returns: 1-D index for square facet array.
;

;     initiates the array for the pixel number -> (x,y) mapping
    COMMON pix2xy, pix2x, pix2y
    sz = SIZE(pix2x)
    IF (sz(sz(0)+1) EQ 0) THEN init_pix2xy ; initiate pix2x and pix2y

    ip_low   = iface MOD 1024   ; content of the last 10 bits
    ip_trunc = iface / 1024     ; truncation of the last 10 bits
    ip_med   = ip_trunc MOD 1024 ; content of the next 10 bits
    ip_hi    = ip_trunc / 1024  ; content of the high weight 10 bits
    ip_trunc = 0                ; free memory
    
;     computes (x,y) coordinates on the face (actually x is in reverse order):
    ix = 1024*pix2x(ip_hi) + 32*pix2x(ip_med) + pix2x(ip_low)
    iy = 1024*pix2y(ip_hi) + 32*pix2y(ip_med) + pix2y(ip_low)
    ip_hi = 0 & ip_med = 0 & ip_low = 0 ; free memory

; get array index for each hp pixel
    RETURN, iy*nside + (nside-1) - ix
END
