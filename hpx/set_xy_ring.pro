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
FUNCTION set_xy_ring, nside, npix, ngrid, start, fmt
;
; Set up look-up table for rings:
;
; Arguments: nside, npix: Usual HEALPix parameters
;            ngrid      : Number of pixels on side of grid 
;            start,fmt  : start time in seconds and message format for
;                         debug/optimisation
; Returns: 1-D index for grid.
;

; Useful constants to avoid in-loop computations:
    ns4 = 4*nside & ns3 = 3*nside & ns2 = 2*nside
    off1 = ns3 - 1 & off2 = 9*nside - 1

    iside = LINDGEN(nside)+1
; 1/4 of the number of pixels in each ring:
    nring = [iside,REPLICATE(nside,2*nside-1),REVERSE(iside)]
; Actual number:
    nr4 = 4*nring
; Get index for start pixel of each ring:
    tot = TOTAL(nr4,/CUMULATIVE)
; TOTAL should be /INTEGER but not available until version 6.1
    tring = LONARR(ns4)
    tring[1] = tot              ; Fills in whole array after [1].

    ix = LONARR(npix,/NOZERO)
    iy = LONARR(npix,/NOZERO)

; Index for pixel number on the equatorial rings, modulo starting offset.
    eq_ring = (15*nside/2) - 1 - LINDGEN(ns4)

    FOR i=1L,nside DO BEGIN
        im1 = i - 1L
        ix[tring[im1]] = REPLICATE(i,nr4[im1])
        ipix = LINDGEN(nring[im1])
        iy[tring[im1]] = $
          [off1-ipix, (off1-nside)-ipix, (off1-ns2)-ipix, (off1+nside)-ipix]
    ENDFOR
    FOR i=nside+1L,ns3-1L DO BEGIN
        im1 = i - 1L
;        iy[tring[im1]] = REPLICATE(i,ns4)
;        ix[tring[im1]] = eq_ring
        ix[tring[im1]] = REPLICATE(i,ns4)
        iy[tring[im1]] = (eq_ring - (i/2L)) MOD ns4
    ENDFOR 
    FOR i=ns3,ns4-1L  DO BEGIN
        im1 = i - 1L
        ix[tring[im1]] = REPLICATE(i,nr4[im1])
        ipix = LINDGEN(nring[im1])
        iy[tring[im1]] = [(off2-i)-ipix, (off2-i-nside)-ipix, $
                          (off2-i-ns2)-ipix, (off2-i-ns3)-ipix] MOD ns4
    ENDFOR
    PRINT, SYSTIME(1)-start, "ix, iy done", FORMAT=fmt

    ix += iy
    PRINT, SYSTIME(1)-start, "ix offset", FORMAT=fmt

    wrap = WHERE(ix LT ns2)
    iy[wrap] += ns4
    ix[wrap] += ns4
    wrap = 0

    RETURN, TEMPORARY(iy)*ngrid + TEMPORARY(ix)  - ns2
END
