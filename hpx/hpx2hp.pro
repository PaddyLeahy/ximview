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
FUNCTION hpx2hp, nside, x, y, order, $
                 RING = ring, NESTED = nest, BOTH = both, R_OFF = tring
;
; J. P. Leahy 2008
;
; Converts HEALPix grid (HPX projection) pixel coords to healpix pixel coords
;
; Inputs:
;    nside: Usual HEALPix parameter
;    x, y : arrays of x and y pixel coordinates on the grid
;    {order, /NESTED, /RING, /BOTH}: specifies required HEALPix sort order
;   
;    tring: Array giving pixel number of first pixel in each ring, for
;           order = 'RING' (calculated internally if not provided).
;
; Returns: Array of HEALPix pixel numbers. For "/BOTH" it is a 2-D
;          array containing [[ring pixels], [nested pixels]]
;
; Revision history:
;  3 August 2009: RING pixel index bug fixed (as per grid2hp_index)
;
gridfacet = [[ 6,  9, -1, -1, -1], $
             [ 1,  5,  8, -1, -1], $ 
             [-1,  0,  4, 11, -1], $
             [-1, -1,  3,  7, 10], $
             [-1, -1, -1,  2,  6]]

nsmin1 = nside - 1

; Check inputs:
xsize = SIZE(x)
ysize = SIZE(y)
IF ARRAY_EQUAL(xsize, ysize) EQ 0B THEN $
  MESSAGE, 'Mis-matched input x and y arrays'

nx = xsize[xsize[0]+2]
ngrid =  5*nside
IF MIN(x) LT 0 OR MIN(y) LT 0 OR MAX(x) GT ngrid - 1 OR MAX(y) GT ngrid - 1 $
  THEN MESSAGE, 'Some image pixel coordinates are out of range'

; Decide which sort order we want:
order_set = N_ELEMENTS(order) GT 0
ring = KEYWORD_SET(ring)  &  nest = KEYWORD_SET(nest)
both = KEYWORD_SET(both)
keyword = ring OR nest OR both
IF order_set + keyword GT 1 THEN $
  MESSAGE, 'Please specify ordering only once!'

IF ~order_set THEN BEGIN
    IF (ring AND nest) OR both THEN order = 'BOTH' ELSE $
      order = ring ? 'RING' : 'NESTED'
ENDIF
test = STRMID(STRUPCASE(order),0,4)
CASE test OF
    'BOTH': BEGIN
        ring = 1B
        nest = 1B
    END
    'RING': ring = 1B
    'NEST': nest = 1B
    ELSE: MESSAGE, 'Unrecognised pixel order'
ENDCASE

; Find facet
xx = x / nside  &  yy = y / nside
facet = gridfacet[xx,yy]

; Reject off-sky pixels:
IF nx GT 1 THEN BEGIN
    goodpix = WHERE(facet GE 0, COMPLEMENT = badpix)
    hpi = LONARR(nx,/NOZERO)

    IF badpix[0] GE 0 THEN hpi[badpix] = -1 ; Pixels are off sky

    IF goodpix[0] EQ -1 THEN RETURN, hpi

; find pixel index in facet
    xx = x[goodpix] MOD nside  &  yy = y[goodpix] MOD nside
    facet = facet[goodpix]
    
ENDIF ELSE BEGIN ; just one pixel: faster code
    IF facet EQ -1 THEN RETURN, -1  
    xx = x MOD nside  &  yy = y MOD nside
    goodpix = 0
    hpi = 0L
ENDELSE

; Translate pixel numbers:
IF ring THEN BEGIN
    IF ~N_ELEMENTS(tring) THEN tring = set_ring_offset(nside)
    row = facet / 4
    IF MIN(row) NE MAX(row) THEN BEGIN
        MESSAGE, "Can't handle pixels in different facet rows"
    ENDIF
    set_hpi_ring, nside, xx, yy, row[0], ringnum, ringpos, offset
    ringpos += (facet MOD 4) * offset
    ringpos MOD= 4*nside
    hpi[goodpix] = tring[ringnum] + ringpos
ENDIF

IF nest THEN BEGIN
    IF ring THEN BEGIN         ; Make second pixel array, we want both
        hpr = lindgen(nx)
        hpr[badpix] = -1L
        hpr[goodpix] = set_hpi_nest(nside, xx, yy) + facet*nside^2
        hpi = [[hpi], [hpr]]
    ENDIF ELSE hpi[goodpix] = set_hpi_nest(nside, xx, yy) + facet*nside^2
ENDIF

RETURN, hpi
END
