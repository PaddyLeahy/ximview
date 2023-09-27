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
FUNCTION xph2hp, nside, x, y, proj, order, NPOLE=npole, SPOLE = spole, $
                    RING = ring, NESTED = nest, BOTH = both, R_OFF = tring
;
; J. P. Leahy 2008
;
; Converts "butterfly" grid pixel coords to healpix pixel coords
;
; Inputs:
;    nside: Usual HEALPix parameter
;    x, y : arrays of x and y pixel coordinates on the grid
;    {proj, /NPOLE, /SPOLE}: specifies "projection" of gridded image
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
npfacet = [[8, 4, 4, 11], $
           [5, 0, 3,  7], $
           [5, 1, 2,  7], $
           [9, 6, 6, 10]]
nprot   = [[3, 3, 0, 0], $
           [3, 3, 0, 0], $
           [2, 2, 1, 1], $
           [2, 2, 1, 1]]
spfacet = [[1, 6,  6, 0], $
           [5, 9, 10, 7], $
           [5, 8, 11, 7], $
           [2, 4,  4, 3]]
sprot   = [[1, 1, 2, 2], $
           [1, 1, 2, 2], $
           [0, 0, 3, 3], $
           [0, 0, 3, 3]]

nsmin1 = nside - 1

; Decide which projection we are in
proj_set = N_ELEMENTS(proj) GT 0
npole = KEYWORD_SET(npole) & spole = KEYWORD_SET(spole)
IF proj_set + npole + spole GT 1 THEN $
  MESSAGE, 'Please specify projection only once!'

IF ~proj_set THEN BEGIN
    proj = 'NPOLE'
    IF npole THEN proj = 'NPOLE'
    IF spole THEN proj = 'SPOLE'
ENDIF

; Check inputs:
xsize = SIZE(x)
ysize = SIZE(y)
IF ARRAY_EQUAL(xsize, ysize) EQ 0B THEN $
  MESSAGE, 'Mis-matched input x and y arrays'

nx = xsize[xsize[0]+2]
ngrid = 4*nside
IF MIN(x) LT 0 OR MIN(y) LT 0 OR MAX(x) GT ngrid - 1 OR MAX(y) GT ngrid - 1 $
  THEN MESSAGE, 'Some image pixel coordinates are out-of-range'

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
CASE STRUPCASE(proj) OF
    'NPOLE': BEGIN
        facet = npfacet[xx,yy]
        rot   = nprot[xx,yy]
    END
    'SPOLE': BEGIN
        facet = spfacet[xx,yy]
        rot   = sprot[xx,yy]
    END
    ELSE: MESSAGE, 'Unrecognised projection'
ENDCASE

IF nx GT 1 THEN BEGIN
    goodpix = WHERE(facet GE 0, COMPLEMENT = badpix)
    hpi = LONARR(nx,/NOZERO)

    IF badpix[0] GE 0 THEN hpi[badpix] = -1 ; Pixels are off sky

    IF goodpix[0] EQ -1 THEN RETURN, hpi

; find pixel index in facet
    xx = x[goodpix] MOD nside  &  yy = y[goodpix] MOD nside
    facet = facet[goodpix]
    
    rot = rot[goodpix]
    zz = xx
    ii = WHERE(rot EQ 1)
    IF ii[0] NE -1 THEN BEGIN
        xx[ii] = yy[ii]
        yy[ii] = nsmin1 - zz[ii]
    ENDIF
    ii = WHERE(rot EQ 2)
    IF ii[0] NE -1 THEN BEGIN
        xx[ii] = nsmin1 - xx[ii]
        yy[ii] = nsmin1 - yy[ii]
    ENDIF
    ii = WHERE(rot EQ 3)
    IF ii[0] NE -1 THEN BEGIN
        xx[ii] = nsmin1 - yy[ii]
        yy[ii] = zz[ii]
    ENDIF
ENDIF ELSE BEGIN ; just one pixel: faster code
    IF facet EQ -1 THEN RETURN, -1  
    xx = x MOD nside  &  yy = y MOD nside
    CASE rot OF
        0:                      ; no change
        1: BEGIN                ; 90 deg anticlockwise
            zz = xx
            xx = yy
            yy = nsmin1 - zz
        END
        2: BEGIN                ; 180 deg
            xx = nsmin1 - xx
            yy = nsmin1 - yy
        END
        3: BEGIN                ; 270 deg anticlockwise
            zz = xx
            xx = nsmin1 - yy
            yy = zz
        END
    ENDCASE
    goodpix = 0
    hpi = 0L
ENDELSE

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
    IF ring THEN BEGIN          ; make second array as we want both
        hpr = lindgen(nx)
        hpr[badpix] = -1L
        hpr[goodpix] = set_hpi_nest(nside,xx,yy)+facet*nside^2
        hpi = [[hpi], [hpr]]
    ENDIF ELSE hpi[goodpix] = set_hpi_nest(nside,xx,yy)+facet*nside^2
ENDIF

RETURN, hpi
END
