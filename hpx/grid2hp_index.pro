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
FUNCTION grid2hp_index, nside, x, y, proj, order, NPOLE=npole, SPOLE = spole, $
                        RING = ring, NESTED = nest, BOTH = both, $
                        R_OFF = tring, SILENT = silent
;+
; NAME:
;       GRID2HP_INDEX
;
; PURPOSE:
;       Converts grid pixel coords to healpix pixel coords
;
; CATEGORY:
;       Mapping
;
; CALLING SEQUENCE:
;
;       Result = GRID2HP_INDEX(Nside, X, Y, Proj, Order)
;
; INPUTS:
;       Nside: Parameter specifying size of HEALPix array
;
;       X:     Array of x pixel coordinates on the grid
;
;       Y:     Array of y pixel coordinates on the grid
;
; OPTIONAL INPUTS:
;       Proj:  'GRID', 'NPOLE', or 'SPOLE': defines input
;              "projection" (= ordering of HEALpix facets). Can
;              also be specified by keywords. If none of these are
;              specified, default is 'GRID'. 
;
;       Order: String describing output HEALPix sort order: 'RING' or
;              'NESTED', or 'BOTH'. These can also be specified by
;              keywords. If none of this are specified, the value
;              from the header is used if available, otherwise the
;              default is 'RING'. 
;   
; KEYWORD PARAMETERS:
;       NPOLE:  Defines input projection as "butterfly", with the
;               north coordinate pole at the centre. 
;
;       SPOLE:  As above, but with the south pole at the centre.
;
;       RING:   Set this flag if the output pixel indices are for Ring
;               order. (Redundant as this is the default).
;
;       NESTED: Set this flag if the output pixels are for Nested order.
;
;       R_OFF:  Scratch array. For multiple calls to GRID2HP_INDEX
;               with the same Nside, efficiency can be slightly
;               increased if this is set to an undefined variable on
;               the first call; don't change between calls.
;
;       SILENT  If set does not issue an error if pixel coords are out
;               of bounds, just returns coords = NaN
; OUTPUTS:
;       Returns an array of HEALPix pixel numbers. For "/BOTH" it is a
;       2-D array containing [[ring pixels], [nested pixels]]
;
; COMMON BLOCKS:
;       XY2PIX: standard HEALPix low-level common.
;
; SIDE EFFECTS:
;       Initialises XY2PIX if not done already.
;
; RESTRICTIONS:
;       Requires HEALPix IDL library.
;       For RING order the pixels should all be in the same facet row.
;
; EXAMPLE:
;       Converting a "GRID" ordered pixel coordinate [1001,2001]:
;  
;           PRINT, GRID2HP_INDEX(1024L, 1001, 2001, /NEST)
;
;       IDL prints:      1745686
;
; MODIFICATION HISTORY:
;      Written by:      J. P. Leahy, January 2008
;      Bug fix: indices in eastern half of facets in row 2 were returned
;      offset by 4*Nside. Fixed by JPL, 3 August 2009 
;
;-
COMPILE_OPT IDL2

gridfacet = [[ 6,  9, -1, -1, -1], $
             [ 1,  5,  8, -1, -1], $ 
             [-1,  0,  4, 11, -1], $
             [-1, -1,  3,  7, 10], $
             [-1, -1, -1,  2,  6]]
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

silent = KEYWORD_SET(silent)

; Decide which projection we are in
proj_set = N_ELEMENTS(proj) GT 0
npole = KEYWORD_SET(npole) & spole = KEYWORD_SET(spole)
IF proj_set + npole + spole GT 1 THEN $
  MESSAGE, 'Please specify projection only once!'

IF proj_set THEN proj = STRUPCASE(proj) ELSE BEGIN
    proj = 'GRID'
    IF npole THEN proj = 'NPOLE'
    IF spole THEN proj = 'SPOLE'
ENDELSE


; Check inputs:
xsize = SIZE(x)
ysize = SIZE(y)
IF ARRAY_EQUAL(xsize, ysize) EQ 0B THEN $
  MESSAGE, 'Mis-matched input x and y arrays'

nx = xsize[xsize[0]+2]
ngrid =  proj EQ 'GRID' ? 5*nside : 4*nside
maxgrid = ngrid - 1
bad = WHERE(x LT 0 OR y LT 0 OR x GT maxgrid OR y GT maxgrid, bcount)
IF bcount GT 0 && ~silent THEN $
  MESSAGE, 'Some image pixel coordinates are out of range'

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

; Quick return for out-of-range pixel
IF bcount EQ nx THEN BEGIN
    nan = !values.F_NAN
    both = ring && nest
    IF nx EQ 1 THEN idx = both ? [nan, nan] : nan $
    ELSE idx = both ? REPLICATE(nan, nx, 2) : REPLICATE(nan, nx)
    RETURN, idx
ENDIF

; Find facet
xx = x / nside  &  yy = y / nside
IF bcount GT 0 THEN BEGIN
    xx[bad] = 0  &  yy[bad] = 0 ; avoids indexing problems below
ENDIF
CASE proj OF
    'GRID': BEGIN
        facet = gridfacet[xx,yy]
        rot   = nx GT 1 ? INTARR(nx) : 0
    END
    'NPOLE': BEGIN
        facet = npfacet[xx,yy]
        rot   = nprot[xx,yy]
    END
    'SPOLE': BEGIN
        facet = spfacet[xx,yy]
        rot   = sprot[xx,yy]
    END
ENDCASE

IF nx GT 1 THEN BEGIN
    hpi = LONARR(nx,/NOZERO)
    IF proj EQ 'GRID' THEN BEGIN
        IF bcount GT 0 THEN facet[bad] = -1
        goodpix = WHERE(facet GE 0, COMPLEMENT = badpix)
        IF badpix[0] GE 0 THEN hpi[badpix] = -1 ; Pixels are off sky
    ENDIF ELSE BEGIN ; find off-sky pixels in butterfly projections
        dg1 = x + y + REPLICATE(1-4*nside,nx)
        dg2 = y - x
        goodpix = WHERE((dg1 GE -nside AND dg1 LE nside) OR $
                        (dg2 GE -nside AND dg2 LE nside), COMPLEMENT = badpix)
        dg1 = 0 & dg2 = 0
    ENDELSE
    IF goodpix[0] EQ -1 THEN BEGIN
        IF ring AND nest THEN hpi = [[hpi],[hpi]]
        RETURN, hpi
    ENDIF

; find pixel index in facet
    xx = x[goodpix] MOD nside  &  yy = y[goodpix] MOD nside
    facet = facet[goodpix]
    
    IF proj NE 'GRID' THEN BEGIN
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
    ENDIF
ENDIF ELSE BEGIN ; just one pixel: faster code
    IF proj EQ 'GRID' THEN bad = facet EQ -1 ELSE BEGIN
      dg1 = x+y-4*nside + 1
      dg2 = y-x
      right = x GE 2*nside
      top   = y GE 2*nside
      bad =    (~right &&  top && (dg1 LT -nside || dg1 GE nside)) $
            || (~right && ~top && (dg2 LT -nside || dg2 GE nside)) $     
            || ( right && ~top && (dg1 LE -nside || dg1 GT nside)) $     
            || ( right &&  top && (dg2 LE -nside || dg2 GT nside))     
    ENDELSE
    IF bad THEN BEGIN
        IF ring AND nest THEN hpi = [-1,-1] ELSE hpi = -1
        RETURN, hpi
    ENDIF
 
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
    badpix = -1
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
        IF badpix[0] NE -1 THEN hpr[badpix] = -1L
        hpr[goodpix] = set_hpi_nest(nside,xx,yy)+facet*nside^2
        hpi = [[hpi], [hpr]]
    ENDIF ELSE hpi[goodpix] = set_hpi_nest(nside,xx,yy)+facet*nside^2
ENDIF

RETURN, hpi
END
