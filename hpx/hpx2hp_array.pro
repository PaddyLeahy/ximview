FUNCTION hpx2hp_array, hpx, order
;+
; NAME: 
; 
;    HPX2HP_ARRAY
;
; PURPOSE:
;
;    Converts healpix grid map to a healpix array
;
; CATEGORY:
;
;    Healpix grid routines
;
; CALLING SEQUENCE:
;
;    hparray = hpx2hp_array(hpx,order)
;
; INPUTS:
;
;    hpx      Healpix grid representation of the sky, i.e. a 5*nside
;             square with the diagonal representing the equator.
;    order    'RING' or 'NEST' to set the order of the output HP array
;
;
; OPTIONAL INPUTS:
;
;
;
; KEYWORD PARAMETERS:
;
;
;
; OUTPUTS:
;
;   returns a healpix array with same nside and data type as input.
;
;
; OPTIONAL OUTPUTS:
;
;
;
; RESTRICTIONS:
;
;    Input should be 2D array with appropriate size for healpix grid map.
;
; PROCEDURE:
;
;   Uses grid2hp_index.pro to convert grid pixel coords to healpix pixel coords.
;
; EXAMPLE:
;
;    > read_fits_map, 'dipole_1024_ring.fits', dp, primary_header,ext_hdr
;    > dpx = hpgrid(dp,ext_hdr) ; convert to grid
;    > dp2 = hpx2hp_array(dpx,'RING') ; convert back
;    >                      ; Check that we got back to the starting point:
;    > PRINT, ARRAY_EQUAL(dp, dp2)
;         1
;    >                      ; Write out result, specifying a few
;    >                      ;  details for the header:
;    > write_fits_map, 'my_map.fits', dp2, UNIT='K', /RING, COORDSYS='C'
;
; MODIFICATION HISTORY:
;
;   Written by J. P. Leahy  Version 1.0 17/06/2020
;
;-
; Coordinates of facet BLC in units of nside for facets 0 through 12,
; where facet 12 is the repeat of facet 6 at the top right of grid.
fx = [1,0,3,2, 2,1,0,3, 2,1,4,3, 4]
fy = [2,1,4,3, 2,1,0,3, 1,0,3,2, 4]

; check input files match in size:
s1= SIZE(hpx) ; Size function returns metadata about its argument, 
              ; eg. # dimensions.

IF S1(0) NE 2 THEN MESSAGE, 'Input array should be 2-D'
; Note: could be upgraded to handle healpix N-D object where each 2D
;       slice is a HPX grid. For the moment allow only 2D.

nx = S1(1)
ny = S1(2)
IF nx NE ny THEN MESSAGE, 'Input array should be square'
nside = nx / 5

IF 5*nside NE nx THEN MESSAGE, 'Input array is not 5*nside squared'

npix = NSIDE2NPIX(nside, ERROR=err)
IF err THEN MESSAGE, 'Input array has improper N_side'

type = S1(3) ; data type code e.g. 4 for floating point.

; Create grids of X and Y coordinates for one facet 
xx = LINDGEN(nside) ; array [0,1,....,nside-1]
yy = xx
yy = REFORM(yy,1,nside) ; turn into 1 x nside array
yy = REBIN(yy,nside,nside) ; replicate columns 
xx = REBIN(xx,nside,nside) ; replicate rows

hparray = MAKE_ARRAY(npix,TYPE=type) ; creates array of same type as input.

; Loop over facets
FOR ifacet = 0, 11 DO BEGIN

; Create look-up table from healpix grid to healpix array:
   i1 = nside*fx[ifacet]
   i2 = i1 + nside - 1L
   x = xx + i1
   j1 = nside*fy[ifacet]
   j2 = j1 + nside - 1L
   y = yy + j1
   ipx = grid2hp_index(nside, x, y,'GRID',order)

; Fill in values
   IF ifacet LT 12 THEN hparray[ipx] = hpx[x,y] ELSE BEGIN

; Facet 12 is a duplicate of facet 6 but may be differently blanked:
      f6 = hparray[ipx]
      f12 = hpx[x,y]
; Find blanked pixels:
      bad = WHERE(~FINITE(f6) OR f6 EQ !healpix.bad_value, nbad) ; NaNs etc
      IF nbad GT 0 THEN f6[bad] = f12[bad]
      bad = WHERE(~FINITE(f12) OR f12 EQ !healpix.bad_value, nbad)
      IF nbad GT 0 THEN f12[bad] = f6[bad]

; Could also be blanked with zeros:
      bad = WHERE(f6 EQ 0, nbad)
      IF nbad GT 0 THEN f6[bad] = f12[bad]
      bad = WHERE(f12 EQ 0, nbad)
      IF nbad GT 0 THEN f12[bad] = f6[bad]

      mismatch = WHERE(f6 NE f12, nmiss)
      IF nmiss GT 0 THEN PRINT, nmiss, FORMAT= $
            "('Warning:',I8,' pixels in overlap facets mismatched')"
      hparray[ipx] = f6
   ENDELSE
ENDFOR

RETURN, hparray

END
