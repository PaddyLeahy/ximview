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
FUNCTION hpgrid, hparray, header, order, proj, RING = ring, NESTED = nest, $
                 NPOLE = npole, SPOLE = spole, UNIT = unit, $
                 HP_INDEX=hp_index, VERBOSE = VERBOSE
;+
; NAME:
;       HPGRID
;
; PURPOSE:
;       Sorts healpix 1D pixel list into a 2D grid, as defined by
;       M. R. Calabretta and B. F. Roukema, 2007: "Mapping on the HEALPix
;       grid", MNRAS, 381, 865--872 
;
;       Output grid pixels correspond exactly to HEALPix pixels or are set
;       to NaN or -1 depending on data type. For the default ordering,
;       HEALPix facet 6 (centred on theta=pi/2, phi = pi) is included twice,
;       at the bottom left and top right of the grid.
;
;       The default 'GRID' projection is the 45-degree rotated pseudo-
;       cylindrical projection illustrated by Calabretta & Roukema,
;       with theta = pi/2, phi = 0 at the centre of an N_grid =
;       5*N_side square grid, with North (theta = 0) towards the upper
;       left and phi increasing towards the bottom left. This has FITS
;       WCS code "HPX".

;       The other projection available is the the polar "butterfly"
;       projection mentioned by Calabretta & Roukema (the two are
;       related by shifts and rotations of the HEALPix facets). This has
;       a pole (theta = 0 or pi) at the centre of an ngrid = (4*nside)
;       square grid, with phi = 0 at the bottom and increasing
;       clockwise as displayed by IDL TV (which puts [0,0] at the
;       bottom left).  Following the "HPXcvt" program by
;       M. Calabretta, the unofficial WCS code "XPH" is used.
;
; CATEGORY:
;       Array manipulation, mapping
;
; CALLING SEQUENCE:
;
;       Result = HPGRID(HParray, Header, Order, Proj)
;
; INPUTS:
;       HParray:  The HEALPix dataset to grid. May be
;                 o  A 1-D array (must be a valid HEALPix length N_pix)
;                 o  An N-D array with the first dimension N_pix
;                 o  An N-D array of pointers, each pointing to a 1-D
;                    array with the same length N_pix (but not
;                    necessarily the same data type; e.g. floats and
;                    integers can be handled together).
;
; OPTIONAL INPUTS:
;       Header:   FITS-format header describing the data. (Modified on
;                 output). If not specified a header is constructed
;                 using the dimensions of HParray 
;
;       Order:    String describing HEALPix sort order: 'RING' or
;                 'NESTED'. These can also be specified by
;                 keywords. If none of this are specified, the value
;                 from the header is used if available, otherwise the
;                 default is 'RING'. 
;
;       Proj:     'GRID', 'NPOLE', or 'SPOLE': defines output
;                 "projection" (= ordering of HEALpix facets). Can
;                 also be specified by keywords. If none of these are
;                 specified, default is 'GRID'. 
;
; KEYWORD PARAMETERS:
;       RING:     Set this flag if the input array is in Ring
;                 order. (Redundant as this is the default).
;
;       NESTED:   Set this flag if the input array is in Nested order.
;
;       NPOLE:    Set to make output map the "butterfly" ordering,
;                 with the north coordinate pole at the centre.
;
;       SPOLE:    As above, but with the south pole at the centre.
;
;       UNIT:     Brightness unit. If specified, overrides header value,
;                 if any. 
;
;       HP_INDEX: Returns a 2-D grid giving the HEALPix pixel numbers
;                 corresponding to each grid pixel, or -1 for grid pixels
;                 not mapped onto the sphere.
;
;       VERBOSE:  Enables diagnostic printout (mostly timing).
;
; OUTPUTS:
;       The returned value is either a 2-D image of size ngrid x ngrid,
;       or (if the input was multidimensional) an array of pointers to
;       such images. The pointer array has the same dimensions as the
;       second and higher dimension of the input. 
;
; OPTIONAL OUTPUTS:
;       Header:  If specified this returns a FITS-like header with
;                astrometry keywords set assuming the data is stored
;                in a FITS primary HDU or image extension
;
;                NB: header does not fully conform to the FITS
;                standard, e.g. TTYPE* and TUNIT* keywords are
;                preserved for reference even though these should not
;                be present in a primary HDU. 
;
; SIDE EFFECTS:
;       For multi-dimensional output, creates heap variables pointed
;       to by the elements of "grid". Remember to free them after use!
;
; RESTRICTIONS:
;       32-bit IDL installations cannot handle arrays with N_side =
;       8192, nominally the largest allowed for HEALPix.
;
; EXAMPLE:
;
;       silly = RANDOMN(seed, 12*32L*32L) ; length for N_side = 32
;       grid  = HPGRID(silly)
;
; MODIFICATION HISTORY:
;       Written by:     J. P. Leahy, November 2007
;       Jan 2008:       Pointer input and output
;       March 2008:     Added HEALPix blanks and enabled mixed-type
;                       input.
;
;-
COMPILE_OPT IDL2
ON_ERROR, 1

start = SYSTIME(1)
verbose = KEYWORD_SET(verbose)

; Check inputs
intype = SIZE(hparray,/TYPE)
S = SIZE(hparray)
nextra = S[0] - 1
ntot = S[S[0]+2]
npix = S[1]
ndim = ntot / npix         ; total size / first dimension
inpointer = intype EQ 10

IF inpointer THEN BEGIN 
    ndim = ntot
; Squash extra dimensions
    extra_dims = S[1:S[0]]
    IF nextra GT 0 THEN hparray = REFORM(hparray,ndim,/OVERWRITE)
    nextra = nextra + 1
    S = SIZE(*hparray[0])
    IF S[0] NE 1 THEN MESSAGE, $
      'Pointers must each reference a single HEALPix map'
    type = INTARR(ndim)
    type[0] = S[2]
    FOR i=1,ndim-1 DO BEGIN
        T = SIZE(*hparray[i])        
        IF ARRAY_EQUAL(S[0:S[0]], T[0:T[0]]) EQ 0B THEN MESSAGE, $
          'All referenced HEALPix arrays must be the same size'
        type[i] = T[2]
    ENDFOR
    npix = S[1]
ENDIF ELSE BEGIN
    extra_dims = nextra GT 0 ? S[2:S[0]] : 1
    IF nextra GT 0 THEN hparray = REFORM(hparray,npix,ndim,/OVERWRITE)
    type = REPLICATE(intype,ndim)
ENDELSE

fpdata = (type GE 4 AND type LE 6) OR type EQ 9
IF MAX(fpdata) EQ 0 THEN nan = INTARR(ndim) ELSE BEGIN
    nan = FLTARR(ndim)          ; Zeros
    nan[WHERE(fpdata)] = !values.F_NAN
ENDELSE

ispointer = ndim GT 1 || inpointer

nside = npix2nside(npix)
IF nside EQ -1 THEN MESSAGE, $
  'First data dimension not a valid Healpix size'

; Initial data:

fmt = "('HPgrid: ',F7.3,' seconds: ',A)"

; coordinate of the lowest corner of each face, in grid, N-pole and
; S-pole projections. (NB: Face 12-15 in butterfly projections is
; second copy of face 4-7).

xcgrd = [1, 0, 3, 2, 2, 1, 4, 3, 2, 1, 4, 3] 
ycgrd = [2, 1, 4, 3, 2, 1, 4, 3, 1, 0, 3, 2] 
xcnp  = [1, 1, 2, 2, 1, 0, 2, 3, 0, 0, 3, 3, 2, 0, 1, 3]
ycnp  = [1, 2, 2, 1, 0, 2, 3, 1, 0, 3, 3, 0, 0, 1, 3, 2]
xcsp  = [0, 0, 3, 3, 2, 0, 1, 3, 1, 1, 2, 2, 1, 0, 2, 3]
ycsp  = [3, 0, 0, 3, 3, 2, 0, 1, 2, 1, 1, 2, 3, 1, 0, 2]

; Facet rotations for butterfly projections:
ronp  = [3, 2, 1, 0, 3, 2, 1, 0, 3, 2, 1, 0, 0, 3, 2, 1]
rosp  = [0, 1, 2, 3, 3, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 3]

lut = ARG_PRESENT(hp_index)

; Do we have a header to update?
ishead = N_ELEMENTS(header) GT 0 

; Define projection:

proj_set = N_ELEMENTS(proj) GT 0
npole = KEYWORD_SET(npole) & spole = KEYWORD_SET(spole)
IF proj_set + npole + spole GT 1 THEN $
  MESSAGE, 'Please specify projection only once!'

IF proj_set THEN BEGIN
    npole = STRCMP(proj, 'NPOLE', 4, /FOLD_CASE)
    spole = STRCMP(proj, 'SPOLE', 4, /FOLD_CASE)
ENDIF

butterfly = npole || spole
IF npole THEN BEGIN
    proj = 'NPOLE'
    xcbut = xcnp
    ycbut = ycnp
    robut = ronp
ENDIF
IF spole THEN BEGIN
    proj = 'SPOLE'
    xcbut = xcsp
    ycbut = ycsp
    robut = rosp
ENDIF

IF ~butterfly THEN proj = 'GRID'

; Set "ring" switch according to ordering:

inhead = 0
IF ishead THEN BEGIN
    horder = SXPAR(header, 'ORDERING', COUNT = count)
    inhead = count GT 0
ENDIF

ring = KEYWORD_SET(ring)
nest = KEYWORD_SET(nest)
order_set = N_ELEMENTS(order) GT 0

spec = order_set + ring + nest
CASE spec OF
    0: IF  inhead THEN BEGIN
            order = horder
            order_set = 1
        ENDIF
    1:   ; Do nothing
    ELSE: MESSAGE, 'Please specify ordering only once!'
ENDCASE

IF order_set THEN nest = STRCMP(order,'NEST',4,/FOLD_CASE)
order = nest ? 'NESTED' : 'RING'
ring  = order EQ 'RING'

IF inhead && ~STRCMP(horder, order, 4, /FOLD_CASE) $
  THEN MESSAGE, /INFORMATIONAL, $
  'Input parameter overriding HEALPix order specified in header' 

IF butterfly THEN BEGIN
    ngrid = 4 * nside
    xcorn = xcbut
    ycorn = ycbut
ENDIF ELSE BEGIN
    ngrid = 5 * nside
    xcorn = xcgrd
    ycorn = ycgrd
ENDELSE

; Make header if there is something to return it to:
IF ARG_PRESENT(header) THEN BEGIN
    IF ~ishead THEN MKHDR, header, type[0], [ngrid, ngrid]
    header = grid_header(header, nside, proj, nextra, extra_dims, $
                         UNIT=unit)
ENDIF

nface = nside*nside
iface = LINDGEN(nface)
ix = FIX(iface MOD nside)  ; Convert to short integers as max is 8191
iy = FIX(iface / nside)    ;
iface = 0

IF ring THEN tring = set_ring_offset(nside)

IF ~ring THEN BEGIN
    hpi = set_hpi_nest(nside, ix, iy) - nface
    ix = 0  &  iy = 0
ENDIF

IF verbose THEN PRINT, SYSTIME(1)-start, "Indices done", FORMAT=fmt

; Make output arrays:
IF ispointer THEN BEGIN
    grid = PTRARR(ndim, /ALLOCATE_HEAP)
    FOR idim = 0L,ndim-1 DO *grid[idim] =   $
      MAKE_ARRAY(ngrid,ngrid,TYPE = type[idim], VALUE=nan[idim])
ENDIF ELSE grid = MAKE_ARRAY(ngrid,ngrid,TYPE = type, VALUE=nan)

IF lut THEN hp_index = MAKE_ARRAY(ngrid,ngrid, TYPE = 3, VALUE = -1L)

IF ispointer THEN face = PTRARR(ndim, /ALLOCATE_HEAP)
; ELSE face = MAKE_ARRAY(nface,ndim,TYPE=type) ; doesn't seem necessary

ns4 = 4*nside

FOR facet = 0,11 DO BEGIN
    x0 = xcorn[facet]*nside & y0 = ycorn[facet]*nside
    x1 = x0 + nside-1 & y1 = y0 + nside-1
    IF ring THEN BEGIN
        row = facet/4
        IF 4*row EQ facet THEN BEGIN
            set_hpi_ring, nside, ix, iy, row, ringnum, ringpos, offset
            ringoff = tring[TEMPORARY(ringnum)]
        ENDIF
        hpi = ringoff + ringpos
        ringpos += offset
        IF facet EQ 4 THEN ringpos MOD= ns4 
    ENDIF ELSE hpi += nface ; For nested just shift indices by one facet.

    IF ispointer THEN FOR idim=0L,ndim-1 DO BEGIN
        fptr = face[idim]
        IF inpointer THEN ptr = hparray[idim]
        *fptr = inpointer ? (*ptr)[hpi] : hparray[hpi,idim]
    ENDFOR ELSE face = hparray[hpi]

    IF ~butterfly THEN BEGIN
        IF lut THEN hp_index[x0:x1,y0:y1] = hpi

        IF ispointer THEN FOR idim=0L,ndim-1 DO BEGIN
            gptr = grid[idim]
            fptr = face[idim]
            (*gptr)[x0:x1,y0:y1] = *fptr
        ENDFOR ELSE grid[x0:x1,y0:y1] = face

    ENDIF ELSE BEGIN    
        IF lut THEN BEGIN
          hpi = REFORM(hpi, nside, nside, /OVERWRITE)
          hp_index[x0:x1,y0:y1] = ROTATE(hpi,robut[facet])
        ENDIF  
        IF facet GE 4 && facet LE 7 THEN BEGIN
            f8 = facet+8
            x2 = xcbut[f8]*nside & y2 = ycbut[f8]*nside
            x3 = x2+nside-1 & y3 = y2 + nside-1
            IF lut THEN hp_index[x2:x3,y2:y3] = ROTATE(hpi,robut[f8])
        ENDIF

        IF ispointer THEN FOR idim=0L,ndim-1 DO BEGIN
            fptr = face[idim]
            gptr = grid[idim]
            dolut = lut && idim EQ 0
            *fptr = REFORM(*fptr, nside, nside, /OVERWRITE)

            (*gptr)[x0:x1,y0:y1] = ROTATE(*fptr,robut[facet])

            IF facet GE 4 && facet LE 7 THEN $
              (*gptr)[x2:x3,y2:y3] = ROTATE(*fptr,robut[f8])

        ENDFOR ELSE BEGIN
            face = REFORM(face, nside, nside, /OVERWRITE)

            grid[x0:x1,y0:y1] = ROTATE(face,robut[facet])

            IF facet GE 4 && facet LE 7 THEN $
              grid[x2:x3,y2:y3] = ROTATE(face,robut[f8])

        ENDELSE
    ENDELSE

ENDFOR                          ; End of loop over facets
hpi = 0  &  ix = 0  &  iy = 0   ; save space
IF ispointer THEN PTR_FREE, face ELSE face = 0
IF verbose THEN PRINT, SYSTIME(1)-start, "Grids filled", FORMAT=fmt
IF verbose THEN HELP, /MEMORY

IF butterfly THEN BEGIN
; Blank V-shaped regions around edge containing duplicate copies of
; facets 4 to 7

    ; Get list of x,y pixel indices for middle bottom of grid:
    ns2 = 2*nside
    list = LINDGEN(ns2*nside)
    y1 = FIX(list / ns2)
    x1 = FIX((list MOD ns2) + nside)
    list = 0
    nss = FIX(nside)
    dx = ABS(x1 - 2S*nss)       ; x offset from grid centre
    dy = nss - y1               ; y offset from top of first row of facets

    bad = WHERE(TEMPORARY(dx) LT TEMPORARY(dy))  ; duplicate pixels
    x1 = x1[bad]
    y1 = y1[bad]
    bad = [x1, ngrid-1-x1,         y1, ngrid-1-y1] + $
          [y1, ngrid-1-y1, ngrid-1-x1,         x1]*ngrid
    x1 = 0 & y1 = 0 

    IF lut THEN hp_index[bad] = -1L
    nbad = 4*nside^2
    IF ispointer THEN FOR idim = 0L,ndim-1 DO BEGIN
        gptr = grid[idim]
        nanlist = REPLICATE(nan[idim], nbad)
        (*gptr)[bad] = nanlist
    ENDFOR ELSE BEGIN
        nanlist = REPLICATE(nan, nbad)
        grid[bad] = nanlist
    ENDELSE
    nanlist = 0  &  bad = 0

    IF verbose THEN PRINT, SYSTIME(1)-start, "Bad pixels blanked", FORMAT=fmt
ENDIF ELSE BEGIN
; Wrap last equatorial longitude facet:
    x0 = 4L*nside  &  y0 = 4L*nside 

    IF lut THEN hp_index[0:nside-1,0:nside-1] = hp_index[x0:*,y0:*]

    IF  ispointer THEN FOR idim = 0L,ndim-1 DO BEGIN
        gptr = grid[idim]
        (*gptr)[0:nside-1,0:nside-1] = (*gptr)[x0:*,y0:*]
    ENDFOR ELSE grid[0:nside-1,0:nside-1] = grid[x0:*,y0:*]
ENDELSE

; Restore squashed dimensions:
hparray = inpointer ? REFORM(hparray, extra_dims, /OVERWRITE) $
                    : REFORM(hparray, S[1:S[0]], /OVERWRITE)

IF ispointer THEN BEGIN
    FOR idim = 0L, ndim-1 DO BEGIN
        gptr = grid[idim]
        *gptr = REFORM(*gptr, ngrid, ngrid, /OVERWRITE)
    ENDFOR
    grid = REFORM(grid, extra_dims, /OVERWRITE) 
ENDIF ELSE grid = REFORM(grid, ngrid, ngrid, /OVERWRITE)

IF verbose THEN PRINT, SYSTIME(1)-start, "All done", FORMAT=fmt
IF verbose THEN HELP, /MEMORY
RETURN, grid
END


