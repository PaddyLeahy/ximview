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
; Module cut4grid
;
; Contents:
;   rotate_coord   Reorders X and Y coordinate lists to achieve rotation
;   checklim       Checks to see if pixels in a facet extend the
;                  limits of the gridded data
;   cut4grid       Main routine.
;
; This module contains a standard IDL documentation header after the
; declaration of cut4grid
;
PRO rotate_coord, x, y, idx, n, rot
;
; Rotates pixel coordinates in a square n by n grid by rotation rot
;
m = n-1
CASE rot OF
    0: ; do nothing
    1: BEGIN
        temp = x[idx]
        x[idx] = m - y[idx]
        y[idx] = temp
    END
    2: BEGIN
        x[idx] = m - x[idx]
        y[idx] = m - y[idx]
    END
    3: BEGIN
        temp = y[idx]
        y[idx] = m - x[idx]
        x[idx] = temp
    END
ENDCASE
END

PRO checklim, ix, iy, xcorn, ycorn, mm, face, idx, mmp
;
; Checks pixel limits in a facet if needed
;
test = (xcorn[face] EQ mm[0]) + 2*(xcorn[face] EQ mm[1])
CASE test OF
    0:                          ; Not a possible extremum
    1: mmp[0] = MIN(ix[idx]) < mmp[0]
    2: mmp[1] = MAX(ix[idx]) > mmp[1]
    3: BEGIN                    ; Max & min in same facet
        m1 = MIN(ix[idx], MAX = m2)
        mmp[0] = m1 < mmp[0]
        mmp[1] = m2 > mmp[1]
    END
ENDCASE
test = (ycorn[face] EQ mm[2]) + 2*(ycorn[face] EQ mm[3])
CASE test OF
    0:                          ; Not a possible extremum
    1: mmp[2] = MIN(iy[idx]) < mmp[2]
    2: mmp[3] = MAX(iy[idx]) > mmp[3]
    3: BEGIN                    ; Max & min in same facet
        m1 = MIN(iy[idx], MAX = m2)
        mmp[2] = m1 < mmp[2]
        mmp[3] = m2 > mmp[3]
    END
ENDCASE

END

FUNCTION cut4grid, pixels, hparray, header, order, proj, $
                   RING = ring, NESTED = nest, NPOLE = npole, SPOLE = spole, $
                   UNIT = unit, HP_INDEX=hp_index, VERBOSE = VERBOSE
;+
; NAME:
;       CUT4GRID
;
;
; PURPOSE:
;       Sorts HEALPix cut4-format 1D pixel list into a 2D grid, as
;       defined by M. R. Calabretta and B. F. Roukema, 2007: "Mapping
;       on the HEALPix grid", MNRAS, 381, 865--872. 
; 
;       Based on hpgrid which deals with unindexed full-sky images.
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
;       Result = CUT4GRID(Pixels, Signal, Header, Order, Proj)
;
; INPUTS:
;       Pixels:   List of HEALPix pixel values
;
;       Signal:  The HEALPix dataset to grid. May be
;                 o  A 1-D array (must be a valid HEALPix length N_pix)
;                 o  An N-D array with the first dimension N_pix
;                 o  An N-D array of pointers, each pointing to a 1-D
;                    array with the same length N_pix, but not
;                    necessarily the same data type. E.g. normal
;                    CUT4 files containg columns "SIGNAL" (float),
;                    "N_OBS" (integer) and "SERROR" (float): these can
;                    all be handled at once as elements of Signal.
;
;       Header:   FITS-format header describing the data. (Modified on
;                 output). 
;
; OPTIONAL INPUTS:
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
;       The returned value is either a 2-D image with the minimal size 
;       needed to contain all the pixels, or (if the input was
;       multidimensional) an array of pointers to such images. The
;       pointer array has the same dimensions as the second and higher
;       dimension of the input. 
;
;       Images in the 'GRID' projection are automatically "rolled" in
;       longitude to avoid splitting the sky region between sections
;       separated by nearly 360 degrees. 
;
; OPTIONAL OUTPUTS:
;       Header:  If specified this returns a FITS-like header with
;                astrometry keywords set assuming the data is stored
;                in a FITS primary HDU.
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
;       Written by:     J. P. Leahy, March 2008
;
;-
COMPILE_OPT IDL2

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
    IF nextra GT 0 THEN hparray = REFORM(hparray,ndim, /OVERWRITE)
    nextra = nextra + 1
    S = SIZE(*hparray[0])
    IF S[0] NE 1 THEN MESSAGE, $
      'Pointers must each reference a single HEALPix array'
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
    IF nextra GT 0 THEN hparray = REFORM(hparray, npix, ndim, /OVERWRITE)
    type = REPLICATE(intype,ndim)
ENDELSE

ispointer = ndim GT 1 || inpointer

; Do we have a header to update?
ishead = N_ELEMENTS(header) GT 0 

lpixarr = N_ELEMENTS(pixels)

IF lpixarr NE npix THEN MESSAGE, $
  'Pixel array different length from data arrays'

IF ishead EQ 0 THEN MESSAGE, 'Pixel list supplied, but no header'
nside = SXPAR(header,'NSIDE')
IF nside EQ 0 THEN MESSAGE, $
  'Pixel list supplied but no N_side value found in header'
nside = LONG(nside)
nss = FIX(nside)

; Initial data:

fmt = "('Cut4grid: ',F7.3,' seconds: ',A)"

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

fpdata = (type GE 4 AND type LE 6) OR type EQ 9
IF MAX(fpdata) EQ 0 THEN nan = INTARR(ndim) ELSE BEGIN
    nan = FLTARR(ndim)          ; Zeros
    nan[WHERE(fpdata)] = !values.F_NAN
ENDELSE

IF butterfly THEN BEGIN
    npanel = 4
    xcorn = xcbut
    ycorn = ycbut
ENDIF ELSE BEGIN
    npanel = 5
    xcorn = xcgrd
    ycorn = ycgrd
ENDELSE
ngrid = npanel * nside
nface = nside * nside

IF ring THEN tring = set_ring_offset(nside)
llshift = 0

                                ; First find weighted centre of region
fnpix = FLOAT(npix)
skyfrac = NSIDE2NPIX(nside) / fnpix
IF ring THEN PIX2VEC_RING, nside, pixels, vec $
        ELSE PIX2VEC_NEST, nside, pixels, vec
    
centroid = TOTAL(vec,1) / fnpix
vec = 0
vec2ang, centroid, bb0, llon0, /ASTRO
PRINT, llon0, bb0, FORMAT= $
  "('Area-weighted centre of FOV is at long, lat: ',2F8.3)"
vec2pix_nest, nside, centroid, cpix
cfacet = cpix / nface
IF verbose THEN PRINT, 'Centre is in HEALPix facet', cfacet

IF ~butterfly THEN BEGIN
    llshift = ROUND(llon0[0]/90d0)
    IF llshift GT 2 THEN llshift -= 4 ; wraps into [-1:2]
    xcorn += llshift
    ycorn += llshift
    idx = WHERE(xcorn GT 4 OR ycorn GT 4)
    IF idx[0] NE -1 THEN BEGIN
        xcorn[idx] -= 4
        ycorn[idx] -= 4
    ENDIF
    idx = WHERE(xcorn LT 0 OR ycorn LT 0)
    IF idx[0] NE -1 THEN BEGIN
        xcorn[idx] += 4
        ycorn[idx] += 4
    ENDIF
    idx = 0
ENDIF 

xcorn *= nside
ycorn *= nside

                                ; Partition pixels into facets
hp2hpx, pixels, order, nside, jface, ix, iy
facets = jface[ UNIQ( jface, SORT(jface) )]
IF verbose THEN PRINT, 'Pixels occur in facets:', facets

nfacets = N_ELEMENTS(facets)
nmin1 = nss - 1S
mmpix = [ngrid, -1, ngrid, -1]
mm = [minmax(xcorn[facets]), minmax(ycorn[facets])]

FOR jj = 0,nfacets-1 DO BEGIN ; Find split facets and rotate accordingly
    face = facets[jj]
    idx = WHERE(jface EQ face)
    IF butterfly THEN BEGIN
        IF face GE 4 && face LE 7 THEN BEGIN
            in = npole ? WHERE(ix[idx] LT nss - iy[idx], COMPLEMENT = out) $
                       : WHERE(ix[idx] GE nmin1 - iy[idx], COMPLEMENT = out)
            IF out[0] GT -1 THEN BEGIN
                f8 = face + 8
                out = idx[out]
                jface[out] = f8 
                facets = [facets, f8]
                rotate_coord, ix, iy, out, nside, robut[f8]
                mm = [minmax(xcorn[facets]), minmax(ycorn[facets])]
                checklim, ix, iy, xcorn, ycorn, mm, f8, out, mmpix
            ENDIF
            idx = in[0] GT -1 ? idx[in] : -1
            in = 0  &  out = 0
        ENDIF 
        IF idx[0] GT -1 THEN rotate_coord, ix, iy, idx, nside, robut[face]
    ENDIF
    checklim, ix, iy, xcorn, ycorn, mm, face, idx, mmpix
ENDFOR       
idx = 0
IF verbose && butterfly THEN PRINT, 'including divided facets:', facets

mmpix += mm
minx = mmpix[0]  &  maxx = mmpix[1]  &  miny = mmpix[2]  &  maxy = mmpix[3]
nx   = maxx - minx + 1  &  ny   = maxy - miny + 1
blc = [minx, miny]  &  trc = [maxx, maxy]
IF verbose THEN PRINT, 'Bottom left and top right corners:', blc, trc

ix  += xcorn[jface]
ix  -= minx
hpi  = TEMPORARY(ix)
iy  += ycorn[jface]
jface = 0
iy  -= miny
iy  *= nx
hpi += TEMPORARY(iy)

IF verbose THEN PRINT, SYSTIME(1)-start, "Indices done", FORMAT=fmt

; Make header if there is something to return it to:
IF ARG_PRESENT(header) THEN BEGIN
    IF ~ishead THEN MKHDR, header, type[0], [nx, ny]
    header = grid_header(header, nside, proj, nextra, extra_dims, $
                         blc, trc, llshift, UNIT = unit)
ENDIF

; Make output arrays:
IF ispointer THEN BEGIN
    grid = PTRARR(ndim, /ALLOCATE_HEAP)
    FOR idim = 0L,ndim-1 DO BEGIN
        fptr = grid[idim]
        *fptr = MAKE_ARRAY(nx, ny, TYPE = type[idim], VALUE = nan[idim])
        IF inpointer THEN ptr = hparray[idim]
        (*fptr)[hpi] = inpointer ? *ptr : hparray[*,idim]
    ENDFOR
ENDIF ELSE BEGIN
    grid = MAKE_ARRAY(nx, ny, TYPE = type, VALUE = nan)
    grid[hpi] = hparray
ENDELSE
IF lut THEN BEGIN
    hp_index = MAKE_ARRAY(nx, ny, TYPE = 3, VALUE = -1L)
    hp_index[hpi] = pixels
ENDIF

IF verbose THEN PRINT, SYSTIME(1)-start, "Grids filled", FORMAT=fmt
IF verbose THEN HELP, /MEMORY

; Restore squashed dimensions:
hparray = inpointer ? REFORM(hparray, extra_dims, /OVERWRITE) $
                    : REFORM(hparray, S[1:S[0]], /OVERWRITE)

IF ispointer THEN grid = REFORM(grid,extra_dims,/OVERWRITE) 

IF verbose THEN PRINT, SYSTIME(1)-start, "All done", FORMAT=fmt
IF verbose THEN HELP, /MEMORY

RETURN, grid

END


