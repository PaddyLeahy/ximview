; -----------------------------------------------------------------------------
;
;  Copyright (C) 2007-2020   J. P. Leahy
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
FUNCTION flux_unit, unit, ftype
;
; Finds unit of flux density given unit of intensity.
; 
; INPUT
;     unit     FITS-compliant unit string for map intensity
;
; OUTPUT
;    
;     ftype   = 2 unit is flux per pixel: no scaling of output
;             = 1 unit is flux per beam: divide output by number of
;                 pixels per beam
;             = 0 unit is flux per sr: multiply output by pixel
;                 area in sr.
;
;     function returns FITS-compliant unit string for flux
;
; ALGORITHM
;
;    I assume unit is relatively simple, e.g. "Jy/beam" "Jy(pix)^(-1)".
;    Arbitrarily complicated unit expressions are allowed in the
;    standard but will not be interpreted correctly.
;   
;    If the unit contains an element referring to solid angle, namely
;    "pix", "pixel", "beam", or "sr", and it is in the denominator, i.e. after
;    a solidus or raised to a -1 power, then that element is
;    removed. If it is the only element in the denominator, so is the
;    solidus or -1 power. Thus "mJy / beam" becomes "mJy", but "W / (Hz sr)"
;    becomes "W / (Hz)". The ftype flag is set as appropriate.
;   
;    If no area unit is recognised, then it is assumed the unit refers
;    to surface brightness (e.g. "K"): " sr" is appended to the
;    unit, and ftype is set to 0. 
;
ftype = [0,2,1,2] ; correspond to 'sr','pix','beam','pixel' 
ulen = STRLEN(unit)

; Do we have an area unit in the unit string?
a1 = STREGEX(unit,'pixel',/FOLD_CASE, LENGTH = alen)
IF a1 LT 0 THEN a1 = STREGEX(unit,'beam',/FOLD_CASE, LENGTH = alen) 
IF a1 LT 0 THEN a1 = STREGEX(unit,'pix',/FOLD_CASE, LENGTH = alen) 
IF a1 LT 0 THEN a1 = STREGEX(unit,'sr', LENGTH = alen)
start = a1

IF start GE 0 THEN BEGIN  ; we have an area unit. 
   aunit = STRMID(unit,a1,alen)
   a2 = a1 + alen - 1
   next  = a2 + 1

   bra = STRPOS(unit,'(',start,/REVERSE_SEARCH)     ; Is it in parentheses?
   bracketed = 0B
   IF bra GE 0 THEN BEGIN
      ket = STRPOS(unit,')',next)
      bracketed = ket GT 0
   ENDIF
   IF bracketed THEN BEGIN
      start = bra
      next = ket + 1
   ENDIF
   blen = next - start
;   PRINT, 'bracket is: '+ STRMID(unit,start,blen)

                                ; Find where denominator starts
   solidus = STRPOS(unit,'/',start, /REVERSE_SEARCH)
   per = solidus GE 0 AND solidus LT start
   IF per THEN BEGIN 
      numerator = STRMID(unit,0,solidus)
      denominator = STRMID(unit,solidus+1,ulen-solidus)
;      PRINT, 'Denominator: "'+ denominator+'"'
   ENDIF ELSE BEGIN ; look for a -1 power
      tail = STRMID(unit,next,ulen-next+1)
      test = STREGEX(tail,'-1|\(-1\)|^-1|^\(-1\)|\*\*-1|\*\*\(-1\)',LENGTH=elen)
      per = test EQ 0
      IF per THEN BEGIN
         numerator = STRMID(unit,0,start)
         denominator = bracketed ? STRMID(unit,bra+1,blen-2) : $
                       STRMID(unit,start,blen)
;         PRINT, 'Denominator: "'+denominator+'"'
      ENDIF
   ENDELSE
   IF per THEN BEGIN
      IF STRTRIM(denominator,2) NE aunit THEN BEGIN ; we have other stuff below.
         ubits = STRSPLIT(unit,' *'+aunit, /EXTRACT, /REGEX)
         fluxunit = STRJOIN(ubits) ; area unit removed from denominator
      ENDIF ELSE fluxunit = STRTRIM(numerator)
   ENDIF ELSE start = -1
ENDIF

IF start LT 0 THEN BEGIN ; Default is intensity unit * sr.
   fluxunit = unit  + ' sr'
   alen = 2
ENDIF

ftype = ftype[alen-2]

RETURN, fluxunit
END
;
PRO box2xy, image, boxsize, xpix, ypix, T, $
            coords, imval, x, y, line, nvalid, ierr
;
; Ancillary routine to return lists of x, y coords given input box
;
; INPUTS
;     image      pointer to 2D image
;     boxsize    size of box in (x,y). If scalar returned as 2-element
;                array
;     xpix, ypix coords of box centre
;     T          SIZE array for image
; 
; OUTPUTS
;     coords     coordinates of box corner
;     imval      Image values in box
;     x, y       arrays of x, y coords for pixels in box
;     line       text string describing box for IMSTATS output
;     nvalid     number of pixels in box
;     ierr       0 for no problem, 1 for box outside image
;
ierr = 0B
B = SIZE(boxsize)
IF B[1] EQ 1 THEN boxsize = REPLICATE(boxsize,2)

xypix = [xpix, ypix]
half = boxsize / 2
xy0 = ((xypix - half) > 0) < (T[1:2] - 1)
xy1 = ((xypix + half) > 0) < (T[1:2] - 1)
dxy = xy1 - xy0 + 1
nvalid = dxy[0]*dxy[1]

coords = [[xy0], [xy0[0], xy1[1]], [xy1], [xy1[0], xy0[1]]] ; Returned

IF nvalid EQ 0 THEN BEGIN
   MESSAGE, /INFORMATIONAL, 'Box outside image, no data to analyse'
   ierr = 1B
   RETURN
ENDIF
imval = (*image)[xy0[0]:xy1[0],xy0[1]:xy1[1]]
    
x = LINDGEN(dxy[0]) + xy0[0]  &  y = LINDGEN(dxy[1]) + xy0[1]
x = REBIN(x,dxy[0],dxy[1])  & y = TRANSPOSE(REBIN(y,dxy[1],dxy[0]))

line = STRING(xy0[0], xy1[0], xy0[1], xy1[1], FORMAT = $
"('IMSTATS: Image statistics for ',I5,' <= x <=',I5,', ',I5,' <= y <=',I5)")
END
;
FUNCTION imstats, image, xpix, ypix, boxsize, ASTROM = astrom, $
                  BADZERO = badzero, LUN = lun, UNIT = unit, BEAMA = beam, $
                  MULT = mult, THRESHOLD = threshold, NT = ncode 
;+
; NAME:
;       IMSTATS
;
; PURPOSE:
;       Prints basic statistics on either (1) an image box surrounding
;       specified position or (2) a region defined by a coordinate
;       list. Blanks (NaNs) are ignored. Optionally, exactly zero
;       pixels are treated as blanks.
;
;       Results are:
;        * Position and values of minimum and maximum pixels (World
;          Coordinates as well as pixel coordinates if ASTROM is set)  
;        * Mean, median, and standard deviation of values
;        * Integrated flux density (if there is enough astrometric
;          information to derive a pixel size, or if units are flux/pixel.)
;        * Number of pixels above a given threshold brightness, and
;          the largest distance between pixels within that set (NB:
;          not necessarily a contiguous island).
;          If astrometry is available, also gives solid angle and LAS
;
; CATEGORY:
;       Image processing, Statistics
;
; CALLING SEQUENCE:
;   
;       Result = IMSTATS(Image, Xpix, Ypix, Boxsize)
;
; INPUTS:
;     Image:  pointer to 2-D image from which statistics are extracted
;
;     Xpix:   EITHER X coordinate of box centre
;             OR     list of X coordinates for all pixels in region
;
;     Ypix:    Y coordinates as above
;
; OPTIONAL INPUTS:
;     Boxsize: 2-element array specifying size of box, or one value if
;              box is square. If specified, Xpix & Ypix must be scalar.
;
; KEYWORD PARAMETERS:
;     UNIT:    intensity unit for output. If not equal to the unit of
;              the image, the factor MULT must be specified to
;              compensate. IMSTATS cannot check that this is done correctly.
;
;     MULT:    factor by which image must be multiplied for
;              consistency with unit. 
;
;     BEAMA:   beam area in square degrees
;
;     BADZERO  Set to treat zero as blank (Strongly deprecated!)
;
;     LUN:     logical unit number for log file
;
;     ASTROM:  astrolib astrometry structure (if any, otherwise 0)
;
;     THRESHOLD: list number of pixels above threshold value 
;
;     NT:      Noise type. 0 = unknown, 1 = independent in each pixel,
;              2 = correlated over beam.
;
; OUTPUTS:
;     Returns structure containing:
;        .coord   coordinates of statistics box BLC and TRC, after
;                 adjustment for edges of image (for pixel list, coords
;                 give minimal box covering ROI).
;        .unit    Unit for intensity values, as specified on input
;        .mean    mean of valid pixels in box
;        .median  median ditto
;        .stddev  standard deviation, ditto
;        .npix    number of valid pixels in box
;
;        .immax     Value of brightest pixel in box
;        .maxpix    = [xmax,ymax] coords of brightest pixel in box
;        .maxpos    Astrometric coordinates of maxpix
;        .immin     Value of faintest pixel in box
;        .minpix    = [xmin,ymin] ditto, faintest pixel
;        .minpos    Astrometric coordinates of minpix
;        .flux      Flux density in box
;        .fluxunit  unit for flux density, derived from unit
;        .flux_b    box flux density, corrected using background mean.
;        .e_flux    estimated error in flux based on background std. dev.
;
;        .threshold Flux density threshold for measuring solid angle
;        .nsource   Number of pixels above flux density threshold.
;        .omega     Solid angle in square degrees equivalent to nsource
;        .size_pix  Maximum distance between pixels above threshold,
;                   in pixels
;        .size_deg  Maximum distance between pixels above threshold in
;                   degrees
;
;        .e1        ellipse structure specifying inner radius of
;                   background annulus
;        .e2        ditto, for outer radius
;        .ebox      corner coordinates of box around background annulus.
;        .b_med     Median of pixels in background annulus
;        .b_mean    Mean, ditto (after 3-sigma outlier clip)
;        .b_stdd    standard deviation, ditto
;        .b_npix    number of valid pixels in background annulus
;                   (excludes >3 sigma outliers as well as blanks)
;
; SIDE EFFECTS:
;     Prints results to screen and to log file if LUN is set.
;
; EXAMPLE:
;
;          test = PTR_NEW(RANDOMN(seed,50,50))
;          boxcoord = IMSTATS(test, 30, 30, 33)
;
;  Prints:  
;
;    IMSTATS: Image statistics for    14 <= x <=   46,    14 <= y <=   46
;             Includes   1089 sky pixels of which  1089 contain valid data
; 
;             Maximum:   3.396E+00 at pixel (   14,   14)
;             Minimum:  -3.212E+00 at pixel (   44,   23)
;
;             Mean:          -1.669E-02  Median: -2.421E-03
;             Standard Dev:   9.850E-01  Unit:   unknown
;     
; MODIFICATION HISTORY:
;       Written by:      J. P. Leahy, 2008
;       Aug 2013:        set coords to -1 for tiny image.
;       Aug 2016:        added threshold parameter
;       May 2020:        Escaped parentheses in STREGEX strings (used
;                        for fixing flux unit as per pixel).
;       Sep 2020:        Added beam area to get flux scaling, and
;                        elliptical background annulus. Image now
;                        passed as pointer, results returned as
;                        structure; removed duplicate output if exact
;                        zeros but added BADZERO keyword. Sorted out
;                        flux density unit, returns LAS.

;-
COMPILE_OPT IDL2
ON_ERROR, 0

r2d = 180d0 / !dpi
line = ['']
T = SIZE(*image)
IF T[1] LT 3 OR T[2] LT 3 THEN BEGIN
    MESSAGE, /INFORMATIONAL, 'Image too small to analyse!'
    coords = -1
    GOTO, QUIT
ENDIF

; check that version is high enough to run SmEll:
order = SORT(['8.3',!version.release])
smell_ok =  order[0] EQ 0 

badzero = KEYWORD_SET(badzero)
IF N_ELEMENTS(ncode) EQ 0 THEN ncode = 0
IF N_ELEMENTS(beam) EQ 0 THEN beam = !values.F_NAN

llbb = KEYWORD_SET(astrom) ; we can work out longs & lats
IF llbb THEN BEGIN
    radec = astrom.coord_sys EQ 'C'
    cdelt = astrom.CDELT
    IF astrom.REVERSE THEN cdelt = REVERSE(cdelt)  

; Find appropriate precision for pixel positions
    IF radec THEN prec = CEIL(-ALOG10(cdelt[1]*3600.0)) ELSE BEGIN 
        prec = CEIL(-ALOG10(ABS(cdelt))) + 1
        nchar = prec + 5
        nchar = STRTRIM(STRING(nchar),2)
        ndp   = STRTRIM(STRING(prec),2)
        ll_fmt = "(' long:',F"+STRJOIN(nchar+'.'+ndp,",'  lat:',F+")+')'
    ENDELSE
ENDIF

no_unit = N_ELEMENTS(unit) EQ 0
IF no_unit THEN unit = '' ELSE unit = ' '+unit

IF N_ELEMENTS(mult) EQ 0 THEN mult = 1

XS = SIZE(xpix)  & YS = SIZE(ypix)
IF ~ARRAY_EQUAL(XS, YS) THEN MESSAGE, 'Mismatched X and Y pixel arrays'
npix = N_ELEMENTS(xpix)

is_box = N_PARAMS() EQ 4

IF is_box THEN BEGIN            ; Find pixel coords in box
   IF npix GT 1 THEN $
      MESSAGE, 'Coordinate list and box size should not both be supplied'

   box2xy, image, boxsize, xpix, ypix, T, coords, box, x, y, line1, nvalid, ierr
   IF ierr THEN GOTO, quit

   line = [' ',line1]
ENDIF ELSE BEGIN                ; Pixel list supplied
   impix = T[1]*T[2]
   idx = xpix + ypix*T[1]
   good = WHERE(idx GE 0 AND idx LT impix, nvalid)
   idx = idx[good]
                                ; Find bounding box of region. NB: may
                                ; go outside image. 
   x01 = MINMAX(xpix)
   y01 = MINMAX(ypix)

   coords = [[x01[0], y01[0]], [x01[0], y01[1]], $
             [x01[1], y01[1]], [x01[1], y01[0]]]

   IF nvalid EQ 0 THEN BEGIN
      MESSAGE, /INFORMATIONAL, 'Region outside image, no data to analyse'
      GOTO, QUIT
   ENDIF

   IF nvalid LT npix THEN BEGIN
      x = xpix[good]  &  y = ypix[good]
   ENDIF ELSE BEGIN
      x = xpix  &  y = ypix
   ENDELSE
   good = 0

   box = (*image)[idx]
   line = [' ',STRING( npix, FORMAT = $
     "('IMSTATS: Image statistics for region with ',I8,' pixels')")]
ENDELSE

; Scale box values as requested in inputs
box *= mult

; Get astronomical coordinates of each point in box, in degrees
IF llbb THEN XY2AD, x, y, astrom, ll, bb

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Exclude invalid pixels from box
;
; check that pixels are on sky and not in a gore
valid = llbb ? WHERE(FINITE(ll) AND FINITE(bb), nvalid) : LINDGEN(nvalid)

; check that pixels are not masked:
IF nvalid GT 0 THEN BEGIN
   good = WHERE(FINITE(box[valid]), ngood) 
   valid = valid[good]
ENDIF ELSE ngood = 0

; Are zeros used for null pixels?
IF badzero THEN BEGIN
   good = WHERE(box[valid] NE 0.0, ngood)
   valid = valid[good]
ENDIF ELSE BEGIN ; check to see if /badzero should have been set
   nnan = nvalid - ngood
   IF nnan EQ 0 THEN BEGIN
      zeroes = WHERE(box EQ 0.0, nzero)
      IF FLOAT(nzero)/nvalid GT 1e-6 THEN BEGIN
         MESSAGE, /INFORMATIONAL, STRING( nzero, FORMAT = $
                                          "('Found',I6,' exact zeroes.')")
         MESSAGE, /INFORMATIONAL,  $
  'If these represent blank pixels replace with NaNs or re-run with /BADZERO'
      ENDIF
   ENDIF
ENDELSE

line = [line, STRING(nvalid, ngood, FORMAT = $
"( '         Includes ',I6,' sky pixels of which',I6,' contain valid data')")]


IF ngood EQ 0 THEN BEGIN
    line = [line, 'IMSTATS: No further analysis possible!']
    GOTO, QUIT
ENDIF

; Find pixel and beam solid angles (in sr) if we have astrometry
IF llbb THEN BEGIN
   beamarea = beam*!dtor^2 ; beam is in square degrees
   wcs = STRMID(astrom.CTYPE[0],5)
   hpx = wcs EQ 'HPX' || wcs EQ 'XPH'
   IF hpx THEN BEGIN
      pixarea = !dpi / (3L*nside^2)
      pixadif = 0d0
   ENDIF ELSE BEGIN
      area = DBLARR(4)
      FOR ipix = 0,3 DO area[ipix] = get_pix_area(coords[*,ipix], astrom)
      pixarea = MEAN(area)
      pixadif = MAX(area) - MIN(area)
   ENDELSE
   pixperbeam = beamarea/pixarea
ENDIF

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Find pixels above threshold and work out angular size
;
size_pix = !values.F_NAN        ; defaults
size_deg = !values.F_NAN
nsource = 0
omega = !values.F_NAN
IF FINITE(threshold) THEN BEGIN
   ids = WHERE(box GE threshold, nsource)
   omega = nsource * pixarea*!radeg^2 ; convert to square degrees
   aunit = 'degrees'
   omegap = omega ; scale for printout
   IF omega LT 0.1 THEN BEGIN
      omegap *= 3600
      aunit = 'arcmin'
      IF omegap LT 0.1 THEN BEGIN
         omegap *= 3600
         aunit = 'arcsec'
      ENDIF
   ENDIF
   line = [line, STRING(nsource, threshold, unit, FORMAT = $
                        "(9X,I6,' pixels above threshold of',G10.4,' ',A)"), $
           STRING(omegap, aunit, FORMAT="(9X,'=',G10.4,' square ',A)")]
   IF nsource GT 0 THEN BEGIN
      xs = x[ids]
      ys = y[ids]
      ; eliminate internal points:
      QHULL, xs, ys, tr
      id = REFORM(tr[0,*])
      nh = N_ELEMENTS(id)
      xh = xs[id]
      yh = ys[id]
      size_pix = 0.0
      FOR jd = 0,nh-2 DO FOR kd = jd+1,nh-1 DO $
         size_pix = ((xh[jd]-xh[kd])^2 + (yh[jd]-yh[kd])^2) > size_pix
      size_pix = SQRT(size_pix)
      size_deg = !values.F_NAN
      line = [line, ' ', STRING(size_pix, FORMAT = $
                   "(9X,'Largest size of thresholded patch:',G12.4,' pix')")]
      IF llbb THEN BEGIN ; Redo in spherical geometry. May give different hull:
         lls = ll[ids]
         bbs = bb[ids]
         QHULL, lls, bbs, tr, BOUNDS=id, SPHERE=sphere 
         nh = N_ELEMENTS(id)
         llh = lls[id]
         bbh = bbs[id]
         size_deg = 0.0
         FOR jd = 0,nh-2 DO FOR kd = jd+1,nh-1 DO BEGIN
            GCIRC, 2, llh[jd], bbh[jd], llh[kd], bbh[kd], dis
            size_deg = (dis/3600.0) > size_deg
         ENDFOR
         line = [line, STRING(size_deg*60, FORMAT = $
                              "(9X,'Largest angular size:',G12.4,' arcmin')")]
      ENDIF
   ENDIF
ENDIF

box = box[valid]
x   = x[valid]
y   = y[valid]
IF llbb THEN BEGIN
   ll  = ll[valid]
   bb  = bb[valid]
ENDIF
line = [line, ' ']
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Find largest and smallest pixel value in box
boxmax = MAX(box)
maxpix = WHERE(box EQ boxmax, nmax)
IF nmax GT 1 THEN BEGIN
    line = [line, 'IMSTATS: Multiple maxima! Quoting the first found']
    maxpix = maxpix[0]
ENDIF
boxmin = MIN(box)
minpix = WHERE(box EQ boxmin, nmin)
IF nmin GT 1 THEN BEGIN
    line = [line, 'IMSTATS: Multiple minima! Quoting the first found']
    minpix = minpix[0]
ENDIF

maxcoor = [x[maxpix], y[maxpix]]
mincoor = [x[minpix], y[minpix]]

normfm = "(9(' '),A,':',G12.4,A,' at pixel (',I5,',',I5,')')"
spacer = STRJOIN(REPLICATE(' ',40))
IF llbb THEN BEGIN
    lmax = ll[maxpix] & bmax = bb[maxpix]
    lmin = ll[minpix] & bmin = bb[minpix]
    linem = STRING('Maximum', boxmax, unit, maxcoor, FORMAT = normfm)
    linen = radec ? '  RA/Dec:' + ADSTRING(lmax, bmax, prec) $
            : STRING(lmax, bmax, FORMAT=ll_fmt)
    line = [line, linem, spacer+linen]

    IF hpx THEN BEGIN
        nside = wcs EQ 'HPX' ? T[1] / 5 : T[1] / 4
        theta = (90d0 - bmax)/r2d  & phi = lmax/r2d
        ang2pix_ring, nside, theta, phi, maxpixr
        ang2pix_nest, nside, theta, phi, maxpixn
        theta = (90d0 - bmin)/r2d  & phi = lmin/r2d
        ang2pix_ring, nside, theta, phi, minpixr
        ang2pix_nest, nside, theta, phi, minpixn

        hpfmh  = "(13(' '),'at HEALPix pixels: ',I9,' (nest); '" + $
                                               ",I9,' (ring)')"

        line = [line, STRING(maxpixn, maxpixr, FORMAT = hpfmh)]
    ENDIF
    linem = STRING('Minimum', boxmin, unit, mincoor, FORMAT = normfm)
    linen = radec? '  RA/Dec:' + ADSTRING(lmin, bmin, prec) $
                   : STRING(lmin, bmin, FORMAT=ll_fmt)
    line = [line, '', linem, spacer+linen]

    IF hpx THEN line = [line, STRING( minpixn, minpixr, FORMAT = hpfmh)]

ENDIF ELSE BEGIN
    line = [line, STRING( 'Maximum', boxmax, unit, maxcoor, FORMAT=normfm)]
    line = [line, STRING( 'Minimum', boxmin, unit, mincoor, FORMAT=normfm)]
    lmax = !values.F_NAN & lmin = !values.F_NAN
    bmax = !values.F_NAN & bmin = !values.F_NAN
ENDELSE

boxmean = MEAN(box)
boxmed  = MEDIAN(box)
boxstd  = STDDEV(box)
boxmad  = MEDIAN(ABS(box-boxmed))

line = [line, ' ', STRING( boxmean, unit, boxmed, unit, FORMAT = $
                "(9(' '),'Mean:        ',G12.4,A,'  Median:',G11.4,A)")]
line = [line, STRING( boxstd, unit, boxmad, unit, FORMAT = $
                "(9(' '),'Standard Dev:',G12.4,A,'  MAD:   ',G11.4,A)")]

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Photometry section
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

fluxunit = no_unit ? flux_unit('unknown', ftype) : flux_unit(unit, ftype)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Backgrounds
;   For rectangular box, we probably don't want a background
;   but if we do it should be rectangular.
;
;   For pixel lists, use elliptical background annulus just outside ROI.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
IF is_box OR ~smell_ok THEN BEGIN 
   e1 = {sup: 0} ; null ellipse
   e2 = e1
   bbcoord = LONARR(2,4)
   npix_bck = 0                 ; no background
   mean_bck = !values.F_NAN
   median_bck = !values.F_NAN
   stddev_bck = !values.F_NAN
ENDIF ELSE BEGIN                ; Define elliptical background annulus

                                ; smallest ellipse around selected
                                ; points
   e1 = smell([[x],[y]])
   param = centrep2param(e1)    ; parameters of ellipse
   axrat = param[3]/param[2]    ; axial ratio
   tau = ATAN(param[5],param[4])*!radeg ; Position angle of ellipse

; Background annulus should have same solid angle as photometry
; aperture, where we define the photometry aperture as the area
; used for integration, i.e ngood pixels, and NOT the area above the
; threshould, i.e. nsource pixels

   nellipse = !dpi/SQRT(DETERM(e1.M/e1.z)) ; ellipse area in pixels
   factor = ngood/nellipse + 1d0 
   e2 = e1
   e2.z *= factor
                                ; find pixels in background ellipse
; first find pixels in box surrounding outer ellipse:
   ebox = ellipse_box(e2)
   ecentre = ROUND(e1.c)        ; nearest pixel to ellipse centre
   dbx = MAX(ABS(ebox[0,*]-ecentre[0])) 
   dby = MAX(ABS(ebox[1,*]-ecentre[1])) 
   bboxsize = 2*FIX([dbx,dby]) + 3 ; make sure box fully encloses ellipse

   box2xy, image, bboxsize, ecentre[0], ecentre[1], T, $
           bbcoord, bckval, x_bck, y_bck

; Scale background values as requested in inputs
   bckval *= mult 
; Now check for inclusion between inner and outer ellipses:
   r = e1.m[0,0]
   s = e1.m[1,1]
   t = e1.m[0,1]
   x_bc = x_bck - e1.c[0]       ; centred coords
   y_bc = y_bck - e1.c[1] 
   conic = r*x_bc^2 + s*y_bc^2 + (2*t)*x_bc*y_bc

   ibck = WHERE(conic GT e1.z AND conic LE e2.z, nbck)
   IF nbck GT 0 THEN bckval = bckval[ibck] ELSE bckval = !values.F_NAN

   bg_good = WHERE(FINITE(bckval), nbg)
   line1 = [STRING('Using elliptical background annulus:', FORMAT="(9X,A)"), $
            STRING(param[0:1], FORMAT = $
                   "(12X,'centre:',F7.1,',',F7.1,' pixels')"), $
            STRING(param[2]*[1,factor], FORMAT = $
        "(12X,'inner and outer semi-major axis:',F6.1,',',F6.1,' pixels')"), $
            STRING(axrat, tau, FORMAT = $ 
        "(12X,'axial ratio:',F6.3,' at',F6.1,' degrees CCW from x-axis')")]
   IF nbg GE 7 THEN BEGIN
      median_bck = MEDIAN(bckval)
      sample = bckval[bg_good]
      nin_old = nbg
      REPEAT BEGIN
         null = MOMENT(sample, MAXMOMENT=2, MEAN = bmean, SDEV = bsdev)
; get robust mean & standard deviation of background by eliminating outliers
         in = WHERE(ABS(sample) LE bmean + 3*bsdev, nin)
         nout = nin_old - nin
         nin_old = nin
         IF nin GE 7 THEN sample = sample[in] ELSE BREAK
      ENDREP UNTIL nout EQ 0
      IF nin GE 7 THEN BEGIN
         mean_bck   = bmean
         stddev_bck = bsdev 
      ENDIF
      npix_bck   = nin
         
      line2 = STRING(mean_bck, unit, median_bck, unit, FORMAT = $
                "(12X,'Background mean:',G12.4,A,'  Median:',G11.4,A)")
      line3 =STRING(stddev_bck, unit, npix_bck, FORMAT = $
                "(12X,'   Standard Dev:',G12.4,A,' over',I8,' pixels')")
      line = [line, ' ', line1, line2, line3]
      if llbb THEN line = [line, STRING(npix_bck/pixperbeam, FORMAT = $
                 "(50X,'(',G11.4,' beam areas)')")]
   ENDIF ELSE line = [line, $
      '          Too few valid points in annulus for background estimate']
ENDELSE 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; Calcuated integrated flux density if possible
flux = !values.F_NAN            ; Default to NaN
flux_bsub = flux
e_flux = flux
IF llbb || ftype EQ 2 THEN BEGIN
   CASE ftype OF
      0: factor = pixarea                     ; Flux per sr
      1: factor = 1/pixperbeam                ; Flux per beam
      2: factor = 1.0                         ; Flux per pixel
      ELSE: MESSAGE, 'Unrecognised flux scaling type'
   ENDCASE
   flux = boxmean * ngood * factor
   fluxstr = numunit(flux, fluxunit, PRECISION = 4)
   line = [line, STRING( fluxstr, FORMAT = $
                         "(9X,'Integrated flux density: ',A)")]

   IF npix_bck GT 7 THEN BEGIN
      flux_bsub = flux - ngood*factor*mean_bck

; Very approximate error estimate: caveat emptor!
      bgxs = 1 + ngood/npix_bck
      CASE ncode OF 
         0: BEGIN
            eflux = !values.F_NAN
            line = [line,'Set noise correlation type to get error estimate', $
                    'using Analysis -> Set image properties']
         END
         1: e_flux = stddev_bck * SQRT(ngood * bgxs)
         2: e_flux = stddev_bck * SQRT(ngood * bgxs/pixperbeam)
         ELSE: MESSAGE, 'unknown noise type'
      ENDCASE
      fluxstr = FINITE(e_flux) ? numunit(flux_bsub, fluxunit, e_flux) $
                : numunit(flux_bsub, fluxunit, PRECISION = 4)
      line1 = STRING(fluxstr, FORMAT = "(12X,'background-corrected: ',A)") 
      line = [line,line1]
   ENDIF

; Report variation of pixel area if it might be significant:
   fracerr = pixadif / pixarea
   minerr = FINITE(e_flux) ? 0.5*e_flux/flux_bsub : 1e-4
   IF fracerr GT minerr THEN line = [line, STRING( 100.*fracerr, FORMAT = $
           "(9X,'Pixel area changes by',F6.2,'% across region;')"), $
                                   '         representative value used.']
ENDIF

QUIT:

IF line[0] NE '' THEN BEGIN
   line = [line,' ']
   PRINT, line, FORMAT = "(A)"  ; Should put each string on its own line
   IF N_ELEMENTS(lun) NE 0 THEN PRINTF, lun, line, FORMAT="(A)"
ENDIF

RETURN, {box: coords, unit: STRTRIM(unit,2), $
         immax: boxmax, maxpix: maxcoor, maxpos: [lmax, bmax], $
         immin: boxmin, minpix: mincoor, minpos: [lmin, bmin], $
         mean: boxmean, stddev: boxstd, median: boxmed, npix: ngood, $
         flux: flux, fluxunit: STRTRIM(fluxunit,2), flux_bsub: flux_bsub, $
         e_flux: e_flux, $
         threshold: threshold, nsource: nsource, omega: omega, $         
         size_pix: size_pix, size_deg: size_deg, $
         e1: e1, e2: e2, ebox: bbcoord, b_med: median_bck, b_mean: mean_bck, $
         b_stdd: stddev_bck, b_npix: npix_bck} 

END

