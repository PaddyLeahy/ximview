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
FUNCTION maxfit, image, xpix, ypix, boxin, ASTROM = astrom, UNIT = unit, $
                 MULT = mult, LUN = lun, EXTREMUM=extra
;+
; NAME:
;       MAXFIT
;
; PURPOSE:
;       Prints maximum or minimum value and position in a section of a
;       2-D image, interpolated to sub-pixel accuracy.
;       Blanks (NaNs) are ignored.
;
;       Results are:
;        * Position and value of minimum or maximum pixel (World
;          coordinates as well as pixel coordinates if ASTROM is set)  
;        * Position and value of fitted peak.
;
; CATEGORY:
;       Image processing, Mathematics - Surface fitting
;
; CALLING SEQUENCE:
;   
;       Result = MAXFIT(Image, Xpix, Ypix, Boxsize)
;
; INPUTS:
;     Image:    2-D image from which statistics are extracted
;
;     Xpix:     X coordinate of box centre
;
;     Ypix:     Y coordinate of box centre
;
; OPTIONAL INPUTS:
;     Boxsize:  Size of (square) search box in pixels.
;
; KEYWORD PARAMETERS:
;     UNIT:     intensity unit of image
;     MULT:     factor by which image values must be multiplied 
;               for consistency with UNIT
;
;     LUN:      logical unit number for log file
;
;     ASTROM:   astrolib astrometry structure (if any, otherwise 0)
;
;     EXTREMUM: = 0 -> Search for maximum
;               = 1 -> Search for largest extremum, independent of
;                     sign
;               = 2 -> Search for minimum.
; OUTPUTS:
;     Returns structure with pixel & astrometric coordinates of peak,
;     and peak value
;
; SIDE EFFECTS:
;     Prints results to screen and to log file if LUN is set.
;
; PROCEDURE:
;     Fits a 2-D quadratic to the 3x3 pixels surrounding the
;     maximum or minimum pixel found in the search box. If an
;     astrometry array is present, the fit is done in world
;     coordinates, which are not assumed to be on a regular grid.
;
;     A warning is issued if the fitted peak falls outside the 3x3
;     grid. If the fit finds a saddle point instead of a peak the
;     returned value is the maximum pixel. If it finds a
;     flat image (all pixels equal), the returned value is the input
;     search box centre.
;     
;     To be done: 
;        Fractional image pixel for HEALPix maps
;
; EXAMPLE:
;
;          test = RANDOMN(seed,50,50)
;          coords = MAXFIT(test, 30, 40, 7, unit = 'K')
;
;  Prints:  
;
;    MAXFIT: Max pixel in search box: 991.5 mK at pixel   27   39
;            Fitted peak 806.6 mK at pixel coord   26.37   39.32
;     
; MODIFICATION HISTORY:
;       Written by:      J. P. Leahy, 2008
;       Updated 2013: output in HMS/DMS for RA/Dec, minor bug fixes.
;       Updated 2020: boxsize default, formatting, return structure
;-
COMPILE_OPT IDL2
ON_ERROR, 0

r2d = 180d0 / !dpi
line = ''
llbb = KEYWORD_SET(astrom)      ; if set, we can find longs and lats 

IF N_ELEMENTS(unit) EQ 0 THEN unit = 'unknown'

IF N_ELEMENTS(mult) EQ 0 THEN mult = 1

IF N_ELEMENTS(extra) EQ 0 THEN extra = 0

coord = [xpix, ypix]    ; Will return these input coords if fit fails.

IF N_ELEMENTS(boxin) EQ 0 THEN boxsize = 7 ELSE boxsize = boxin

T = SIZE(image)
xsize = T[1] & ysize = T[2]
IF xsize LT 3 || ysize LT 3 THEN BEGIN
    MESSAGE, /INFORMATIONAL, 'Image too small to fit!'
    GOTO, QUIT
ENDIF

half = boxsize / 2
x0 = ((xpix - half )> 0) < (T[1]-1)
y0 = ((ypix - half )> 0) < (T[2]-1)
x1 = ((xpix + half )> 0) < (T[1]-1)
y1 = ((ypix + half )> 0) < (T[2]-1)
npix = (x1-x0+1)*(y1-y0+1)

IF npix LE 6 THEN BEGIN
    MESSAGE, /INFORMATIONAL, $
      'Box too small for fitting. Increase box size or avoid edge of image'
    GOTO, QUIT
ENDIF
box = mult*image[x0:x1,y0:y1]

CASE extra OF
    0: BEGIN
        boxmax = MAX(box, /NAN)
        maxpix = WHERE(box EQ boxmax, count)
        peakstr = 'Max'
    END
    1: BEGIN
        boxmax = MAX(ABS(box), /NAN)
        maxpix = WHERE(ABS(box) EQ boxmax, count)
        peakstr = 'Peak'
    END
    2: BEGIN
        boxmax = MIN(box, /NAN)
        maxpix = WHERE(box EQ boxmax, count)
        peakstr = 'Min'
    END
ENDCASE

CASE count OF 
    0: BEGIN
        MESSAGE, /INFORMATIONAL, 'No valid pixels found to fit!'
        GOTO, QUIT
    END
    1: ; OK
    ELSE: BEGIN
        cstr = STRTRIM(STRING(count),2)
        MESSAGE, /INFORMATIONAL, 'WARNING: ' + cstr+ $
          ' pixels with brightness equal to maximum in box'
        MESSAGE, /INFORMATIONAL, 'Fitting around first one found'
        maxpix = maxpix[0]
    END
ENDCASE
maxcoor = ARRAY_INDICES(box,maxpix) + [x0,y0]

; fancy formatting:
digits = FIX(ALOG10(MAX([xsize,ysize]))) + 1
intf = 'I'+STRTRIM(STRING(digits),2)
flf  = 'F'+STRTRIM(STRING(digits+3),2)+'.2'

fmt = "('MAXFIT: ',A,' pixel in search box: ',A,' at pixel '," + $
      intf+",', ',"+intf+")"
flux = numunit(boxmax, unit, PRECISION = 4)
IF STRLEN(flux) GT 27 THEN flux = STRMID(flux,0,27)
line = [' ', STRING(peakstr, flux, maxcoor, FORMAT = fmt)]

; Re-centre 3x3 box on peak:
x0 = (maxcoor[0] - 1) > 0
y0 = (maxcoor[1] - 1) > 0
x1 = x0 + 2
y1 = y0 + 2
IF x1 GT T[1] - 1 THEN BEGIN
    x1 = T[1] - 1
    x0 = T[1] - 4
ENDIF
IF y1 GT T[2] - 1 THEN BEGIN
    y1 = T[2] - 1
    y0 = T[2] - 4
ENDIF
box = mult*image[x0:x1,y0:y1]
npix = 9
nvalid = npix

IF llbb THEN BEGIN 
    x = REBIN(LINDGEN(3),3,3)
    y = TRANSPOSE(x) + y0
    x = x + x0
    XY2AD, x, y, astrom, ll, bb
    valid = WHERE(FINITE(ll) AND FINITE(bb), nvalid)
ENDIF ELSE valid = LINDGEN(nvalid)

IF nvalid GT 0 THEN good = WHERE(FINITE(box[valid]), ngood) ELSE ngood = 0

temps = box[valid]
temps = temps[good]
IF llbb THEN BEGIN
    ll = ll[valid]  &  bb = bb[valid]
    ll = ll[good]   &  bb = bb[good]
                                ; Check for wrap at 360deg
    IF MAX(ll) - MIN(ll) GT 180d THEN ll[WHERE(ll GT 180d)] -= 360d
                                ; Subtract means to avoid round-off errors  
    llbar = MEAN(ll)  & bbbar = MEAN(bb) 
    data = TRANSPOSE([[ll-llbar], [bb-bbbar], [temps]])

; Parabolic fit on actual pixel sky coords:
    model = SFIT(data,2,/IRREGULAR,KX=A,/MAX_DEGREE)
ENDIF ELSE BEGIN
    IF ngood EQ npix THEN $
      model = SFIT(box, 2, KX=A, /MAX_DEGREE) $ ; fit on pixel coords
    ELSE BEGIN
        x = INDGEN(3) # REPLICATE(1,3)
        y = TRANSPOSE(x)
        data = TRANSPOSE([[x[good]],[y[good]],[temps]])
        model = SFIT(data,2,/IRREGULAR,KX=A,/MAX_DEGREE)
    ENDELSE
ENDELSE

; Find peak coords (x2,y2) and value, by analysis:
; Check for degeneracies and saddlepoints
kgauss = 4*A[2]*A[5] - A[4]^2
IF kgauss GT 0. THEN BEGIN
    kx = A[4]/(2.*A[2])  &  ky = A[4]/(2.*A[5])
    x2 = (kx*A[1] - A[3])/(2.*A[5] - kx*A[4])
    y2 = (ky*A[3] - A[1])/(2.*A[2] - ky*A[4])
    peak = A[0] + y2*(A[1] +y2*A[2]) + x2*(A[3] + y2*A[4] + x2*A[5])
ENDIF ELSE IF kgauss LT 0. THEN BEGIN
    line = [line, '        Fit gives saddle point not peak.']
    coord = maxcoor
    GOTO, QUIT
ENDIF ELSE BEGIN                ; Zero curvature
    line = [line, '        Fit degenerate.']
    GOTO, QUIT
ENDELSE

flux = numunit(peak, unit, PRECISION = 4)
IF STRLEN(flux) GT 34 THEN flux = STRMID(flux,0,34)

IF llbb THEN BEGIN
    cdelt = astrom.CDELT
    IF astrom.REVERSE THEN cdelt = REVERSE(cdelt)
    IF astrom.COORD_SYS EQ 'C' THEN BEGIN
        prec = CEIL(-ALOG10(cdelt[1]*3600.0)) > 0
    ENDIF ELSE BEGIN 
        prec = CEIL(-ALOG10(ABS(cdelt))) + 1
        nchar = prec + 6
        ns1 = 9-nchar[0] 
        ns2 = 9-nchar[1]
        mid  = ns1 GT 0 ? ",'"+STRJOIN(REPLICATE(' ',ns1),'')+"'" : ''
        tail = ns2 GT 0 ? ",'"+STRJOIN(REPLICATE(' ',ns2),'')+"'" : ''
        nchar = STRTRIM(STRING(nchar),2)
        ndp   = STRTRIM(STRING(prec),2)
        lls = ["('   Long:',",",',  Lat:',"]
        ll_fmt = STRJOIN(lls+'F'+nchar+'.'+ndp)+')'
    ENDELSE
    
    x2 = x2 + llbar  &  y2 = y2 + bbbar
    IF x2 LT 0d0 THEN x2 += 360d0
    AD2XY, x2, y2, astrom, xpix, ypix
    coord = [xpix, ypix]
                                ; Check that position is in fitting box:
    maxpix = ROUND(coord)
    xy = WHERE(x EQ maxpix[0] AND y EQ maxpix[1])
    warn = xy EQ -1
    IF ASTROM.COORD_SYS EQ 'C' THEN BEGIN
       llstr =  adstring(x2, y2, prec)
       rastr = STRMID(llstr,0,11+prec)
       decstr = STRMID(llstr,12+prec)
       llstr = '  RA:' + rastr + ', Dec:' + decstr
    ENDIF ELSE llstr = STRING(x2, y2, FORMAT = ll_fmt)
    line = [line, STRING(flux, FORMAT = "('        Fitted peak ',A,';')")+ $
            llstr]
    fmt = "('        Image pixel coordinates '," + flf+",', ',"+flf+")"
    line = [line, STRING(coord, FORMAT = fmt)]
    hpx = STRCMP(STRMID(astrom.CTYPE[0],5), 'HPX')
    IF (hpx) THEN BEGIN
        nside = T[1] / 5
        theta = (90d0 - y2)/r2d  & phi = x2/r2d
        ang2pix_ring, nside, theta, phi, maxpixr
        ang2pix_nest, nside, theta, phi, maxpixn
        line = [line, STRING(maxpixn, maxpixr, FORMAT = $ 
               "( '        in HEALPix pixels: ',I9,' (nest); ',I9,' (ring)')")]
    ENDIF
ENDIF ELSE BEGIN ; no astrometry
    coord = [x2, y2] + [x0, y0]
    fmt = "('        Fitted peak ',A,' at pixel coord ',"+flf+",', ',"+flf+")"
    line = [line, STRING(flux, coord, FORMAT = fmt)]
    ; Reality check
    warn = (x2 LT 0 || x2 GT 2 || y2 LT 0 || y2 GT 2)
ENDELSE
IF warn THEN line = [line, $
                     '        WARNING: Fitted peak outside 3x3 fitting box']

QUIT:

IF line[0] NE '' THEN BEGIN
    line = [line,' ']
    PRINT, line, FORMAT = "(A)"  ; Should put each string on its own line
    IF N_ELEMENTS(lun) NE 0 THEN PRINTF, lun, line, FORMAT="(A)"
ENDIF

IF N_ELEMENTS(x2) EQ 0 THEN BEGIN ; fit failed, use maximum in box.
   peak = boxmax
   IF llbb THEN XY2AD, coord[0], coord[1], astrom, x2, y2
   warn = 1B
ENDIF

pos = llbb ? [x2,y2] : REPLICATE(!values.F_NAN,2)
result = {maxpix: coord, pos: pos, flux: peak, unit: unit, warn: warn}
RETURN, result
END

