; -----------------------------------------------------------------------------
;
;  Copyright (C) 2007-2013   J. P. Leahy
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
FUNCTION scale_image, image, range, wrap, trfunc, zero, beta, BAD = bad, $
                      BOTCOL = botcol, TOPCOL = topcol, $,
                      ABS_RANGE = abs_range, VERBOSE = verbose
; J. P. Leahy 2008
;
; Inputs revised 2013: BAD, BOT, TOP moved to keywords; defaults given in
; GR_GLOBAL common.
;
; Scales an image to bytes (intended for TV display)
;
; Inputs:
;   image:   The array to scale (not necessarily a 2-D image).
;            Assumed single precision floating point.
;   wrap:>0: values outside range are mapped to x MOD (top-bot), where
;            x is the value they would have if the scaling continued
;            smoothly outside [bot,top].
;         0: values are truncated to [bot,topcol].
;        <0: values below range[0] are truncated, values above
;            range[1] are wrapped.
;   trfunc:  Scaling transfer function to use (default 0 = linear)
;   zero:    estimated true zero level        (default 0)
;   beta:    Beta parameter for Asinh scaling (default 1)
;   verbose: If set, details of ranges are printed.
;
; KEYWORD INPUTS:
;            Not needed if GR_GLOBAL common is set up.
;   bad:     Special byte value for invalid (NaN etc) elements of image.
;            Default 0
;   botcol:  Byte assigned to range[0]. Default: 0
;   topcol:  Byte assigned to range[1]. Default: Max in colour table.
;
; Input/Output
;   range:   Values mapped to bytes botcol and topcol.
;            If only one value supplied, it is used as the top of the
;            range and the image minimum is used as the bottom.
;            If not set on input, is set internally to the full image
;            range.
;
; Keyword Outputs:
;   abs_range:  min,max on image.
;
; Returns:
;   Byte image
;
COMPILE_OPT IDL2, HIDDEN
ON_ERROR, 2
COMMON gr_global, windev, redraw_req, colmap, badcol, syscol, ibot, itop

start = SYSTIME(1)

IF N_ELEMENTS(badcol) EQ 0 THEN badcol = 0B
IF N_ELEMENTS(syscol) EQ 0 THEN syscol = -1
IF N_ELEMENTS(ibot) EQ 0 THEN ibot = 0B
IF N_ELEMENTS(itop) EQ 0 THEN itop =  (!D.table_size-1) < 255

verbose = KEYWORD_SET(verbose)

; Check for integer data types
is_float = WHERE(SIZE(image, /TYPE) EQ [1,2,3,12,13,14,15]) EQ -1

r2 = MAX(image, /NAN, MIN = r1)
abs_range = [r1, r2]
IF verbose THEN MESSAGE, /INFORMATIONAL, STRING(abs_range, FORMAT = $
                "('Minimum & Maximum on image:   ', 2(1X,E11.3))")

; Check that range is defined. It may have been undefined by being set
; to a string:
is_str = WHERE(SIZE(range, /TYPE) EQ 7) NE -1
IF ~is_str THEN BEGIN
    r2 = MAX(range)
    IF N_ELEMENTS(range) GT 1 THEN r1 = MIN(range)
ENDIF

IF r2 EQ r1 THEN BEGIN ; Stretch null range.
    IF abs_range[0] LT r1 THEN r1 = abs_range[0] ELSE $
      IF abs_range[1] GT r2 THEN r2 = abs_range[1] ELSE r2 += 1
    IF verbose THEN MESSAGE, /INFORMATIONAL, $
      'Min = Max! Stretching range for display'
ENDIF
range = [r1, r2]
IF verbose THEN MESSAGE,/INFORMATIONAL,  STRING(range, FORMAT = $
                "('Minimum & Maximum for scaling:', 2(1X,E11.3))")

IF N_ELEMENTS(bad)    EQ 0 THEN bad    = badcol
IF N_ELEMENTS(botcol) EQ 0 THEN botcol = ibot
IF N_ELEMENTS(topcol) EQ 0 THEN topcol = itop
IF N_ELEMENTS(trfunc) EQ 0 THEN trfunc = 0
IF N_ELEMENTS(zero)   EQ 0 THEN zero   = 0.0
IF N_ELEMENTS(beta)   EQ 0 THEN beta   = 1.0

absent = !P.background
IF syscol[0] NE -1 THEN BEGIN
    greys = syscol
    test = WHERE(absent EQ greys, hit)
    IF hit EQ 0 THEN greys = [greys,absent]
    test = WHERE(bad EQ greys, hit)
    IF hit EQ 0 THEN greys = [greys, bad]
ENDIF ELSE greys = [absent, bad]

greys  = greys[ UNIQ( greys, BSORT(greys) ) ]

idx = WHERE(greys GE botcol AND greys LE topcol, ngrey)
IF ngrey GT 0 THEN greys = greys[idx]
IF botcol GT 0 THEN BEGIN
    bots  = LINDGEN(botcol)
    greys = [bots,greys]
    ngrey = ngrey+botcol
ENDIF
gnext =  ngrey GT 1 ? [greys[1:*],topcol+1] : [topcol+1]
topcol = topcol - ngrey

ntop = topcol + 1

npix = N_ELEMENTS(image)
dims = npix GT 1 ? SIZE(image, /DIMENSIONS) : 1

dbyte  = topcol
drange = FLOAT(r2 - r1)
; Scale to appropriate intensity range (and squash NANs):
CASE trfunc OF
    0: BEGIN  ; Linear
        scale = dbyte / drange
    END
    1: BEGIN  ; Asinh
        asr   = ASINH((range - zero)/beta)
        scale = dbyte / (asr[1] - asr[0])
        r1 = zero
    END
    2: BEGIN  ; Sqrt
        scale = dbyte / SQRT(drange)
    END
    3: BEGIN  ; Histogram equalization: no memory shortcut!
        IF npix EQ 1 THEN RETURN, botcol ; Can't hist-equal one pixel!

        test = FINITE(image)
        good = WHERE(test, ngood)
        test = 0
        IF ngood EQ npix THEN BEGIN
            good = 0
            ngood = -1
        ENDIF
        CASE ngood OF
            -1: bchunk = HIST_EQUAL(image, MAXV = r2, MINV = r1, $
                                     TOP = dbyte)
            0: bchunk = 0       ; Do nothing
            ELSE: BEGIN
                bchunk = HIST_EQUAL(image[good], MAXV = r2, MINV = r1, $
                                    TOP = dbyte)
            END
        ENDCASE
        IF ngrey GT 0 THEN BEGIN
            null = HISTOGRAM(bchunk, BINSIZE = 1B,MIN = greys[0], $
                             OMIN = omin, OMAX = omax, REVERSE_INDICES = ri)
            FOR ig = 0, ngrey-1 DO BEGIN
                i0 = greys[ig] - omin - ig
                i1 = (gnext[ig] < (omax+1)) - omin - ig - 2
                p0 = ri[i0]  & p1 = ri[i1+1]
                IF p0 LT p1 THEN bchunk[ ri[p0 : p1-1] ] += (ig+1)
;                idx = WHERE(bchunk GE greys[ig])
;                IF idx[0] NE -1 THEN bchunk[idx] += 1
            ENDFOR
        ENDIF
        idx = 0
        CASE ngood OF
            -1: byte_image = TEMPORARY(bchunk)
             0: byte_image = REPLICATE(bad, dims)
             ELSE: BEGIN
                 byte_image = REPLICATE(bad, dims)
                 byte_image[good] = bchunk
             END
         ENDCASE
                                ; Skip wrapping etc: N/A in this case
         GOTO, DONE
    END
    ELSE: MESSAGE, 'Unknown scaling function:' + STRING(trfunc)
ENDCASE

lchunk = 1024L^2 ; Loop overheads should be negligible for chunks this big
IF npix LT lchunk THEN lchunk = npix
nchunk = DIVUP(npix, lchunk)
idc = LINDGEN(lchunk)
last = npix - (nchunk-1L)*lchunk
byte_image = REPLICATE(bad, dims)

FOR i = 0, nchunk-1 DO BEGIN
    chunk = is_float ? image[idc] : FLOAT(image[idc])
    goodpix = WHERE(FINITE(chunk))
    IF goodpix[0] EQ -1 THEN byte_image[idc] = bad ELSE BEGIN
        chunk = chunk[goodpix]
                                ; Do these ops with compound
                                ; assignment to avoid memory
                                ; overheads:
        chunk -= r1
        CASE trfunc OF
            0:                  ; Linear
            1: BEGIN            ; Asinh
                chunk /= beta
                chunk = ASINH(TEMPORARY(chunk))
                chunk -= asr[0]
            END
            2: BEGIN            ; Sqrt
                chunk >= 0.0
                chunk = SQRT(TEMPORARY(chunk))
            END
        ENDCASE
        chunk *= scale

        IF wrap LE 0 THEN lows  = WHERE(chunk LE 0.)           ELSE lows = -1
        IF wrap EQ 0 THEN highs = WHERE(chunk GE FLOAT(dbyte)) ELSE highs = -1

        chunk MOD= ntop
        negs = WHERE(chunk LT -0.5)
        IF negs[0] NE -1 THEN chunk[negs] += ntop
        bchunk = BYTE(ROUND(chunk))

        IF lows[0]  NE -1 THEN bchunk[lows]  = botcol
        IF highs[0] NE -1 THEN bchunk[highs] = topcol

        IF ngrey GT 0 THEN BEGIN
            null = HISTOGRAM(bchunk, BINSIZE = 1B,MIN = greys[0], $
                             OMIN = omin, OMAX = omax, REVERSE_INDICES = ri)
            FOR ig = 0, ngrey-1 DO BEGIN
                i0 = greys[ig] - omin - ig
                i1 = (gnext[ig] < (omax+1)) - omin - ig - 2
                p0 = ri[i0]  & p1 = ri[i1+1]
                IF p0 LT p1 THEN bchunk[ ri[p0 : p1-1] ] += (ig+1)
;                idx = WHERE(bchunk GE greys[ig])
;                IF idx[0] NE -1 THEN bchunk[idx] += 1
            ENDFOR
        ENDIF

        byte_image[idc[goodpix]] = bchunk
    ENDELSE

    IF i LT nchunk - 2 THEN idc += lchunk ELSE $
      idc = LINDGEN(last) + (nchunk-1L)*lchunk
ENDFOR

DONE:

topcol = topcol + ngrey

IF verbose THEN $
   MESSAGE, /INFORMATIONAL, STRING(minmax(byte_image),FORMAT="('output range:',2I4)")
; IF verbose THEN PRINT, 'Scaled in ',SYSTIME(1)-start,' seconds'

RETURN, byte_image

END
