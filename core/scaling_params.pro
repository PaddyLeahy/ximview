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
; Module SCALING_PARAMS
;
; J. P. Leahy February 2008
;
; Makes a robust estimate of mode and standard deviation of an
; image. Part of the Ximview Auto-scale system.
;
;
FUNCTION good_data, map, nmap
;
; Finds how many pixels in a floating-point dataset (size nmap) are NANs,
; infinite, or suspiciously zero.
; If the latter, sets zeros = NAN.
;
COMPILE_OPT IDL2, HIDDEN
ON_ERROR, 2

bad = ~FINITE(map)
bad = WHERE(bad, nbad)
; Temporarily mask bad values as windows complains about
;  comparing NANs to anything
IF nbad GT 0L THEN map[bad] = -1.0
noughts = map EQ 0.0
IF nbad GT 0L THEN map[bad] = !values.F_NAN
bad = 0
noughts = WHERE(noughts, nn)
IF (nn GT 1E-5 * nmap) THEN BEGIN 
;    MESSAGE, /INFORMATIONAL, 'Seems to be using zero as blank'
    nbad += nn
    map[noughts] = !values.F_NAN
ENDIF

RETURN, nmap - nbad
END

FUNCTION get_chan, thing, data, imap, oned, howto, col
;
; Returns one channel of the data, as controlled by "howto" and "oned"
;
; Inputs:
;    thing or data: image or pointer to it
;    imap:          the requested channel (of those "wanted")
;    oned:          "image" is 1-D array
;    howto:         How the data is stored
;    col:           Array of wanted columns
COMPILE_OPT IDL2, HIDDEN

IF oned THEN BEGIN
    CASE howto OF
        0: map = thing[*,col[imap]-1]
        1: map = data[*,imap]
        2: map = (*data)[*,imap]
        3: map = *data[imap]
    ENDCASE
ENDIF ELSE BEGIN
    CASE howto OF
        0: map = thing[*,*,col[imap]-1]
        1: map = data[*,*,imap]
        2: map = (*data)[*,*,imap]
        3: map = *data[imap]
    ENDCASE
ENDELSE

IF WHERE(SIZE(map,/TYPE) EQ [4,5]) EQ -1 THEN map = FLOAT(map)

RETURN, map
END

PRO scaling_params, thing, data, imap, oned, howto, col, ar, mode, sdev, $
                    PLOT = do_plot
;
; Extracts robust estimate of zero level and standard deviation from map.
;
; Inputs:
;   thing, data, imap, oned, howto, col: see inputs for get_chan.
;
; Outputs:
;   ar:   full range of image
;   mode: Mode of image (robustly estimated)
;   sdev: negative-side RMS of image relative to zero.
;
COMPILE_OPT IDL2, HIDDEN
ON_ERROR, 2

do_plot = KEYWORD_SET(do_plot)
IF do_plot THEN BEGIN
    col0 = 240
    linecols, col0, ci
    !P.color = 255 ; white
    !P.background = 0B ; black
ENDIF
map = get_chan(thing, data, imap, oned, howto, col)
nmap = N_ELEMENTS(map)

mode = 0.
sdev = -1

CATCH, on_error
IF on_error NE 0 THEN BEGIN
    CATCH, /CANCEL
    MESSAGE, /REISSUE_LAST
    MESSAGE, /INFORMATIONAL, '# pixels (total, valid): '+ $
      STRJOIN(STRTRIM(STRING([nmap, ngood]),2),' ')
    IF N_ELEMENTS(bmax) NE 0 THEN $
      MESSAGE, /INFORMATIONAL, 'Histogram min, max, #bins:' + $
      STRJOIN(STRTRIM(STRING([bmin, bmax, nbin]),2),' ')
    IF do_plot THEN $
      XYOUTS, 0.5, 0.5, 'No scaling possible', ALIGNMENT = 0.5, /NORMAL
    RETURN
ENDIF

ngood = good_data(map, nmap)
ar = MINMAX(map, /NAN)

IF ngood EQ 0 THEN MESSAGE, 'Map all blank or zero'

IF ar[0] EQ ar[1] THEN BEGIN ; flat map
    mode = ar[0]
    sdev = 0.0
    RETURN
ENDIF  



mu = TOTAL(map, /NAN) / ngood
IF ~FINITE(mu^2) THEN BEGIN
   PRINT, 'Mean of map seems to be:', mu
   PRINT, 'Number of good pixels:', ngood
   MESSAGE, 'Mean is too close to machine maximum value', /CONTINUE
   MESSAGE, 'likely due to unmasked "bad pixel" values.', /CONTINUE
   MESSAGE, 'Convert bad to NaN before scaling.'
ENDIF

; find standard deviation with minimal use of memory
map -= mu
map ^= 2
sdev = SQRT( TOTAL(map, /NAN) / ngood )
map = 0

; First guess:
mode = mu 
; Choose binsize to make random fluctuation of bins next to max unlikely
; to change order (2 sigma cut)
nbin = ROUND((16.*nmap)^0.2)

map = get_chan(thing, data, imap, oned, howto, col)
ngood = good_data(map, nmap)
idx = INDGEN(5) - 2
FOR i=0,10 DO BEGIN
    bmin = (mode - 2*sdev) > ar[0]
    bmax = (mode + 2*sdev) < ar[1]
    hgm = HISTOGRAM(map, MIN = bmin, MAX = bmax, NBIN = nbin, $
                    /NAN, LOCATIONS = locs)
    null = MAX(hgm, imode)
;    imode = WHERE(hgm EQ MAX(hgm))
;    imode = imode[0]
    binsize = locs[1] - locs[0]
    locs += binsize / 2

                                ; Next guess in case we have to repeat:
    mode = locs[imode]
    sdev = MIN(ABS([bmin-mode, bmax-mode]))/2.
                                ; Check for very lopsided distribution
    IF imode GE nbin/5 AND imode LE (4*nbin)/5 THEN BREAK
ENDFOR

x = idx + imode
coeff = POLY_FIT(idx, hgm[x], 2, COVAR = covar, $
                 MEASURE_ERRORS = SQRT(hgm[x]), CHISQ=chisq)

mode = bmin + (imode + 0.5 - coeff[1]/(2.0*coeff[2]) )*binsize


IF do_plot THEN BEGIN
    PLOT, locs, hgm, PSYM = 10, XSTYLE = 1, TITLE = 'Mode analysis', $
      XTITLE = 'Brightness', YTITLE = 'Pixel count'

    OPLOT,locs[x], (coeff[0] + idx*(coeff[1] + idx*coeff[2])), COLOR = ci[1]
    ARROW, mode, MAX(hgm)/4, mode, MAX(hgm)/2, /DATA, COLOR = 2+col0
    XYOUTS, mode, MAX(hgm)/8, 'Adopted Mode', ALIGNMENT = 0.5, COLOR = ci[2]
ENDIF

; Find the one-sided SD (from points below the mode, unless the mean
; is below the mode, suggesting negative skewness, in which case from points
; above the mode.

idx = mu GT mode ? WHERE(map LE mode) : WHERE(map GE mode)
IF idx[0] NE -1 THEN low = map[idx] ELSE low = mode

map = 0
nmap = N_ELEMENTS(low)
ngood = good_data(low, nmap)

low -= mode
low ^= 2
sdev =SQRT( TOTAL(low, /NAN) / ngood)

END
