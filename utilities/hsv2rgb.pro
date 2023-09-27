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
PRO hsv2rgb, hue, sat, value,  badcol, rgbptr, satmin, satmax, HVONLY = hvonly
;+
; NAME:
;       HSV2RGB
;
; PURPOSE:
;       Converts images representing hue, saturation and value to
;       byte images representing Red, Green and Blue.
;
; CATEGORY:
;       Color Table Manipulation, Image Processing
;
; CALLING SEQUENCE:
;
;       HSV2RGB, Hue, Sat, Value,  Badcol, Rgbptr, Satmin, Satmax
;
; INPUTS:
;       Hue:    real-valued array: should be scaled 0 to 360 on input
;       Sat:    real-valued array: will be scaled 0 to 1
;       Value:  either byte-valued array: scaled 0 to 255 on input,
;               or real-valued array: scaled internally.
;       Badcol: byte value for blanks
;
; KEYWORD PARAMETERS:
;       HVONLY: Set saturation to 1 throughout (in which case Sat
;               should be a dummy value)
;
; OUTPUTS:
;      Rgbptr:  Array of 3 pointers to byte-scaled R, G, B arrays
;               NB: new heap variables are allocated.
;      Satmin:  The input sat value that maps to greys (always zero)
;      Satmax:  The input sat value that maps to fully-saturated
;               colours (equal to the maximum of sat on input).
;
; PROCEDURE:
;       Algorithm from A. R. Smith (1978, Computer Graphics, 12, pp 12-19),
;       via Wikipedia! Tuned to minimize memory overheads.
;
; EXAMPLE:
;       To display an image, first using pseudocolour and then using
;       hue/saturation to code the phase and relative amplitude of a
;       complex image:
;
;          vals = BYTSCL(image)
;          DEVICE, DECOMPOSED = 0
;          TV, vals                               ; Monochrome/Pseudocolour
;
;          hue = !radeg * ATAN(comp_image, /PHASE)
;          negs = WHERE(hue LT 0., nneg)
;          IF nneg GT 0 THEN hue[negs] += 360.0
;          sat = ABS(comp_image) / image          ; NB no scaling required
;          hsv2rgb, hue, sat, vals, 0B, rgb       
;          DEVICE, /DECOMPOSED
;          TV, *rgb[0], *rgb[1], *rgb[2], TRUE=3
;
; MODIFICATION HISTORY:
;       Written by:      J. P. Leahy March 2008
;       5/Jul/2014: Bug fixed: max in histogram needs setting.
;-
COMPILE_OPT IDL2

do_sat = ~KEYWORD_SET(hvonly)
hsize = SIZE(hue)
vsize = SIZE(value)
ssize = SIZE(sat)
good1 = ARRAY_EQUAL(hsize, ssize) || ~do_sat
good2 = ARRAY_EQUAL(hsize[0:hsize[0]], vsize[0:vsize[0]])
IF ~(good1 && good2) THEN MESSAGE, 'Input arrays are mismatched in size'

scale_v = vsize[vsize[0]+1] NE 1
npix = hsize[hsize[0]+2]

IF do_sat THEN satmin = MIN(sat, MAX = satmax, /NAN) ELSE BEGIN
    sat = 1.0
    satmin = 0.0
    satmax = 1.0
ENDELSE
satmin = 0.0

IF scale_v THEN BEGIN
    valmin = MIN(value, MAX=valmax, /NAN)
    scale  = 255/(valmax - valmin)
ENDIF
rgbptr = PTRARR(3,/ALLOCATE_HEAP)
FOR i = 0,2 DO *rgbptr[i] = BYTARR(hsize[1:hsize[0]],/NOZERO)
rchan = [0,2,1,1,3,0]
gchan = SHIFT(rchan,2)
bchan = SHIFT(rchan,-2)

lchunk = 1024L^2 ; Loop overheads should be negligible for chunks this big
IF npix LT lchunk THEN lchunk = npix
nchunk = DIVUP(npix, lchunk)
idc = LINDGEN(lchunk)
last = npix - (nchunk-1L)*lchunk
FOR ichunk=0,nchunk-1 DO BEGIN   
    hfin = FINITE(hue[idc])
    vfin = scale_v ? FINITE(value[idc]) : value[idc] NE badcol
    hvfin = hfin OR vfin
    bad = do_sat ? WHERE(~hvfin OR ~FINITE(sat[idc])) : WHERE(~hvfin)

    hfin = 0  &  vfin = 0

    f = hue[idc]    
    f /= 60
    hi = BYTE(f)
    hi MOD= 6
    f -= FLOAT(hi)
    IF do_sat THEN BEGIN
        p = sat[idc]                         
        p -= satmin
        p /= (satmin - satmax)  ; = -s (now s scaled from 0 to 1)
    ENDIF ELSE p = REPLICATE(-1.0, N_ELEMENTS(idc))

    q = p                       ; = -s
    p += 1.0                    ; = (1-s)
    q *= f                      ; = -f*s
    f = 0
    t = p                       ; = (1-s)
    t -= q                      ; = (1-s + f*s) = (1 - (1-f)*s)
    q += 1.0                    ; = (1-f*s)
    v = value[idc]
    IF scale_v THEN BEGIN
        v -= valmin
        v *= scale
    ENDIF

    hash = PTRARR(4, /ALLOCATE_HEAP)
    p = BYTE(p*v)
    *hash[1] = TEMPORARY(p)
    q = BYTE(q*v)
    *hash[2] = TEMPORARY(q)
    t = BYTE(t*v)
    *hash[3] = TEMPORARY(t)
    v = BYTE(v)
    *hash[0] = TEMPORARY(v)
    void = HISTOGRAM(hi, BINSIZE = 1, MIN = 0, MAX = 5, REVERSE_INDICES = ri)
    FOR i=0,5 DO BEGIN
        IF ri[i] NE ri[i+1] THEN BEGIN
           idx = ri[ri[i]:ri[i+1]-1]
           idcx = idc[idx]
           (*rgbptr[0])[idcx] = (*hash[rchan[i]])[idx]
           (*rgbptr[1])[idcx] = (*hash[gchan[i]])[idx]
           (*rgbptr[2])[idcx] = (*hash[bchan[i]])[idx]
        ENDIF
    ENDFOR
    ri = 0 
    PTR_FREE, hash

    IF bad[0] GE 0 THEN BEGIN
        idcx = idc[bad]
        FOR j=0,2 DO (*rgbptr[j])[idcx] = badcol
        idcx = 0
    ENDIF

    IF ichunk LT nchunk - 2 THEN idc += lchunk ELSE $
      idc = LINDGEN(last) + (nchunk-1L)*lchunk
ENDFOR

END 
