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
pro linecols, start, colindex
;+
; NAME:
;       LINECOLS
;
; PURPOSE:
;       Sets up simple colour table for plotting graphs.
;
; CATEGORY:
;       Color table manipulation
;
; CALLING SEQUENCE:
;
;       LINECOLS, Start, Colindex
;
; INPUTS:
;       Start:    Offset from 0 for starting colour
;
; OUTPUTS:
;       Colindex: Colour indices for:
;                    0 Black
;                  1:3 primaries:       red   green   blue
;                  4:6 complementaries: cyan  magenta yellow
;                  7:9 primary pastels,
;                10:12 complementary pastels 
;                   13 orange
;                 last white  (if start is too high, some colours may
;                 be missed, but the last index is always for white).
;
; SIDE EFFECTS:
;       If the current device does not use decomposed colour, LINECOLS
;       replaces colour table values from start to start+13. It also
;       always sets colour table index 0 to black and 255 to white.
;
;       For decomposed devices, the returned Colindex values are 24-bit
;       colour codes.
;
; EXAMPLE:
;       
;       LINECOLS, 0, cols
;       phi = FINDGEN(100)*0.1
;       PLOT, SIN(phi)
;       OPLOT, COS(phi), COLOR = col[1]  ; For red line.
;
; MODIFICATION HISTORY:
;       Written by:      J. P. Leahy March 2008 
;-
DEVICE, GET_DECOMPOSED = decomp

red   = [0, 255, 0, 0, 0, 255, 255, 255, 127, 127, 127, 255, 255, 255, 225]
green = [0, 0, 255, 0, 255, 0, 255, 127, 255, 127, 255, 127, 255, 127, 255]
blue  = [0, 0, 0, 255, 255, 255, 0, 127, 127, 255, 255, 255, 127,   0, 255]

IF decomp EQ 1 THEN BEGIN
    colindex = [red + 256L*(green + 256*blue)]
ENDIF ELSE BEGIN
    ncolours = (255 - start) < N_ELEMENTS(red)
    idx = LINDGEN(ncolours)
    TVLCT, red[idx], green[idx], blue[idx], start
    TVLCT,  0B,   0B,   0B,   0B 
    TVLCT, 255B, 255B, 255B, 255B
    colindex = idx + start
ENDELSE

END
