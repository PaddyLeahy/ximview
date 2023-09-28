; SIMPLE_COLOUR
;
; Copyright (C) 2007, 2012, 2023 J. P. Leahy 
;
;    This program is free software: you can redistribute it and/or modify
;    it under the terms of the GNU General Public License as published by
;    the Free Software Foundation, either version 3 of the License, or
;    (at your option) any later version.
;
;    This program is distributed in the hope that it will be useful,
;    but WITHOUT ANY WARRANTY; without even the implied warranty of
;    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
;    GNU General Public License for more details.
;
;    You should have received a copy of the GNU General Public License
;    along with this program.  If not, see <https://www.gnu.org/licenses/>.
;
;+
; NAME:
;      SIMPLE_COLOUR
;
; PURPOSE:
;      Sets up the plot device and colour table for simple line 
;      diagrams, and optionally opens an output EPS file.
;
; CATEGORY:
;      DIRECT GRAPHICS
;
; CALLING SEQUENCE:
;      simple_colour, c1, cols, PS=psfile     
; INPUTS:
;      None
; OPTIONAL INPUTS:
;      None
; KEYWORD PARAMETERS:
;      psfile    Name of file for encapsulated postscript output.
; OUTPUTS:
;      None
; OPTIONAL OUTPUTS:
;      c1    Original value of !P.color
;      cols  Original colour table (to allow later restoration)
; COMMON BLOCKS:
;      None
; SIDE EFFECTS:
;      If PS is unset, sets the plot device to WIN or X as appropriate,
;      with RETAIN=2 and DECOMPOSE=0. Otherwise sets the device
;      to PS and opens the specified plot file as /COLOR, /ENCAPSULATED.
;
;      Sets colour table as follows:
;      1 red 2 green 3 blue 4 cyan 5 magenta 6 yellow 7 pink 8 apple 9
;      sky blue 10 Pale cyan 11 Pale magenta 12 pale yellow 13 white/black
;      14 orange.
;
;      Colour 13 is set white unless the postscript device is selected
;      in which case it is black. The default background colour is
;      unchanged, also colours 15 and higher.
; RESTRICTIONS:
;      None. Don't forget to close the postscript device.
; PROCEDURE:
;      Sets colours 1-14 explicitly. Changes device to PS if PS is set.
; EXAMPLE:
;
;      simple_colour, c1, oldlut
;      PLOT, x, fun1             ; white on black
;      OPLOT, x, fun2, col=1     ; red
;      OPLOT, x, fun2, col=3     ; blue
;
;      simple_colour, PS='myplot.eps'
;      PLOT, x, fun1             ; black on white
;      PLOT, x, fun2, col=1      ; red
;      DEVICE, /CLOSE
;      simple_colour ; back to X or WIN
;
;      !P.color = c1
;      TVLCT, oldlut  ; back to whatever colour table we had before.
;
; MODIFICATION HISTORY:
;      J. P. Leahy, ages ago.
;      Version 1.0  Doc head added 12 Oct 2012
;      Version 2.0  Darker colours for PS printing
;-
pro simple_colour, c1, cols, PS=psfile
;

; Remember previous colour setting just in case:
TVLCT, r, g, b, /GET
cols = [[r],[g],[b]]
c1 = !P.color

IF (KEYWORD_SET(psfile)) THEN BEGIN

    SET_PLOT,'PS'
    DEVICE, /COLOR, FILENAME=psfile, /ENCAPSULATED
    !P.THICK = 3
    !P.CHARTHICK = 2
    !P.CHARSIZE = 1.1
    !X.THICK = 3
    !Y.THICK = 3
;    Set colour table to primaries, complementaries, and dark pastels:
    TVLCT, [200, 0, 0, 0, 255, 127, 127,  64,  64,  64, 127, 127, 0, 200], $
           [0, 127, 0, 200, 0, 127,  64, 127,  64, 127,  64, 127, 0, 113], $
           [0, 0, 255, 200, 255, 0,  64,  65, 127, 127, 127,  64, 0, 0], 1

ENDIF ELSE BEGIN
    system = STRUPCASE(!VERSION.os_family)
    windev = STRCMP(system, 'WINDOW', 6) ? 'WIN' : 'X'
    SET_PLOT,windev
    DEVICE, RETAIN=2, DECOMPOSE=0
    !P.THICK = 1
    !P.CHARTHICK = 1
    !P.CHARSIZE = 1
    !X.THICK = 1
    !Y.THICK = 1
;    Set colour table to primaries, complementaries, and pastels:
    TVLCT, [255, 0, 0, 0, 255, 255, 255, 127, 127, 127, 255, 255, 255, 225], $
      [0, 255, 0, 255, 0, 255, 127, 255, 127, 255, 127, 255, 255, 127], $
      [0, 0, 255, 255, 255, 0, 127, 127, 255, 255, 255, 127, 255, 0], 1
ENDELSE

END
