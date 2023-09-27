; -----------------------------------------------------------------------------
;
;  Copyright (C) 2016-2020   J. P. Leahy
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
PRO ximgo, a0, a1, a2, a3, a4, a5, PIXELS = pixels, MARK = domark
;+
; NAME: 
;       XIMGO
;
; PURPOSE: 
;       Re-centres Ximview display at given position. Equivalent to running
;       Display -> Set View Centre (i.e. ximview_goto).
;
; CATEGORY:
;       Widget helpers
;
; CALLING SEQUENCE: 
;
;        XIMGO, a0, a1, [a2, a3, a4, a5,] [/PIXELS,] [/MARK]
;
; INPUTS:
;       a0 to a5  ra (a0 to a2) & dec (a3 to a5) in sexagecimal or
;                 ra (a0, a1) & dec (a2, a3) om sexagecimal (h,m; d,m)
;                 longitude (a0) & latitude (a1) in degrees or
;                 pixel coords (a0, a1)
;
; KEYWORD PARAMETERS:
;       PIXELS:  If set, longitude and latitude are x, y pixel
;                coordinates. Ignored with a warning if coords are hexadecimal.
;       MARK:    If set, put the mark on the new centre
;
; OUTPUTS:
;       None
;       
; SIDE EFFECTS:
;       Specified position becomes centre of Ximview display
;
; RESTRICTIONS:
;       None
;
; EXAMPLE:
;       Display some WMAP map data with XIMVIEW and copy the total
;       intensity image to the command line environment:
;
;            XIMVIEW, 'wmap_band_iqumap_r9_3yr_K_v2', '*', COL=[1,2,3]
;            XIMGO,  309.516,  19.417  ; Centres view on Centaurus A
;
; MODIFICATION HISTORY:
;       Written by:      J. P. Leahy, June 2016
;                        (cloned from ximget).
;                        Sept 2020: XIM_CATCH made external call, sexagesimal
;                        input options added
;-
COMPILE_OPT IDL2
ON_ERROR, 2
start = SYSTIME(1)

pixels = KEYWORD_SET(pixels)
domark = KEYWORD_SET(domark)
nargs = N_PARAMS()
CASE nargs OF
   2: BEGIN
      longitude = a0
      latitude  = a1
   END
   4: BEGIN ; Assume long = RA if given in sexagesimal
      longitude = TEN(a0,a1,0d0)*15d0
      latitude  = TEN(a2,a3,0d0)
      IF pixels THEN MESSAGE, /INFORMATIONAL, 'Ignoring /PIXELS'
      pixels = 0B
   END
   6: BEGIN   
      longitude = TEN(a0,a1,a2)*15d0
      latitude  = TEN(a3,a4,a5)
      IF pixels THEN MESSAGE, /INFORMATIONAL, 'Ignoring /PIXELS'
      pixels = 0B
   END  
   ELSE: MESSAGE, 'Two position coordinates must be specified'
ENDCASE

; Find the ximview widget and its control structures
xim_catch, top, state, mode
WIDGET_CONTROL, state.TABS, SENSITIVE = 0

swap_lut, *state.XIM_GRAPH, (*state.TABARR)[0], old_graph

pixvals = STRCOMPRESS(STRING(mode.X_CENTRE, mode.Y_CENTRE))

xhalf = mode.XHALF  &  yhalf = mode.YHALF
TVCRS, xhalf, yhalf

pix = [longitude, latitude]
r2d = 180d0 / !dpi

IF pixels THEN BEGIN
    xpix = pix[0]  &  ypix = pix[1]
ENDIF ELSE BEGIN
    IF ~state.IS_ASTROM THEN MESSAGE, "No astrometry available: can't convert coordinates to pixels"
    ll = pix[0]  &  bb = pix[1]
    AD2XY, ll, bb, *state.ASTROM, xpix, ypix
ENDELSE
  
xpix = (ROUND(xpix) > 0) < (state.IMSIZE[1]-1)
ypix = (ROUND(ypix) > 0) < (state.IMSIZE[2]-1)

; Now we've found the point, update the screen

IF mode.OVERVIEW THEN BEGIN
  mode.OVERVIEW = 0B
  WIDGET_CONTROL, state.ZOOMCOL, /SENSITIVE
  WIDGET_CONTROL, state.READOUT, /SENSITIVE
ENDIF

ierr = 0
coord = gscroll(void, xpix, ypix, xhalf, yhalf, $
  mode.ZOOM_FACTOR, ierr, 1, done, DO_WRAP = mode.ROLL)
IF ierr NE 0 THEN MESSAGE, 'GSCROLL error '+STRING(ierr)

mode.OXTV     = xhalf  &  mode.OYTV     = yhalf
mode.X_CENTRE = xpix   &  mode.Y_CENTRE = ypix
mode.XPIX     = xpix   &  mode.YPIX     = ypix
mode.DONE     = done   &  mode.NEW_VIEW = 0
mode.DRAG     = 0B

IF domark THEN BEGIN
   mode.XPT2 = mode.XPT1 & mode.YPT2 = mode.YPT1
   mode.XPT1 = mode.XPT  & mode.YPT1 = mode.YPT
   mode.XPT  = xpix      & mode.YPT = ypix
ENDIF
                                ; Re-plot graphics ovelays
overlay, mode, state.ASTROM

; If panels remain to be loaded, send timer event to DRAW window
; requesting re-draw:
IF done EQ 0 THEN WIDGET_CONTROL, (*state.TABARR)[0].DRAW, TIMER = 0.5

WIDGET_CONTROL, state.TABS, SET_UVALUE=mode
pix_print, state, 0, start

restore_lut, dummy, old_graph

WIDGET_CONTROL, state.TABS, /SENSITIVE

END

