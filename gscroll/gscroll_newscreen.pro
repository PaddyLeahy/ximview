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
PRO gscroll_newscreen, iscreen, top_array, zoom_factor, xpix, ypix, $
                       xhalf, yhalf, done, nocopy
;
; Changes the current screen.
;
; Inputs:
;    iscreen:      The screen to make current
;    top_array:    Structure array from the high-level program with
;                  contents describing each screen. Returned shifted
;                  so that the current screen is in index 0.
;                  must include tags .BYTE_PTR with pointer to byte
;                  images and .SCREEN with permanent index of screen
;    zoom_factor:  What it says
;    xpix, ypix:   Image coords of central pixel on view window
;    xhalf, yhalf: Coords of central pixel on view window.
;    nocopy:       Update pixmap but don't copy to view window
; Output:
;    done:         True if all pixmaps fully updated
;
COMPILE_OPT IDL2, HIDDEN
COMMON GRID_GSCROLL

IF ~N_ELEMENTS(nocopy) THEN nocopy = 0

; Save current pixmap panel info
pixmap_panel = screens[0].PIXMAP_PANEL

; Give it a reasonable time to complete operations:
maxtime = 0.3

IF started EQ 0 THEN MESSAGE, $
  'Call GSCROLL_SETUP before calling GSCROLL_NEWSCREEN.'

; Check for deleted screens (assumes at most one):
dead = (WHERE(top_array.SCREEN EQ -1))[0]
IF dead NE -1 THEN BEGIN
    null = WHERE(screens.viewwin EQ old_window, old)

    IF ~old THEN WDELETE, screens[dead].viewwin ELSE BEGIN
        IF screens[dead].viewwin EQ old_window THEN BEGIN
            dscreen = WHERE(top_array.SCREEN EQ 0)
            old_window = screens[dscreen].viewwin
        ENDIF
    ENDELSE

    WDELETE, screens[dead].pixwin

    ishift = dead + 1           ; puts deleted element last
    top_array = SHIFT(top_array, -ishift)
    screens = SHIFT(screens, -ishift)

                                ; Remove last element
    ntab = N_ELEMENTS(screens)
    top_array = top_array[0:ntab-2]
    screens = screens[0:ntab-2]
ENDIF

; Find wanted screen
ishift = WHERE(top_array.SCREEN EQ iscreen)
IF ishift EQ -1 THEN MESSAGE, "Can't find requested screen:" + $
  STRING(iscreen,FORMAT="(I2)")

; Shift high- and low-level screen arrays
top_array = SHIFT(top_array, -ishift)
screens   = SHIFT(screens, -ishift)

IF started EQ 1 THEN IF nocopy THEN BEGIN
    done = 1
    RETURN
ENDIF ELSE MESSAGE, 'GRID_GSCROLL common not properly initialised'

; Update new current screen without shifting virtual_panel, so that
; the different screens stay synchronised

; Find the image panel containing the current central pixel
; NB: can't go via virtual panel as we may have loaded the wrapped
; equivalent of the wanted panel:
cip2 = [xpix/xp_image + nx_border, ypix/yp_image + ny_border]
cip =  cip2[0] + cip2[1]*nx_in
cur_ip = image_panel[cip]

; Check
iip = pixmap_panel[virtual_panel[focus]].idx
IF iip EQ -1 THEN MESSAGE, 'Previous screen not initialised'
IF iip NE cip AND iip NE cur_ip.WRAP THEN MESSAGE, $
  'Mismatched image panel indices!'

; Decode zoom factor to integers:
zoom = zoom_factor GE 1.0
zfac = zoom ? FIX(zoom_factor) : FIX(1./zoom_factor)

; Find the position of our pixel within the panel
IF zoom THEN BEGIN
    xoff = (xpix - cur_ip.x0)*zfac + zfac/2
    yoff = (ypix - cur_ip.y0)*zfac + zfac/2
ENDIF ELSE BEGIN
    xoff = divup(xpix - cur_ip.x0,zfac)
    yoff = divup(ypix - cur_ip.y0,zfac)
ENDELSE

iip2 = ARRAY_INDICES(image_panel,iip)
image = top_array.BYTE_PTR

gscroll_load, image, zoom, zfac, xoff, yoff, xhalf, yhalf, iip2, $
  maxtime, done, nocopy

END
