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
PRO gscroll_tidy, CIAO=ciao
;
; J. P. Leahy 2008
;
; Tidies up and restores graphics state.
;
; Inputs:
;   ciao:  Does not delete pixmaps
;
; Changes:
;   Original version Nov 2007-Jan 2008
;   Version 2: Handles multiple screens: Feb 2008
;   
;
COMPILE_OPT IDL2
COMMON GRID_GSCROLL

IF started EQ 0 THEN RETURN

IF KEYWORD_SET(ciao) THEN BEGIN
    PRINT, pixwin, blankwin, FORMAT= $
      "(/'GSCROLL: Leaving pixmap (hidden) windows numbers',I3,' and',I3)"
    started = 1
ENDIF ELSE BEGIN
    IF N_ELEMENTS(blankwin) NE 0 THEN WDELETE, blankwin
    nscreen = N_ELEMENTS(screens)
    null = WHERE(screens.viewwin EQ old_window, old)
    IF old EQ 0 THEN BEGIN
        FOR i=0,nscreen-1 DO WDELETE, screens[i].viewwin
        SET_PLOT, old_device
        WSET, old_window
    ENDIF

    screens.started = 0
    FOR i=0,nscreen-1 DO BEGIN
        WDELETE, screens[i].pixwin
        ; Cancel LUT pointer which might have been locally generated.
        IF PTR_VALID(screens[i].LUT) THEN PTR_FREE, screens[i].LUT
        ; All other pointers are generated in calling program.
    ENDFOR
    DEVICE, DECOMPOSED = old_decomposed
    started = 0
ENDELSE

END

