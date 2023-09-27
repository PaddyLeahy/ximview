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
;
; Module gscroll_load:
;
;   Business end of gscroll that actually updates the pixmap(s) and
;   view window(s).
;
; Routines:
;   visible_panels:   Finds which panels are visible in the view window
;   load_panel:       Loads a panel onto a pixmap
;   gscroll_pix2view: Copies from pixmap to view window
;   gscroll_load:     Calls above routines as needed.
;
PRO visible_panels, xoff, yoff, xhalf, yhalf, iip2, w_ipp, w_iip
;
; Finds which panels are visible when view is centred at xoff, yoff in
; the focus panel.
;
; Inputs: 
;    xoff, yoff:    coordinates of target ("central") pixel in focus
;                   panel.
;    xhalf, yhalf:  coordinates of central pixel in view window
;    iip2:          coordinates of the focus panel in image_panel
;
; Outputs:
;    w_ipp:         List of indices for pixmap_panel array
;    w_iip:         List of indices for image_panel array
;                   (sorted on distance from focus panel).
;
COMPILE_OPT IDL2, HIDDEN
COMMON GRID_GSCROLL

x0 = xoff - xhalf & y0 = yoff - yhalf
x1 = x0 + !D.x_vsize - 1   & y1 = y0 + !D.y_vsize - 1

x0 = divup(x0,xpanel) & y0 = divup(y0,ypanel) 
x1 = x1/xpanel & y1 = y1/ypanel

focus2 = [focus MOD nx, focus/nx]
wanted = cellindex(x0, x1, y0, y1, nx, focus2[0], focus2[1], /SORT)

nw = N_ELEMENTS(wanted)
                                ; Find which image panels we want:
w_ipp = virtual_panel[wanted]
w_vp2 = ARRAY_INDICES(virtual_panel, wanted) 
w_iip2 = w_vp2 + REBIN(iip2-focus2, 2, nw, /SAMPLE)
w_iip = REFORM(w_iip2[0,*] + w_iip2[1,*]*nx_in)

                                ; Test for panels off the edge of the
                                ; panel array: 
bad = WHERE(w_iip2[0,*] LT 0 OR w_iip2[0,*] GT nx_in-1 OR $
            w_iip2[1,*] LT 0 OR w_iip2[1,*] GT ny_in-1)
IF bad[0] NE -1 THEN w_iip[bad] = -2

RETURN

END

PRO load_panel, ipp, iip, pixmap_panel, zoom, zfac, image
;
; Loads a given image panel into a given pixmap panel
; The pixmap must be the current window on entry.
;
; ipp:          index to pixmap panel
; iip:          index to image panel 
; pixmap_panel: associated  with current screen
; zoom:         true if zoomed in
; zfac:         integer factor for zoom or rebin.
; image:        pointer to the image to load (bytes)
;
COMPILE_OPT IDL2, HIDDEN
COMMON GRID_GSCROLL

pp = pixmap_panel[ipp]

; Check that wanted panel exists, and, if not, use blank panel:
blank =  iip LT 0 || iip GT nx_in*ny_in-1 
IF ~blank THEN BEGIN
    ip = image_panel[iip]
    blank = ip.dummy
    IF do_wrap && blank && ip.wrap GT -1 THEN BEGIN
        ; Wanted pixel does not exist but it has a wrapped counterpart
        iip=ip.wrap 
        ip = image_panel[iip]
        blank = 0
    ENDIF 
ENDIF 

IF blank THEN BEGIN
  DEVICE, COPY=[0, 0, xpanel, ypanel, pp.xc, pp.yc, blankwin] 
  GOTO, QUIT
ENDIF

; Try to find panel on existing pixmap:
pan = gscroll_find_panel(iip, pixmap_panel, rotation)
IF pan NE -1 THEN BEGIN
    pan = pixmap_panel[pan]
    DEVICE, COPY=[pan.xc, pan.yc, xpanel, ypanel, pp.xc, pp.yc]
ENDIF ELSE BEGIN
; RGB or monochrome?
    nchan = N_ELEMENTS(image)
    FOR i=0,nchan-1 DO IF PTR_VALID(image[i]) THEN T = SIZE(*image[i])

; Make sure we don't try to bin data off the edge of the image:
; NB: only happens at top & right as we start at (0,0)
    x0  = ip.x0               &  y0  = ip.y0 
    x1c = ip.x1 < (T[1] - 1)  &  y1c = ip.y1 < (T[2] - 1)
    dx  = x1c - x0 + 1        &  dy  = y1c - y0 + 1
    IF zoom THEN BEGIN
        dxout = dx*zfac  &  dyout = dy*zfac
    ENDIF ELSE BEGIN
        dxout = dx/zfac              &  dyout = dy/zfac
        x1c   = dxout*zfac + x0 - 1  &  y1c   = dyout*zfac + y0 - 1
    ENDELSE

; Bin the part of the image that we have:
    test = T[T[0]+1] ; data type
    IF test NE 1 THEN MESSAGE, STRING(test, FORMAT = $
                                      "('Data type =',I3,' not bytes.')") 

    IF nchan EQ 1 THEN BEGIN
        b1 = REBIN((*image)[x0:x1c,y0:y1c], dxout, dyout, /SAMPLE)
        IF rotation NE 0 THEN b1 = ROTATE(b1, rotation) 
                                ; worry about borders!
        TV, b1, pp.xc, pp.yc
    ENDIF ELSE BEGIN
        block = dxout*dyout
        b1 = BYTARR(nchan*block, /NOZERO)
        FOR i=0,nchan-1 DO BEGIN
            IF PTR_VALID(image[i]) THEN BEGIN
                b0 = REBIN((*image[i])[x0:x1c,y0:y1c], dxout, dyout, /SAMPLE)
                IF rotation NE 0 THEN b0 = ROTATE(b0, rotation) 
                b0 = REFORM(b0, block, /OVERWRITE)
            ENDIF ELSE b0 = BYTARR(block)
            b1[i*block] = b0
        ENDFOR
        b1 = REFORM(b1, dxout, dyout, nchan, /OVERWRITE)
        TV, b1, pp.xc, pp.yc, TRUE = 3
    ENDELSE

    IF dx*dy NE xp_image*yp_image THEN BEGIN
        T1 = SIZE(b1)
        rb = xpanel - T1[1]  &  tb = ypanel - T1[2]
; Pad it with blanks to make up to the right size:
        IF tb GT 0 THEN DEVICE, COPY=[0,0, xpanel, tb, pp.xc, pp.yc+T1[2], $
                                      blankwin]
        IF rb GT 0 THEN DEVICE, COPY=[0,0, rb, ypanel, pp.xc+T1[1], pp.yc, $
                                      blankwin]
    ENDIF
ENDELSE
; Update pixmap_panel index:

QUIT:

pixmap_panel[ipp].idx = iip
pixmap_panel[ipp].rot = rotation

END

PRO gscroll_pix2view, xhalf, yhalf, zoom, zfac, pixwin
;
; Copies from pixmap to view window (assumed set).
;
COMPILE_OPT IDL2, HIDDEN

COMMON GRID_GSCROLL

; Check that the screen window isn't wrapped in the pixmap:
x0 = xcen - xhalf & y0 = ycen - yhalf
x1 = x0+!D.x_vsize-1 & y1 = y0+!D.y_vsize-1

IF zoom && ((x0/zfac)*zfac NE x0) THEN $
  PRINT, 'Coordinate trouble:', x0, xcen, xhalf

x0c = x0 > 0    &  y0c = y0 > 0 
xst = x0c - x0  &  yst = y0c - y0 
xwc = (x1 < (xwidth-1)) - x0c + 1 
ywc = (y1 < (ywidth-1)) - y0c + 1
DEVICE, COPY = [x0c, y0c, xwc, ywc, xst, yst, pixwin]
;IF xst NE 0 || xwc NE !D.x_vsize THEN $
;    PRINT, 'GSCROLL partial X copy', xst, xwc

; Fill in edges
IF x0 LT 0 THEN DEVICE, COPY=[xwidth+x0, y0c, -x0, ywc, 0, yst, pixwin] $
           ELSE IF x1 GE xwidth THEN $
                DEVICE, COPY=[0, y0c, !D.x_vsize-xwc, ywc, xwc, yst, pixwin]

IF y0 LT 0 THEN DEVICE, COPY=[x0c, ywidth+y0, xwc, -y0, xst, 0, pixwin] $
           ELSE IF y1 GE ywidth THEN $
                DEVICE, COPY=[x0c, 0, xwc, !D.y_vsize-ywc, xst, ywc, pixwin]

; fill in corners
IF x0 LT 0 && y0 LT 0 THEN $
    DEVICE, COPY=[xwidth+x0,ywidth+y0,-x0,-y0,0,0,pixwin] $
ELSE IF x0 LT 0 && y1 GE ywidth THEN $
    DEVICE, COPY=[xwidth+x0,0,-x0,ywidth-ywc,0,ywc,pixwin] $
ELSE IF x1 GE xwidth && y0 LT 0 THEN $
    DEVICE, COPY=[0,ywidth+y0,xwidth-xwc,-y0,xwc,0,pixwin] $
ELSE IF x1 GE xwidth && y1 GE ywidth THEN $
    DEVICE, COPY=[0,0,xwidth-xwc,ywidth-ywc,xwc,ywc,pixwin]

END

PRO gscroll_load, image, zoom, zfac, xoff, yoff, xhalf, yhalf, iip2, maxtime, $
                  done, nocopy
;
; Check that the visible panels surrounding the current panel are
; loaded and if not, load them
;
; Inputs:
;  image:        Array of pointers to byte images 
;  (no longer used: data passed in screens array in GRID_GSCROLL instead)
;  zoom:         True if zoomed in
;  zfac:         integer zoom factor
;  xoff, yoff:   Coords of central pixel in focus panel
;  xhalf, yhalf: Coords of central pixel on screen
;  iip2:         Coords of central panel in image_panel
;  maxtime:      Time allowed for loading before returning to caller
;  nocopy:       If set, don't copy to view window.
;
; Output:
;  done:         True if all pixmaps fully updated.
;
COMPILE_OPT IDL2, HIDDEN
COMMON GRID_GSCROLL

IF ~N_ELEMENTS(nocopy) THEN nocopy = 0

visible_panels, xoff, yoff, xhalf, yhalf, iip2, w_ipp, w_iip
IF N_ELEMENTS(w_ipp) EQ 0 THEN MESSAGE, 'Apparently no panels in FOV!!'

done = 1
start = SYSTIME(1)
nscreen = N_ELEMENTS(screens)
FOR isc = 0, nscreen-1 DO BEGIN
    screen = 0
    screen = screens[isc]
    lut = screen.lut
    pixwin = screen.pixwin
    pixmap_panel = screen.pixmap_panel

    IF do_wrap THEN BEGIN
                                ; re-label any wrapped panels (assumes rot = 0)
        wraps = WHERE(w_iip GE 0 AND pixmap_panel[w_ipp].idx GE 0 AND $
                      pixmap_panel[w_ipp].idx NE w_iip AND $
                      pixmap_panel[w_ipp].idx EQ image_panel[w_iip].wrap)

;        IF wraps[0] NE -1 THEN pixmap_panel[w_ipp[wraps]].idx = w_iip[wraps]
        IF wraps[0] NE -1 THEN w_iip[wraps] = pixmap_panel[w_ipp[wraps]].idx
    ENDIF

    j = WHERE(pixmap_panel[w_ipp].idx NE w_iip, jcount)
    
    WSET, pixwin

    IF screens[isc].IS_RGB THEN BEGIN
        DEVICE, /DECOMPOSED
        im_scr = screens[isc].RGB
    ENDIF ELSE im_scr = screens[isc].IMPTR

    IF redraw || isc EQ 0 THEN TVLCT, (*lut).R, (*lut).G, (*lut).B

    FOR i=0,jcount-1 DO BEGIN
        ji = j[i]
        load_panel, w_ipp[ji], w_iip[ji], pixmap_panel, zoom, zfac, im_scr
        IF SYSTIME(1)-start GT maxtime THEN BEGIN
            done = 0
            BREAK
        ENDIF
    ENDFOR
    screens[isc].pixmap_panel = pixmap_panel

                                ; Finally, update view window
    IF (isc*top_only + nocopy) EQ 0 THEN BEGIN 
        WSET, screens[isc].viewwin
        gscroll_pix2view, xhalf, yhalf, zoom, zfac, pixwin
    ENDIF

    IF screens[isc].IS_RGB THEN DEVICE, DECOMPOSED = 0
    IF done EQ 0 THEN BREAK     ; abandon later screens
ENDFOR

IF isc GT 0 THEN WSET, screens[0].viewwin

END
