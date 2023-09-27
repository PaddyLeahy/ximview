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
PRO gscroll_roi, xpix, ypix, image_size, x_centre, y_centre, xhalf, yhalf, $
                 zoom_factor
;
; Draws a Region of Interest (ROI). Works best with pseudo-colour.
;
; Inputs:
;    xpix, ypix:  Lists of pixel coords in ROI
;     ... the rest as for main gscroll function
;
COMPILE_OPT IDL2, HIDDEN

COMMON GRID_GSCROLL

screen = screens[0]

zoom = zoom_factor GE 1.0
zfac = zoom ? ROUND(zoom_factor) : ROUND(1./zoom_factor)

x0 = MIN(xpix)  &  x1 = MAX(xpix) 
y0 = MIN(ypix)  &  y1 = MAX(ypix)

xmax = image_size[0]-1
ymax = image_size[1]-1
IF x0 LT 0 OR x1 GT xmax OR y0 LT 0 OR y1 GT ymax THEN $
    MESSAGE, 'ROI indices fall off edge of map'

; Find blc panel...
blc_ip2 = [x0/xp_image + nx_border, y0/yp_image + ny_border]
blc_ip =  blc_ip2[0] + blc_ip2[1]*nx_in
x0 = (image_panel[blc_ip]).x0  &  y0 = (image_panel[blc_ip]).y0
; ...and trc
trc_ip2 = [x1/xp_image + nx_border, y1/yp_image + ny_border]
trc_ip =  trc_ip2[0] + trc_ip2[1]*nx_in
x1 = (image_panel[trc_ip]).x1  &  y1 = (image_panel[trc_ip]).y1

ncolx = x1 - x0 + 1  &  ncoly = y1 - y0 + 1

idx = (xpix - x0) + ncolx*(ypix - y0)

; Find current central pixmap panel and image panel
cur_pp = screen.pixmap_panel[virtual_panel[focus]]
cur_ip = cur_pp.idx
cur_ip2 = ARRAY_INDICES(image_panel,cur_ip)
; First find the position of centre pixel within the focus panel
cenpan = image_panel[cur_ip]
IF zoom THEN BEGIN
    xoff = (x_centre - cenpan.x0)*zfac + zfac/2
    yoff = (y_centre - cenpan.y0)*zfac + zfac/2
ENDIF ELSE BEGIN
    xoff = divup(x_centre - cenpan.x0, zfac)
    yoff = divup(y_centre - cenpan.y0, zfac)
ENDELSE
    
; Find the visible panels
visible_panels, xoff, yoff, xhalf, yhalf, cur_ip2, w_ipp, w_iip

IF do_wrap THEN BEGIN ; Doesn't work
                                ; re-label any wrapped panels (assumes rot = 0)
    wraps = WHERE(w_iip GE 0 AND screen.pixmap_panel[w_ipp].idx GE 0 AND $
                  screen.pixmap_panel[w_ipp].idx NE w_iip AND $
                  screen.pixmap_panel[w_ipp].idx EQ image_panel[w_iip].wrap)

    IF wraps[0] NE -1 THEN screen.pixmap_panel[w_ipp[wraps]].idx = w_iip[wraps]
ENDIF

; Find panels covered by ROI
w_iip2 = ARRAY_INDICES(image_panel,w_iip)
good = WHERE(w_iip2[0,*] GE blc_ip2[0] AND w_iip2[0,*] LE trc_ip2[0]-1 AND $
             w_iip2[1,*] GE blc_ip2[1] AND w_iip2[1,*] LE trc_ip2[1]-1,ngood)
IF ngood GT 0 THEN BEGIN 
    w_iip = w_iip[good] 
    w_ipp = w_ipp[good]
ENDIF ELSE MESSAGE, "Can't find ROI on image"

overlay = bytarr(ncolx, ncoly)
overlay[idx] = 255B
PRINT, ncolx, ncoly, N_ELEMENTS(xpix)
PRINT, MIN(idx), MAX(idx)
WSET, screen.pixwin
DEVICE, SET_GRAPHICS_FUNCTION = 7 ; bitwise OR operation

dx  = xp_image  &  dy  = yp_image
FOR i=0,ngood-1 DO BEGIN
    ip = image_panel[w_iip[i]]
    pp = screen.pixmap_panel[w_ipp[i]]
    x0p = ip.x0 - x0 &  y0p = ip.y0 - y0
    x1p = ip.x1 - x0 &  y1p = ip.y1 - y0
    PRINT, x0p, x1p, y0p, y1p
    IF MAX(overlay[x0p:x1p,y0p:y1p]) EQ 0B THEN CONTINUE ; Don't copy empties.
    PRINT, 'got data'
    IF zoom THEN BEGIN
        b1 = REBIN(overlay[x0p:x1p,y0p:y1p], dx*zfac, dy*zfac, /SAMPLE)
    ENDIF ELSE BEGIN
        dxout = dx/zfac               &  dyout = dy/zfac
        x1p   = dxout*zfac + x0p - 1  &  y1p   = dyout*zfac + y0p - 1
        b1 = REBIN(overlay[x0p:x1p,y0p:y1p], dxout, dyout, /SAMPLE)
    ENDELSE
; Bin the part of the image that we have:
    IF rotation NE 0 THEN b1 = ROTATE(b1,rotation) ; worry about borders!
    TV, b1, pp.xc, pp.yc
ENDFOR
DEVICE, SET_GRAPHICS_FUNCTION = 3 ; back to normal copy from source
WSET, screen.viewwin

                                ; Finally, update view window
xcen = cur_pp.xc + xoff  &  ycen = cur_pp.yc + yoff

; gscroll_pix2view, xcen, ycen, xhalf, yhalf, zoom, zfac

END
