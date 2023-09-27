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
; Module GSCROLL
;
; J. P. Leahy 2008
;
; Scrolling subroutine inspired by Google Earth
;
;
; Method:   We maintain a number of coord systems:
;           * Original image(s) (they must all be the same size).
;             Notionally subdivided into a number of panels. Each
;             panel contains xpanel x ypanel *output* pixels.
;             The images passed to gscroll are assumed to be scaled to bytes.
;
;             Information about each panel is contained in the
;             structure array image_panel, which contains corner
;             coords and information about wrapping. It includes a
;             border of "virtual" panels which can  be filled in to
;             include wrapping information. Only one image_panel array
;             is maintained as it is common to all input maps.
;
;           * A pixmap, divided into a number of panels.
;             This is about the size of the full graphics screen
;             of the display device, but consists of an integral
;             number of panels. It will be slightly larger than the
;             graphics screen if allowed, otherwise slightly smaller.
;
;             There are two indices the pixmap:
;                pixmap_panel (structure) indexes the corresponding
;                image_panel, its rotation, and the corner coords on pixmap.
;                virtual_panel (INT array) lists the pixmap panels
;                sorted so that the "focus" panel, i.e. the one
;                containing the demanded pixel, is always the centre one.
;
;                The key feature of gscroll is that the central
;                panels of the virtual pixmap are kept up to date by
;                either copying panels from elsewhere in the pixmap if
;                they are loaded (with the correct rotation), or
;                loading them explicitly via the TV command if not.
;
;             There is a separate pixmap_panel array for each input
;             image, stored in the screens structure for that image,
;             since the different images may not be updated
;             synchronously. But there is only one virtual_panel
;             array, since the *wanted* arrangement of panels on the
;             pixmap is the same for all images.
;
;           * The current windows are updated by copying the relevant
;             section of the pixmaps. If the field of view extends off
;             the pixmap, the remaining region is filled in from a
;             special "blank" pixmap.
;
;  "housekeeping" parameters + the index arrays are stored in the
;  common block GRID_GSCROLL. Some meanings:
;
;  Started: not set or zero: all arrays & windows should be
;                            initialised.
;           1:               Partly initialised, some parameters
;                            (e.g. wrap control) may change.
;           2:               Fully initialised.
;
;  There is also a started field in the structure for each field, with
;  the same meanings
;
; Changes:
;   Original version Nov 2007-Jan 2008
;   Version 2: Handles multiple screens: Feb 2008
;              Minor update July 2013
;              Suppress "rolling will be disabled" message August 2013
;
FUNCTION shorter, vec1, vec2
COMPILE_OPT IDL2, HIDDEN

RETURN, TOTAL(vec1*vec1) < TOTAL(vec2*vec2)
END

PRO set_wrap, option, ierr
;
; set up wraps.
; Watch out for edges of image not coinciding with edges of panels.
;
COMPILE_OPT IDL2, HIDDEN

COMMON GRID_GSCROLL

CASE option OF
    "healpix_grid": BEGIN
        nside = (nx_in-2*nx_border)*xp_image/5
        npps_x = nside/xp_image
        npps_y = nside/yp_image
        IF npps_x*xp_image NE nside || npps_y*yp_image NE NSIDE THEN BEGIN
;            PRINT, 'SET_WRAP: non-integral number of panels in each facet'
;            PRINT, '          rolling will be disabled.'
            ierr = 8
            GOTO, QUIT
        ENDIF
        x1 = npps_x + nx_border-1  &  y1 = npps_y + ny_border-1
        image_panel[0:x1,0:y1].wrap = $
          cellindex(0, x1, 0, y1, nx_in, 4*npps_x, 4*npps_y)
        image_panel[0:x1,0:y1].wrap_rot = 0
        x0 = 4*npps_x + nx_border & y0 = 4*npps_y + ny_border
        x1 = nx_in - 1 & y1 = ny_in - 1
        image_panel[x0:x1,y0:y1].wrap = $
          cellindex(x0,x1,y0,y1,nx_in,-4*npps_x,-4*npps_y)
        image_panel[x0:x1,y0:y1].wrap_rot = 0
    END
    "torus": BEGIN
        image_panel[0,*].wrap = cellindex(0,nx_in-1,ny_in-1,ny_in-1,nx_in)
        image_panel[*,0].wrap = cellindex(nx_in-1,nx_in-1,0,ny_in-1,nx_in)
        image_panel[nx_in-1,*].wrap = cellindex(0,nx_in-1,1,1,nx_in)
        image_panel[*,ny_in-1].wrap = cellindex(1,1,0,ny_in-1,nx_in)
        image_panel[0,0].wrap = (ny_in-1)*nx_in - 2
        image_panel[nx_in,ny_in].wrap = nx_in + 1
        image_panel[0,ny_in].wrap = 2*nx_in - 2
        image_panel[nx_in,0] = (ny_in-2)*nx_in + 1
    END
    "off": BEGIN
        image_panel.wrap = -1L
        image_panel.wrap_rot = 0
    END
    ELSE: ; Do nothing
ENDCASE

QUIT:

END

PRO new_zoom, zoom, zfac, xhalf, yhalf, image_size, ierr
;
; Zoom or view size has changed; reset lookups
;
COMPILE_OPT IDL2, HIDDEN

COMMON GRID_GSCROLL

; Remember old data if it exists... we may need to reset.
IF started EQ 2 THEN BEGIN
    old_xpi = xp_image  &  old_ypi = yp_image
ENDIF
IF zoom THEN BEGIN
    xp_image = xpanel/zfac
    yp_image = ypanel/zfac
    IF xp_image*zfac NE xpanel || yp_image*zfac NE ypanel THEN BEGIN
        PRINT, 'GSCROLL: Incompatible zoom and panel size'
        ierr = 3
        GOTO, QUIT
    ENDIF
    IF xp_image LE 4 && yp_image LE 4 THEN BEGIN
;       PRINT, 'GSCROLL: Maximum zoom in reached'
        ierr = 6
    ENDIF
    IF xp_image LT 4 && yp_image LT 4 THEN BEGIN
        IF started EQ 2 THEN BEGIN
            xp_image = old_xpi  &  yp_image = old_ypi
        ENDIF
        ierr = 7
        GOTO, QUIT
    ENDIF
ENDIF ELSE BEGIN
    xp_image = xpanel*zfac
    yp_image = ypanel*zfac
    IF xp_image GE image_size[1] && $
      yp_image GE image_size[2] THEN BEGIN
;       PRINT, 'GSCROLL: Maximum zoom out reached'
        ierr = 6
    ENDIF
    IF started EQ 2 && xp_image GT 2*image_size[1] && $
      yp_image GT 2*image_size[2] THEN BEGIN
        IF started EQ 2 THEN BEGIN
            xp_image = old_xpi  &  yp_image = old_ypi
        ENDIF
        ierr = 7
        GOTO, QUIT
    ENDIF
; Note that very small maps are allowed on initialization (started < 2).
ENDELSE

; Analyse input map into panels
IF image_size[0] NE 2 THEN BEGIN
    ierr = 4
    PRINT, 'GSCROLL: image should be 2-D!'
    GOTO, QUIT
ENDIF
; Border should fill half a screen.
nx_border = divup(!D.x_vsize, 2*xpanel)
ny_border = divup(!D.y_vsize, 2*ypanel)
nx_in = divup(image_size[1], xp_image) + 2*nx_border
ny_in = divup(image_size[2], yp_image) + 2*ny_border

; Remove any border nasties from the template image panel structure,
; and turn into an array:
image_panel[0,0].wrap = -1L & image_panel[0,0].wrap_rot = 0
image_panel = REPLICATE(image_panel[0,0], nx_in, ny_in)

; Fill in corner coords of each panel
x0 = (LINDGEN(nx_in)-nx_border)*xp_image ; NB border has negative x0 and y0
x1 = x0 + xp_image - 1
FOR j=0,ny_in-1 DO BEGIN
    image_panel[*,j].x0 = x0
    image_panel[*,j].x1 = x1
    image_panel[*,j].y0 = REPLICATE((j-ny_border)*yp_image,nx_in)
    image_panel[*,j].y1 = REPLICATE((j-ny_border+1)*yp_image-1,nx_in)
ENDFOR
x0 = 0 & x1 = 0
; Mark dummy border panels (only for use with wrapping)

image_panel.dummy = image_panel.x0 LT 0 OR image_panel.y0 LT 0 OR $
                       image_panel.x0 GT image_size[1]-1 OR $
                       image_panel.y0 GT image_size[2]-1
indices = 0

; Try to set wraps if requested
ierr2 = 0
IF do_wrap THEN set_wrap, 'healpix_grid', ierr2
ierr = ierr + ierr2

; Reset indices for pixmap and virtual map:
screens.pixmap_panel.idx = -1
virtual_panel = LINDGEN(nx,ny)

; Reset effective central pixel on view window:
xhalf = !D.x_vsize / 2
yhalf = !D.y_vsize / 2
IF zoom THEN BEGIN          ; make xhalf, yhalf land on a pixel centre
    nbigpix = xhalf/zfac
    xhalf = zfac*nbigpix + zfac/2
    nbigpix = yhalf/zfac
    yhalf = zfac*nbigpix + zfac/2
ENDIF

QUIT:

END

FUNCTION gscroll, image, xpix, ypix, xhalf, yhalf, zoom_factor, $
                  ierr, new_view, done, DO_WRAP=wrap
;
; Main program: Updates pixmap as needed, and returns pixmap
; coordinates of image pixel x,y.
;
; gscroll_setup must be called before gscroll is called.
;
; Input:
;     image:       Pointer array, each to a large byte image
;                  (no longer used: data passed via screens array in COMMON).
;     xpix, ypix:  image pixel coords for pixel to be displayed at
;                  centre of screen
;     xhalf,yhalf: Offset wanted for xpix, ypix from BLC of screen
;     zoom_factor: >1 means each input pixel is mapped to a
;                   zoom_factor x zoom_factor grid of output pixels.
;                  <1 means that output pixels sample every
;                  1/zoom_factor input pixel (no averaging).
;     new_view:    = 1 Indicates new screen view: program allows itself up
;                  to a second to fill in the pixels. (Normally allows
;                  30 msec for pixmap loading operations).
;                  = 2 The window has been re-sized. Call new_zoom to
;                  re-set indices.
;     do_wrap:     Set if you want to wrap. Wrap must be specified in
;                  the set_wrap subroutine
; Outputs:
;     xhalf, yhalf are updated on change of zoom.
;     ierr:        returns with non-zero value if it detects a
;                  problem.
;     done:        = 1 if all wanted panels loaded, = 0 if it ran out
;                  of time.
; Returns: coordinates of wanted pixel in pixmap
;
; Output:    The pixmap is updated to contain as much as possible of
; the required field of view at the requested zoom factor, within the
; time allowed (maxtime). This is then copied to the view window.
;
COMPILE_OPT IDL2, HIDDEN

COMMON GRID_GSCROLL

ierr = 0
maxtime = 0.03
; in seconds. Adjust for efficient processing. NB, automatically reset to
; 1 second for the call after a change of zoom.

IF  N_ELEMENTS(started) EQ 0 THEN started = 0

IF started EQ 0 THEN BEGIN
    PRINT, 'GSCROLL: Call GSCROLL_SETUP before calling GSCROLL'
    ierr = 5
    GOTO, QUIT
ENDIF

wrap = KEYWORD_SET(wrap)
; Check to see if do_wrap is initialised:
IF N_ELEMENTS(do_wrap) EQ 0 THEN do_wrap = 0
; This will be set to true later if wrapping turns out to be possible.

T = size_image
IF xpix LT 0 || xpix GE T[1] || ypix LT 0 || ypix GE T[2] THEN BEGIN
    PRINT, xpix, ypix, $
      FORMAT="(/'GSCROLL: Requested pixel (',I6,',',I6,') is not on image.')"
    ierr = 1
    GOTO, QUIT
ENDIF
IF N_ELEMENTS(new_view) EQ 0 THEN new_view = 0

screen = screens[0]

WSET, screen.viewwin

; Decode zoom factor to integers:
zoom = zoom_factor GE 1.0
zfac = zoom ? FIX(zoom_factor) : FIX(1./zoom_factor)

IF (zoom_factor NE old_zoomf) || (wrap && ~do_wrap) $
  || (new_view EQ 2) THEN BEGIN
    maxtime = 1.0 ; Give it time to fill the screen
    do_wrap = wrap
    new_zoom, zoom, zfac, xhalf, yhalf, T, ierr
    IF ierr EQ 8 THEN BEGIN
        ; Problem with wrapping
        do_wrap = 0
    ENDIF
    IF ierr EQ 7 THEN BEGIN
                                ; exceeded max/min zoom.
        zoom_factor = old_zoomf
        zoom = zoom_factor GE 1.0
        zfac = zoom ? FIX(zoom_factor) : FIX(1./zoom_factor)
        IF new_view NE 2 THEN maxtime = 0.03
    ENDIF
    IF ierr GT 0 && ierr LT 6 THEN GOTO, QUIT

    screen = screens[0]  ; pixmap_panel has been updated
    old_zoomf = zoom_factor
ENDIF
started = 2

IF new_view THEN maxtime = 1.0

; End of setup. Now locate requested FOV, update indices, update
; pixmap, update view window, and we're done!

; Find the image panel containing the wanted pixel:
iip2 = [xpix/xp_image + nx_border, ypix/yp_image + ny_border]
iip =  iip2[0] + iip2[1]*nx_in
cur_ip = image_panel[iip]

; Find the position of our pixel within the panel
IF zoom THEN BEGIN
    xoff = (xpix - cur_ip.x0)*zfac + zfac/2
    yoff = (ypix - cur_ip.y0)*zfac + zfac/2
ENDIF ELSE BEGIN
    xoff = divup(xpix - cur_ip.x0,zfac)
    yoff = divup(ypix - cur_ip.y0,zfac)
ENDELSE

; Is the wanted panel loaded? Is it the current focus panel?
cur_ipp = gscroll_find_panel(iip, screen.pixmap_panel, rotation)
old_ipp = virtual_panel[focus]
old_pp  = screen.pixmap_panel[old_ipp]

; If the wanted pixel is not in the current panel, move the current
; panel to where the pixel is.
IF cur_ipp NE -1 THEN BEGIN
    IF cur_ipp NE old_ipp THEN BEGIN
        ; find shift
        temp = ARRAY_INDICES(virtual_panel,[cur_ipp,old_ipp])
        shift = temp[*,1] - temp[*,0]
        virtual_panel = SHIFT(virtual_panel, shift)
    ENDIF
ENDIF ELSE IF old_pp.idx NE -1 THEN BEGIN
; Current panel is not loaded, but old one is, so find where we are
; relative to that.
    temp = ARRAY_INDICES(image_panel,old_pp.idx)
    shift = temp - iip2

; Check for shortcut on wrapping
    IF do_wrap EQ 2 THEN BEGIN
        iwrap = image_panel[iip].wrap
        owrap = image_panel[old_pp.idx].wrap
        IF iwrap GE 0 THEN BEGIN
            temp1 = ARRAY_INDICES(image_panel,iwrap)
            shift1 = temp - temp1
            shift = shorter(shift,shift1)
        ENDIF
        IF owrap GE 0 THEN BEGIN
            temp1 = ARRAY_INDICES(image_panel,owrap)
            shift1 = temp1 - iip2
            shift = shorter(shift,shift1)
        ENDIF
    ENDIF
    virtual_panel = SHIFT(virtual_panel,shift)
ENDIF
; NB no shift is done if we are in current panel or we can't locate
; it.

; Reset current pixmap panel allowing for shift of focus, to get the
; pixmap coordinates of the central pixel:
cur_pp = screen.pixmap_panel[virtual_panel[focus]]
xcen = cur_pp.xc + xoff  &  ycen = cur_pp.yc + yoff

; Now load any data needed into the pixmap and then update the view
; window(s):
gscroll_load, image, zoom, zfac, xoff, yoff, xhalf, yhalf, iip2, maxtime, $
  done, 0B

RETURN, [xcen, ycen]

QUIT:

END
