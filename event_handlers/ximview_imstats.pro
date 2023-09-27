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
;  Module XIMVIEW_IMSTATS
;
;  J. P. Leahy February 2008
;              Updated 2013 new LUT handling.
;              Updated 2020 passes beam area to imstats, passes image
;                           pointer instead of image, passes beam area
;                           picks up structure returned by imstats.
;
;  Contains the event handlers that manage calling the IMSTATS
;  routine.
;
FUNCTION ximview_imstats, event, imstat_out
;
;  Widget front end for imstats
;
;  INPUT 
;       event  calling event
;
;  RETURNS
;       ierr   0 no problem
;              1 no centre set for box
;              2 user abandoned ROI setting
;              3 illegal zoom mode for ROI
;
; OUTPUT
;      imstat_out  Structure with results from IMSTATS
;
COMPILE_OPT IDL2, HIDDEN
ON_ERROR, 1

ierr = 0

tagname_full = TAG_NAMES(event,/STRUCTURE_NAME)
WIDGET_CONTROL, event.TOP,   GET_UVALUE = state

ximview_resize, event.TOP, state
 
tagname = STRMID(tagname_full,0,12)
CASE tagname OF
   'WIDGET_BUTTO' : BEGIN  ; Actually this is the only tagname recognised here.
; get statistics
      log = state.LOGLUN
      astrom = state.IS_ASTROM ? *state.ASTROM : 0

      tabarr = state.TABARR
      str = (*tabarr)[0]
      im_ptr = str.IM_PTR  &  unit = str.UNIT & mult = str.MULT
      barea = str.BEAM.BEAM_AREA
      ncode = str.NOISE_TYPE
      threshold = str.THRESHOLD
      is_ptr = PTR_VALID(im_ptr)
      IF ~is_ptr && PTR_VALID(str.RGB[0]) THEN BEGIN
         MESSAGE, /INFORMATIONAL, 'Cannot measure statistics on RGB image'
         RETURN, 3
              ; Image is stored directly: move temoprarily to stack:
         im_ptr = PTR_NEW(/ALLOCATE_HEAP)
         WIDGET_CONTROL, state.LABEL, GET_UVALUE = *im_ptr, /NO_COPY
      ENDIF

      swap_lut, *state.XIM_GRAPH, str, oldgraph
       
      WIDGET_CONTROL, state.TABS,  GET_UVALUE = mode
      zoom = mode.ZOOM          &  zfac = mode.ZFAC
      xhalf = mode.xhalf        &  yhalf = mode.YHALF
      x_centre = mode.X_CENTRE  &  y_centre = mode.Y_CENTRE
      
      IF str.SCREEN NE state.LASTTAB THEN BEGIN
         WIDGET_CONTROL, str.BASE, GET_UVALUE = tname
         tname = 'Now on tab: ' + tname
         PRINT, tname
         PRINTF, log, tname
         state.LASTTAB = str.SCREEN
         WIDGET_CONTROL, event.TOP, SET_UVALUE = state
      ENDIF

      IF state.ROI THEN BEGIN
         IF zoom LT 0 THEN BEGIN ; Would give ROI with only subset of pixels
            ierr = 3
            msg = ['ROIs not enabled when zoom < 1', 'Zoom in and try again']
         ENDIF
         IF mode.OVERVIEW THEN BEGIN
            ierr = 3
            msg = ['ROIs not enabled in overview mode', $
                   'Enter zoom mode and try again']
         ENDIF
         IF ierr GT 0 THEN BEGIN
            ok = DIALOG_MESSAGE(msg, /ERROR, TITLE = 'XIMVIEW_IMSTATS', $
                                DIALOG_PARENT=event.top) 
            RETURN, ierr
         ENDIF
         zoomin = zoom GE 0 ? [zfac, zfac] : [1,1]
         offset = tv2im(0L,0L, mode)
         
         oldpan = mode.PAN      ; Disable panning
         mode.PAN = 0
                                ; Need to disable any button 3 actions
                                ; as well

         WIDGET_CONTROL, state.TABS, SET_UVALUE = mode
            
         imsize = state.IMSIZE[1:2]
         nx = imsize[0]  &  ny = imsize[1]

; Call CW_DEFROI to get the ROI, but don't use its OFFSET & IMAGE_SIZE
; keywords as this restricts usage with ROLL (and tends to start the
; cursor in a silly place due to a CW_DEFROI bug).
;
         REPEAT BEGIN
            region = CW_DEFROI(str.DRAW, ZOOM = zoomin)
            IF region[0] EQ -1 THEN BEGIN
               ans = DIALOG_MESSAGE('Try again?', /QUESTION, $
                        TITLE = 'XIMVIEW_IMSTATS', DIALOG_PARENT=event.top)
               IF ans EQ 'No' THEN BEGIN
                  ierr = 2
                  BREAK
               ENDIF
             ENDIF
         ENDREP UNTIL region[0] NE -1
         mode.PAN = oldpan      ; Re-activate normal pan/zoom mode
         WIDGET_CONTROL, state.TABS, SET_UVALUE = mode
         
         IF region[0] EQ -1 THEN RETURN, ierr ; operation cancelled

                                ; find real pixel values
         nxscreen = !D.x_vsize / zoomin[0]
         IF zoom LT 0 THEN BEGIN
            xpix = (region MOD nxscreen)*zfac + offset[0]
            ypix = (region / nxscreen)*zfac   + offset[1]
         ENDIF ELSE BEGIN
            xpix = (region MOD nxscreen) + offset[0]
            ypix = (region / nxscreen)   + offset[1]
         ENDELSE

         IF mode.ROLL THEN BEGIN
            ns4 = state.NS4
            top= nx+ns4
            idx = WHERE(xpix + ypix GT top)
            is_top = idx[0] GT -1
            IF is_top THEN BEGIN
               xpix[idx] -= ns4  &  ypix[idx] -= ns4
            ENDIF
            bot = state.NSIDE
            idx = WHERE(xpix + ypix LT bot)
            is_bot = idx[0] GT -1
            IF is_bot THEN BEGIN
               xpix[idx] += ns4  &  ypix[idx] += ns4
            ENDIF
            idx = 0
         ENDIF

;            gscroll_roi, xpix, ypix, imsize, x_centre, y_centre, $
;              xhalf, yhalf, mode.ZOOM_FACTOR

         imstat_out = imstats(im_ptr, xpix, ypix, ASTROM = astrom, $
                              LUN = log, UNIT = unit, BEAMA = barea, $
                              MULT = mult, THRESHOLD=threshold, NT = ncode) 

; Draw background ellipses if available:
         IF imstat_out.e1.sup GT 0 THEN BEGIN
            epar = centrep2param(imstat_out.e1)
            xy = param2points(epar)
            txy = im2tv(xy[*,0], xy[*,1], mode)
            PLOTS, txy[*,0], txy[*,1], /DEVICE
            epar = centrep2param(imstat_out.e2)
            xy = param2points(epar)
            txy = im2tv(xy[*,0], xy[*,1], mode)
            PLOTS, txy[*,0], txy[*,1], /DEVICE
            EMPTY
         ENDIF
      ENDIF ELSE BEGIN          ; Box around marked point

         IF mode.XPT LT 0 THEN BEGIN
            ierr = 1
            MESSAGE, /INFORMATIONAL, $
                     'Mark centre with middle mouse button first!'
            GOTO, QUIT
         ENDIF

         xpt = FIX(mode.XPT)  &  ypt = FIX(mode.YPT)

         imstat_out = imstats(im_ptr, xpt, ypt, state.STATBOX, $
                              MULT = mult, ASTROM = astrom, LUN = log, $
                              UNIT = unit, BEAMA = barea, NT = ncode, $
                              THRESHOLD = threshold) 

; Draw box:
         box = imstat_out.box
         IF ~mode.OVERVIEW && box[0] NE -1 THEN BEGIN
            blc = im2tv(box[0,0]-0.5, box[1,0]-0.5, mode)
            trc = im2tv(box[0,2]+0.5, box[1,2]+0.5, mode)
            PLOTS, [blc[0], blc[0], trc[0], trc[0], blc[0]], $
                   [blc[1], trc[1], trc[1], blc[1], blc[1]], /DEVICE
            EMPTY
         ENDIF
      ENDELSE

QUIT:
      restore_lut, dummy, oldgraph
      IF ~is_ptr THEN BEGIN
         WIDGET_CONTROL, state.LABEL, SET_UVALUE = *im_ptr, /NO_COPY
         PTR_FREE, im_ptr
      ENDIF
   END
   ELSE: MESSAGE, 'Unrecognised event type: ' + tagname_full
ENDCASE

RETURN, ierr
END
