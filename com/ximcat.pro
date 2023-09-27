; -----------------------------------------------------------------------------
;
;  Copyright (C) 2016   J. P. Leahy
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
PRO ximcat, catalog, coord_sys, hits
;+
; NAME: 
;       XIMCAT
;
; PURPOSE: 
;       Loads a catalog into ximview, which will be plotted.
;
; CATEGORY:
;       Widget helpers
;
; CALLING SEQUENCE: 
;
;        XIMCAT, catalog, coord_sys, hits
; INPUTS:
;       catalog: structure containing catalog data
;
; OPTIONAL INPUTS:
;       coord_sys: 'C' = RA/Dec, 'G' Galactic, 'E' = ecliptic
;                 (default: 'C') 
;       
; KEYWORD PARAMETERS:
;       None
;
; OUTPUTS:
;       hits   Byte array, same length as catalog. =1 (true) if source
;              is on current image and not at a blank pixel, = 0 otherwise
;       
; SIDE EFFECTS:
;       Specified positions are plotted by ximview
;
; RESTRICTIONS:
;       None
;
; EXAMPLE:
;       Display SUMSS map with XIMVIEW, extract sources in this
;       map from the SUMSS catalog, and plot these positions:
;
;            XIMVIEW, 'SUMSS\J0040M40.FITS'
;            this = WHERE(sumsscat.mosaic EQ 'J0040M40')
;            thiscat = {ra: sumsscat.ra[this], dec: sumsscat.dec[this]}
;            XIMCAT, thiscat, 'C', hits
;
; MODIFICATION HISTORY:
;       Written by:      J. P. Leahy, June 2016
;                        (cloned from ximget).
;                        Sept 2020: calls xim_catch instead
;                        of doing it all here.
;-
COMPILE_OPT IDL2
ON_ERROR, 0

IF N_ELEMENTS(coord_sys) EQ 0 THEN coord_sys = 'C'
ll = catalog.ra
bb = catalog.dec
symb   = 1        ; cross
colour = !P.COLOR ; Fixed for the moment.
do_label = 0B     ; no thanks

; Find Ximview
xim_catch, top, state, mode

; Make sure coords are consistent with astrometry of image
; and if not, convert.
co_sys = (*state.ASTROM).COORD_SYS
IF co_sys NE coord_sys THEN BEGIN
  IF co_sys NE 'E' and co_sys NE 'G' THEN MESSAGE, /INFORMATIONAL, $
    "I'm not going to try to convert catalogue coords!" $
  ELSE BEGIN
    ang2vec, bb, ll, vec, /ASTRO
    outvec = ROTATE_COORD(vec, INCO=coord_sys, OUTCO=co_sys)
    vec2ang, outvec, bb, ll, /ASTRO
    vec = 0
  ENDELSE
ENDIF

; Calculate image pixel coordinates:

AD2XY, ll, bb, *STATE.astrom, xpix, ypix

; only save positions actually on image
xp = ROUND(xpix) & yp = ROUND(ypix)
ax = (*state.astrom).naxis
on_map = WHERE(xp GE 0 AND yp GE 0 AND xp LT ax[0] AND yp LT ax[1], nin)

IF nin GT 0 THEN BEGIN 
   ll = ll[on_map]
   bb = bb[on_map]
   xpix = xpix[on_map]
   ypix = ypix[on_map]

   IF ARG_PRESENT(hits) THEN BEGIN
                                ; check that at valid pixels
                                ; NB only for displayed map, so keep
                                ; all sources on the grid.
      im_ptr = (*state.TABARR)[0].IM_PTR
      hit_id = WHERE(FINITE((*im_ptr)[xp[on_map],yp[on_map]]))
      hit_id = on_map[hit_id]
      hits = BYTARR(N_ELEMENTS(ll))
      hits[hit_id] = 1B
   ENDIF

   newcat = {SYMBOL: symb, COLOUR: colour, DO_LABEL: do_label, $
             LON: ll[on_map], LAT: bb[on_map], XPIX: xpix, YPIX: ypix}
  
; Stash the catalog. Note that mode.catalog is a pointer to an array
; of pointers.
   pcat = PTR_NEW(newcat)
   cats = *MODE.catalog
   IF N_ELEMENTS(cats) EQ 0 THEN cats = [pcat] ELSE cats = [cats,pcat]

   *mode.catalog = cats
   mode.CATPLOT = 1B

   WIDGET_CONTROL, state.TABS, SET_UVALUE=mode
ENDIF ELSE MESSAGE, /INFORMATIONAL, $
                    'No source positions in this catalogue on the image'

;ncat = N_ELEMENTS(cats)
;PRINT, 'Number of catalogs: ', ncat
END

