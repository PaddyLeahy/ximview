; -----------------------------------------------------------------------------
;
;  Copyright (C) 2007-2020   J. P. Leahy
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
; Module ximview: widget-based image inspection tool
;
; J. P. Leahy 2008 - 2020
;
; This file contains an IDL documentation block after the ximview
; procedure is declared. As well as ximview itself the file contains the
; event handlers for the overall widget and the image screen, the
; shorter menu-item event handlers, and an number of ancilliary
; routines. On start-up the main program calls the routines in
; make_ximview.pro to create the main widget. 
;
; Where is the data?
;
;  As a rule, bulk data is stored on the heap and accessed via IDL pointers.
;  Pointers, along with other metadata, are stored in structures, accessed
;  via the uservalues (UVALUE) of high-level widgets. 
;
;  The displayed image is copied from a hidden pixmap window maintained by 
;  the low-level GSCROLL software.
;
;  .TOP uservalue = "state" structure: contains widget IDs for key
;  widgets, and various parameters which are fixed or change rarely
;  during execution: values of input parameters (except the actual image),
;  widget geometry, SIZE output for input map (.IMSIZE),
;  colour indices, etc. The state structure is created in the main
;  Ximview routine.
;
;  state.TABARR is a pointer to a structure array, each member of
;  which contains the data for each tabbed screen, including the
;  pointer to the byte-scaled image, widget IDs and scale information,
;  and usually the pointer to the original image. The template for the
;  tab structures is also defined in the main Ximview routine.
;
;  state.TABS uservalue is a structure ("mode") containing all the
;  variable geometry parameters, switches for different operation modes,
;  zoom details. The mode structure is created in make_ximview, which
;  also creates the main ximview widget.
;
;  The mode structure contains a pointer to a list of catalogues that
;  may be overplotted.
;
;  state.LABEL uservalue used to be used for the image itself in some 
;  circumstances. This should no longer occur, but the code to use it
;  remains in place just in case.
;
;  state.XIM_GRAPH is a pointer to a structure containing the
;  graphics state variables required for ximview. Event handlers swap this
;  for the current graphics state while they run, and swap back at the end.
;
;  the .LUT tag in the tab structures is a pointer to a structure
;  containing the R,G,B arrays for the colour table asssociated with
;  the tab, along with the colour indices for absent pixels and line
;  graphics.
;--------------------------------------------------------------------------
;
; Utility routines for colour and scaling:
;
FUNCTION colour_schemes, index
; Returns the name of the colour scheme associated with index.
;
COMPILE_OPT IDL2, HIDDEN

schemes = ['Rainbow', 'Heat', 'Blue-yellow-white', 'Greyscale', $
           'Red-black-blue', 'Cyclic']

IF index EQ -1 THEN RETURN, schemes ELSE RETURN, schemes[index]
END

FUNCTION scale_funs, index
; Returns the name of the scaling function with index.
;
COMPILE_OPT IDL2, HIDDEN

trfunc = ['Linear','Asinh','Sqrt', 'Hist eq']

IF index EQ -1 THEN RETURN, trfunc ELSE RETURN, trfunc[index]
END

PRO set_scale, scale_pars, str
;
; Sets RANGE in tab structure str based on scale_pars and the colour table
; in effect.
;
COMPILE_OPT IDL2, HIDDEN
ON_ERROR, 2
; Global graphics parameters: 
;  Set in main Ximview procedure.
;     windev:     'X' (unix) or 'WIN' (MSWIN)
;     colmap:     true if a color map is available
;     redraw_req: graphics must be re-drawn to use a new color map
;     badcol:     index in LUT used for bad values
;     syscol:     indices in LUT used for grey levels indicating off-sky
;                 and bad data. Ideally, system standard colours.
;     ibot, itop  index in LUT for maximum and minimum of range.
; To be Set in ximview_lut (not yet implemented
;     m2_col      Colour for second-last marked point
;     m3_col      Colour for third-last marked point
;     grid_col    Colour for coordinate grid lines
;     lab_col     Colour for coordinate labels
;     geo_col     Colour for geometry lines 

COMMON gr_global, windev, redraw_req, colmap, badcol, syscol, ibot, itop, $
   m2_col, m3_col, grid_col, lab_col, geo_col

mindat = scale_pars[0]
maxdat = scale_pars[1]
zero   = scale_pars[2]
sdev   = scale_pars[3]

old_izero = str.izero

; Find which way up the intensity scale is:
ratio = (maxdat-zero) / ABS(mindat - zero)
top = 0
IF ratio GT 5 THEN top = 1 ELSE IF ratio LT 0.2 THEN top = -1
CASE str.TRFUNC OF
    1: range = [mindat, maxdat] ; Asinh scale
    3: BEGIN                    ; Histogram equalization
                                ; Avoid under-digitization on initial
                                ; analysis: of 3000 bins we want at
                                ; least 15 to cover +/- 1 sigma.
        range = ratio GT 1 ? mindat + [0, 400*sdev] : maxdat - [0, 400*sdev]
    END
    ELSE: BEGIN                 ; Func 0 or 2: linear or sqrt
        offsymm = [25., 12., 25., 12., 9., 25., 50.]
        IF str.TRFUNC EQ 2 THEN offsymm = (2./3.)*offsymm^2
        offasym = offsymm - 3.
        range = top*offasym[str.COLTAB] + [-1, 1]*offsymm[str.COLTAB]
        range = range*sdev + zero
    END
ENDCASE
range = (range > mindat) < maxdat
beta = [2., 2., 2., 2., 1.5, 2, 2]

str.ABSRANGE = [mindat, maxdat]
str.RANGE    = range
str.BETA     = sdev*beta[str.COLTAB]
str.ZERO     = zero
str.MODE     = zero
str.SDEV     = sdev
str.IZERO    = scale_image(zero, range, str.WRAP, str.TRFUNC, zero, str.BETA)

END

PRO fill_gores, nside, imsize, astrptr, byteptrs
;
; Set off-sky pixel values to "absent". In future should deal with
; projections other than HPX.
;
; Inputs:
;   nside: HEALpix parameter
;   imsize: SIZE array for image
;   astrptr: Pointer to astrometry structure
;   byteptrs: Array of pointers to byte image arrays.
;
COMPILE_OPT IDL2, HIDDEN

ntab = N_ELEMENTS(byteptrs)
astrom = *astrptr

proj = STRMID(astrom.CTYPE[0],5,3)

CASE proj OF
    'HPX': BEGIN
        nsmin1 = nside - 1L
        nx = imsize[1]
        ny = imsize[2]
        blc = 2.5*nside + 0.5 - astrom.CRPIX
        trc = blc + [nx, ny] - 1
        cropped = ~ARRAY_EQUAL(blc, [0,0]) || ~ARRAY_EQUAL(trc, 5*nside*[1,1])

        IF ~cropped THEN BEGIN
            ix  = REBIN( INDGEN(nside), nside, nside)
            iy  = TRANSPOSE(ix)
            idx = TEMPORARY(ix) + nx*TEMPORARY(iy)
            idx = REFORM(idx, nside*nside, /OVERWRITE)
            absar = REPLICATE(!P.background, nside*nside)
        ENDIF
                                ; BLC coordinates for empty panels:
        xc = nside*[2L, 3L, 4L, 3L, 4L, 0L, 4L, 0L, 1L, 0L, 1L, 2L]
        yc = nside*[0L, 0L, 0L, 1L, 1L, 2L, 2L, 3L, 3L, 4L, 4L, 4L]
        x0 = (xc - blc[0]) > 0
        y0 = (yc - blc[1]) > 0
        x1 = (xc+nside - 1L - blc[0]) < (nx - 1)
        y1 = (yc+nside - 1L - blc[1]) < (ny - 1)
        offsets = x0 + nx*y0
        dx = x1 - x0 + 1  &  dy = y1 - y0 + 1
                                ; Locate each missing facet and blank
        FOR i=0,11 DO BEGIN
            IF cropped THEN BEGIN
                nsx = dx[i]  &  nsy = dy[i]
                IF nsx LT 0 || nsy LT 0 THEN CONTINUE
                ix  = REBIN( INDGEN(nsx), nsx, nsy)
                iy  = TRANSPOSE( nsx EQ nsy ? ix : REBIN(INDGEN(nsy),nsy,nsx))
                idx = TEMPORARY(ix) + nx*TEMPORARY(iy)
                npix = nsx*nsy
                idx = REFORM(idx, npix, /OVERWRITE)
                absar = REPLICATE(!P.background, npix)
            ENDIF
            indices = idx + offsets[i]
            FOR j=0,ntab-1 DO BEGIN
                tptr = byteptrs[j]
                IF PTR_VALID(tptr) THEN (*tptr)[indices] = absar
            ENDFOR
        ENDFOR
    END
    'XPH': BEGIN
                ; Uses same code as hpgrid does to set pixel values to nan
 
        ns2 = 2L*nside
        ngrid = 4L*nside
        list = LINDGEN(ns2*nside)
        y1 = FIX(list / ns2)
        x1 = FIX((list MOD ns2) + nside)
        list = 0
        nss = FIX(nside)
        dx = ABS(x1 - 2S*nss)       ; x offset from grid centre
        dy = nss - y1               ; y offset from top of first row of facets

        bad = WHERE(TEMPORARY(dx) LT TEMPORARY(dy))  ; duplicate pixels
        x1 = x1[bad]
        y1 = y1[bad]
        bad = [x1, ngrid-1-x1,         y1, ngrid-1-y1] + $
              [y1, ngrid-1-y1, ngrid-1-x1,         x1]*ngrid
        x1 = 0 & y1 = 0 

        nbad = 4L*nside^2
        gore = REPLICATE(!P.background, nbad)
        FOR idim = 0L,ntab-1 DO BEGIN
            gptr = byteptrs[idim]
            (*gptr)[bad] = gore
        ENDFOR
    END
    ELSE: ; Don't know what to do with other projections
ENDCASE
END

FUNCTION ramp, nlevel, DOWN = down
; Calculates a linear ramp from 0B to 255B over nlevel levels
;
COMPILE_OPT IDL2, HIDDEN

ramp = BYTE( (255L * LINDGEN(nlevel)) / (nlevel-1) )
IF KEYWORD_SET(down) THEN ramp = REVERSE(ramp, /OVERWRITE)
RETURN, ramp
END

PRO ximview_lut, coltab, izero, decomp
;
; Set up colour table for XIMVIEW, with special values.
;
; Inputs:
;   coltab:     requested colour table
;   izero:      colour level corresponding to zero
;               (for blue-black-red scale only)
;
; Outputs:
;   decomp:     required value of device 'decomposed' state
;
; Set in gr_global: m2_col, m3_col etc
;
COMPILE_OPT IDL2, HIDDEN
COMMON gr_global

; Revise colour table
black_colour  = [  0B,   0B,   0B]
white_colour  = [255B, 255B, 255B]

black = 0B
grey = syscol
absent = syscol[0]

; Greys are special levels that must be left as gaps in the colour table
; to avoid flashing & make RGB images work well.
IF ~redraw_req THEN grey = [0B, grey]
ngrey = N_ELEMENTS(grey)
                                 ; temporarily shift zero index level
                                 ; if there are reserved greylevels below it
null = WHERE(grey LT izero, zgcount)
izero = izero - zgcount

bad_colour    = REPLICATE(badcol,3) ;  Neutral grey, as in MOLLVIEW et al.
absent_colour = REPLICATE(absent,3) ; Darker grey for off-sky pixels

ncol = !D.table_size
white = ncol - 1
line  = white - 1
bot = ibot
top = itop - ngrey

decomp = 0B

ntop = top + 1
tail = BYTARR(ncol-ntop)

CASE coltab OF
   0: BEGIN                  ; Rainbow pseudo colour
      basic = 39             ; Standard "Rainbow-white". Saturates at level 235
      line_colour   = [255B, 175B, 175B] ;  pink.
      LOADCT, basic, /SILENT             ; Expensive because reads from file!
      TVLCT, rnew,gnew,bnew, /GET
      rnew[top] = white_colour[0]
      gnew[top] = white_colour[1]
      bnew[top] = white_colour[2]
   END
   1: BEGIN                                ; Heat (black-red-white)
      line_colour   = [175B, 255B, 175B]   ; pale green
      rnew = [ramp(172), REPLICATE(255B, ntop-172), tail]
      gnew = [BYTARR(116), ramp(ntop-116), tail]
      bnew = [BYTARR(186), ramp(ntop-186), tail]
   END
   2: BEGIN                     ; Blue-yellow-white (colourblind equivalent)
      line_colour   = [255B, 175B, 175B] ;  pink.
      parab = BYTE(255*(1. - ((3.0/ntop)*FINDGEN(2*ntop/3) - 1.0)^2))
      n2 = ntop - (2*ntop/3) - 1
      bnew = [parab, parab[0:n2], tail]
      rnew = [BYTARR(ntop/3), ramp(ntop/4), $
              REPLICATE(255B, ntop - (ntop/3) - (ntop/4)), tail]
      gnew = rnew
   END
   3: BEGIN                                ; Greyscale
      absent_colour = [100B, 100B, 200B]   ; blue-grey
      bad_colour    = [200B, 255B, 200B]
      line_colour   = [100B, 255B, 100B] ; green
      rnew = [ramp(ntop), tail]
      gnew = rnew  &  bnew = rnew
   END
   4: BEGIN                                ; Blue-black-red
      line_colour   = [175B, 255B, 175B]   ; pale green
      pos_range = top - izero
      neg_range = izero
      peak = pos_range > neg_range
      p2  = peak/2
      p2b = peak - p2
      
      ramp0 = ramp(p2)
      ramp1 = [ramp0, REPLICATE(255B, p2b)]
      ramp2 = [BYTARR(p2b), TEMPORARY(ramp0)]
      IF neg_range GT 0 THEN BEGIN
         rnew =  [REVERSE(ramp1[0:neg_range-1]), 0B]
         gnew =  [REVERSE(ramp2[0:neg_range-1]), 0B]
         bnew =  [REVERSE(ramp2[0:neg_range-1]), 0B]
      ENDIF ELSE BEGIN
         rnew = [0B]  &  gnew = [0B]  &  bnew = [0B]
      ENDELSE
      IF pos_range GT 0 THEN BEGIN
         rnew = [rnew, ramp2[0:pos_range-1], tail]
         gnew = [gnew, ramp2[0:pos_range-1], tail]
         bnew = [bnew, ramp1[0:pos_range-1], tail]
      ENDIF ELSE BEGIN
         rnew = [rnew, tail]
         gnew = [gnew, tail]
         bnew = [bnew, tail]
      ENDELSE
   END
   5: BEGIN                                ; Cyclic
      line_colour   = [255B, 255B, 255B]   ; white
      null = FLTARR(ncol - ntop)
      one = [REPLICATE(1., ntop), null]
      hue = [360.* FINDGEN(ntop) / top, null] ; hue is in degrees
      
      TVLCT, hue, one, one, /HSV
      TVLCT, rnew, gnew, bnew, /GET
   END
   ELSE: MESSAGE, 'Unknown colour table'
ENDCASE

; Set colours for older marker points as faded versions of original:
m2_col = BYTE(line_colour*0.8)
m3_col = BYTE(line_colour*0.6)

DEVICE, DECOMPOSED = decomp

r = rnew  &  g = gnew  &  b = bnew

; Leave gaps in colour table for special greys, including absent and bad:
i1 = line-1  &  i2 = line-2
FOR ii = 0,ngrey-1 DO BEGIN
    i0 = grey[ii]
    r[i0+1:i1] = r[i0:i2]
    g[i0+1:i1] = g[i0:i2]
    b[i0+1:i1] = b[i0:i2]
    r[i0] = i0  &  b[i0] = i0  &  g[i0] = i0
ENDFOR

izero += zgcount

new_colours = TRANSPOSE([[absent_colour], [bad_colour], [line_colour], $
                         [white_colour]])
cols = [absent, badcol, line, white]
r[cols] = new_colours[*,0]
g[cols] = new_colours[*,1]
b[cols] = new_colours[*,2]

TVLCT, r, g, b

!P.color      = line
!P.background = absent

END

FUNCTION invert_scale, ivalue, str
;
; Finds the image value corresponding to given byte value (if possible).
;
COMPILE_OPT IDL2, HIDDEN
COMMON gr_global

value = ivalue
bad = WHERE(ivalue EQ syscol, count)
IF bad[0] GT -1 THEN value[bad] =  !values.F_NAN
IF count EQ N_ELEMENTS(ivalue) THEN RETURN, value

bot = ibot
top = itop
ngrey = N_ELEMENTS(syscol)
FOR ig=0,ngrey-1 DO BEGIN
    idx = WHERE(ivalue GT syscol[ig])
    IF idx[0] NE -1 THEN value[idx] -= 1
ENDFOR
IF bot GT 0 THEN value -= bot
top = top - ngrey - bot

r1 = str.RANGE[0]  &  r2 = str.RANGE[1]
CASE str.TRFUNC OF
    0: scale = top / (str.RANGE[1] - str.RANGE[0])
    1: BEGIN
        asr = ASINH((str.RANGE - str.ZERO)/str.BETA)
        scale = top / (asr[1] - asr[0])
        r1 = str.ZERO
    END
    2:  scale = top / SQRT(r2 - r1)
    3: MESSAGE, /INFORMATIONAL, 'Cannot invert Histogram equalization'
ENDCASE

value /= scale
CASE str.TRFUNC OF
    0: ; No action for linear
    1: BEGIN  ; Asinh
        expval = EXP(value + asr[0])
        value = str.BETA*(expval - 1d0/expval)/2d0
    END
    2: value = value^2
    3: value = !values.F_NAN
ENDCASE
value += r1

RETURN, value

END
;
PRO set_colour_bar, str
;
; Draws an intensity scale bar in the current graphics window.
;
; Input:
;    str: structure describing tab to label
;
COMPILE_OPT IDL2, HIDDEN
COMMON gr_global

CASE str.COLLAB[0] OF
    'mono': GOTO, PSEUDOCOLOUR  ; carry on
    'HSV':  hsv_label, str
    ELSE:   rgb_label, str
ENDCASE
RETURN

PSEUDOCOLOUR:

WSET, str.SCALE_INDEX
ERASE

nxpix = !D.x_vsize - 40  &  nypix = 10

absr  = str.ABSRANGE     &  range  = str.RANGE
wrap  = str.WRAP         &  trfunc = str.TRFUNC
zero  = str.ZERO         &  beta   = str.BETA

absent = !P.background

temp = [[absr], [range]]
r1 = MIN(temp[0,*])
r2 = MAX(temp[1,*])
IF wrap LE 0 THEN r1 = range[0]
IF wrap EQ 0 THEN r2 = range[1]

scale = (r2 - r1) * FINDGEN(nxpix) / (nxpix-1) + r1
scale = scale_image(scale, range, wrap, 0)

x0 = 20               &  x1 = x0 + nxpix - 1
y0 = !D.y_vsize - 13  &  y1 = y0 + nypix - 1

TV, REBIN(scale, nxpix, nypix), x0, y0

                                ; Now the hard part: label the scale
position = CONVERT_COORD([x0,x1], [y0,y1], /DEVICE, /TO_NORMAL)
xold = !X  &  yold = !Y
!Y.window = [position[1,0], position[1,1]]
!X.window = [position[0,0], position[0,1]]

dbyte = itop - ibot
IF trfunc EQ 1 || trfunc EQ 2 THEN BEGIN
                                ; Choose tick values for non-linear scales
    iz2 = scale_image(0.0, range, wrap, trfunc, zero, beta)
    rough = iz2[0] + (LINDGEN(13) - 6)*(dbyte / 5)
                                ; Purge out-of-range values:
    good = WHERE(rough GE ibot AND rough LE itop, ngood)
    rough = rough[good]
    idzero = WHERE(rough EQ iz2[0])
                                ; convert rough values to round numbers
    imrough = invert_scale(rough, str)
    IF idzero GE 0 && idzero LE ngood-1 THEN BEGIN
        imrough[idzero] = 2.0   ; (avoid zero and unities)
        ntoav = idzero < (ngood - 1 - idzero)
        IF ntoav GT 0 THEN BEGIN
                                ; symmetrize range around zero
            negs = -REVERSE(imrough[idzero-ntoav:idzero-1])
            pos  =  imrough[idzero+1:idzero+ntoav]
            avs  =  0.5*(pos + negs)
            imrough[idzero-ntoav:idzero-1] = -REVERSE(avs)
            imrough[idzero+1:idzero+ntoav] = avs
        ENDIF
    ENDIF
    logs    = FLOOR( ALOG10( ABS(imrough) ) )
    imscale = imrough / 10d0^logs
    leads   = 1.0*ROUND(imscale) ; leading digit
                                 ; Use two digits roughly halfway between 1 & 2
    unities = WHERE(ABS(imscale) GT 1.2 AND ABS(imscale) LT 1.7)
    IF unities[0] NE -1 THEN $
      leads[unities] = ROUND(imrough[unities] / 10d0^(logs[unities]-1)) / 10d0

    ticks = leads * 10d0^logs
                                ; Restore zero to actual zero:
    IF idzero GE 0 && idzero LE ngood - 1 THEN ticks[idzero] = 0d0

                                ; Purge again:
    good = WHERE(ticks GE range[0] AND ticks LE range[1], nticks)
    ticks   = ticks[good]
                                ; purge duplicates
    ticks = ticks[UNIQ(ticks)]
    
                                ; Find position on linear (ie. colour
                                ; index) scale:
    iticks = scale_image(ticks, range, wrap, trfunc, zero, beta)
    tnames = STRARR(N_ELEMENTS(ticks))
    aticks = ABS(ticks)
    low = WHERE(aticks GE 1.7d-3 AND aticks LT 10d0, nlow)
    IF nlow GT 0 THEN tnames[low] = $
      STRTRIM(STRING(ticks[low],FORMAT="(G6.2)"),2)
    mid = WHERE(aticks GE 10d0 AND aticks LT 1d4, nmid)
    IF nmid GT 0 THEN tnames[mid] = $
      STRTRIM(STRING(ticks[mid],FORMAT="(F5.0)"),2)
    high = WHERE(aticks LT 1.7d-3 OR aticks GE 1d4,nhi)
    IF nhi GT 0 THEN tnames[high] = $
      STRTRIM(STRING(ticks[high], FORMAT="(E9.1)"),2)
    idzero = WHERE(ticks EQ 0d0)
    IF idzero NE -1 THEN tnames[idzero] = '0.0'
    tnames = '!3'+tnames
ENDIF

                                ; Restore original unit
absmax = 0.1d0*MAX(ABS(absr)) * str.MULT^2
test = numunit(absmax, str.UNIT, OUT_UNIT = ounit, /FORCE)

!X.title = '!3'+ounit
CASE trfunc OF
   0: AXIS, XRANGE = [r1, r2], XSTYLE = 1, XTICKLEN = 0.4
   1: BEGIN                     ; Asinh
      AXIS, XRANGE = [ibot, itop], XSTYLE = 1, XTICKLEN = 0.4, $
            XTICKS = nticks-1, XTICKV = iticks, XTICKNAME = tnames
   END
   2: BEGIN                     ; Sqrt
      AXIS, XRANGE = [ibot, itop], XSTYLE = 1, XTICKLEN = 0.4, $
            XTICKS = nticks-1, XTICKV = iticks, XTICKNAME = tnames
   END
   3: BEGIN                     ; Histogram equalization.
                                ;  God knows what the LUT is, just mark
                                ;  the beginning and end
      ticks = [r1, r2]
      AXIS, XRANGE = [r1, r2], XSTYLE = 1, XTICKLEN = 0.4, $
            XTICKS = 1, XTICKV = ticks
   END
ENDCASE
AXIS, XAXIS = 1, XTICKS = 1, XSTYLE = 0
!X = xold  &  !Y = yold
WSET, str.WINDOW

END
;
PRO rgb_label, str
;
; Draws labels for the R, G, B channels in the space usually used for
; the pseudo-colour scale bar.
;
; Inputs:
;     str:        Structure describing tab
;
COMPILE_OPT IDL2, HIDDEN

WSET, str.SCALE_INDEX
xsize = !D.x_vsize

DEVICE, GET_CURRENT_FONT = oldfont
oldfontcode = !P.font
oldxchar = !D.x_ch_size  &  oldychar = !D.y_ch_size
!P.font = 1
absent = !P.background
DEVICE,/DECOMPOSED
!P.background = 0
ERASE

DEVICE, SET_FONT = "Helvetica Bold", /TT_FONT
DEVICE, SET_CHARACTER_SIZE = [10,16]

strspace = (xsize - 40) / 3
y0 = 0.5*(!D.y_vsize - !D.y_ch_size)
FOR i=0,2 DO XYOUTS, 20 +(i+0.5)*strspace, y0, str.COLLAB[i], $
  COLOR = 255L*(256L^i), /DEVICE, ALIGNMENT = 0.5

                                ; Restore graphics/font state
DEVICE, DECOMPOSED = 0
!P.background = absent    &  !P.font = oldfontcode
DEVICE, SET_CHARACTER_SIZE = [oldxchar, oldychar]
CASE !P.font OF
    -1:                         ; Hershey don't need explicit setting
    0: DEVICE, SET_FONT = oldfont ; Device font
    1: DEVICE, SET_FONT = oldfont, /TT_FONT ; True type
ENDCASE

WSET, str.WINDOW

END
;------------------------------------------------------------------------------
; Utility routines for plotting
;
PRO get_centre, zoom_factor, xhalf, yhalf
;
; Finds effective central pixel on view window
;
COMPILE_OPT IDL2, HIDDEN

xhalf = !D.x_vsize / 2
yhalf = !D.y_vsize / 2

zfac = ROUND(zoom_factor)
IF zfac GT 1 THEN BEGIN ; make xhalf, yhalf land on a pixel centre
    nbigpix = xhalf/zfac
    xhalf = zfac*nbigpix + zfac/2
    nbigpix = yhalf/zoom_factor
    yhalf = zfac*nbigpix + zfac/2
ENDIF

END

FUNCTION im2tv, xpix, ypix, mode
; x_centre, y_centre, xhalf, yhalf, zoom, zfac from mode
; Converts from image pixel to screen pixel displayed by XIMVIEW
;
COMPILE_OPT IDL2, HIDDEN

IF mode.zoom LE 0 THEN BEGIN
    xtv = (xpix - mode.x_centre)/mode.zfac + mode.xhalf
    ytv = (ypix - mode.y_centre)/mode.zfac + mode.yhalf
ENDIF ELSE BEGIN
    xtv = (xpix - mode.x_centre)*mode.zfac + mode.xhalf
    ytv = (ypix - mode.y_centre)*mode.zfac + mode.yhalf
ENDELSE

RETURN, [[xtv], [ytv]]
END

FUNCTION tv2im, xtv, ytv, mode
; x_centre, y_centre, xhalf, yhalf, zoom, zfac
; Converts from screen pixel displayed by XIMVIEW to image pixel
;
COMPILE_OPT IDL2
ON_ERROR, 0

IF mode.zoom LE 0 THEN BEGIN
    xshift = xtv - mode.xhalf  &  yshift = ytv - mode.yhalf
    xpix = xshift*mode.zfac + mode.x_centre
    ypix = yshift*mode.zfac + mode.y_centre
ENDIF ELSE BEGIN
    xpix = (xtv/mode.zfac) - (mode.xhalf/mode.zfac) + mode.x_centre
    ypix = (ytv/mode.zfac) - (mode.yhalf/mode.zfac) + mode.y_centre
ENDELSE

RETURN, [[xpix], [ypix]]

END
;
FUNCTION in_view, xy, x0, y0, x1, y1
;
; Returns true if point (x,y) is inside window, false otherwise
;
COMPILE_OPT IDL2, HIDDEN

IF N_ELEMENTS(x0) EQ 0 THEN BEGIN ; check if on device window
   x0 = 0
   y0 = 0
   x1 = !d.X_VSIZE
   y1 = !d.Y_VSIZE
ENDIF
RETURN, xy[0] GT x0 && xy[0] LT x1 && xy[1] GT y0 && xy[1] LT y1
END
;
PRO marker, mode  
; x_centre, y_centre, xhalf, yhalf, zoom, zfac from mode
; Draws a marker at *image* pixel xpix, ypix
;
COMPILE_OPT IDL2, HIDDEN
line = [0,2,1]
xpix = [mode.xpt, mode.xpt1, mode.xpt2]
ypix = [mode.ypt, mode.ypt1, mode.ypt2]
coord = im2tv(xpix, ypix, mode) 
xcirc = [ 8., 7.608, 6.472, 4.702, 2.472, 0.,-2.472,-4.702,-6.472,-7.608,-8., $
          -7.608,-6.472,-4.702,-2.472, 0., 2.472, 4.702, 6.472, 7.608, 8.]
ycirc = [ 0., 2.472, 4.702, 6.472, 7.608, 8., 7.608, 6.472, 4.702, 2.472, 0.,$
          -2.472,-4.702,-6.472,-7.608,-8.,-7.608,-6.472,-4.702,-2.472, 0.] 
;phase = 2*!dpi*findgen(21)/20d0
; x = 8*COS(phase) & y = 8*SIN(phase)

FOR i=0,2 DO BEGIN
   IF xpix[i] GE 0 THEN BEGIN
      xtv = coord[i,0]  &  ytv = coord[i,1]
      PLOTS, [-12,-4]+xtv, [  0, 0]+ytv, /DEVICE, LINE=line[i]
      PLOTS, [ 12, 4]+xtv, [  0, 0]+ytv, /DEVICE, LINE=line[i]
      PLOTS, [  0, 0]+xtv, [-12,-4]+ytv, /DEVICE, LINE=line[i]
      PLOTS, [  0, 0]+xtv, [ 12, 4]+ytv, /DEVICE, LINE=line[i]
      PLOTS, xcirc + xtv, ycirc + ytv, /DEVICE, line=line[i]
    ENDIF
ENDFOR
END
;
PRO plotcat, catalog,  mode
; x_centre, y_centre, xhalf, yhalf, zoom, zfac from mode
;  Plots data in a catalog structure
;
ON_ERROR, 0
label = catalog.do_label
symb  = catalog.symbol  ; currently ignored
;oldcol = !P.color
;!P.color = catalog.colour

; coordinates of points on each cross
cx = [-5,5,0,0,0]
cy = [0,0,0,5,-5]

coord = im2tv(catalog.xpix, catalog.ypix, mode) 

xtv = REFORM(coord[*,0])  &  ytv = REFORM(coord[*,1])

nsource = N_ELEMENTS(xtv)
FOR is=0,nsource-1 DO PLOTS, xtv[is]+cx, ytv[is]+cy, /DEVICE
IF label THEN XYOUTS, xtv, ytv, catalog.labels, ALIGNMENT=-0.2, /DEVICE

;!P.color = oldcol
END
;
PRO grid_lines, range, inc_in, sex, lat_lon, lines, nl, label, lablen
;
; Sets values of grid lines to plot given increment. Aim is to make
; sure there are lines at round numbers
;
; INPUTS
;     range   coordinate range of visible window
;     inc_in  increment between grid lines
;     sex     if true, interpret as sexagesimal coordinate in degrees
;     lat_lon = 1 for longitude axis, = 2 for latitude, = 0 for neither.
;
; OUTPUTS
;     lines   list of coordinate values at which to plot coordinate
;             lines
;     nl      number of lines
;     label   suitably-formatted strings for line labels
;     lablen  Approximate length of string in device units
;
; TODO: avoid rounding up to 24h / 360d ?
;
IF sex && lat_lon EQ 1 THEN BEGIN
   factor = 15d0
   ifactor = 1 / factor
ENDIF ELSE ifactor = 1
inc   = inc_in   * ifactor
lower = range[0] * ifactor
upper = range[1] * ifactor

IF sex && inc LT 1d0 THEN BEGIN
   inv = 1/inc
   IF inv LT 60d0 THEN inc0 = 1d0 ELSE $
      IF inv LT 3600d0 THEN inc0 = 1/60d0 ELSE inc0 = 1/3600d0
ENDIF ELSE inc0 = 10^(FLOOR(ALOG10(inc)) + 1) 
      
lmin = FLOOR(lower/inc0)*inc0

; tweak start point if coordinates are really degrees
IF lat_lon EQ 1 THEN lmin >= -180d0 ELSE IF lat_lon EQ 2 THEN lmin >= -90d0

nl   = FIX((upper-lmin)/inc) + 2 ; add one line to accomodate tweaks
lines  = lmin - inc + DINDGEN(nl)*inc

; tweak start point to ensure zero coord has a line
IF lower*upper LE 0 && MIN(ABS(lines),imin) NE 0 THEN lines -= lines[imin]

in = WHERE(lines GE lower AND lines LE upper, nl)
lines = lines[in]

IF lat_lon EQ 1 THEN BEGIN ; put in range 0 to 360 or 0 to 24
   neg = WHERE(lines LT 0,nn)
   turn = sex ? 24 : 360
   IF nn GT 0 THEN lines[neg] += turn
   big = WHERE(lines GT turn, nn)
   IF nn GT 0 THEN lines[big] -= turn
ENDIF

grid_labels, lines, inc, lat_lon, sex, label, lablen

IF lat_lon EQ 1 THEN lines *= factor

;HELP, inc_in, inc, inc0, lmin
;PRINT, lines
;print, label
END
;
PRO label_line, tvxy, npt, lablen, label, cs
;
; Labels a grid line provided that it gets close to left edge and has
; a slope of <= 45 deg, or gets close to bottom edge and has a slope of >
; 45 deg
;
; INPUTS
;   tvxy   array of coordinates for line points (device coords)
;   npt    number of points in line
;   lablen Length of label in device coords, used to decide offset
;          from edge
;   label  Text string to label line with
;   cs     character size scaling factor
;
; TODO: offset to centre strings
; TODO: colour
;  
margin = 0.6*lablen
loff = MIN(ABS(tvxy[*,0] - margin),il)
IF loff LT 0.3*lablen THEN BEGIN ; line gets close enough to left edge
   il1 = il LT npt - 2 ? il+1 : il-1
   dxy = tvxy[il1,*] - tvxy[il,*]
   IF ABS(dxy[1]) LE ABS(dxy[0]) && in_view(tvxy[il,*]) THEN BEGIN
      angle = !radeg*ATAN(dxy[1],dxy[0])
      IF ABS(angle) GT 90 THEN angle += 180
      xyouts, tvxy[il,0],tvxy[il,1], label, /DEVICE, ALIGNMENT = 0.5, $
                ORIENTATION = angle, CHARSIZE=cs
      RETURN                    ; don't label the same line twice
   ENDIF
ENDIF

loff = MIN(ABS(tvxy[*,1] - margin),il)
IF loff LT 0.3*lablen THEN BEGIN
   il1 = il LT npt - 2 ? il+1 : il-1
   dxy = tvxy[il1,*] - tvxy[il,*]
   IF ABS(dxy[1]) GE ABS(dxy[0]) && in_view(tvxy[il,*]) THEN BEGIN
      angle = !radeg*ATAN(dxy[1],dxy[0])
      IF ABS(angle) GT 90 THEN angle += 180
      xyouts, tvxy[il,0],tvxy[il,1], label, /DEVICE, ALIGNMENT = 0.5, $
              ORIENTATION = angle, CHARSIZE=cs
   ENDIF
ENDIF

END
;
PRO overlay, mode, astrom, pole_pix, nside
;
; Plots graphics overlays if switched on. NB Call swap_lut first to
; make sure we have the right graphics state.
;
COMPILE_OPT IDL2
ON_ERROR, 0
                                ; First find image pixel range
                                ; corresponding to view window 
nx = !D.X_VSIZE
ny = !D.Y_VSIZE

xy = tv2im([0,nx],[0,ny], mode)
x0 = xy[0,0]
x1 = xy[1,0]
y0 = xy[0,1]
y1 = xy[1,1]

; Plots coordinate grid after display changes
dogrid = mode.GRID_PLOT AND N_ELEMENTS(pole_pix) GT 0
IF dogrid THEN BEGIN 
   npt = nx > ny           ; number of points to plot on each line
   cs = (0.8*npt/512) < 1.25   ; character size for line labels

                                ; Estimate range of longitudes and
                                ; latitudes currently displayed. 

                                ; estimate solid angle in units of 4 pi
   npix = (x1-x0+1)*(y1-y0+1)
   IF nside GT 0 THEN omega = FLOAT(npix) / NSIDE2NPIX(nside) ELSE BEGIN
      cdelt = (*astrom).cdelt
      omega = npix*cdelt[0]*cdelt[1]/360L^2/!dpi
   ENDELSE
   
   IF omega GT 0.4d0 THEN BEGIN ; assume all-sky image
      lrange = [0d0,360d0]
      brange = [-90d0,90d0]
      IF mode.GRID_CALC THEN mode.grid_step = [30d0,30d0]
   ENDIF ELSE BEGIN
                                ; Is there a pole inside the view window?
      npole = in_view(pole_pix[0:1], x0, y0, x1, y1) 
      spole = in_view(pole_pix[2:3], x0, y0, x1, y1) 
      hpx = (*astrom).projection EQ 'HPX'
      IF hpx THEN BEGIN ; check for alternate poles
         FOR ip=1,3 DO BEGIN
            IF ~npole THEN npole = in_view(pole_pix[0:1]+ip*nside, $
                                           x0, y0, x1,y1)
            IF ~spole THEN spole = in_view(pole_pix[2:3]+ip*nside, $
                                           x0, y0, x1,y1)
         ENDFOR
      ENDIF   
      all_lon = (npole || spole) && ~hpx 
      IF all_lon THEN lrange = [0d0,360d0] 

      IF npole && spole && ~hpx THEN brange =  [-90d0, 90d0] ELSE BEGIN
                                ; Get coords along sides
         xpix = LINDGEN(nx)
         ypix = LINDGEN(ny)
         side1 = tv2im(LONARR(ny), ypix, mode)       ; left
         side2 = tv2im(xpix, REPLICATE(ny,nx), mode) ; top
         side3 = tv2im(xpix, LONARR(nx), mode)       ; bottom
         side4 = tv2im(REPLICATE(nx,ny), ypix, mode) ; right
         xtest = [side1[*,0],side2[*,0],side3[*,0],side4[*,0]]
         ytest = [side1[*,1],side2[*,1],side3[*,1],side4[*,1]]
         xy2ad, xtest, ytest, *astrom, ll,bb
         brange = MINMAX(bb)
         IF npole THEN brange[1] =  90d0
         IF spole THEN brange[0] = -90d0
         IF ~all_lon THEN BEGIN
            lrange = MINMAX(ll)
                                ; check for longitude wrap
            dlong = lrange[1] - lrange[0]
            IF dlong GT 180d0 THEN BEGIN ; Probably longitude wrap
               id = WHERE(ll GT 180d0,nn)
               IF nn GT 0 THEN ll[id] -= 360d0
               lrange = MINMAX(ll)
            ENDIF
         ENDIF
                                ; check for gores 
         bad = ~FINITE(bb)
         nn = TOTAL(bad)
         IF nn EQ 2*(nx+ny) THEN RETURN ; we are entirely in a gore,
                                        ; nothing to plot 
      ENDELSE
   ENDELSE                      ; coordinate range set
   
   sex = (*astrom).COORD_SYS EQ 'C'

   IF mode.GRID_STEP[0] LT 0d0 THEN BEGIN
                                   ; choose new grid step
      set_grid_interval, lrange, 1, linc, lmin, nl, SEXAGESIMAL = sex
      set_grid_interval, brange, 2, binc, bmin, nb, SEXAGESIMAL = sex
      mode.GRID_STEP = [linc, binc]
   ENDIF
                                ; choose grid lines on each axis:
   linc = mode.GRID_STEP[0]
   grid_lines, lrange, linc, sex, 1, ll0, nl, l_label, l_len

   binc = mode.GRID_STEP[1]
   grid_lines, brange, binc, sex, 2, bb0, nb, b_label, b_len

   lablen = (l_len > b_len)*cs

                                ; Plot lines of constant latitude:
   ll = lrange[0] + DINDGEN(npt)*(lrange[1]-lrange[0])/(npt-1)
   FOR ib = 0L, nb-1 DO BEGIN
      bb = REPLICATE(bb0[ib], npt)
      ad2xy, ll, bb, *astrom, xx, yy
      get_jumps, xx, yy, bad, njump
      FOR j=0,njump DO BEGIN
         lensegment = bad[j+1] - bad[j]
         IF lensegment GT 1 THEN BEGIN
            u1 = xx[bad[j]:bad[j+1]-1]
            v1 = yy[bad[j]:bad[j+1]-1]
            tvxy = im2tv(u1, v1, mode)
            PLOTS, tvxy[*,0], tvxy[*,1], /DEVICE
            label_line, tvxy, lensegment, lablen, b_label[ib], cs
         ENDIF 
      ENDFOR
   ENDFOR
                                ; plot longitude lines:
   bb = brange[0] + DINDGEN(npt)*(brange[1]-brange[0])/(npt-1)
   lablen = STRLEN(b_label[0])*!d.X_CH_SIZE*cs
   FOR il = 0L, nl-1 DO BEGIN
      ll = REPLICATE(ll0[il], npt)
      ad2xy, ll, bb, *astrom, xx, yy
      get_jumps, xx, yy, bad, njump
      FOR j=0,njump DO BEGIN
         lensegment = bad[j+1] - bad[j]
         IF lensegment GT 1 THEN BEGIN
            u1 = xx[bad[j]:bad[j+1]-1]
            v1 = yy[bad[j]:bad[j+1]-1]
            tvxy = im2tv(u1, v1, mode)
            PLOTS, tvxy[*,0], tvxy[*,1], /DEVICE
            label_line, tvxy, lensegment, lablen, l_label[il], cs
         ENDIF
      ENDFOR
   ENDFOR
ENDIF                           ; End of grid plotting
   
                                ; Mark current point(s)
marker, mode

                                ; Plot any catalogs                           
IF mode.CATPLOT THEN BEGIN
  cats = *mode.CATALOG
  ncat = N_ELEMENTS(cats)
  FOR icat = 0, ncat-1 DO plotcat, *cats[icat], mode
ENDIF

EMPTY

END
;
PRO overview, byteptr, mode, resamp, corner, NOCENTRE = no_centre
;
; Plots overview with FOV of zoomed-in field outlined.
;
; Inputs:
;     byteptr:      pointer to Byte array with image data
;                   (in RGB mode, array of 3 pointers)
;     mode:         structure containing
;       zoom_factor:  what it says
;       xpix, ypix:   Image coordinates at which to centre cursor
;       xhalf, yhalf: Display coords for display centre
;       x_centre, y_centre:  Image coords for display centre
;     no_centre:    Don't move the cursor to the centre pixel.
;
; Outputs:
;     resamp:       Resample factors [x,y]
;     corner:       Coordinates of BLC on display
;
COMPILE_OPT IDL2, HIDDEN

do_centre = ~KEYWORD_SET(no_centre)

IF ~is_gdl() THEN DEVICE, /CURSOR_CROSSHAIR

nchan = N_ELEMENTS(byteptr)
FOR i = 0, nchan-1 DO IF PTR_VALID(byteptr[i]) THEN T = SIZE(*(byteptr[i]))

IF N_ELEMENTS(T) EQ 0 THEN MESSAGE, $
  'Internal error: no valid pointers received'

                                ; Box showing current FOV
xlo  = mode.x_centre - mode.xhalf/mode.zoom_factor
xhi  = xlo  + !D.x_vsize/mode.zoom_factor
ylo  = mode.y_centre - mode.yhalf/mode.zoom_factor
yhi  = ylo  + !D.y_vsize/mode.zoom_factor
xbox = [xlo,xhi,xhi,xlo,xlo]
ybox = [ylo,ylo,yhi,yhi,ylo]

ERASE
resamp = REPLICATE( divup(T[1],!D.x_vsize) > divup(T[2],!D.y_vsize), 2)
IF resamp[0] EQ 1 THEN BEGIN  ; Try zooming in
    resamp = REPLICATE( (!D.x_vsize/T[1]) < (!D.y_vsize/T[2]) , 2)

    xysiz   = T[1:2] * resamp
    xcorner = (!D.x_vsize - xysiz[0] )/2
    ycorner = (!D.y_vsize - xysiz[1] )/2
    xbox    = xbox*resamp[0] + xcorner
    ybox    = ybox*resamp[1] + ycorner
    xtv     = mode.xpix*resamp[0] + xcorner
    ytv     = mode.ypix*resamp[1] + ycorner

    block = xysiz[0]*xysiz[1]
    image = BYTARR(nchan*block)

    FOR i=0,nchan-1 DO IF PTR_VALID(byteptr[i]) THEN image[i*block] = $
      REFORM(REBIN(*byteptr[i], xysiz[0], xysiz[1], /SAMPLE), block)

    resamp = 1.0 / resamp
ENDIF ELSE BEGIN
    xcorner = (!D.x_vsize - T[1] / resamp[0])/2
    ycorner = (!D.y_vsize - T[2] / resamp[1])/2
    xbox    = xbox/resamp[0] + xcorner
    ybox    = ybox/resamp[1] + ycorner
    xtv     = mode.xpix/resamp[0] + xcorner
    ytv     = mode.ypix/resamp[1] + ycorner

    xysiz   = DIVUP(T[1:2], resamp)
    block = xysiz[0]*xysiz[1]
    image = BYTARR(nchan*block)

    FOR i=0,nchan-1 DO IF PTR_VALID(byteptr[i]) THEN image[i*block] = $
      REFORM((*byteptr[i])[0:*:resamp[0], 0:*:resamp[1]], block)
ENDELSE

image = REFORM(image, xysiz[0], xysiz[1], nchan, /OVERWRITE)

IF nchan EQ 3 THEN BEGIN
    DEVICE, /DECOMPOSED
    TV, image, xcorner, ycorner, TRUE = 3
    DEVICE, DECOMPOSED = 0
ENDIF ELSE TV, image, xcorner, ycorner

                                ; Plot box
PLOTS, xbox, ybox, /DEVICE
                                ; Put cursor on current point
IF do_centre THEN TVCRS, xtv, ytv

corner = [xcorner, ycorner]

END
;
PRO update_screen, tabarr, itab, mode, done
;  Redraws one tab
;
;  Inputs:
;     tabarr: Pointer to structure array describing tabs
;     itab:   element of array required
;     mode:   Usual mode structure
;     done:   usual flag.
;
COMPILE_OPT IDL2, HIDDEN
ON_ERROR, 0
str = (*tabarr)[itab]

mono = str.COLLAB[0] EQ 'mono'
set_colour_bar, str

xhalf = mode.XHALF  &  yhalf = mode.YHALF
IF mode.OVERVIEW THEN BEGIN
    done = 1
    tptr = mono ? str.BYTE_PTR : str.RGB
    overview, tptr, mode, resamp, corner
    mode.RESAMP = resamp  &  mode.CORNER = corner
    mode.NEW_VIEW = 2
ENDIF ELSE BEGIN
    ierr = 0
    coord = gscroll(dummy, mode.X_CENTRE, mode.Y_CENTRE, xhalf, yhalf, $
                    mode.ZOOM_FACTOR, ierr, 2, done, DO_WRAP = mode.ROLL)
    IF ierr NE 0 && ierr NE 6 THEN MESSAGE, $
      'GSCROLL error '+STRING(ierr)
ENDELSE

END
;
PRO prep_screen, state, mode, tabarr, oldgraph
;
; Prepares for creation of new tab. Cancels blinking, stashes old
; graphics state, shifts tab arrays to put screen 0 at front,
;
; Inputs:
;     state:    usual state structure
;     mode:     usual mode structure
; Outputs:
;     tabarr:   array of structures for each tab
;     oldgraph: structure with old graphics state from swap_lut
;
COMPILE_OPT IDL2, HIDDEN

WIDGET_CONTROL, state.TABS, SENSITIVE = 0

IF mode.BLINK THEN BEGIN ; Switch off blinking
    WIDGET_CONTROL, mode.BBASE, /DESTROY
    mode.BLINK = 0B
    mode.BWIN = -1  &  mode.BSWIN = -1   & mode.BBASE = -1
    WIDGET_CONTROL, state.TABS, SET_UVALUE = mode
    gscroll_setpar, /HIDDEN
ENDIF

tabarr = *state.TABARR
swap_lut, *state.XIM_GRAPH, tabarr[0], oldgraph

IF state.FIRST THEN RETURN

gscroll_newscreen, 0, tabarr, mode.ZOOM_FACTOR, mode.X_CENTRE, mode.Y_CENTRE, $
  mode.XHALF, mode.YHALF, done, 1B

*state.TABARR = tabarr

END

PRO ximview_resize, top, state
;
; Resizes widget. Called when the program notices
; that the window has been re-sized.
;
COMPILE_OPT IDL2, HIDDEN
COMMON gr_global
ON_ERROR, 2
WIDGET_CONTROL, state.TABS,  GET_UVALUE = mode

dx0 = state.NEWWIDTH  - state.XSIZE ; NB really these are SCR_sizes for unix
dy0 = state.NEWHEIGHT - state.YSIZE ;

IF dx0 EQ 0 && dy0 EQ 0 THEN RETURN

tabarr = state.TABARR
ntab = N_ELEMENTS(*tabarr)
str = (*tabarr)[0]
;
; Increase graphics size by change in widget size (up to a point):
;
geom = WIDGET_INFO(str.DRAW, /GEOMETRY)
xdraw = geom.DRAW_XSIZE  &  ydraw = geom.DRAW_YSIZE

geom = WIDGET_INFO(state.PAD1, /GEOMETRY)
xp1 = geom.XSIZE
geom = WIDGET_INFO(state.PAD2, /GEOMETRY)
xp2 = geom.XSIZE

IF state.NEWWIDTH EQ -1 THEN BEGIN  ; Code to reset to default size
    newx = 512  &  newy = 512
    dx = newx - xdraw
    xtra = 0  &  ytra = 0
ENDIF ELSE BEGIN
    geom = WIDGET_INFO(state.READLAB, /GEOMETRY)
    len  = geom.XSIZE > 256 ; string length or min allowed view region
    dx = MAX([dx0, len-xdraw, -(xp1+xp2)])
    dy = dy0 > (256-ydraw)
    maxdx = state.MAXWIN[0] - xdraw
    maxdy = state.MAXWIN[1] - ydraw
    dx = dx < maxdx
    dy = dy < maxdy

    newx = xdraw+dx  &  newy = ydraw + dy
ENDELSE

dx1 = dx/2 > (-xp1)

is_unix = STRCMP(!version.OS_FAMILY,'UNIX', 4,/ FOLD_CASE)
IF dx NE 0 AND is_unix THEN WIDGET_CONTROL, top, UPDATE = 0

FOR itab = 0,ntab-1 DO WIDGET_CONTROL, (*tabarr)[itab].DRAW,  $
  DRAW_XSIZE = newx, DRAW_YSIZE = newy

IF dx NE 0 THEN BEGIN
    FOR itab = 0,ntab-1 DO $
      WIDGET_CONTROL, (*tabarr)[itab].SCALE, DRAW_XSIZE = xdraw+dx
    WIDGET_CONTROL, state.PAD1, XSIZE = xp1 + dx1
    WIDGET_CONTROL, state.PAD2, XSIZE = xp2 + dx - dx1
ENDIF

IF mode.BLINK THEN BEGIN
    draw = WIDGET_INFO(mode.BBASE, /CHILD)
    scale = WIDGET_INFO(draw, /SIBLING)
    WIDGET_CONTROL, draw, DRAW_XSIZE = newx, DRAW_YSIZE = newy
    IF dx NE 0 THEN WIDGET_CONTROL, scale, DRAW_XSIZE = newx
ENDIF

IF dx NE 0 AND is_unix THEN WIDGET_CONTROL, top, /UPDATE;
; Record revised geometry:
;
geom = WIDGET_INFO(top,/GEOMETRY)
IF is_unix THEN BEGIN
    state.XSIZE     = geom.XSIZE
    state.NEWWIDTH  = geom.XSIZE
                                ; MBAR can be innaccurate if widget
                                ; did not fit on screen.
    state.YSIZE = state.MBAR GT 20 ? geom.YSIZE + state.MBAR : $
                                     geom.SCR_YSIZE
    state.NEWHEIGHT = state.YSIZE
ENDIF ELSE BEGIN
    state.XSIZE    = geom.XSIZE  &  state.YSIZE     = geom.YSIZE
    state.NEWWIDTH = geom.XSIZE  &  state.NEWHEIGHT = geom.YSIZE
ENDELSE

WIDGET_CONTROL, top, SET_UVALUE = state

IF state.FIRST THEN RETURN  ; nothing to re-draw.
;
; Re-draw screens
;
swap_lut, *state.XIM_GRAPH, str, old_graph

get_centre, mode.ZOOM_FACTOR, xhalf, yhalf
mode.XHALF = xhalf  &  mode.YHALF = yhalf

xcent = mode.X_CENTRE  &  ycent = mode.Y_CENTRE
zoomf = mode.ZOOM_FACTOR

null3 = REPLICATE(PTR_NEW(),3)
IF dx NE 0 || mode.OVERVIEW THEN FOR itab = 0,ntab-1 DO BEGIN
    stri = (*tabarr)[itab]
    mono = stri.COLLAB[0] EQ 'mono'

    lutptr = stri.LUT
    IF redraw_req THEN TVLCT, (*lutptr).R, (*lutptr).G, (*lutptr).B
    !P.background = (*lutptr).ABSENT
    !P.color      = (*lutptr).LINE

    IF mode.OVERVIEW THEN BEGIN
        WSET, stri.WINDOW
        tptr = mono ? stri.BYTE_PTR : stri.RGB
                                ; Unlike normal, don't centre cursor
                                ; in middle of window... a bad idea if
                                ; you are currently dragging the edge
                                ; of the window!
        overview, tptr, mode, $
; zoomf, mode.XPIX, mode.YPIX, xhalf, yhalf, xcent, ycent, 
                  resamp, corner, /NOCENTRE
    ENDIF
    IF dx NE 0 THEN set_colour_bar, stri
ENDFOR
;
; Re-draw zoomed-in image:
;
IF mode.OVERVIEW THEN BEGIN
   mode.RESAMP = resamp
   mode.CORNER = corner
   mode.NEW_VIEW = 2
ENDIF ELSE BEGIN
   ierr = 0
   zoom = mode.ZOOM       &  zfac = mode.ZFAC
   mode.XPIX = xcent      &  mode.YPIX = ycent
   
;    tptr = (*tabarr).BYTE_PTR
   coord = gscroll(tptr, xcent, ycent, xhalf, yhalf, $
                   zoomf, ierr, 2, done, DO_WRAP = mode.ROLL)
   IF ierr NE 0 && ierr NE 6 THEN MESSAGE, 'GSCROLL error '+STRING(ierr)
   
                                ; Disable panning if image fits in screen
   bltv = tv2im(0, 0, mode)
; xcent, ycent, xhalf, yhalf, zoom, zfac)
   trtv = tv2im(!D.x_vsize, !D.y_vsize, mode)

   mode.PAN =  MAX(bltv) GT 0 OR MAX(state.IMSIZE[1:2] - trtv) GT 1
   
   mode.DONE = done
                                ; tell overlay to recalculate grid step
                                ; for new window size 
   IF mode.GRID_CALC THEN mode.GRID_STEP = [-1,-1]

                                ; Plot any graphics overlays
    overlay, mode, state.astrom, state.pole_pix, state.nside

                                ; Request another go if loading not finished
    IF done EQ 0 THEN WIDGET_CONTROL, (*tabarr)[0].DRAW, TIMER = 0.5
ENDELSE

WIDGET_CONTROL, state.TABS,  SET_UVALUE = mode
WIDGET_CONTROL, state.TABS, /SENSITIVE

restore_lut, *state.XIM_GRAPH, old_graph

END

;------------------------------------------------------------------------------
; Utility routine for processing metadata, called by make_tabs and
; also ximview_setcoord
;
PRO parse_header, T, header, column, roll, verbose, astrom, is_astrom, $
                  pole_pix, csystem, proj, ounit, beam, title, $
                  nside, ns4, ll_fmt, prec, boxsize
 
;
; Extracts astrometry & unit info from fits header read in via parse_input
;
; Inputs:
;   T:       Size array for the input or data
;   header:  FITS header
;   column:  List of columns/slices extracted from original file
;   roll:    If true, ry to force interpretation as HPX grid
;   verbose: Pixel printout will include timing info
;
; Outputs:
;   astrom:    WCS astrometry structure from header
;   is_astrom: Use the astrom structure
;   csystem:   String descriping coordinate system
;   proj:      HEALPix projection: 'GRID', 'NPOLE', 'SPOLE', or '' (N/A)
;   ounit:     Intensity unit from header, possibly with revised prefix
;   beam:      Structure with BMAJ, BMIN, BPA (all deg), beam_area (sq dg)
;   title:     structure containing strings to use for pixel printout title
;   nside:     HEALPix nside, if any
;   ns4:       4*nside (offset for rolling)
;   ll_fmt:    Format string to write out long & lat
;   prec:      Number of decimals for arcsec readout (RA/Dec only)
;   boxsize:   suggested size of imstats box (12 x 12 beams, roughly)
;
COMPILE_OPT IDL2, HIDDEN

unit = SXPAR(header,'BUNIT', COUNT=got_unit)
IF ~got_unit THEN BEGIN ; Look for TUNITi cards:
    unit = SXPAR(header,'TUNIT*', COUNT=got_unit)
    IF got_unit GE MAX(column) THEN unit = unit[column-1] ELSE $
      IF got_unit GE 1 THEN unit = unit[0] ELSE unit = 'unknown'
ENDIF
unit = STRTRIM(unit,2)

; get beam area
bmaj = SXPAR(header,'BMAJ', COUNT=got_bmaj)
bmin = SXPAR(header,'BMIN', COUNT=got_bmin)
bpa = SXPAR(header,'BPA', COUNT=got_bpa)
; Calculate beam area in square degrees
bmaj = got_bmaj EQ 0 ? !values.F_NAN : FLOAT(bmaj)
bmin = got_bmin EQ 0 ? !values.F_NAN : FLOAT(bmin)
bpa  = got_bpa  EQ 0 ? !values.F_NAN : FLOAT(bpa)
beam_area = bmaj*bmin * !pi/(4*ALOG(2.0)) ; in sq deg
beam = {BEAM, bmaj: bmaj, bmin: bmin, bpa: bpa, beam_area: beam_area}

oldq = !quiet ; Turn off astrolib chat
!quiet = 1
EXTAST, header, astrom, noparam
!quiet = oldq

is_astrom = noparam GT 0
nside = 0B  &  ns4 = 0B  & proj = ''

IF roll THEN BEGIN  ; check format is right:
    error = 0
    error = error && T[1] NE T[2]
    nside = T[1]/5
    npix = NSIDE2NPIX(nside)    ; Is this a valid NSIDE?
    error = error || npix EQ -1 ;
    IF error THEN BEGIN
        MESSAGE, /INFORMATIONAL, 'Roll requested but input map format is wrong'
        roll = 0B
        nside = 0B
    ENDIF ELSE IF ~is_astrom THEN BEGIN ; Make HPX astrom structure
        proj = 'GRID'
        cdelt = [-1,1]*(90d0/nside)
        make_astr, astrom, CRPIX = 0.5d0*(T[1:2] + 1), CRVAL = [0d0, 0d0], $
                   CD=0.5*[[1,-1],[1,1]], DELTA = cdelt, NAXIS = T[1:2], $
                   CTYPE = ['XLON-HPX', 'XLAT-HPX'], PV2 =  [4d0, 3d0]
        is_astrom = 1B
    ENDIF
    ns4 = 4*nside
ENDIF

IF is_astrom THEN BEGIN
                                ; Turn off astrom if WCS type not recognised
    is_astrom = astrom.KNOWN
                                ; Get overall astrometry details

    ii = WHERE( astrom.COORD_SYS EQ ['G', 'E', 'C', 'S', 'T', 'H', 'X'] )
    IF ii EQ -1 THEN ii = 6
    csystem = (['Galactic', 'Ecliptic', 'Equatorial', 'Supergalactic', $
                'Terrestrial', 'Helioecliptic', 'unknown'])[ii]
                                  
    equinox = astrom.EQUINOX
    eqset = FINITE(equinox)
    radesys = astrom.RADECSYS
    
    IF ii LE 3 || ii EQ 5 THEN csystem = csystem+' '+radesys
    CASE STRMID(radesys,0,3) OF
      'FK4': borj = 'B'
      'FK5': borj = 'J'
      'ICR': borj = 'J'
      ELSE:  borj = ''
    ENDCASE
    IF eqset THEN csystem += $
       STRING(borj, equinox, FORMAT = "('  Equinox ',A,F6.1)")

                                ; Find out about specific WCS systems
                                ; (a) for rolling through +/- 180 deg
                                ; (b) for getting pixel areas
    wcs = astrom.PROJECTION
    hpx = wcs EQ 'HPX'

    cylindrical = WHERE(['CYP', 'CEA', 'CAR', 'MER'] EQ wcs)
    equiareal = WHERE(['ZEA', 'CEA', 'SFL', 'GLS', 'PAR', $
                       'MOL', 'AIT', 'COE', 'BON', 'QSC'] EQ wcs)
    equiareal = equiareal GT -1
    IF equiareal THEN pixarea = get_pix_area(astrom.CRPIX, astrom) ELSE BEGIN
                                ; maybe pixel area is independent of
                                ; native longitude:
        lon_eq = WHERE(['TAN','STG', 'ARC', 'ZPN', 'AIR', 'CYP', 'CAR', $
                        'MER', 'COP', 'COD', 'COO'] EQ wcs)
        IF ~lon_eq THEN BEGIN
            CASE wcs OF
                'AZP': lon_eq = astrom.PV2[1] EQ 0
                'SZP': lon_eq = ABS(astrom.PV2[2]) EQ 90d0
                'SIN': lon_eq = astrom.PV2[0] EQ 0d0 AND astrom.PV2[1] EQ 0d0
                ELSE: ; it's a mess.
            ENDCASE
        ENDIF
    ENDELSE

; Cylindrical, equiareal, pixarea and lon_eq are not returned at
; present... for future expansion.

    IF hpx THEN BEGIN ; Find nside. Don't assume that image is not cropped.
        ad2xy, [90d0, 0d0], [0d0, 0d0], astrom, x, y
        nside = ROUND(x[1] - x[0])
        proj = 'GRID'
                      ; Enable roll if we have full sky coverage:
        roll = T[1] EQ 5*nside
    ENDIF ELSE IF wcs EQ 'XPH' THEN BEGIN
        ad2xy, [315d0, 0d0] , [0d0, 90d0], astrom, x, y
        nside = ABS(ROUND(x[1] - x[0]))
        roll = 0B
        proj = astrom.CRVAL[1] GT 0. ? 'NPOLE' : 'SPOLE'
    ENDIF
    IF nside GT 0 THEN BEGIN
;        PRINT, 'Healpix N_side:', nside
        pixarea = !dpi / (3L*nside^2)
        equiareal = 1B
        ns4 = 4*nside
    ENDIF

    IF FINITE(beam_area) THEN BEGIN   ; estimate boxsize for imstats
       IF N_ELEMENTS(pixarea) EQ 0 $ ; get rough estimate of pixel size
       THEN pixarea = ABS(astrom.CDELT[0]*astrom.CDELT[1]) $ ; square degrees
       ELSE pixarea *= (180d0/!dpi)^2 ; convert from radians
       boxsize = 12*ROUND(SQRT(beam_area/pixarea))
                                ; Make sure boxsize is odd, so there
                                ; is a central pixel:
       IF (boxsize/2)*2 EQ boxsize THEN boxsize+= 1
    ENDIF ELSE boxsize = 33

    IF is_astrom THEN BEGIN     ; Find pole pixel coordinates 
       IF roll && hpx THEN pole_pix = [0,2*nside,2*nside,0] - 0.5d0 ELSE BEGIN
          AD2XY, [0d0,0d0],[90d0,-90d0], astrom, xpole, ypole
          lastx = T[1] - 1
          lasty = T[2] - 1
          poles = xpole GT 0 AND xpole LT lastx AND ypole GT 0 AND ypole LT lasty
          
          pole_pix = [xpole[0],ypole[0],xpole[1],ypole[1]]
          IF ~poles[0] THEN pole_pix[0:1] = !values.D_NAN ; north pole off map
          IF ~poles[1] THEN pole_pix[2:3] = !values.D_NAN ; south pole off map
       ENDELSE
    ENDIF ELSE pole_pix = REPLICATE(!values.D_NAN,4)

                                 ; Work out precision needed for lat/long
    cdelt = astrom.CDELT
    IF astrom.REVERSE THEN cdelt = REVERSE(cdelt)
    IF astrom.COORD_SYS EQ 'C' THEN BEGIN
        prec = CEIL(-ALOG10(cdelt[1]*3600.0))
        len = prec GT 0 ? 23+2*prec : 22
        nchar = INTARR(2)
        nchar[0] = len / 2
        nchar[1] = len - nchar[0] 
        ll_fmt = ''
    ENDIF ELSE BEGIN 
        prec = CEIL(-ALOG10(ABS(cdelt))) + 1
        nchar = prec + 6
        ns1 = 9-nchar[0] 
        ns2 = 9-nchar[1]
        mid  = ns1 GT 0 ? ",'"+STRJOIN(REPLICATE(' ',ns1),'')+"'" : ''
        tail = ns2 GT 0 ? ",'"+STRJOIN(REPLICATE(' ',ns2),'')+"'" : ''
        nchar = STRTRIM(STRING(nchar),2)
        ndp   = STRTRIM(STRING(prec),2)
        ll_fmt = '(1X,'+STRJOIN('F'+nchar+'.'+ndp,mid+',')+tail+')'
    ENDELSE
    prec = MAX(prec)
ENDIF ELSE BEGIN
    hpx = 0B
    astrom = 0B
    csystem = ''
    boxsize = 33
ENDELSE


ounit    = STRARR(N_ELEMENTS(column))
ounit[*] = unit

; Set title for intensity readout: units or "Brightness"
mid = form_unit(ounit[0])
tail = '   Zoom  '
IF verbose THEN tail = tail +'  dt (sec)'

thead = 'X pix Y pix '
IF nside GT 0 THEN thead = thead + ' NEST pix RING pix' ; HEALPix grid.
IF is_astrom THEN BEGIN
  IF astrom.COORD_SYS EQ 'C' THEN BEGIN
      thead += padstring('   RA  ',(nchar[0]-2) > 7,'R')
      thead += eqset ? STRING(borj,equinox,FORMAT="(A1,I4)") : '     '
      thead += padstring('  Dec  ',(nchar[1]-3) > 7,'R')
  ENDIF ELSE BEGIN
      thead += padstring(' longitude',1+(nchar[0]> 9),'R')
      thead += padstring(' latitude',nchar[1] > 9,'R')
  ENDELSE
ENDIF
  
title = {head: thead, unit: mid, tail: tail}

END
;------------------------------------------------------------------------------
; Miscellaneous utility routines
;
FUNCTION get_log_name
; Finds a unique name for the log file.
;
COMPILE_OPT IDL2, HIDDEN

existing = FILE_SEARCH('ximview_*.log')
IF existing[0] EQ '' THEN logfile = 'ximview_1.log' ELSE BEGIN
                                ; Find maximum index among existing files:
    existing = STRMID(existing, 8) ; removes 'ximview_' prefix
                                ; now remove '.log' suffix:
    nums    = STRSPLIT( STRJOIN(existing), '.log', /REGEX, /EXTRACT)
    newnum  = STRTRIM( STRING( MAX(FIX(nums)) + 1 ), 2)
    logfile = 'ximview_'+newnum+'.log'
ENDELSE

RETURN, logfile
END

FUNCTION form_unit, newunit
; Formats unit into 12-letter string that (ideally) will be centred
; when two spaces are added to the end, unless "unknown"
;
COMPILE_OPT IDL2, HIDDEN

ounit = STRTRIM(newunit,2)
IF STRCMP(ounit, 'unknown', 6, /FOLD_CASE) THEN mid = '  Brightness' $
ELSE BEGIN
    lunit = STRLEN(ounit)
    IF lunit GE 10 THEN mid = '  '+STRMID(ounit,0,10) $
                   ELSE BEGIN
        mid = padstring(ounit,14,'C')
        mid = STRMID(mid,0,12)
    ENDELSE
ENDELSE 

RETURN, mid
END

PRO pix_print, state, log, start
;
; Prints pixel details to readout and possibly terminal and logfile
;
; Inputs:
;         state:       structures with global parameters
;         log:         Write out to terminal and logfile
;         start:       start time of this cycle, for diagnostic
;                      printing.
;
COMPILE_OPT IDL2, HIDDEN
ON_ERROR, 2

r2d = 180d / !dpi
WIDGET_CONTROL, state.TABS, GET_UVALUE = mode
xpix        = mode.XPIX
ypix        = mode.YPIX
zoom_factor = mode.ZOOM_FACTOR

current = WIDGET_INFO(state.TABS, /TAB_CURRENT)
ntab    = WIDGET_INFO(state.TABS, /TAB_NUMBER)
blinking = mode.BLINK && current EQ ntab-1
iscreen = blinking ? WHERE((*state.TABARR).SCREEN EQ (*state.BLINK_SEQ)[0]) $
                   : 0

nside = state.NSIDE
                                ;  Monochrome or RGB?
mono = (*state.TABARR)[iscreen].COLLAB[0] EQ 'mono'
IF mono THEN BEGIN
    im_ptr = (*state.TABARR)[iscreen].IM_PTR
    mult   = (*state.TABARR)[iscreen].MULT
    fmt = (*state.TABARR)[iscreen].FORMAT
    IF ~im_ptr THEN BEGIN
        WIDGET_CONTROL, state.LABEL, GET_UVALUE = image, /NO_COPY
        imval = image[xpix,ypix]
        WIDGET_CONTROL, state.LABEL, SET_UVALUE = image, /NO_COPY
    ENDIF ELSE imval = (*im_ptr)[xpix,ypix]
    imval = imval * mult
ENDIF ELSE BEGIN
    rgbptr = (*state.TABARR)[iscreen].RGB
    rval = 0 & gval = 0 & bval = 0
    IF PTR_VALID(rgbptr[0]) THEN rval = (*rgbptr[0])[xpix,ypix]
    IF PTR_VALID(rgbptr[1]) THEN gval = (*rgbptr[1])[xpix,ypix]
    IF PTR_VALID(rgbptr[2]) THEN bval = (*rgbptr[2])[xpix,ypix]
ENDELSE

IF state.IS_ASTROM THEN BEGIN ; Get long and lat if available from astrom.
    XY2AD, xpix, ypix, *state.ASTROM, ll, bb
    goodpix = FINITE(ll) && FINITE(bb)
ENDIF ELSE goodpix = 1B


IF nside GT 0 THEN BEGIN        ; Get HEALPix pixel numbers
                                ; This whole section can probably be replaced
                                ; by calls to ang2pix_nest & ang2pix_ring.
                                ;
    is_grid = state.PROJ EQ 'GRID'
    np2 = is_grid ? 2.5 : 2.0
    blc = np2*nside + 0.5 - (*state.ASTROM).CRPIX

    nfull = is_grid ? 5*nside : 4*nside

    IF is_grid THEN BEGIN
        llshift = ROUND((*state.ASTROM).CRVAL[0] / 90d0)
        blc = blc - llshift*nside
    ENDIF

    xshift = (xpix + blc[0])
    yshift = (ypix + blc[1])

    IF is_grid THEN BEGIN
        IF xshift GE nfull || yshift GE nfull THEN BEGIN
            xshift -= state.NS4
            yshift -= state.NS4
        ENDIF ELSE IF xshift LT 0 || yshift LT 0 THEN BEGIN
            xshift += state.NS4
            yshift += state.NS4
        ENDIF
    ENDIF
    idx = grid2hp_index(nside, ROUND(xshift), ROUND(yshift), state.PROJ, $
                        /BOTH, R_OFF = *state.RING0, /SILENT)
    rix = idx[0]  &  nix = idx[1]
    idx = 0
    ngood = nix[0] GE 0
    rgood = rix[0] GE 0
    IF state.PROJ EQ 'GRID' && (goodpix NE ngood || goodpix NE rgood) THEN BEGIN
        MESSAGE, /INFORMATIONAL, 'Astrometry mismatch'
        MESSAGE, /INFORMATIONAL, STRING(xpix, ypix, ll, bb, nix, rix, $
                             FORMAT = "('Params:',2I6,2F12.6,2I10)")
    ENDIF
ENDIF ELSE nix = -1


; Now construct the output string:

pos = STRING(xpix, ypix, FORMAT = "(I5,1X,I5,' ')")

IF nside GT 0 THEN pos += goodpix ? STRING(nix, rix, FORMAT = "(2I9)") $
                                  : '                  '

IF state.IS_ASTROM THEN BEGIN
   IF goodpix THEN BEGIN
      pos +=  (*state.ASTROM).COORD_SYS EQ 'C' ?  $
              adstring(ll, bb, state.PREC)     :  $
              STRING(ll, bb, FORMAT = state.LL_FMT)
   ENDIF ELSE pos += STRING(' No sky pixel', FORMAT = "(A,6(' '))")
ENDIF

IF mono THEN pos += goodpix ? STRING(imval, FORMAT = fmt) $
                            : '            ' $
        ELSE pos += STRING(rval,gval,bval,FORMAT="(1X,2(I3,','),I3)")
                           
pos += STRING(zoom_factor, FORMAT = "(1X,F6.3,' ')")

dt = SYSTIME(1) - start

IF state.VERBOSE THEN BEGIN
   pos += STRING(dt, FORMAT = "(F10.5)")
;  nest2ring, nside, nix, rix2
;  pos +=STRING(rix-rix2, FORMAT = "(I10)") 
ENDIF

IF log THEN BEGIN
    PRINT, pos
    PRINTF, state.LOGLUN, pos
ENDIF

WIDGET_CONTROL, state.READOUT, SET_VALUE = pos

END

;--------------------------------------------------------------------------
;
; Cleanup routines:
;
PRO ximview_cleanup, id
; Ensures widget dies tidily
;
COMPILE_OPT IDL2, HIDDEN
ON_ERROR, 1

WIDGET_CONTROL, id, GET_UVALUE = state

ximview_tidy, state

PRINT, ''
PRINT, 'Ximview finished'

END
PRO ximview_tidy, state
; Does the real clean up. Called by ximview_cleanup and also by the
; catch block in the main program.
;
COMPILE_OPT IDL2, HIDDEN
ON_ERROR, 2
IF N_ELEMENTS(state) EQ 0 THEN RETURN

IF state.LOGLUN GT 0 THEN FREE_LUN, state.LOGLUN, /FORCE
;
; Free pointers:
;
tabarr = *state.TABARR
ntab = N_ELEMENTS(tabarr)
FOR itab = 0, ntab-1 DO BEGIN
    str = tabarr[itab]
    IF PTR_VALID(str.BYTE_PTR) THEN PTR_FREE, str.BYTE_PTR
    IF PTR_VALID(str.LUT)      THEN PTR_FREE, str.LUT
    FOR i=0,2 DO IF PTR_VALID(str.RGB[i]) THEN PTR_FREE, str.RGB[i]
    IF str.TEMPORARY AND PTR_VALID(str.IM_PTR) THEN PTR_FREE, str.IM_PTR
ENDFOR

IF PTR_VALID(state.BLINK_SEQ) THEN PTR_FREE, state.BLINK_SEQ
IF PTR_VALID(state.RING0)     THEN PTR_FREE, state.RING0
IF PTR_VALID(state.ASTROM)    THEN PTR_FREE, state.ASTROM
IF PTR_VALID(state.TABARR)    THEN PTR_FREE, state.TABARR
IF PTR_VALID(state.XIM_GRAPH) THEN PTR_FREE, state.XIM_GRAPH
PTR_FREE, state.TAB_TEMPLATE
IF PTR_VALID(state.FILES) THEN BEGIN
    IF N_ELEMENTS(*state.FILES) GT 0 THEN BEGIN
        file_str = *state.FILES
        PTR_FREE, file_str.HEADER
    ENDIF
    PTR_FREE, state.FILES
ENDIF
IF PTR_VALID(state.control) THEN BEGIN
  IF N_ELEMENTS(*state.control) GT 0 THEN BEGIN
     nt = N_TAGS(*state.control)
     FOR itag=0,nt-1 DO $ ; find pointers and free them
        IF SIZE((*state.control).(itag),/TYPE) EQ 10 THEN $
           PTR_FREE, (*state.control).(itag)
  ENDIF
  PTR_FREE, state.CONTROL
ENDIF

gscroll_set = ~state.FIRST
IF gscroll_set THEN gscroll_tidy

END
;
PRO clear_state, state
;
; Resets state structure in preparation for loading a new image with
; no astrometry match to previous (/OVERWRITE).
;
state.LASTTAB = 1L
state.TITLE = {head: ' ', unit: ' ', tail: ' '}
state.FIRST = 1B
state.PROJ = ''
state.ROLL = 0B
state.IN_COUNT = 0S
PTR_FREE, state.FILES
state.FILES = PTR_NEW(/ALLOCATE_HEAP)
*state.BLINK_SEQ = -1
state.NSIDE = 0L
state.NS4 = 0L
PTR_FREE, state.RING0 
state.RING0 = PTR_NEW(/ALLOCATE_HEAP)
state.IS_ASTROM = 0B
state.pole_pix = !values.D_NAN
state.LL_FMT =  "(' ',2F9.3)"
state.PREC = 0S
PTR_FREE, state.ASTROM
state.csystem = ''
;*state.TABARR = [*state.TAB_TEMPLATE]
END
;--------------------------------------------------------------------------
;
; Event handlers for the menubar options
;
; NB: the following large event handlers are in separate files:
; ximview_fitsload, ximview_rdimage, ximview_scale, ximview_blink, 
; ximview_rgb, ximview_pol, ximview_maxfit, ximview_imstats.
;
PRO ximview_newlog, event
; Starts a new log file
;
COMPILE_OPT IDL2, HIDDEN
ON_ERROR, 1

WIDGET_CONTROL, event.ID, GET_VALUE = label
WIDGET_CONTROL, event.TOP,  GET_UVALUE = state

CASE label OF
    'Overwrite old file': logfile = 'ximview.log'
    'New sequence #': logfile = get_log_name()
    'Named...': BEGIN
        logfile = get_user_datum(event, 'File name', 30)
        IF STRCMP(logfile, 'CANCEL', 6, /FOLD_CASE) THEN RETURN
    END
    ELSE: MESSAGE, /INFORMATIONAL, 'Option ' + label + ' not yet available.'
ENDCASE

FREE_LUN, state.LOGLUN
OPENW, loglun, logfile, /GET_LUN

; Write header lines to logfile here..
WIDGET_CONTROL, state.LABEL,   GET_VALUE = name
WIDGET_CONTROL, state.READLAB, GET_VALUE = title
nside = state.NSIDE
PRINTF,loglun,  state.VERSION, SYSTIME(), name, FORMAT = $
  "('XIMVIEW Version ',A,' Restarted log at ',A,//'Dataset: ',A)"
line = ['']
IF state.IS_ASTROM THEN line = [line, 'Coordinate system: '+state.CSYSTEM, '']
IF nside GT 0 THEN line = [line, 'Seems to be HEALPix with N_side' + $
                                 STRING(nside,FORMAT="(I5)"), '']
line = [line,title]
PRINTF, loglun, line, FORMAT="(A)"

state.LOGLUN = loglun
WIDGET_CONTROL, event.TOP, SET_UVALUE = state

END

PRO ximview_deltab, event
; Deletes a tab
;
COMPILE_OPT IDL2, HIDDEN
COMMON gr_global
ON_ERROR, 1

label = WIDGET_INFO(event.ID, /NAME)
IF label EQ 'BUTTON' THEN WIDGET_CONTROL, event.ID,  GET_VALUE = label

WIDGET_CONTROL, event.TOP, GET_UVALUE = state
tabarr = *state.TABARR

WIDGET_CONTROL, state.TABS, SENSITIVE = 0

current = WIDGET_INFO(state.TABS, /TAB_CURRENT)
ntab = N_ELEMENTS(tabarr)
funny = current GE ntab  ; non-standard tab (blink)

swap_lut, *state.XIM_GRAPH, tabarr[0], oldgraph

IF ntab LT 1 THEN BEGIN
    ok = DIALOG_MESSAGE("No tabs to delete", TITLE='XIMVIEW', $
                        DIALOG_PARENT = event.TOP)
    RETURN
ENDIF

tabs = get_tab_uvals(state.TABS, tabid)
CASE label OF
    'Current tab': IF funny THEN BEGIN
        deadtab = WIDGET_INFO(event.TOP, FIND_BY_UNAME = tabs[current])
        id = -1
        index = current
    ENDIF ELSE BEGIN
        id = 0
        deadtab = tabarr[0].BASE
        index   = tabarr[0].SCREEN
    ENDELSE
    'Specify': BEGIN
        tabs = [tabs,'Cancel']
        index = get_user_item(event, 'Select tab to delete', tabs)
        IF index EQ ntab THEN RETURN ELSE deadtab = tabid[index]
        id = WHERE(tabarr.SCREEN EQ index)
    END
    'BASE' : BEGIN ; Directed via PAD2
        index = event.TAB
        deadtab = tabid[index]
        id = WHERE(tabarr.SCREEN EQ index)
    END
    ELSE: MESSAGE, /INFORMATIONAL, 'Option ' + label + ' not yet available.'
ENDCASE

IF ntab GT 1 THEN WIDGET_CONTROL, deadtab, /DESTROY

WIDGET_CONTROL, state.TABS, GET_UVALUE = mode
current = WIDGET_INFO(state.TABS, /TAB_CURRENT)

; Update screen numbers for tabs which fall after the one we just deleted:
idx = WHERE(tabarr.SCREEN GT index)
IF idx[0] NE -1 THEN tabarr[idx].SCREEN -= 1

IF id GE 0 THEN BEGIN
   str = tabarr[id]
   IF ntab GT 1 THEN tabarr[id].SCREEN = -1

; delete old pixmaps, structures and pointers, and/or bulk data;
; update new current screen:
   IF PTR_VALID(str.LUT) THEN BEGIN
      lut = *str.LUT
      PTR_FREE, str.LUT
   ENDIF ELSE lut = *state.XIM_GRAPH
   IF PTR_VALID(str.BYTE_PTR) THEN BEGIN ; not RGB tab either
      IF str.TEMPORARY THEN IF PTR_VALID(str.IM_PTR) THEN $
         PTR_FREE, str.IM_PTR ELSE BEGIN
         WIDGET_CONTROL, state.LABEL, GET_UVALUE = image, /NO_COPY
         image = 0
      ENDELSE
      PTR_FREE, str.BYTE_PTR
   ENDIF ELSE IF str.TEMPORARY THEN FOR i=0,2 DO IF PTR_VALID(str.RGB[i]) $
         THEN PTR_FREE, str.RGB[i]
ENDIF
; Delete low-level pixmap & window structures, update view, and also
; remove structure describing deleted tab fro tabarr and shift other elements
; in tabarr if necessary: 

IF ntab GT 1 THEN BEGIN
   gscroll_newscreen, current, tabarr, mode.ZOOM_FACTOR, mode.X_CENTRE, $
                    mode.Y_CENTRE, mode.XHALF, mode.YHALF, done, mode.OVERVIEW

   IF redraw_REQ THEN DEVICE, DECOMPOSED = 0B $
                 ELSE DEVICE, DECOMPOSED = tabarr[0].DECOMPOSED

   IF ~funny THEN ntab = ntab - 1
ENDIF ELSE BEGIN
   ; shut down scrolling
   window =  tabarr[id].WINDOW
   gscroll_tidy
   WSET, window
   ERASE
   ; replace LUT structure (cleared by gscroll_tidy)
   tabarr[id].LUT = PTR_NEW(lut)
   done = 1
   ntab = 0
   ; clear list of files on Headers pull-down
   buttons = WIDGET_INFO(state.DATASETS, /ALL_CHILDREN)
   nbut = N_ELEMENTS(buttons)
   FOR ib=0L,nbut-1 DO WIDGET_CONTROL, buttons[ib], /DESTROY
   ; Clear mark
   ; reset state structure
   clear_state, state
   WIDGET_CONTROL, event.TOP, SET_UVALUE = state
   ; Turn off overlays
   mode.XPT = -1
   mode.YPT = -1
   mode.XPT1 = -1
   mode.YPT1 = -1
   mode.XPT2 = -1
   mode.YPT2 = -1
   IF PTR_VALID(mode.CATALOG) THEN BEGIN
      PTR_FREE, mode.CATALOG
      mode.CATALOG = PTR_NEW(/ALL)
   ENDIF
   mode.CATPLOT = 0B
   mode.GRID_PLOT = 0B
ENDELSE

*state.TABARR = tabarr

; disable blinking if now only one tab
IF ntab LT 2 THEN BEGIN
    WIDGET_CONTROL, state.BLINK,  SENSITIVE = 0
    WIDGET_CONTROL, state.FRAMES, SENSITIVE = ntab
ENDIF

mode.DONE = done
WIDGET_CONTROL, state.TABS, SET_UVALUE = mode
                                ; Request another go if loading not finished
IF done EQ 0 THEN WIDGET_CONTROL, tabarr[0].DRAW, TIMER = 0.5

restore_lut, dummy, oldgraph

IF ntab GT 0 THEN WIDGET_CONTROL, state.TABS, /SENSITIVE

END

PRO ximview_2png, event
; Dumps current screen as a PNG file
;
COMPILE_OPT IDL2, HIDDEN
ON_ERROR, 1

WIDGET_CONTROL, event.ID, GET_VALUE = label
WIDGET_CONTROL, event.TOP,  GET_UVALUE = state

CASE label OF
    'Image and scale bar': both = 1B 
    'Image only': both = 0B
    ELSE: MESSAGE, /INFORMATIONAL, 'Option ' + label + ' not yet available.'
ENDCASE

name = GET_USER_DATUM(event, 'Name for PNG file', 30, 'ximview.png')
IF name EQ 'Cancel' THEN RETURN

WIDGET_CONTROL, event.TOP,  GET_UVALUE = state
str = (*state.TABARR)[0]

swap_lut, *state.XIM_GRAPH, str, oldgraph

IF ~state.FOCUS THEN WIDGET_CONTROL, str.DRAW, /INPUT_FOCUS

image = TVRD(TRUE = 1)
WSET, str.SCALE_INDEX
scalebar = TVRD(TRUE = 1)

IF both THEN BEGIN
   WSET, str.WINDOW
   image = [[[scalebar]],[[image]]]
ENDIF

WRITE_PNG, name, image

restore_lut, dummy, oldgraph

END

PRO ximview_reset, event
; Restores possibly damaged program state: draw panel to sensitive,
; pan mode to true.
;
COMPILE_OPT IDL2, HIDDEN

WIDGET_CONTROL, event.TOP,  GET_UVALUE = state
WIDGET_CONTROL, state.TABS, GET_UVALUE = mode
mode.PAN = 1B
WIDGET_CONTROL, state.TABS, SET_UVALUE = mode
WIDGET_CONTROL, state.TABS, /SENSITIVE

END

PRO ximview_dump, event
;
; Prints debug info
;
COMPILE_OPT IDL2, HIDDEN
WIDGET_CONTROL, event.TOP,  GET_UVALUE = state
WIDGET_CONTROL, event.ID, GET_VALUE = option
CASE option OF
   'State info': BEGIN
      HELP, state, /STRUCTURE
      PRINT, state
      ntab = N_ELEMENTS(*state.TABARR)
      PRINT, 'Length of TABARR array is', ntab
   END
   'Tab info': BEGIN 
      str = (*state.TABARR)[0]
      PRINT, 'Tab structure:'
      HELP, str, /STRUCTURE
      PRINT, str
   END
   'Mode info': BEGIN
      WIDGET_CONTROL, state.TABS, GET_UVALUE=mode
      PRINT, 'mode structure'
      HELP, mode, /STRUCTURE
      PRINT, mode
   END
   'LUT info': BEGIN
      str = (*state.TABARR)[0]
      PRINT, 'LUT structure:'
      HELP, *str.LUT, /STRUCTURE
      PRINT, 'R:', (*str.LUT).R
      PRINT, 'G:', (*str.LUT).G
      PRINT, 'B:', (*str.LUT).B
   END
   'Astrometry': BEGIN
      PRINT, 'Astrometry structure:'
      HELP, *state.astrom, /STRUCTURE
      PRINT, *state.astrom
   END
ENDCASE
END

PRO ximview_exit, event
COMPILE_OPT IDL2, HIDDEN

WIDGET_CONTROL, event.TOP, /DESTROY

END

PRO ximview_clear_mark, event
; Unsets the marked point and re-draws screen to remove it if in view.
;
COMPILE_OPT IDL2, HIDDEN
ON_ERROR, 1

WIDGET_CONTROL, event.TOP,  GET_UVALUE = state

IF state.FIRST THEN RETURN ; havn't started yet.
WIDGET_CONTROL, state.TABS, GET_UVALUE = mode
WIDGET_CONTROL, event.ID, GET_VALUE = option
CASE option OF
   'Last': BEGIN
      mode.XPT = mode.XPT1  & mode.YPT = mode.YPT1
      mode.XPT1 = mode.XPT2 & mode.YPT1 = mode.YPT2
      mode.XPT2 = -1 & mode.YPT2 = -1
   END
   'All': BEGIN
      mode.XPT = -1 & mode.YPT = -1
      mode.XPT1 = -1 & mode.YPT1 = -1
      mode.XPT2 = -1 & mode.YPT2 = -1
   END
   ELSE: MESSAGE, /INFORMATIONAL, 'Option '+option+' not recognised'
ENDCASE
WIDGET_CONTROL, state.TABS, SET_UVALUE = mode

swap_lut, *state.XIM_GRAPH, (*state.TABARR)[0], oldgraph
gscroll_refresh, mode.XHALF, mode.YHALF
restore_lut, dummy, oldgraph

END
;
PRO ximview_newlut, event
; Sets new colour table
;
COMPILE_OPT IDL2, HIDDEN
COMMON gr_global
ON_ERROR, 1

WIDGET_CONTROL, event.ID,  GET_VALUE = label
WIDGET_CONTROL, event.TOP, GET_UVALUE = state

IF state.FIRST THEN RETURN
IF ~colmap THEN BEGIN
    MESSAGE, /INFORMATIONAL, $
      'Colour table manipulation is not possible on this system'
    RETURN
ENDIF

str = (*state.TABARR)[0]
schemes = colour_schemes(-1)
index = (WHERE(STRCMP(label, schemes, /FOLD_CASE)))[0]
IF index NE str.COLTAB THEN BEGIN
    WIDGET_CONTROL, state.TABS, SENSITIVE = 0
    done = 1

    swap_lut, *state.XIM_GRAPH, dummy, old_graph

    ximview_lut, index, str.IZERO,decomp
    str.COLTAB     = index
    str.DECOMPOSED = decomp
    TVLCT, r, g, b, /GET
    *str.LUT = {r:r, g:g, b:b, line: !P.color, absent: !P.background}
                                  ; Check for multiple uses of the old LUT
    (*state.TABARR)[0] = str
    gscroll_setpar, /BLANK ; sets the blank pixmap to the new !P.background

    WIDGET_CONTROL, state.TABS, GET_UVALUE = mode

    IF state.GLOBAL_COLOUR THEN BEGIN
        ntab = N_ELEMENTS(*state.TABARR)
        (*state.TABARR).COLTAB     = index
        (*state.TABARR).DECOMPOSED = decomp
        FOR i=0,ntab-1 DO *(*state.TABARR)[i].LUT = *str.LUT

        IF redraw_req THEN FOR i = 0,ntab-1 DO $
          update_screen, state.TABARR, i, mode, done

    ENDIF ELSE IF redraw_req THEN $
      update_screen, state.TABARR, 0, mode, done

    IF done EQ 0 THEN WIDGET_CONTROL, str.DRAW, TIMER=0.5
    mode.DONE = done

    restore_lut, *state.XIM_GRAPH, old_graph
    
    WIDGET_CONTROL, state.TABS, SET_UVALUE = mode
    WIDGET_CONTROL, state.TABS, /SENSITIVE
ENDIF

END
;
PRO ximview_goto, event
;
; Sets view centre to a specified position
;
COMPILE_OPT IDL2, HIDDEN
ON_ERROR, 1

start = SYSTIME(1)

WIDGET_CONTROL, event.TOP, GET_UVALUE = state
WIDGET_CONTROL, state.TABS, GET_UVALUE = mode
WIDGET_CONTROL, state.TABS, SENSITIVE = 0

swap_lut, *state.XIM_GRAPH, (*state.TABARR)[0], old_graph

pixvals = STRCOMPRESS(STRING(mode.X_CENTRE, mode.Y_CENTRE))
bnames = 'image pixel'
coindex = [0]
nside = state.NSIDE
IF nside NE 0 THEN BEGIN
    bnames = bnames+'|HP ring pixel|HP nest pixel'
    coindex = [coindex, 1, 2]
ENDIF
IF state.IS_ASTROM THEN  BEGIN  ; || state.NSIDE NE 0
   bnames = bnames+'|(long\, lat)'
   coindex = [coindex, 3]
ENDIF
buttons = N_ELEMENTS(coindex) EQ 1 ? '0, BUTTON,'+bnames+', ROW, TAG=coord' : $
          '0, BUTTON,' +bnames+ ', ROW, EXCLUSIVE, SET_VALUE=0, TAG=coord'
          
desc = ['0, TEXT, ' + pixvals + $
        ', LABEL_LEFT=Enter position:, WIDTH=20, TAG=text', buttons, $
        '1, BASE,, ROW', '0, BUTTON, Done, QUIT', $
        '2, BUTTON, Cancel, QUIT, TAG=cancel']
form = CW_FORM(desc, /COLUMN, TITLE = 'Set centre of view')

IF form.CANCEL THEN RETURN

xhalf = mode.XHALF  &  yhalf = mode.YHALF
TVCRS, xhalf, yhalf
pix = STRSPLIT(form.TEXT,', ',/EXTRACT, COUNT = count)
r2d = 180d0 / !dpi
IF coindex[form.COORD] LT 3 THEN pix = DOUBLE(pix)
CASE coindex[form.COORD] OF
   0: BEGIN                     ; pixel coordinates
      IF count NE 2 THEN BEGIN
         ok = DIALOG_MESSAGE('Enter two numbers for pixel coordinates', $
                             DIALOG_PARENT=event.TOP, /ERROR)
         RETURN
      ENDIF
      xpix = pix[0]  &  ypix = pix[1]
   END
   1: BEGIN                     ; HEALPix ring pixel number
      IF count NE 1 THEN BEGIN
         ok = DIALOG_MESSAGE('Enter only one number for HP pixel', $
                             DIALOG_PARENT=event.TOP, /ERROR)
         RETURN
      ENDIF
      pix = ROUND(pix[0])
      order = 'RING'
      PIX2ANG_RING, nside, pix, theta, phi
      ll = phi[0]*r2d  &  bb = 90d - theta[0]*r2d
   END
   2: BEGIN                     ; HEALPix nest pixel number
      IF count NE 1 THEN BEGIN
         ok = DIALOG_MESSAGE('Enter only one number for HP pixel', $
                             DIALOG_PARENT=event.TOP, /ERROR)
         RETURN
      ENDIF
      pix = ROUND(pix[0])
      order = 'NESTED'
      PIX2ANG_NEST, nside, pix, theta, phi
      ll = phi[0]*r2d  &  bb = 90d - theta[0]*r2d
   END
   3: BEGIN                     ; latitude and longitude
      CASE count OF
         2: BEGIN
            msg = 'Cannot parse position'
            ll = string2coord(pix[0],1,err)
            IF err EQ 2 THEN  ok = DIALOG_MESSAGE(msg, /ERROR, $
                                                  DIALOG_PARENT=event.TOP)
            bb = string2coord(pix[1],0,err)
            IF err EQ 2 THEN ok = DIALOG_MESSAGE(msg, /ERROR, $
                                                  DIALOG_PARENT=event.TOP)
            IF err EQ 1 THEN ok = DIALOG_MESSAGE( $
                   'angle in hours not allowed for latitude', $
                                 /ERROR, DIALOG_PARENT=event.TOP)
            IF err NE 0 THEN RETURN
         END         
         4: BEGIN               ; Assume long = RA if given in sexagesimal
            pix = FLOAT(pix)
            ll = TEN(pix[0],pix[1],0d0)*15d0
            bb = TEN(pix[2],pix[3],0d0)
         END
         6: BEGIN   
            pix = FLOAT(pix)
            ll = TEN(pix[0],pix[1],pix[2])*15d0
            bb = TEN(pix[3],pix[4],pix[5])
         END  
         ELSE: ok = DIALOG_MESSAGE('Enter two numbers for coordinates', $
                                   DIALOG_PARENT=event.TOP, /ERROR)
      ENDCASE
   END
ENDCASE

IF form.COORD GT 0 && state.IS_ASTROM THEN $
    AD2XY, ll, bb, *state.ASTROM, xpix, ypix 

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

                                ; Re-plot any graphics overlays:
overlay, mode, state.astrom, state.pole_pix, state.nside

; If panels remain to be loaded, send timer event to DRAW window
; requesting re-draw:
IF done EQ 0 THEN WIDGET_CONTROL, (*state.TABARR)[0].DRAW, TIMER = 0.5

WIDGET_CONTROL, state.TABS, SET_UVALUE=mode
pix_print, state, 0, start

restore_lut, dummy, old_graph

WIDGET_CONTROL, state.TABS, /SENSITIVE
END

PRO ximview_autoscale_all, event
; Autoscales all tabs
;
COMPILE_OPT IDL2, HIDDEN
COMMON gr_global
ON_ERROR, 1

WIDGET_CONTROL, event.TOP,  GET_UVALUE = state
WIDGET_CONTROL, state.TABS,  GET_UVALUE = mode

ntab = N_ELEMENTS(*state.TABARR)
IF ntab EQ -1 THEN BEGIN
    MESSAGE, /INFORMATIONAL, 'No image to scale'
    RETURN
ENDIF

col = INDGEN(ntab)
tablab = get_tab_uvals(state.TABS)

WIDGET_CONTROL, state.TABS, SENSITIVE = 0

do_plot = 0
acted = 0B

FOR itab=0,ntab-1 DO BEGIN
    str = (*state.TABARR)[itab]
    IF ~PTR_VALID(str.BYTE_PTR) THEN CONTINUE ; funny tab

    name = tablab[str.SCREEN]
    MESSAGE, /INFORMATIONAL, STRING(name, FORMAT = "('Tab: ',A)")

    newscale = 0B
    IF str.MODE EQ 0.0 THEN BEGIN
        IF ~str.IM_PTR THEN BEGIN
            howto = 1
            imap = str.SCREEN
            WIDGET_CONTROL, state.LABEL, GET_UVALUE = data, /NO_COPY
        ENDIF ELSE BEGIN
            howto = 3
            data = str.IM_PTR
            imap = 0
        ENDELSE
                                ;Get mode and rms estimate if not done already
        scaling_params, dummy, data, imap, 0, howto, col, $
          ar, zero, sdev, PLOT = do_plot
        newscale = 1B
        IF howto EQ 1 THEN $
          WIDGET_CONTROL, state.LABEL, SET_UVALUE = data, /NO_COPY

        MESSAGE, /INFORMATIONAL, 'Found mode: ' + $
          numunit(zero*str.MULT, str.UNIT) + $
          ' and estimated rms: ' + numunit(sdev*str.MULT, str.UNIT)
    ENDIF ELSE BEGIN
        ar   = str.ABSRANGE
        zero = str.MODE
        sdev = str.SDEV
    ENDELSE
    izero = str.IZERO
                                ; Set min and max appropriate for pseudo-col:
    set_scale, [ar, zero, sdev], str
    range = str.RANGE
    IF newscale THEN set_print_fmt, str
    
    (*state.TABARR)[itab] = str
     
    IF str.COLTAB EQ 4 && str.IZERO NE izero THEN BEGIN
        ; Red-Black-Blue colour scale: re-set with correct zero level.
        ximview_lut, str.COLTAB, str.izero, decomp
        TVLCT, r, g, b, /GET
        *str.LUT = {r:r, g:g, b:b, line: !P.color, absent: !P.background}
    ENDIF 

                                ; Rescale byte images:
    IF ~str.IM_PTR THEN BEGIN
        WIDGET_CONTROL, state.LABEL, GET_UVALUE = data, /NO_COPY
        *str.BYTE_PTR = scale_image(data, str.RANGE, str.WRAP, str.TRFUNC, $
                                  str.ZERO, str.BETA, ABS_RANGE = str.ABSRANGE)
        WIDGET_CONTROL, state.LABEL, SET_UVALUE = data, /NO_COPY
    ENDIF ELSE *str.BYTE_PTR = scale_image(*str.IM_PTR, str.RANGE, str.WRAP, $
                                   str.TRFUNC, $
                                   str.ZERO, str.BETA, ABS_RANGE= str.ABSRANGE)
    acted = 1B
ENDFOR

IF ~acted THEN RETURN ; no actual scaling done.
nside = state.NSIDE
IF nside NE 0 && state.IS_ASTROM THEN $
  fill_gores, nside, state.IMSIZE, state.ASTROM, (*state.TABARR).BYTE_PTR
  
swap_lut, *state.XIM_GRAPH, dummy, old_graph

FOR itab=0,ntab-1 DO BEGIN
    IF redraw_req THEN BEGIN
        lutptr = (*state.TABARR)[itab].LUT
        TVLCT, (*lutptr).R, (*lutptr).G, (*lutptr).B
        !P.background = (*lutptr).ABSENT
        !P.color      = (*lutptr).LINE
    ENDIF
    update_screen, state.TABARR, itab, mode, done
ENDFOR

restore_lut, dummy, old_graph

mode.DONE = done
                                ; Request another go if loading not finished
IF done EQ 0 THEN WIDGET_CONTROL, (*state.TABARR)[0].DRAW, TIMER = 0.5

WIDGET_CONTROL, state.TABS, /SENSITIVE
WIDGET_CONTROL, state.TABS,  SET_UVALUE = mode, /NO_COPY

END
;
PRO ximview_grid, event
; Displays coordinate grids
;
COMPILE_OPT IDL2, HIDDEN
ON_ERROR, 0

WIDGET_CONTROL, event.ID,   GET_VALUE = label
WIDGET_CONTROL, event.TOP,  GET_UVALUE = state
WIDGET_CONTROL, state.TABS, GET_UVALUE = mode

astrom = state.ASTROM
IF SIZE(*astrom,/TYPE) NE 8 THEN BEGIN
   ok = DIALOG_MESSAGE('No astrometry information to plot coordinates', $
                       /ERROR, DIALOG_PARENT = top)
   RETURN
ENDIF

CASE label OF
   'On/Off': BEGIN              ; Toggle coord line plotting
      mode.GRID_PLOT = ~mode.GRID_PLOT
   END 
  'Grid interval': BEGIN
     IF mode.GRID_STEP[0] GT 0d0 THEN BEGIN
        fmt = "(F7.4,',',F8.4)"
        step = STRING(mode.GRID_STEP, FORMAT = fmt)
     ENDIF ELSE step = '30.0000, 30.0000'
     grid = get_user_datum(event,'lon/lat grid interval [degrees]', 20, step)
     IF STRCMP(grid, 'CANCEL', 6, /FOLD_CASE) THEN RETURN
     grid = DOUBLE(STRSPLIT(grid,', ',/EXTRACT, COUNT = count))
     IF count EQ 1 THEN grid = [grid, grid]
     IF grid[0] LT 0d0 THEN mode.GRID_CALC = 1B ELSE  BEGIN
        mode.GRID_STEP = grid
        mode.GRID_CALC = 0B
     ENDELSE
  END
  ELSE: MESSAGE, /INFORMATIONAL, 'Option ' + label + ' not yet available.'
ENDCASE

; If requested, plot or re-plot grid right away:
IF mode.grid_plot && ~mode.overview THEN BEGIN
   catplot = mode.CATPLOT
   mode.CATPLOT = 0B ; no neet to replot catalogues
   swap_lut, *state.XIM_GRAPH, (*state.TABARR)[0], old_graph
   overlay, mode, astrom, state.pole_pix, state.nside
   restore_lut, dummy, old_graph
   mode.CATPLOT = catplot
ENDIF
WIDGET_CONTROL, state.TABS, SET_UVALUE = mode, /NO_COPY

END
;
PRO ximview_vo, event
; Plots catalogs, sets off catalog processing loop, downloads
; catalogues from the Virtual observatory.
;
COMPILE_OPT IDL2, HIDDEN
ON_ERROR, 0

WIDGET_CONTROL, event.ID, GET_VALUE = label
WIDGET_CONTROL, event.TOP, GET_UVALUE = state

doprog = 0B
CASE label OF
  'On/Off': BEGIN
     WIDGET_CONTROL, state.TABS, GET_UVALUE = mode, /NO_COPY
     mode.CATPLOT = ~mode.CATPLOT
     WIDGET_CONTROL, state.TABS, SET_UVALUE=mode, /NO_COPY
  END
  'Process catalogue': BEGIN
     cats = ['NVSS','SUMSS','EMU_giants']
     incat = get_user_item(event,'Choose catalogue',cats)
     incat = cats[incat]
     prog = 'gencat_proc'
     state.prog = incat
     WIDGET_CONTROL, event.TOP, SET_UVALUE = state
     WIDGET_CONTROL, state.progcol, EVENT_PRO = prog
     WIDGET_CONTROL, state.progcol, SENSITIVE = 1B
     CALL_PROCEDURE, prog, event
  END
  'Load new catalogue': MESSAGE, /INFORMATIONAL, 'VO access not yet enabled'

  ELSE: MESSAGE, /INFORMATIONAL, 'Option ' + label + ' not yet available.'
ENDCASE

RETURN

; TODO Launch widget to select:
;   - (Vizier or local file)
;   - Catalog (provide list of options in droplist)
;   - Search constraints (string)
;   - add or replace to existing catalogues?

info = queryvizier(cat_name, 'NONE', CONSTRAINT = constraint, $
                   /SILENT, /ALLCOLUMNS)
; Find appropriate columns via UCDs
ucd = SXPAR(info.hdr,'UCD__*', COUNT=nucd)

; Parse UCDs
;
IF nucd GT 0 THEN  BEGIN
   main = BYTARR(nucd)
   ucd = STRLOWCASE(ucd)        ; just in case, should already be lower case
                                ; Split into words
   FOR iu=0,nucd-1 DO BEGIN
      word = STRSPLIT(ucd[iu],';',/EXTRACT,COUNT=nw)
                                ; Strip off any domain prefix:
      FOR iw=0,nw-1 DO BEGIN
         ic = STREGEX(word[iw],':')
         IF ic GE 0 THEN word[iw] = STRMID(word[iw],ic+1)
      ENDFOR
      ucd[iu] = word[0]         ; primary UCD
      null = WHERE(word EQ 'meta.main',nm)
      main = nm EQ 1 
   ENDFOR

   racol = WHERE(ucd EQ 'pos.eq.ra',nra)
   IF nra GE 1 THEN BEGIN 
      deccol = WHERE(ucd EQ 'pos.eq.dec') 
      If nra GT 1 THEN BEGIN    ; use main column
         id = WHERE(main[racol])
         racol = racol[id]
         id = WHERE(main[deccol])
         deccol = deccol[id]
      ENDIF
      racol = racol[0] + 1     ; make a scalar and offset as header is info.(0)
      deccol = deccol[0] + 1
      catsys = 'C'
      lon = info.(racol)
      lat = info.(deccol)
   ENDIF ELSE BEGIN
      llcol = WHERE(ucd EQ 'pos.galactic.lon', nll)
      IF nll GE 1 THEN BEGIN
         bbcol = WHERE(ucd EQ 'pos.galactic.lat')
         IF nll GT 1 THEN BEGIN
            id = WHERE(main[llcol])
            llcol = llcol[id]
            id = WHERE(main[bbcol])
            bbcol = bbcol[id]
         ENDIF
         llcol = llcol[0] + 1
         bbcol = bbcol[0] + 1
         catsys = 'G'
         lon = info.(llcol)
         lat = info.(bbcol)
      ENDIF 
   ENDELSE
   nsource = N_ELEMENTS(lon)
   
   namecol = WHERE(ucd EQ 'meta.id',nid)
   If nid GT 1 THEN BEGIN       ; use main column
      id = WHERE(main[namecol])
      namecol = namecol[id]
   ENDIF
   IF nid GE 0 THEN BEGIN
      namecol = namecol[0] + 1
      labels = info.(namecol)
   ENDIF ELSE labels = STRARR(nsource)
ENDIF ELSE BEGIN
   MESSAGE, /INFORMATIONAL, 'No UCDs: cannot find positions in catalogue'
   RETURN
ENDELSE

; Make sure coords are consistent with astrometry of image
; and if not, convert.
co_sys = *state.ASTROM.COORD_SYS
IF co_sys NE cat_sys THEN BEGIN
  IF co_sys NE 'E' and co_sys NE 'G' THEN MESSAGE, /INFORMATIONAL, $
    "I'm not going to try to convert catalogue coords!" $
  ELSE BEGIN
    ang2vec, lon, lat, vec, /ASTRO
    outvec = ROTATE_COORD(vec, INCO = catsys, OUTCO = co_sys)
    vec2ang, outvec, lon, lat, /astro
    vec = 0
  ENDELSE
ENDIF
symb   = 1        ; cross  
colour = !P.COLOR ; Fixed for the moment.
do_label = 0      ; no thanks

; Calculate image pixel coordinates:

AD2XY, lon, lat, *STATE.astrom, xpix, ypix
catalog = {SYMBOL: symb, COLOUR: colour, DO_LABEL: do_label, $
           LABELS: labels, LON: lon, LAT: lat, XPIX: xpix, YPIX: ypix}

; Stash the catalog:
pcat = PTR_NEW(catalog)
cats = *MODE.catalog
IF N_ELEMENTS(cats) EQ 0 THEN cats = [pcat] ELSE cats = [cats, pcat]

*mode.catalog = cats

WIDGET_CONTROL, state.TABS, SET_UVALUE=mode

; Plot catalog right away
IF ~mode.overview THEN BEGIN
   mode.grid_plot = 0B          ; no need to re-plot grid.
   swap_lut, *state.XIM_GRAPH, (*state.TABARR)[0], old_graph
   overlay, mode, state.ASTROM, state.pole_pix, state.nside
   restore_lut, dummy, old_graph   
ENDIF

END
;
PRO ximview_prog, event
;
; Routine to handle events from programmable buttons A, B,...
; Obsolete, as active program is set as EVENT_PRO for progcol
; Now just used as a placeholder when no active program
;
ON_ERROR, 1

start = SYSTIME(1)

WIDGET_CONTROL, event.ID, GET_VALUE = label
WIDGET_CONTROL, event.TOP, GET_UVALUE = state
prog = ''
IF prog EQ '' THEN BEGIN
   MESSAGE, /INFORMATIONAL, 'No program set'
   RETURN
ENDIF
CASE label OF
   'A': BEGIN
      a = a
   END
   'B': BEGIN
      b = b
   END
   'C': BEGIN
      c = c
   END
  ELSE: 
ENDCASE
swap_lut, *state.XIM_GRAPH, (*state.TABARR)[0], old_graph
restore_lut, dummy, old_graph
END
;
PRO ximview_profile_event, event
;
;  Processes events from profile window
;
COMPILE_OPT IDL2, HIDDEN
ON_ERROR, 1

WIDGET_CONTROL, event.TOP, GET_UVALUE = dat

tag = TAG_NAMES(event, /STRUCTURE_NAME)
CASE tag OF
   'WIDGET_BUTTON': BEGIN
      WIDGET_CONTROL, event.ID, GET_VALUE = label
      CASE label OF
         'Write Postscript': BEGIN 
            name = GET_USER_DATUM(event, 'Name for PS file', 30, $
                                  'ximview_profile.eps')
            IF name EQ 'Cancel' THEN RETURN
            old_dev = !D.name
            SIMPLE_COLOUR, PS = name
            xtitle = 'offset from start of slice ['+dat.l_unit+']'
            ytitle = 'data value [' + dat.unit +']'

            PLOT, dat.lp, dat.pp, XTITLE=xtitle, XSTYLE=1, $
                  YTITLE=ytitle, TITLE=dat.title
            DEVICE, /CLOSE
            SIMPLE_COLOUR ; restore graphics settings to TV mode
            SET_PLOT, old_dev
         END
         'Save Data': BEGIN ; export data in some format or other.
            name = GET_USER_DATUM(event, 'Name for save file', 30, $
                                  'ximview_profile.sav')
            IF name EQ 'Cancel' THEN RETURN
            SAVE, dat, FILE = file
         END
         'Done': WIDGET_CONTROL, event.TOP, /DESTROY
      ENDCASE
   END
   'WIDGET_DRAW': BEGIN ; cursor moved, update readout
      plot_x = FLOAT(event.x)
                                ; Now translate plot x position into
                                ; position on x axis and hence image pos.
      lval = (plot_x/dat.vx - dat.x0)/dat.x1 > 0.
      id = WHERE(dat.length GT lval, nid)
      IF nid EQ 0 THEN id = N_ELEMENTS(dat.pix_label)
      id = id[0] - 1
      WIDGET_CONTROL, dat.READOUT, SET_VALUE = dat.pix_label[id]
   END
   ELSE: MESSAGE, 'unexpected event'
ENDCASE
END
;
PRO boundary_crossings, xa, ya, xx, yy, dx, dy
;
; Finds intersection of line with pixel boundaries
;
; INPUTS
;     xa       Start and end x coordinates of line
;     ya       Start and end y coordinates of line
;
; INPUT/OUTPUT
;     xx       x coordinates of crossings (new points appended)
;     yy       y coordinates of crossings (new points appended)
;
; OUTPUTS
;     dx       x-offset between start and end
;     dy       y-offset between start and end
;
dx = xa[1] - xa[0]
dy = ya[1] - ya[0]

ix  = ROUND(xa)
iy  = ROUND(ya)
nx = ix[1] - ix[0]
ny = iy[1] - iy[0]

xsign = nx GE 0 ? 1 : -1
IF nx NE 0 THEN BEGIN ; intercepts with x boundaries
   xx1 = ix[0] + xsign*(FINDGEN(ABS(nx)) + 0.5)
   slope = dy / dx
   intercept = ya[1] - slope * xa[1]
                                ; Get intersections with pixel
                                ; boundaries between pt & pt1 
   yy1 = slope*xx1 + intercept
   xx = [xx,xx1]
   yy = [yy,yy1]
ENDIF
ysign = ny GE 0 ? 1 : -1
IF ny NE 0 THEN BEGIN ; intercepts with y boundaries
   yy1 = iy[0] + ysign*(FINDGEN(ABS(ny)) + 0.5)
   slope = dx / dy
   intercept = xa[1] - slope * ya[1]
   xx1 = slope*yy1 + intercept
   xx = [xx,xx1]
   yy = [yy,yy1]
ENDIF

END
;
PRO ximview_profile, event, plot_data
;
; Draws a profile along a specified line through the image.
; NB Start of profile is PT1, end is PT, i.e. start is second-to-last
; point marked and end is last point.
;
; Cursor readout gives image position corresponding to any position
; along the slice.
;
; TODO: save in fits format?
; TODO: user selects great circle vs straight line
; TODO: Extend length
; TODO; rubber line
;
WIDGET_CONTROL, event.TOP, GET_UVALUE = state
WIDGET_CONTROL, state.TABS, GET_UVALUE = mode

IF mode.XPT*mode.XPT1 LT 0 THEN BEGIN
   ok = DIALOG_MESSAGE('Mark two successive points to define profile', $
                  DIALOG_PARENT=event.TOP, /ERROR)
   RETURN
ENDIF

xa = [mode.xpt1, mode.xpt]
ya = [mode.ypt1, mode.ypt]

; Find intersections of slice path with pixel boundaries
xx = [mode.xpt1]
yy = [mode.ypt1]
IF state.is_astrom THEN astrom = *state.astrom
IF state.is_astrom THEN BEGIN ; great circle path
   xy2ad, xa, ya, astrom, lon, lat
   interval = TOTAL(ABS(astrom.cdelt))
   arc = gcirc_arc(lon[0], lat[0], lon[1], lat[1], interval, length)
   posang,  1, lon[0]/15, lat[0], lon[1]/15,lat[1], angle

   ad2xy, arc[*,0], arc[*,1], astrom, xa, ya
   np = N_ELEMENTS(xa)
   FOR i=0L,np-2 DO boundary_crossings, xa[i:i+1], ya[i:i+1], xx, yy
   xy2ad, xx, yy, astrom, lonc, latc
   good = WHERE(FINITE(lonc) AND FINITE(latc),ng)
   IF ng EQ 0 THEN BEGIN
      ok = DIALOG_MESSAGE('No good data on slice', $
                          DIALOG_PARENT=event.TOP, /ERROR)
      plot_data = {title: 'no data'}
      RETURN
   ENDIF
   xx = xx[good]
   yy = yy[good]
   lonc = lonc[good]
   latc = latc[good]
   gcirc, 2, lon[0], lat[0], lonc, latc, ll
   l_unit = 'arcsec'
ENDIF ELSE BEGIN ; Simple case: straight line in the projection
   boundary_crossings, xa, ya, xx, yy, dx, dy
   angle = !radeg* ATAN(-dx,dy)
   length = SQRT(dx^2 + dy^2)
   ll = SQRT((xx-xa[0])^2 + (yy-ya[0])^2)
   l_unit = 'pixel'
ENDELSE
xx = [xx,mode.xpt]
yy = [yy,mode.ypt]
ll = [ll,length]

; sort and eliminate duplicate points:
id = UNIQ(ll,SORT(ll))
npt = N_ELEMENTS(id)
ll = ll[id]
xx = xx[id]
yy = yy[id]
; find pixels for each segment:
xp = ROUND((xx[1:*] + xx[0:npt-2])/2)
yp = ROUND((yy[1:*] + yy[0:npt-2])/2)
IF state.is_astrom THEN xy2ad, xp, yp, astrom, lonp, latp 

str = (*state.TABARR)[0]

; Plot slice over image:
get_jumps, xa, ya, bad, njump
xy =im2tv(xa, ya, mode)
swap_lut, *state.XIM_GRAPH, str, old_graph
FOR j=0,njump DO BEGIN
   lensegment = bad[j+1] - bad[j]
   IF lensegment GT 1 THEN BEGIN
      u1 = xy[bad[j]:bad[j+1]-1,0]
      v1 = xy[bad[j]:bad[j+1]-1,1]
      PLOTS, u1, v1, /DEVICE
   ENDIF
ENDFOR
EMPTY
restore_lut, dummy, old_graph

; Get pixel values (in output units)
imptr = str.IM_PTR
pixval = str.mult*(*imptr)[xp,yp]

; Format text label describing each pixel on the profile:
plabel = STRARR(npt-1)
fmt = str.FORMAT
hpx = state.NSIDE GT 0
IF hpx THEN BEGIN
   r2d = 180d0/!dpi
   phi = lonp/r2d
   theta = (90d0 - latp)/r2d
   ANG2PIX_NEST, state.nside, theta, phi, nix
   ANG2PIX_RING, state.nside, theta, phi, rix
ENDIF
FOR ip = 0L, npt-2 DO BEGIN
   plabel[ip] = STRING(xp[ip],yp[ip], FORMAT="(I5,1X,I5,' ')")
   IF hpx THEN plabel[ip] += STRING(nix[ip], rix[ip], FORMAT = "(2I9)")
   IF state.is_astrom THEN BEGIN
      plabel[ip] += astrom.COORD_SYS EQ 'C' ?  $
              adstring(lonp[ip], latp[ip], state.PREC)     :  $
              STRING(lonp[ip], latp[ip], FORMAT = state.LL_FMT)
   ENDIF 
   plabel[ip] += STRING(pixval[ip], FORMAT=fmt)
ENDFOR

; Create arrays to plot stepped-line graph:
; duplicate all points in slice length list except first and last
lp = REFORM(TRANSPOSE(REBIN(ll,npt,2)),2*npt)
lp = lp[1:-2]
; Duplicate all points in pixval (NB one element shorter so no
; trimming needed)
pp = REFORM(TRANSPOSE(REBIN(pixval,npt-1,2)),2*(npt-1))

; Title plot with location of slice:
title = STRING(mode.xpt1,mode.ypt1,angle, FORMAT = $
               "('Starting at pixel (',2F7.1,') at PA =',F6.1,' deg')")


swap_lut, dummy1, dummy, old_graph
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Launch widget for plot
ttit = 'Ximview Profile'
profl  = WIDGET_BASE(/COLUMN, TITLE = ttit)
base   = WIDGET_BASE(profl, /COLUMN)
draw   = WIDGET_DRAW(profl, XSIZE=640, YSIZE=512, RETAIN = 2, /MOTION_EVENTS)

; Set up cursor readout area below image tabs:
                                ; find a font
temp    = WIDGET_BUTTON(profl, VALUE = 'Redshirt')
devfont = WIDGET_INFO(temp, /FONTNAME)
devfont = WIDGET_INFO(temp, /FONTNAME)
WIDGET_CONTROL, temp, /DESTROY
is_unix = STRCMP(!version.OS_FAMILY,'UNIX', 4,/ FOLD_CASE)
fixedfont = is_unix ? devfont : 'lucida console*10'

                                ; make label for readout
label = state.TITLE
title_string = label.HEAD + form_unit(str.UNIT)
readlab  = WIDGET_LABEL(profl, value = title_string, /ALIGN_LEFT, $
                       FONT = fixedfont, /DYNAMIC_RESIZE)
readout  = WIDGET_LABEL(profl, value = 'No pixel assigned', $
                       FONT = fixedfont, /ALIGN_LEFT, /DYNAMIC_RESIZE)

controls = WIDGET_BASE(profl, /ROW)
write_ps = WIDGET_BUTTON(controls, VALUE = 'Write Postscript')
save_dat = WIDGET_BUTTON(controls, VALUE = 'Save Data')
pad      = WIDGET_BASE(controls)
done     = WIDGET_BUTTON(controls, VALUE = 'Done')

bgeom  = WIDGET_INFO(profl,/GEOMETRY)
dgeom  = WIDGET_INFO(done,/GEOMETRY)
cgeom  = WIDGET_INFO(controls, /GEOMETRY)
length = dgeom.XOFFSET + dgeom.XSIZE + cgeom.xpad
dx  = FIX(bgeom.XSIZE - bgeom.XPAD - length)
WIDGET_CONTROL, pad, XSIZE = dx

WIDGET_CONTROL, profl, /REALIZE
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; Now draw profile plot
WIDGET_CONTROL, draw, GET_VALUE = window
WSET, window
DEVICE, decomposed = 0
col0 = 240
linecols, col0, ci
!P.color = 255                  ; white
!P.background = 0B              ; black

xtitle = 'offset from start of slice ['+l_unit+']'
ytitle = 'data value [' + str.unit +']'
PLOT, lp, pp, XTITLE=xtitle, XSTYLE=1, YTITLE=ytitle, TITLE=title
IF 0 THEN BEGIN ; debug option
   OPLOT, ll, pixval, LINE=1, COLOR=col0+2, PSYM=10
ENDIF

; Supply data to widget for readout, PDF output, and save output:
plot_data = {title: title, lp: lp, pp: pp, length: ll, l_unit: l_unit, $
             pixval: pixval, unit: str.UNIT, pix_label: plabel, $
             readout: readout, vx: !D.x_vsize, x0: !X.S[0], x1: !X.s[1]}

restore_lut, dummy, old_graph

WIDGET_CONTROL, profl, SET_UVALUE = plot_data

                                ; Allow user to get rid of plot!
XMANAGER, 'ximview_profile', profl, /NO_BLOCK
 
END
;
PRO ximview_profile_opts, event
;
; Sets options for profiles. Launch widget to allow user to select
;  * Great circle or linear slice (if astrometry is set)
;  * Run slice between last two marked points, or run from last marked
;  point to cursor position until clicked, with slice line
;  continuously updated on screen ("rubber splice")
;  Extend slice beyond end points. Default: 2 beam areas. Specify in
;  beams (default, if available), or pixels.
;
WIDGET_CONTROL, event.TOP, GET_UVALUE = state
WIDGET_CONTROL, state.TABS, GET_UVALUE = mode
;

MESSAGE, /INFORMATIONAL, 'No options to set yet'
END
PRO ximview_geom, event, result
; Ruler / Protractor
;
; TODO: elastic ruler with popup to peak up end point?
; 
ON_ERROR, 0

r2d = 180d0/!dpi
r2as = 3600d0*r2d
WIDGET_CONTROL, event.TOP, GET_UVALUE = state
WIDGET_CONTROL, state.TABS, GET_UVALUE = mode

label = WIDGET_INFO(event.ID, /NAME)
IF label EQ 'BUTTON' THEN WIDGET_CONTROL, event.ID,  GET_VALUE = label


badpt = mode.XPT*mode.XPT1 LT 0
IF badpt && label EQ 'Ruler' THEN BEGIN
   ok = DIALOG_MESSAGE('Mark two successive points to use ruler', $
                  DIALOG_PARENT=event.TOP, /ERROR)
   RETURN
ENDIF
xx = [mode.XPT, mode.XPT1]
yy = [mode.YPT, mode.YPT1]

doangle = label EQ 'Protractor'
IF doangle THEN BEGIN
   badpt = badpt OR mode.xpt2 LT 0
   IF badpt THEN BEGIN
      ok = DIALOG_MESSAGE('Mark three successive points to use protractor', $
                          DIALOG_PARENT=event.TOP, /ERROR)
      RETURN
   ENDIF
   xx = [xx, mode.XPT2]
   yy = [yy, mode.YPT2]
ENDIF

; Calculate distances between last two marked points (XPT, XPT1 etc):
astro = state.IS_ASTROM
fmt1 = "('Distance between last two marked points: ',"
IF astro THEN BEGIN
   XY2AD, xx, yy, *state.astrom, ra, dc
   pixsize = TOTAL(ABS((*state.astrom).cdelt)) ; Actually dx+dy in degrees
   interval = pixsize*10 < 1d0
   arc = gcirc_arc(ra[0], dc[0], ra[1], dc[1], interval, dis)
   AD2XY, arc[*,0], arc[*,1], *state.astrom, xa,ya
   
   unit = 'arcsec'
   PRINT, dis, unit, dis/60.0, FORMAT = fmt1+"G0.5,1X,A,2X,G0.5,' arcmin')"
   posang, 1, ra[1]/15, dc[1], ra[0]/15,dc[0], pa_deg
   PRINT, pa_deg, FORMAT = $
          "('New point at position angle',F7.1,' degrees wrt old point')"
ENDIF ELSE BEGIN
   dx1 = xx[0] - xx[1]
   dy1 = yy[0] - yy[1]
   dis = SQRT(dx1^2 + dy1^2)
   unit = 'pixel'
   PRINT, dis, unit, FORMAT = fmt1+"G0.5,1X,A)"
   xa = xx[0:1]
   ya = yy[0:1]
   pa_deg = !radeg*ATAN(dy1,dx1)
ENDELSE
xy =im2tv(xa, ya, mode)

IF doangle THEN BEGIN
; Find angle 2-1-0

   IF astro THEN BEGIN
      arc2 = gcirc_arc(ra[1], dc[1], ra[2], dc[2], interval, dis2)
      AD2XY, arc2[*,0], arc2[*,1], *state.astrom, xa2,ya2
      posang, 1, ra[1]/15, dc[1], ra[2]/15,dc[2], pa2_deg
      angle = (pa2_deg - pa_deg) ; degrees
                                ; Fold into +/- 180 Deg
      angle MOD= 360d0
      IF angle GT  180d0 THEN angle -= 360d0
      IF angle LT -180d0 THEN angle += 360d0
   ENDIF ELSE BEGIN
      dx2 = xx[2] - xx[1]
      dy2 = yy[2] - yy[1]
      dis2 = SQRT(dx2^2 + dy2^2)
      angle = r2d*ATAN(dx1*dy2-dy1*dx2, dx2*dx1 + dy2*dy1)
      xa2 = xx[1:2]
      ya2 = yy[1:2]
   ENDELSE
   PRINT, angle, FORMAT= $
          "('Angle between last three marked points:',F7.1,' degrees')"

; plot angle marker. NB angle on projection is not angle on sky
   radius = ((0.5 * MAX([dis,dis2])) > 8) < MIN([dis,dis2]) ; arcsec
   nseg = FIX(radius*ABS(angle) / (3600*pixsize)) > 3
   phi = (INDGEN(nseg+1)*(angle/nseg) + pa_deg)/r2d
   IF astro THEN BEGIN
      radius /= 3600            ; degrees
      theta = REPLICATE(radius/r2d,nseg+1)
      ang2vec, theta, !dpi-phi, vec
      mat = euler_matrix_new(ra[1],90d0-dc[1],0d0,/Y,/DEG)
      rotvec = mat ## vec
      vec2ang, rotvec, theta, phi
      rac = phi*r2d
      dcc = 90 - theta*r2d
      ad2xy, rac, dcc, *state.astrom, xc, yc
   ENDIF ELSE BEGIN
      xc = xx[1] + radius*COS(phi)
      yc = yy[1] + radius*SIN(phi)
   ENDELSE
ENDIF ELSE angle = !values.F_NAN
swap_lut, *state.XIM_GRAPH, (*state.TABARR)[0], old_graph

PLOTS, xy[*,0], xy[*,1], /DEVICE
IF doangle THEN BEGIN
   xy2 = im2tv(xa2,ya2,mode)
   PLOTS, xy2[*,0],xy2[*,1], /DEVICE, LINE=2
   xy2 = im2tv(xc,yc,mode)
   PLOTS, xy2[*,0],xy2[*,1], /DEVICE
ENDIF
EMPTY
restore_lut, dummy, old_graph

result = {distance: dis, dunit: unit, posang: pa_deg, $
          angle: angle, aunit: 'degree'}
END
;
PRO ximview_imstat_opts, event
; Sets options for imstats
;
COMPILE_OPT IDL2, HIDDEN
ON_ERROR, 1

WIDGET_CONTROL, event.ID, GET_VALUE = label
WIDGET_CONTROL, event.TOP,  GET_UVALUE = state
CASE label OF
    'Region of Interest': state.ROI = 1
    'Box':                state.ROI = 0
    'Set box size': BEGIN
; adjust imstats window
       default = STRJOIN(STRTRIM(STRING(state.STATBOX),2),', ')
       text = get_user_datum(event, $
                             'Box size in pixels (x,y; or x for square)', $
                             12, default)
       IF STRCMP(text, 'CANCEL', 6, /FOLD_CASE) THEN RETURN
       strings = STRSPLIT(text,' ,',/EXTRACT)
       IF N_ELEMENTS(strings) EQ 1 THEN strings = REPLICATE(strings,2)
       state.STATBOX = FIX(strings)
    END
    'Set threshold': BEGIN
       str = (*state.TABARR)[0]
       default = STRTRIM(STRING(str.THRESHOLD),2)
       unit = str.UNIT
       text = get_user_datum(event, 'threshold for pixel count (' $
                             +unit+')',12, default)
       IF STRCMP(text, 'CANCEL', 6, /FOLD_CASE) THEN RETURN
       (*state.TABARR)[0].THRESHOLD = FLOAT(text)
    END
    ELSE: MESSAGE, /INFORMATIONAL, 'Option ' + label + ' not yet available.'
ENDCASE

WIDGET_CONTROL, event.TOP, SET_UVALUE=state
END

PRO ximview_maxfit_opts, event
; Sets options for maxfit
;
COMPILE_OPT IDL2, HIDDEN
ON_ERROR, 1

WIDGET_CONTROL, event.ID, GET_VALUE = label
WIDGET_CONTROL, event.TOP,  GET_UVALUE = state
CASE label OF
    'Find extremum': state.PEAK = 1
    'Find maximum':  state.PEAK = 0
    'Find minimum':  state.PEAK = 2
    'Set box size': BEGIN
; adjust imstats window
        text = get_user_datum(event, $
                              'Peakfit search box size (pixels)', $
                              3, STRING(state.MAXBOX,FORMAT = "(I3)"))
        IF STRCMP(text, 'CANCEL', 6, /FOLD_CASE) THEN RETURN
        state.MAXBOX = FIX(text)
    END
    ELSE: MESSAGE, /INFORMATIONAL, 'Option ' + label + ' not yet available.'
ENDCASE

WIDGET_CONTROL, event.TOP, SET_UVALUE=state
END
;
PRO ximview_setprop_event, event
;  Processes events from set property dialog
;
COMPILE_OPT IDL2, HIDDEN
ON_ERROR, 1

WIDGET_CONTROL, event.TOP, GET_UVALUE = info

tag = TAG_NAMES(event, /STRUCTURE_NAME)
CASE tag OF
    'WIDGET_DROPLIST': BEGIN
        CASE event.ID OF
            info.TABLIST: WIDGET_CONTROL, info.CURRENT, $
              SET_VALUE = info.PROPLIST[event.INDEX]
            info.CHOOSE: WIDGET_CONTROL, event.ID, SET_UVALUE = event.INDEX
        ENDCASE
    END
    'WIDGET_COMBOBOX': ;   WIDGET_CONTROL, event.ID, SET_UVALUE = event.STR
    'WIDGET_BUTTON': BEGIN
        WIDGET_CONTROL, event.ID, GET_VALUE = label
        CASE label OF
            'Accept': BEGIN     ; Extract results from widget
                tab = info.TABLIST EQ -1 ? -1 : WIDGET_INFO(info.TABLIST, /DROPLIST_SELECT)
;                WIDGET_CONTROL, info.CHOOSE, GET_UVALUE = property
;                IF property EQ -1 THEN $
                type = WIDGET_INFO(info.CHOOSE, /NAME)
                IF type EQ 'DROPLIST' THEN $
                  property = WIDGET_INFO(info.CHOOSE, /DROPLIST_SELECT) ELSE $
                  property = WIDGET_INFO(info.CHOOSE, /COMBOBOX_GETTEXT)
                IF info.CHOOSE2 NE -1 THEN BEGIN
                    WIDGET_CONTROL, info.CHOOSE2, GET_VALUE = prop2
                    IF ~STREGEX(prop2, '^ *$', /BOOLEAN) THEN property = prop2
                ENDIF
                result = {cancel: 0B, tab: tab, property: property}
            END
            'Cancel': result = {cancel: 1B}
        ENDCASE
        WIDGET_CONTROL, info.RETID, SET_UVALUE = result
        WIDGET_CONTROL, event.TOP, /DESTROY
    END
    ELSE: MESSAGE, /INFORMATIONAL, 'Unrecognised event type: ' + tag
ENDCASE

END
PRO ximview_setprop, event, tablab, scr0, property, proplist, options, $
             cancel, tab, new, EDITABLE = editable
; Launches a dialog to choose a property for a specified tab from a
; list (optionally editable, ie you can write your own value if not in
; the list) Mostly redundant as replaced by setprop2 except for coordsys.
;
; Inputs:
;     event:    Used to hold results
;     tablab:   Labels for the tabs available to be set, or 0 for global properties
;     scr0:     Index of tablab corresponding to the current tab
;     property: Name of property to set
;     proplist: Currently-assigned properties for each tab
;     options:  List of options for the property
;     editable: True if the user can enter a new value as well as
;               choosing from a list
; Outputs:
;     cancel:   True if operation abandoned
;     tab:      tab to set (index of tablab)
;     new:      New value of property.
;
COMPILE_OPT IDL2, HIDDEN

editable = KEYWORD_SET(editable)
; if tablab is not a list of strings, property is global:
global = SIZE(tablab,/TYPE) NE 7  
widget_title = global ? 'Set ' + property : 'Set tab ' + property
 
is_unix = STRCMP(!version.OS_FAMILY, 'UNIX', 4, /FOLD_CASE)


query = WIDGET_BASE(GROUP_LEADER = event.TOP, /MODAL, /COLUMN, $
                    TITLE = widget_title)
base    = WIDGET_BASE(query, /ROW)

IF ~global THEN BEGIN
   void    = WIDGET_LABEL(base, VALUE ='Tab to set')
   tablist = WIDGET_DROPLIST(base, VALUE = tablab, UVALUE = scr0)
   WIDGET_CONTROL, tablist, SET_DROPLIST_SELECT = scr0
ENDIF ELSE tablist = -1

void    = WIDGET_LABEL(base, VALUE = ' Current ' + property + ':')
current = WIDGET_LABEL(base, VALUE = proplist[scr0], /DYNAMIC_RESIZE)

now = WHERE(proplist[scr0] EQ options)
propset = now[0] NE -1
base = WIDGET_BASE(query, /ROW)
void = WIDGET_LABEL(base, VALUE = ' New '+ property + ':')
IF editable THEN BEGIN
    IF is_unix && !VERSION.release EQ 6.0 THEN BEGIN
        choose = WIDGET_DROPLIST(base, VALUE = options, UVALUE = now[0])
        IF propset THEN WIDGET_CONTROL, choose, SET_DROPLIST_SELECT = now[0]
        void = WIDGET_LABEL(base, VALUE = 'or:')
        choose2 = WIDGET_TEXT(base, VALUE = '', xsize = 10, /EDITABLE)
    ENDIF ELSE BEGIN
        choose = WIDGET_COMBOBOX(base, VALUE = options, $
                                 UVALUE = proplist[scr0], /EDITABLE)
        IF propset THEN WIDGET_CONTROL, choose, SET_COMBOBOX_SELECT = now[0]
        choose2 = -1
    ENDELSE
ENDIF ELSE BEGIN
    choose = WIDGET_DROPLIST(base, VALUE = options, UVALUE = now[0])
    IF propset THEN WIDGET_CONTROL, choose, SET_DROPLIST_SELECT = now[0]
    choose2 = -1
ENDELSE

bbase   = WIDGET_BASE(query, /ROW)
void   = WIDGET_BUTTON(bbase, VALUE = 'Accept')
void   = WIDGET_BUTTON(bbase, VALUE = 'Cancel')

info = {retID: event.ID, tablist: tablist, proplist: proplist, $
        current: current, choose: choose, choose2: choose2}

WIDGET_CONTROL, event.ID, GET_UVALUE = save
WIDGET_CONTROL, query, SET_UVALUE = info
WIDGET_CONTROL, query, /REALIZE
XMANAGER, 'ximview_setprop', query
WIDGET_CONTROL, event.ID, GET_UVALUE = result
WIDGET_CONTROL, event.ID, SET_UVALUE = save

IF SIZE(result,/TYPE) NE 8 THEN cancel = 1B ELSE cancel = result.CANCEL
IF cancel THEN RETURN

tab = result.TAB  &  new = result.PROPERTY

END
PRO ximview_setprop2_event, event
;  Processes events from set property dialog
;
COMPILE_OPT IDL2, HIDDEN
ON_ERROR, 0

WIDGET_CONTROL, event.TOP, GET_UVALUE = info

tag = TAG_NAMES(event, /STRUCTURE_NAME)

CASE tag OF
   'WIDGET_DROPLIST': BEGIN
      CASE event.ID OF
         info.TABLIST: BEGIN    ; Reset all current values
            WIDGET_CONTROL, info.CURRENT_POL, SET_VALUE = $
                            info.POLLIST[event.INDEX]               
            WIDGET_CONTROL, info.CURRENT_FREQ, SET_VALUE = $
                            info.FREQCODE[event.INDEX]               
            WIDGET_CONTROL, info.CURRENT_BMAJ, SET_VALUE = info.BMAJ[event.INDEX]
            WIDGET_CONTROL, info.CURRENT_BMIN, SET_VALUE = info.BMIN[event.INDEX]
            WIDGET_CONTROL, info.CURRENT_BPA,  SET_VALUE = info.BPA[event.INDEX]
            WIDGET_CONTROL, info.CURRENT_NOISE, SET_VALUE = $
                            info.ncode[event.INDEX]
            WIDGET_CONTROL, info.CURRENT_UNIT, SET_VALUE = $
                            info.units[event.INDEX]
         END
         ELSE: ; do nothing, read final value at end
      ENDCASE
   END
   'WIDGET_COMBOBOX':       ;   WIDGET_CONTROL, event.ID, SET_UVALUE = event.STR
   'WIDGET_BUTTON': BEGIN
      WIDGET_CONTROL, event.ID, GET_VALUE = label
      CASE label OF
         'Accept': BEGIN        ; Extract results from widget
            tab     = WIDGET_INFO(info.TABLIST, /DROPLIST_SELECT)
            polcode = WIDGET_INFO(info.CHOOSE_POL, /DROPLIST_SELECT) 
            ncode   = WIDGET_INFO(info.CHOOSE_NOISE, /DROPLIST_SELECT)
            WIDGET_CONTROL, info.CHOOSE_BMAJ, GET_VALUE = bmaj
            WIDGET_CONTROL, info.CHOOSE_BMIN, GET_VALUE = bmin
            WIDGET_CONTROL, info.CHOOSE_BPA,  GET_VALUE = bpa
            ; Calculate beam area in square degrees
            bmaj = bmaj EQ '' ? !values.F_NAN : FLOAT(bmaj[0])/3600.0
            bmin = bmin EQ '' ? !values.F_NAN : FLOAT(bmin[0])/3600.0
            bpa  = bpa  EQ '' ? !values.F_NAN : FLOAT(bpa[0])
            beam_area = bmaj*bmin * !pi/(4*ALOG(2.0))
            beam = {BEAM, bmaj: bmaj, bmin: bmin, bpa: bpa, beam_area: beam_area}

            WIDGET_CONTROL, info.CHOOSE_UNIT, GET_VALUE = unit
                                ; Do we have a combo box for freq?
            type  = WIDGET_INFO(info.CHOOSE_FREQ, /NAME)
            IF type EQ 'DROPLIST' THEN BEGIN
               id = WIDGET_INFO(info.CHOOSE_FREQ, /DROPLIST_SELECT) 
               freq = info.freqcode[id]
               WIDGET_CONTROL, info.CHOOSE_F2, GET_VALUE = f2
               IF ~STREGEX(f2, '^ *$', /BOOLEAN) THEN freq = f2
            ENDIF ELSE freq = WIDGET_INFO(info.CHOOSE_FREQ, /COMBOBOX_GETTEXT)
            
            result = {cancel: 0B, tab: tab, polcode: polcode, freq: freq, $
                      beam: beam, ncode: ncode, unit: unit}

         END
         'Cancel': result = {cancel: 1B}
      ENDCASE
      WIDGET_CONTROL, info.RETID, SET_UVALUE = result
      WIDGET_CONTROL, event.TOP, /DESTROY
   END
   ELSE: MESSAGE, /INFORMATIONAL, 'Unrecognised event type: ' + tag
ENDCASE

END
;
PRO ximview_setprop2, event, tablab, scr0, cancel, tab, new
; Launches a dialog to set various properties of the image displayed
; on each tab, namely the ones that may not be properly set in the
; image header.
;
; Current properties that can be set are
;   o polarization code
;   o frequency
;   o coordinate system (most likely for HEALPix map
;   o Beam properties (BMAJ, BMIN, BPA)
;   o Noise type: unknown (NT = 0)
;                 noise independent in each pixel (NT = 1)
;                 noise correlated over beam area (NT = 2)
;
; Inputs:
;     event:    Used to hold results
;     tablab:   Labels for the tabs available to be set 
;               (NB some properties are global and will be set for all tabs).
;     scr0:     Index of tablab corresponding to the current tab
;
; Outputs:
;     cancel:   True if operation abandoned
;     tab:      tab to set (index of tablab)
;     new:      New value of property.
;
WIDGET_CONTROL, event.TOP, GET_UVALUE = state

IF state.FIRST THEN RETURN ; haven't started yet.

widget_title = 'Set image properties' 
is_unix = STRCMP(!version.OS_FAMILY, 'UNIX', 4, /FOLD_CASE)

query = WIDGET_BASE(GROUP_LEADER = event.TOP, /MODAL, /COLUMN, $
                    TITLE = widget_title)
base    = WIDGET_BASE(query, /ROW)
;
; Set list of tabs with actual images (as opposed to blink, RGB etc)
;
screens  = (*state.TABARR).SCREEN
good     = WHERE( PTR_VALID( (*state.TABARR).BYTE_PTR ))
scr0     = screens[good[0]]
iscreen  = BSORT(screens)
good     = WHERE( PTR_VALID( (*state.TABARR)[iscreen].BYTE_PTR ))
iscreen  = iscreen[good]
tablab   = get_tab_uvals(state.TABS)
tablab   = tablab[screens[iscreen]]
scr0     = (WHERE(screens[iscreen] EQ scr0))[0]

void    = WIDGET_LABEL(base, VALUE ='Tab to set')
tablist = WIDGET_DROPLIST(base, VALUE = tablab, UVALUE = scr0)
WIDGET_CONTROL, tablist, SET_DROPLIST_SELECT = scr0

; Setting polarization
;
scodes = ['YX', 'XY', 'YY', 'XX', 'LR', 'RL', 'LL', 'RR', 'unknown', $
          'I', 'Q', 'U', 'V', 'Pol Intensity', 'Fractional Pol', 'Pol Angle', $
          'Spectral index', 'Optical depth', 'RM']
polcodes = (*state.TABARR).POLCODE
polstr   = scodes[polcodes+8]
tabcode  = polstr[iscreen]

base = WIDGET_BASE(query, /ROW)
void    = WIDGET_LABEL(base, VALUE = ' Polarization  current:')
current_pol = WIDGET_LABEL(base, VALUE = tabcode[scr0], /DYNAMIC_RESIZE)
now = WHERE(tabcode[scr0] EQ scodes)
propset = now[0] NE -1
void = WIDGET_LABEL(base, VALUE = ' New:')
choose_pol = WIDGET_DROPLIST(base, VALUE = scodes, UVALUE = now[0])
IF propset THEN WIDGET_CONTROL, choose_pol, SET_DROPLIST_SELECT = now[0]
;
; Setting frequency
;
freqcode = (*state.TABARR).FREQCODE
freqcode  = freqcode[iscreen]
nulls = WHERE(freqcode EQ '')
IF nulls[0] NE -1 THEN freqcode[nulls] = '          '

base = WIDGET_BASE(query, /ROW)
void    = WIDGET_LABEL(base, VALUE = ' Frequency   Current:')
current_freq = WIDGET_LABEL(base, VALUE = freqcode[scr0], /DYNAMIC_RESIZE)
now = WHERE(freqcode[scr0] EQ scodes)
propset = now[0] NE -1
void = WIDGET_LABEL(base, VALUE = ' New:')
IF is_unix && !VERSION.release EQ 6.0 THEN BEGIN ; combobox not available.
   choose_freq = WIDGET_DROPLIST(base, VALUE = freqcode, UVALUE = now[0])
   IF propset THEN WIDGET_CONTROL, choose_freq, SET_DROPLIST_SELECT = now[0]
   void = WIDGET_LABEL(base, VALUE = 'or:')
   choose_f2 = WIDGET_TEXT(base, VALUE = '', xsize = 10, /EDITABLE)
ENDIF ELSE BEGIN
   choose_freq = WIDGET_COMBOBOX(base, VALUE = options, $
                            UVALUE = freqcode[scr0], /EDITABLE)
   IF propset THEN WIDGET_CONTROL, choose, SET_COMBOBOX_SELECT = now[0]
   choose_f2 = -1
ENDELSE
;
; Setting beam
;
beam = (*state.TABARR).beam 
beam = beam[iscreen]
bmaj = beam.bmaj
sbmaj = STRTRIM(string(bmaj),2)
bmin = beam.bmin
sbmin = STRTRIM(string(bmin),2)
bpa  = beam.bpa
sbpa = STRTRIM(string(bpa),2)

base = WIDGET_BASE(query, /ROW)
void    = WIDGET_LABEL(base, VALUE = ' BMAJ (arcsec)  Current:')
current_bmaj = WIDGET_LABEL(base, VALUE = sbmaj[scr0], /DYNAMIC_RESIZE)
void = WIDGET_LABEL(base, VALUE = ' New:')
choose_bmaj =  WIDGET_TEXT(base, VALUE = '', xsize = 10, /EDITABLE)
void    = WIDGET_LABEL(base, VALUE = ' BMAJ (arcsec)  Current:')
current_bmin = WIDGET_LABEL(base, VALUE = sbmin[scr0], /DYNAMIC_RESIZE)
void = WIDGET_LABEL(base, VALUE = ' New:')
choose_bmin =  WIDGET_TEXT(base, VALUE = '', xsize = 10, /EDITABLE)
void    = WIDGET_LABEL(base, VALUE = ' BPA (degrees)  Current:')
current_bpa = WIDGET_LABEL(base, VALUE = sbpa[scr0], /DYNAMIC_RESIZE)
void = WIDGET_LABEL(base, VALUE = ' New:')
choose_bpa =  WIDGET_TEXT(base, VALUE = '', xsize = 10, /EDITABLE)
;
; Setting noise properties
;
inoise = (*state.TABARR).noise_type
inoise = inoise[iscreen]
ncode = ['unknown','Independent in each pixel','Correlated over beam']

base = WIDGET_BASE(query, /ROW)
void = WIDGET_LABEL(base, VALUE = ' Noise type  Current:')
now = inoise[scr0]
current_noise = WIDGET_LABEL(base, VALUE =ncode[now],/DYNAMIC_RESIZE)
void = WIDGET_LABEL(base, VALUE = ' New:')
choose_noise = WIDGET_DROPLIST(base, VALUE = ncode, UVALUE = now)
WIDGET_CONTROL, choose_noise, SET_DROPLIST_SELECT = now
;
; Setting unit
;
units = (*state.TABARR).unit
units = units[iscreen]

base = WIDGET_BASE(query, /ROW)
void = WIDGET_LABEL(base, VALUE = ' Brightness unit  Current:')
current_unit = WIDGET_LABEL(base, VALUE = units[scr0], /DYNAMIC_RESIZE)
void = WIDGET_LABEL(base, VALUE = ' New:')
choose_unit =  WIDGET_TEXT(base, VALUE = '', xsize = 10, /EDITABLE)
;
; Finish off widget
;
bbase   = WIDGET_BASE(query, /ROW)
void   = WIDGET_BUTTON(bbase, VALUE = 'Accept')
void   = WIDGET_BUTTON(bbase, VALUE = 'Cancel')

info = {retID: event.ID, tablist: tablist, $
        pollist: tabcode, current_pol: current_pol, choose_pol: choose_pol, $
        freqcode: freqcode, current_freq: current_freq, $
        choose_freq: choose_freq, bmaj: sbmaj, bmin: sbmin, bpa: sbpa, $
        current_bmaj: current_bmaj, current_bmin: current_bmin, $
        current_bpa: current_bpa, choose_bmaj: choose_bmaj, $
        choose_bmin: choose_bmin, choose_bpa: choose_bpa, $
        ncode: ncode, current_noise: current_noise, choose_noise: choose_noise, $
        units: units, current_unit: current_unit, choose_unit: choose_unit}

WIDGET_CONTROL, event.ID, GET_UVALUE = save
WIDGET_CONTROL, query, SET_UVALUE = info
WIDGET_CONTROL, query, /REALIZE
XMANAGER, 'ximview_setprop2', query
WIDGET_CONTROL, event.ID, GET_UVALUE = result
WIDGET_CONTROL, event.ID, SET_UVALUE = save

IF SIZE(result,/TYPE) NE 8 THEN cancel = 1B ELSE cancel = result.CANCEL
IF cancel THEN RETURN

; Interpret results here. Only change values if set.
tab = result.TAB

polcode = result.polcode - 8
oldcode = (*state.TABARR)[iscreen[tab]].POLCODE
IF polcode NE oldcode THEN BEGIN
   (*state.TABARR)[iscreen[tab]].POLCODE = polcode
   MESSAGE, /INFORMATIONAL, 'Setting tab '+tablab[tab]+' polarization to ' $
            + scodes[polcode+8]
ENDIF

oldcode = (*state.TABARR)[iscreen[tab]].noise_type
IF result.ncode NE oldcode THEN $
   (*state.TABARR)[iscreen[tab]].noise_type = result.ncode

freq = result.freq
IF freq NE '' THEN BEGIN
   (*state.TABARR)[iscreen[tab]].FREQCODE = freq
   MESSAGE, /INFORMATIONAL, 'Setting tab ' + tablab[tab] + ' frequency to ' $
  + freq
ENDIF

beam = result.beam
good = FINITE(beam.bmaj) || FINITE(beam.bmin) || FINITE(beam.bpa)
IF good THEN (*state.TABARR)[iscreen[tab]].beam = result.beam

IF result.unit NE '' THEN (*state.TABARR)[iscreen[tab]].unit = result.unit

END
;
PRO ximview_setcoord, event
; Updates coordinate system
; this has knock-on effects for
;   * state.astrom
;   * state.csystem
;   * state.ll_fmt
;   * state.prec
;   * state.title.head
;
; we re-run parse_header to account for all this.
;
COMPILE_OPT IDL2, HIDDEN
ON_ERROR, 1

WIDGET_CONTROL, event.TOP, GET_UVALUE = state
IF state.FIRST THEN RETURN ; havn't started yet.
IF ~state.is_astrom THEN BEGIN
  ok = DIALOG_MESSAGE("Can't set coordinate system without any astrometry", $
                      DIALOG_PARENT = event.TOP, /ERROR)
  RETURN
ENDIF

astrom   = *state.ASTROM
codes = ['G', 'E', 'C', 'S', 'T', 'H', 'X']
ii = WHERE( astrom.COORD_SYS EQ  codes)
IF ii EQ -1 THEN ii = 6

systems = (['Galactic', 'Ecliptic', 'Equatorial', 'Supergalactic', $
  'Terrestrial', 'Helioecliptic', 'unknown'])

ximview_setprop, event, -1, 0, 'Coordinate System', systems[ii], systems, $
  cancel, tab, newsys

IF cancel THEN RETURN

IF newsys NE ii THEN BEGIN
  MESSAGE, /INFORMATIONAL, 'Setting coordinate system to ' + systems[newsys]
  MESSAGE, /INFORMATIONAL, 'Stored headers are NOT updated.'
  lontype = ['GLON','ELON','RA--','SLON','TLON','HLON','XLON']
  lattype = ['GLAT','ELAT','DEC-','SLAT','TLAT','HLAT','XLAT']

  ; Grab a header: doesn't matter which as they all match in projection etc.
  file_str = (*state.FILES)[0]
  header = *file_str.header
  longkey = STRING(astrom.AXES[0], FORMAT="('CTYPE',I1)")
  SXADDPAR, header, longkey, lontype[newsys]+'-'+astrom.PROJECTION
  latkey =   STRING(astrom.AXES[1], FORMAT="('CTYPE',I1)")
  SXADDPAR, header,  latkey, lattype[newsys]+'-'+astrom.PROJECTION

  parse_header, state.IMSIZE, header, 1, 0, state.VERBOSE, $
                astrom, is_astrom, pole_pix, csystem, proj, unit, beam, title, $
                nside, ns4, ll_fmt, prec, boxsize
    
 ; update_astrom, csystem, ll_fmt, prec, title, state
  state.CSYSTEM = csystem 
  state.ASTROM = PTR_NEW(astrom)
  state.POLE_PIX = pole_pix
  state.LL_FMT = ll_fmt
  state.PREC = prec
  state.TITLE = title
  state.STATBOX = boxsize > 33
  astrom = 0
  
  WIDGET_CONTROL, event.TOP, SET_UVALUE = state
  ;   Update units on readout label
  mid = form_unit((*state.TABARR)[0].UNIT)
  
  title_string = state.TITLE.HEAD + mid + state.TITLE.TAIL
  WIDGET_CONTROL, state.READLAB, SET_VALUE=title_string

ENDIF

END

PRO ximview_header, event
; Prints FITS header
;
COMPILE_OPT IDL2, HIDDEN
ON_ERROR, 1

WIDGET_CONTROL, event.ID, GET_UVALUE = incount
WIDGET_CONTROL, event.TOP,  GET_UVALUE = state

file_str = (*state.FILES)[incount - 1]

is_mswin = STRCMP(!version.OS_FAMILY,'Windows', 4,/ FOLD_CASE)
IF is_mswin THEN BEGIN ; Set a fixed-width font:
    old_font = WIDGET_INFO(state.LABEL, /FONTNAME)
    font = 'lucida console*10'
    WIDGET_CONTROL, DEFAULT_FONT = font
ENDIF

title = 'FITS header for '+file_str.FILE
IF file_str.EXTENSION GT 0 THEN title = 'Merged ' + title

XDISPLAYFILE, '', TEXT = *file_str.HEADER, $
  GROUP = event.TOP, DONE_BUTTON = 'Done with FITS header', height = 24, $
  TITLE = title
  
; Restore system default (?)
IF is_mswin THEN WIDGET_CONTROL, DEFAULT_FONT = old_font

END

PRO ximview_help, event
; Launches help windows
;
COMPILE_OPT IDL2, HIDDEN
ON_ERROR, 1

; Find help directory
info = ROUTINE_INFO('ximview', /SOURCE)
dir = STRSPLIT(info.PATH, 'ximview.pro', /REGEX, /EXTRACT) + 'docs/'
unhelp = ['Abandon hope all ye who enter here.', $
          "Lasciate ogne speranza, voi ch'intrate", $
          'Send requests for help to /dev/null', $
          'ZOMG!!! wtf?', 'Please write the help page you need.']
nun = N_ELEMENTS(unhelp)
WIDGET_CONTROL, event.ID, GET_VALUE = label
WIDGET_CONTROL, event.TOP,  GET_UVALUE = state

dir = dir[0]
assistant = 0B
CASE label OF
    'Help': BEGIN
        file = 'help.txt'  &  title = "Don't Panic!"
        done = 'Done with Ximview HELP'
        height = 50
        ver = FLOAT(!version.RELEASE)
        assistant = ver GE 6.2
        book = ver GT 8 ? 'ximview.html' : 'ximview.adp'
;        ok = DIALOG_MESSAGE(unhelp[ FIX(RANDOMU(seed,1)*nun) ], $
;                            DIALOG_PARENT = event.TOP, TITLE = title)
    END
    'Release Notes': BEGIN
        file = 'release_notes.txt' &  title = 'Release Notes'
        done = 'Done with Release notes'  &  height = 24
;        ok = DIALOG_MESSAGE('E E F G G F E D C C D E Ee di D', $
;                            DIALOG_PARENT = event.TOP, TITLE = title)
    END
    'About': BEGIN
        file = 'about.txt' &  title = 'About Ximview v'+ state.version
        done = 'Done with About Ximview'  &  height = 24
;        ok = DIALOG_MESSAGE('Ximview is about 10 cm across', $
;                            DIALOG_PARENT = event.TOP, $
;                            TITLE = title)
    END
    ELSE: MESSAGE, /INFORMATIONAL, 'Option ' + label + ' not yet available.'
ENDCASE

IF assistant THEN BEGIN
    ONLINE_HELP, BOOK=book
ENDIF ELSE BEGIN
    XDISPLAYFILE, dir+file, GROUP = event.TOP, DONE_BUTTON = done, $
      HEIGHT = height, TITLE = title
ENDELSE
END

;------------------------------------------------------------------------------
; Routines handling the most common events, driven by either the widget buttons
; or mouse actions in the draw window, or manipulation of the widget itself:
;
;    pan, zoom, overview mode, tab switching, resizing, etc.
;
PRO ximview_event, event
;
; Processes events from the top level base: resize, timer, keyboard focus.
; Also swallows events from the draw window which are actually
; processed by ximview_scroll and possibly CW_DEFROI.
;
COMPILE_OPT IDL2, HIDDEN
ON_ERROR, 1

WIDGET_CONTROL, event.TOP, GET_UVALUE = state

name = TAG_NAMES(event, /STRUCTURE_NAME)
CASE name OF
   'WIDGET_BUTTON': BEGIN
      WIDGET_CONTROL, event.ID, GET_VALUE = label
      CASE label OF
         'Restore default screen size': BEGIN
            state.NEWWIDTH  = -1
            state.NEWHEIGHT = -1
                                ; Disable events from the draw window
                                ; during resize:
            WIDGET_CONTROL, state.TABS, SENSITIVE = 0
            WIDGET_CONTROL, event.TOP, SET_UVALUE = state
            ximview_resize, event.TOP, state
         END
         ELSE:
      ENDCASE
   END
   '': BEGIN                    ; Probably menu event
      WIDGET_CONTROL, event.ID, GET_VALUE = label
      option = event.VALUE
      CASE option OF
         'Mark point.Middle mouse button': BEGIN
            state.mark_opt = 0B
         END
         'Mark point.Right mouse button': BEGIN
            state.mark_opt = 1B
         END
         'Mark point.Shift-click': BEGIN
            state.mark_opt = 2B
         END
         'Zoom in/out.Mouse wheel': BEGIN
            state.zoom_opt = 0B
         END
         'Zoom in/out.Ctrl-/ctrl-shift- click': BEGIN
            state.zoom_opt = 1B
         END
         ELSE: MESSAGE,/INFORMATIONAL, 'Option ' + option + ' from menu ' $
                       + label +  ' is not yet implemented.'
      ENDCASE
      WIDGET_CONTROL, event.TOP, SET_UVALUE = state
   END
   'WIDGET_BASE': BEGIN         ; Resize event.

; Record current size of tlb. Don't update immediately as many resize
; events are produced when moving the edge of a window.
      state.NEWWIDTH  = event.X
      state.NEWHEIGHT = event.Y
                                ; Disable events from the draw window
                                ; during resize:
      WIDGET_CONTROL, state.TABS, SENSITIVE = 0
      WIDGET_CONTROL, event.TOP, SET_UVALUE = state
                                ; Send off timer event. Last one to
                                ; arrive is acted on (!).
      WIDGET_CONTROL, event.TOP, TIMER = 0.5
   END
   'WIDGET_TIMER': BEGIN        ;... and act here.
      ximview_resize, event.TOP, state
   END
   'WIDGET_KBRD_FOCUS' : IF event.ENTER THEN BEGIN
                                ; When window becomes active, set
                                ; graphics state to what we want and
                                ; save previous values.

; Return immediately if somehow we already have focus:
      IF state.FOCUS THEN RETURN
      state.FOCUS = 1
      WIDGET_CONTROL, event.TOP, SET_UVALUE = state
   ENDIF ELSE BEGIN             ; Lost keyboard focus
      IF ~state.FOCUS THEN RETURN
      state.FOCUS = 0
      WIDGET_CONTROL, event.TOP, SET_UVALUE = state
   ENDELSE
   ELSE: MESSAGE, /INFORMATIONAL, 'Unknown event type received: '+name
ENDCASE

END

PRO ximview_overview, event
;
; Makes overview plot on request from "overview" button.
;
COMPILE_OPT IDL2, HIDDEN
COMMON gr_global
ON_ERROR, 1

WIDGET_CONTROL, event.TOP,   GET_UVALUE = state

ximview_resize, event.top, state ; Checks the window has not been re-sized

tabarr = *state.tabarr
ntab = N_ELEMENTS(tabarr)
WIDGET_CONTROL, state.TABS,  GET_UVALUE = mode

; Switch to overview mode from zoom mode:
mode.OVERVIEW = 1
WIDGET_CONTROL, state.ZOOMCOL, SENSITIVE = 0
WIDGET_CONTROL, state.READOUT, SENSITIVE = 0

; Set cursor on marked point if available, otherwise centre of FOV.
IF mode.XPT GE 0 THEN BEGIN
    mode.XPIX = mode.XPT  &  mode.YPIX = mode.YPT
ENDIF ELSE BEGIN
    mode.XPIX = mode.X_CENTRE  &  mode.YPIX = mode.Y_CENTRE
ENDELSE

swap_lut, *state.XIM_GRAPH, dummy, old_graph

FOR itab = 0,ntab-1 DO BEGIN
    WSET, tabarr[itab].WINDOW

    IF redraw_req THEN BEGIN
        DEVICE, DECOMPOSED=0
        lutptr = tabarr[itab].LUT
        TVLCT, (*lutptr).R, (*lutptr).G, (*lutptr).B
        !P.background = (*lutptr).ABSENT
        !P.color      = (*lutptr).LINE
    ENDIF
    null3 = REPLICATE(PTR_NEW(),3)
    IF ARRAY_EQUAL(tabarr[itab].RGB, null3) THEN tptr = tabarr[itab].BYTE_PTR $
    ELSE tptr = tabarr[itab].RGB
    overview, tptr, mode, resamp, corner
ENDFOR
mode.RESAMP = resamp  &  mode.CORNER = corner
mode.NEW_VIEW = 1

restore_lut, dummy, old_graph

WIDGET_CONTROL, state.TABS,  SET_UVALUE = mode

END

FUNCTION ximview_zoom, event
;
; Changes zoom factor & redraws view as appropriate
;
COMPILE_OPT IDL2, HIDDEN
ON_ERROR, 1

start = SYSTIME(1)

WIDGET_CONTROL, event.TOP,   GET_UVALUE = state

ximview_resize, event.TOP, state ; Just checking (usually)

WIDGET_CONTROL, state.TABS,  GET_UVALUE = mode

tabarr = state.TABARR
zoom = mode.ZOOM
WIDGET_CONTROL, event.ID, GET_UVALUE = uval
CASE uval OF
    'IN' :  zoom += 1          ; zoom in
    '1:1' : zoom  = 0          ; zoom factor 1
    'OUT' : zoom -= 1          ; zoom out
    ELSE : BEGIN ; event sent here from draw widget: mousewheel zoom
       IF event.VALUE EQ 'IN'  THEN zoom += 1 ELSE $
       IF event.VALUE EQ 'OUT' THEN zoom -= 1 
    END
ENDCASE
zoom_factor = 2.^zoom

swap_lut, *state.XIM_GRAPH, (*tabarr)[0], old_graph

xhalf = mode.XHALF     &  yhalf = mode.YHALF
xcent = mode.X_CENTRE  &  ycent = mode.Y_CENTRE
ierr = 0
coord = gscroll(temp, xcent, ycent, xhalf, yhalf, $
                zoom_factor, ierr, 1, done, DO_WRAP = state.ROLL)

lowerr = ierr MOD 8
IF lowerr EQ 6 || lowerr EQ 7 THEN BEGIN
; Can't zoom any more.
    button = zoom LT 0 ? state.zoomout : state.zoomin
    WIDGET_CONTROL, button, SENSITIVE = 0
    state.MAXZOOM = button
    WIDGET_CONTROL, event.TOP, SET_UVALUE = state
    IF lowerr EQ 7 THEN BEGIN
; Max zoom exceeded: may happen from mousewheel zoom 
; Ignore & return
        RETURN, ierr
    ENDIF
ENDIF ELSE IF state.MAXZOOM NE 0 THEN BEGIN ; restore zoom option
    WIDGET_CONTROL, state.MAXZOOM, /SENSITIVE
    state.MAXZOOM = 0
    WIDGET_CONTROL, event.TOP, SET_UVALUE = state
ENDIF

IF lowerr NE 0 && lowerr LT 6 THEN $
  MESSAGE, 'GSCROLL error '+STRING(ierr)

; rolling disabled if ierr = 8
mode.ROLL = ierr-lowerr NE 8 ? state.ROLL : 0

mode.ZOOM = zoom
zfac = 2^ABS(zoom)
mode.ZFAC = zfac
mode.ZOOM_FACTOR = zoom_factor
mode.XHALF = xhalf  &  mode.YHALF = yhalf

                                ; tell overlay to recalculate grid step
                                ; for new window size 
IF mode.GRID_CALC THEN mode.GRID_STEP = [-1,-1]

                                ; Re-plot any graphics overlays:
overlay, mode, state.astrom, state.pole_pix, state.nside

restore_lut, dummy, old_graph

                                ; Disable panning if image fits in screen
bltv = tv2im(0, 0, mode)
trtv = tv2im(!D.x_vsize, !D.y_vsize, mode)

mode.PAN =  MAX(bltv) GT 0 OR MAX(state.IMSIZE[1:2] - trtv) GT 1

mode.DONE = done
; Request another go if loading not finished
IF done EQ 0 THEN WIDGET_CONTROL, (*tabarr)[0].DRAW, TIMER = 0.5

WIDGET_CONTROL, state.TABS,  SET_UVALUE = mode

                                ; Record updated zoom state
pix_print, state, 0, start

RETURN, ierr
END

PRO ximview_tab, event
; Processes tab change events and also swallows events from draw window
;
COMPILE_OPT IDL2, HIDDEN
COMMON gr_global
ON_ERROR, 1

start = SYSTIME(1)

WIDGET_CONTROL, event.TOP, GET_UVALUE = state
WIDGET_CONTROL, state.TABS,  GET_UVALUE = mode
xpix = mode.XPIX  &  ypix = mode.YPIX
zoomf = mode.ZOOM_FACTOR

name = TAG_NAMES(event, /STRUCTURE_NAME)
SWITCH name OF
    'WIDGET_TIMER':          ; Do nothing, processed by ximview_scroll
    'WIDGET_DRAW': BREAK     ; Do nothing, processed by ximview_scroll
    'WIDGET_TAB':            ; Actual tab event or..
    '' :  BEGIN              ; Anonymous event from ximview_blink
        tabarr = state.TABARR
                                ; Check we don't have a temporary tab:
        IF event.TAB GT MAX((*tabarr).SCREEN) THEN RETURN

                                ;   Check for bizarre counting bug:
        current = event.TAB EQ 1 ? WIDGET_INFO(state.TABS, /TAB_CURRENT) $
                                 : event.TAB

        swap_lut, *state.XIM_GRAPH, (*tabarr)[0], old_graph

        gscroll_newscreen, current, (*tabarr), zoomf, $
          mode.X_CENTRE, mode.Y_CENTRE, mode.XHALF, mode.YHALF, done, $
          mode.OVERVIEW

        draw = (*tabarr)[0].DRAW
        mode.DONE = done
        IF done EQ 0B THEN WIDGET_CONTROL, draw, TIMER = 0.5

                                ; Disable analysis tools for RGB
                                ; window
        sens = PTR_VALID((*tabarr)[0].BYTE_PTR)
        WIDGET_CONTROL, state.IMSTATS, SENSITIVE = sens
        WIDGET_CONTROL, state.PEAKFIT, SENSITIVE = sens
        WIDGET_CONTROL, state.PROFILE, SENSITIVE = sens

        IF redraw_req THEN DEVICE, DECOMPOSED = 0 $
                      ELSE DEVICE, DECOMPOSED = (*tabarr)[0].DECOMPOSED

                                ;   Update units on readout label
        mid = form_unit( (*tabarr)[0].UNIT)

        title_string = state.TITLE.HEAD + mid + state.TITLE.TAIL
        WIDGET_CONTROL, state.READLAB, SET_VALUE=title_string

                                ;  Mark current point
        IF ~mode.OVERVIEW THEN BEGIN
           overlay, mode, state.ASTROM, state.pole_pix, state.NSIDE

            pix_print, state, 0, start
        ENDIF

        restore_lut, dummy, old_graph

        BREAK
    END
    ELSE: MESSAGE, /INFORMATIONAL, 'Unknown event type received: '+name
ENDSWITCH
WIDGET_CONTROL, state.TABS,  SET_UVALUE = mode

END

FUNCTION ximview_scroll, event
; Processes events from main draw widget
; Jobs:
;  (0) Overview or zoom mode?
;      If overview, set xpix, ypix on button press and enter zoom
;      mode.
;  Otherwise:
;  (1) Is button 1 down or not?
;      Yes: drag image using gscroll subroutines
;      No:  Record current cursor position
;  (2) Has button 2 been clicked (act on button down)
;      If so, mark spot and record in log
;
COMPILE_OPT IDL2, HIDDEN
COMMON gr_global
ON_ERROR, 1

start = SYSTIME(1)

returnable = event

WIDGET_CONTROL, event.TOP,  GET_UVALUE = state
WIDGET_CONTROL, state.TABS, GET_UVALUE = mode
tabarr = state.TABARR
str = (*tabarr)[0]

swap_lut, *state.XIM_GRAPH, str, oldgraph

xhalf = mode.XHALF  &  yhalf = mode.YHALF

name = TAG_NAMES(event, /STRUCTURE_NAME)
CASE name OF
    'WIDGET_TIMER': BEGIN     ; gscroll didn't finish loading panel
        IF ~mode.DONE THEN BEGIN
;
; Check that we have graphics focus and if not, grab it for a bit
;
            IF ~state.FOCUS THEN BEGIN
                new_view = 1 ; We might as well take our time...
            ENDIF ELSE new_view = 0

            coord = gscroll(void, mode.X_CENTRE, mode.Y_CENTRE, xhalf, yhalf, $
                            mode.ZOOM_FACTOR, ierr, new_view, done, $
                            DO_WRAP = mode.ROLL)

            IF ierr NE 0 THEN MESSAGE, 'GSCROLL error '+STRING(ierr)
            mode.DONE = done
        ENDIF ELSE done = 1
        returnable = 0 ; don't pass this up the chain.
    END
    'WIDGET_DRAW': BEGIN
;
        xtv = event.X  & ytv = event.Y

        CASE state.mark_opt OF
           0: mark = event.press EQ 2 ; middle mouse button
           1: mark = event.press EQ 4 ; right mouse button
           2: mark = event.press && event.modifiers EQ 1 ; shift-click
        ENDCASE
        CASE state.zoom_opt OF
           0: do_zoom = event.type EQ 7  ; mouse-wheel event
           1: do_zoom = event.press && (event.modifiers AND 2B) ; ctrl-click
        ENDCASE
        xpmax = state.IMSIZE[1]-1
        ypmax = state.IMSIZE[2]-1
        zfac  = mode.ZFAC  &  zoom = mode.ZOOM


        IF ~mode.OVERVIEW THEN BEGIN ; We are in pan/zoom mode
            IF do_zoom THEN BEGIN 
                                     ; mousewheel zoom: pass to ximview_zoom
               restore_lut, dummy, oldgraph
               in_or_out = state.zoom_opt ? event.modifiers EQ 2 $
                                          : event.clicks > 0 
               val = in_or_out ? 'IN' : 'OUT'
               newevent = {ID: 0L, TOP: 0L, HANDLER: 0L, VALUE: val}
               ; Why not call ximview_zoom directly here?
               ; WIDGET_CONTROL fills in ID, TOP, HANDLER.
               WIDGET_CONTROL, state.zoomcol, SEND_EVENT=newevent
               RETURN, 0
            ENDIF

            x_centre = mode.X_CENTRE  &  y_centre = mode.Y_CENTRE

            drag = mode.PAN AND (event.PRESS OR (mode.DRAG && ~event.RELEASE))
            mode.DRAG = drag

            IF drag THEN BEGIN
                               ; Set cursor to gripping hand
                if ~is_gdl() THEN cursor_grip
;               DEVICE, CURSOR_STANDARD=52

; Calculate new centre pixel = shifted by opposite amount from cursor
;
                oxtv = mode.OXTV  &  oytv = mode.OYTV
;  Get position of pixel under cursor
                coord = tv2im(oxtv, oytv, mode)

                IF zoom LE 0 THEN BEGIN
                    xshift = xtv - oxtv
                    yshift = ytv - oytv
                    xpix = ((x_centre - xshift*zfac) > 0) < xpmax
                    ypix = ((y_centre - yshift*zfac) > 0) < ypmax
                ENDIF ELSE BEGIN
                    xpix = ((x_centre - (xtv/zfac) + (oxtv/zfac)) > 0) < xpmax
                    ypix = ((y_centre - (ytv/zfac) + (oytv/zfac)) > 0) < ypmax
                ENDELSE
                IF mode.ROLL THEN BEGIN
                    ns4 = state.NS4
                    IF xpix + ypix GT (state.IMSIZE[1] + ns4) THEN BEGIN
                        xpix -= ns4 & ypix -= ns4
                    ENDIF ELSE IF xpix + ypix LT state.NSIDE THEN BEGIN
                        xpix += ns4 & ypix += ns4
                    ENDIF
                ENDIF

                ierr = 0
                null = gscroll(void, xpix, ypix, xhalf, yhalf, $
                        mode.ZOOM_FACTOR, ierr, 0, done, DO_WRAP=mode.ROLL)
                IF ierr NE 0 THEN MESSAGE, 'GSCROLL error '+STRING(ierr)
                mode.DONE = done
                mode.X_CENTRE = xpix
                mode.Y_CENTRE = ypix

                xpix = coord[0]  &  ypix = coord[1]
                IF mode.ROLL THEN BEGIN
                    IF xpix+ypix GT (state.IMSIZE[1] + ns4) THEN BEGIN
                        xpix -= ns4 & ypix -= ns4
                    ENDIF ELSE IF xpix + ypix LT state.NSIDE THEN BEGIN
                        xpix += ns4 & ypix += ns4
                    ENDIF
                ENDIF

                xpix = (xpix > 0) < xpmax
                ypix = (ypix > 0) < ypmax
            ENDIF ELSE BEGIN
                IF ~is_gdl() THEN DEVICE, /CURSOR_CROSSHAIR

                coord = tv2im(xtv, ytv, mode)
                xpix = coord[0]  &  ypix = coord[1]
                IF mode.ROLL THEN BEGIN
                    ns4 = state.NS4
                    IF xpix+ypix GT (state.IMSIZE[1] + ns4) THEN BEGIN
                        xpix -= ns4 & ypix -= ns4
                    ENDIF ELSE IF xpix + ypix LT state.NSIDE THEN BEGIN
                        xpix += ns4 & ypix += ns4
                    ENDIF
                ENDIF
                xpix = (xpix > 0) < xpmax
                ypix = (ypix > 0) < ypmax

                done = 1 ; No gscroll update therefore done by definition!
            ENDELSE

            IF mark THEN BEGIN  ; mark spot
               mode.XPT2 = mode.XPT1 & mode.YPT2 = mode.YPT1
               mode.XPT1 = mode.XPT  & mode.YPT1 = mode.YPT
               mode.XPT = xpix 
               mode.YPT = ypix
                                ; Notify any program running about click
;              IF state.prog NE '' THEN WIDGET_CONTROL, state.progcol, $
;                 SEND_EVENT = {PRESS, ID: 0, TOP: 0, HANDLER: 0}  
            ENDIF

            mode.OXTV = xtv   &  mode.OYTV = ytv
            mode.XPIX = xpix  &  mode.YPIX = ypix

; End of code for zoom mode
        ENDIF ELSE IF event.PRESS THEN BEGIN ; Click in Overview mode...
            mode.OVERVIEW = 0   ; Switch to zoom mode...
            WIDGET_CONTROL, state.ZOOMCOL, /SENSITIVE
            WIDGET_CONTROL, state.READOUT, /SENSITIVE

            IF mode.RESAMP[0] GT 1 THEN BEGIN
                resamp = ROUND(mode.RESAMP)
                xpix = ((xtv - mode.CORNER[0])*resamp[0] > 0) < xpmax
                ypix = ((ytv - mode.CORNER[1])*resamp[1] > 0) < ypmax
            ENDIF ELSE BEGIN
                resamp = ROUND(1.0/mode.RESAMP)
                xpix = ((xtv - mode.CORNER[0])/resamp[0] > 0) < xpmax
                ypix = ((ytv - mode.CORNER[1])/resamp[1] > 0) < ypmax
            ENDELSE

            ierr = 0
            coord = gscroll(void, xpix, ypix, xhalf, yhalf, mode.ZOOM_FACTOR, $
                            ierr, mode.NEW_VIEW, done, DO_WRAP=state.ROLL)

            lowerr = ierr MOD 8
            IF lowerr EQ 6 THEN BEGIN ; Tiny map: max zoom already.
                in = WIDGET_INFO(state.ZOOMCOL, FIND_BY_UNAME = 'IN')
                WIDGET_CONTROL, in, SENSITIVE = 0
                state.MAXZOOM = in
                WIDGET_CONTROL, event.TOP, SET_UVALUE = state
            ENDIF

            IF lowerr NE 0 && lowerr NE 6 THEN $
              MESSAGE, 'GSCROLL error '+STRING(ierr)

; rolling disabled if ierr = 8
            mode.ROLL = ierr NE 8 ? state.ROLL : 0
            mode.DONE = done
            mode.OXTV  = xhalf    &  mode.OYTV  = yhalf
            mode.XHALF = xhalf    &  mode.YHALF = yhalf
            mode.X_CENTRE = xpix  &  mode.Y_CENTRE = ypix
            mode.XPIX = xpix      &  mode.YPIX = ypix
            mode.NEW_VIEW = 0
                                ; Disable panning if image fits in screen
            bltv = tv2im(0, 0, mode)
; xpix, ypix, xhalf, yhalf, zoom, zfac)
            trtv = tv2im(!D.x_vsize, !D.y_vsize, mode)

            mode.PAN =  MAX(bltv) GT 0 OR MAX(state.IMSIZE[1:2] - trtv) GT 1

            mark = 0

            TVCRS, mode.oxtv, mode.oytv

        ENDIF ELSE BEGIN
            restore_lut, dummy, oldgraph
            RETURN, 0  ; overview mode, no click yet.
        ENDELSE
        
        IF mark && str.SCREEN NE state.LASTTAB THEN BEGIN
            WIDGET_CONTROL, str.BASE, GET_UVALUE = tname
            tname = 'Now on tab: ' + tname
            PRINT, tname
            PRINTF, state.LOGLUN, tname
            state.LASTTAB = str.SCREEN
            WIDGET_CONTROL, event.TOP, SET_UVALUE = state
        ENDIF

        pix_print, state, mark, start

    END                         ; End of "WIDGET_DRAW" case
    ELSE: MESSAGE, /INFORMATIONAL, 'Unknown event type received: '+name
ENDCASE
                                ; Plot graphics overlay including
                                ; grid, catalogs, and marked points
overlay, mode, state.astrom, state.pole_pix, state.NSIDE

restore_lut, dummy, oldgraph

; If panels remain to be loaded, send timer event to DRAW window
; requesting re-draw:
IF done EQ 0 THEN WIDGET_CONTROL, event.ID, TIMER = 0.5

WIDGET_CONTROL, state.TABS, SET_UVALUE = mode

RETURN, returnable
END
;
PRO start_gscroll,  ntab, state, imsize, noldtab, mode, start
;
; Enables the gscroll system or adds more screens to it.
; Sets up mode if this is the first load.
;
COMPILE_OPT IDL2, HIDDEN
COMMON gr_global

verbose = state.VERBOSE
first = state.FIRST
tab_arr = *state.TABARR

IF first THEN BEGIN
    xpanel = 128L               ; default size of panel on pixmap.
    maxzoom = xpanel / 4L
                                ; Initialize gscroll common and pixmap windows:
    ierr = 0
    gscroll_setup, ntab, imsize, maxwin, ierr, WINDOW = tab_arr.WINDOW, $
      IMAGE = tab_arr.BYTE_PTR, LUT = tab_arr.LUT, RGB = tab_arr.RGB, $
       REDRAW = redraw_req, /HIDDEN
    IF ierr NE 0 THEN MESSAGE, 'GSCROLL_SETUP error ' + STRING(ierr)
    state.MAXWIN = maxwin

    IF verbose THEN MESSAGE,/INFORMATIONAL, $
      STRING(SYSTIME(1)-start, "Scrolling set up", $
             FORMAT = "(F7.3,' seconds: ', A)")
    IF verbose THEN HELP, /MEMORY

; Set up for initial overview:
    xpix = imsize[1]/2  &  ypix = imsize[2]/2
    x_centre = xpix  &  y_centre = ypix

; Set initial zoom: largest that gets whole map on screen, or 1,
; whichever is larger.
    WSET, tab_arr[0].WINDOW
    ratio = (!D.x_vsize/imsize[1]) < (!D.y_vsize/imsize[2])
    ratio = ratio < maxzoom
    zoom = ratio GT 0 ? FIX(ALOG(ratio) / ALOG(2.0)) : 0
    zoom_factor = 2.^zoom

; Find effective central pixel on view window:
    get_centre, zoom_factor, xhalf, yhalf

    WIDGET_CONTROL, state.TABS, GET_UVALUE = mode

    mode.ZOOM     = zoom      &  mode.ZOOM_FACTOR = zoom_factor
    mode.ZFAC     = 2^ABS(zoom)
    mode.OXTV     = xhalf     &  mode.OYTV     = yhalf
    mode.XPIX     = xpix      &  mode.YPIX     = ypix
    mode.XHALF    = xhalf     &  mode.YHALF    = yhalf
    mode.X_CENTRE = x_centre  &  mode.Y_CENTRE = y_centre

ENDIF ELSE FOR icol = 0,ntab-1 DO BEGIN
                                ; Update new low-level structures:
    itab = icol + noldtab
    str = tab_arr[itab]
    WSET, str.WINDOW
    gscroll_addscreen, str.BYTE_PTR, str.LUT, str.RGB 
ENDFOR

*state.TABARR = tab_arr

END

PRO fill_screens, ntab, noldtab, tabarr, mode, state, start
; Fills DRAW windows of newly-created tabs
;
COMPILE_OPT IDL2, HIDDEN
COMMON gr_global
ON_ERROR, 2
first = state.FIRST

FOR i=0,ntab-1 DO BEGIN
    itab = i + noldtab
    IF redraw_req THEN BEGIN
        lutptr = (*tabarr)[itab].LUT
        TVLCT, (*lutptr).R, (*lutptr).G, (*lutptr).B
        !P.background = (*lutptr).ABSENT
        !P.color      = (*lutptr).LINE
    ENDIF
    update_screen, tabarr, itab, mode, done
ENDFOR

IF ~first THEN BEGIN
                                ; Switch to newly-loaded screen
    WIDGET_CONTROL, state.TABS, SET_TAB_CURRENT = noldtab

    gscroll_newscreen, noldtab, *tabarr, mode.ZOOM_FACTOR, $
      mode.X_CENTRE, mode.Y_CENTRE, mode.XHALF, mode.YHALF, done, mode.OVERVIEW
    
    mode.DONE = done
                                ; Request another go if loading not finished
    IF done EQ 0B THEN WIDGET_CONTROL, (*tabarr)[0].DRAW, TIMER = 0.5

                                ;  Update overlays etc
    IF ~mode.OVERVIEW THEN BEGIN
       overlay, mode, state.ASTROM, state.pole_pix, state.NSIDE
       pix_print, state, 0, start
    ENDIF
ENDIF

END
;
PRO ximview, input, range, proj, order, COLUMN = column, EXTENSION = exten, $
             WRAP = wrap, ROLL = roll, NAME = name, MACSTYLE = macstyle, $
             NPOLE = npole, SPOLE = spole, RING = ring, NESTED=nest, $
             LOG = log, TEMPORARY = temporary, VERBOSE = verbose, $
             HELP = help, OVERWRITE = Overwrite, EXPERT = expert
;+
; NAME:
;       XIMVIEW
;
; PURPOSE:
;       Inspection and basic analysis tool for large images including
;       those stored in HEALPix format. Supports FITS WCS coordinate
;       and unit descriptions. Runs transparently over X-connections
;       (but requires a fast link to operate effectively).
;
;       Multi-dimensional images are displayed on a set of tabbed
;       screens, which support blinking etc.
;
; CATEGORY:
;       Widgets.
;
; CALLING SEQUENCE:
;
;       XIMVIEW  [/VERBOSE]
;                (data can be loaded from files via menu options)
;   or
;       XIMVIEW, Input, [Range, COL=columns, EXT=ext, /OVER, /MAC, $
;                        /WRAP, /ROLL, NAME=name, /TEMP, /VERBOSE] $
;                        [proj | /NPOLE | /SPOLE ] [order | /RING | /NESTED]
;                (options on the last line apply to HEALPix data only).
;   or  XIMVIEW, /HELP
;                (prints list of command-line options)
;
; INPUTS:
;       Input:  Any of the following:
;               o Image array (2 or more dimensions)
;               o HEALPix array (1 or more dimension with the first
;                 being a valid HEALPix size)
;               o Structure containing a FITS header as tag (0),
;                 plus a 2+D image as tag (1), or a set of 1D
;                 HEALPix arrays as tags (1), (2),..., or a set of
;                 pointers to 1D or 2D arrays as tags (1), (2),...
;                 Arrays on all tags must be the same size.
;                 CUT4 format arrays are accepted.
;               o Name of a FITS file containing an image (in the
;                 primary HDU or in an IMAGE extension) or a set of
;                 HEALPix arrays including CUT4 format (in a binary
;                 table extension). The ".fits" or ".FITS" extensions
;                 may be omitted.
;               o A single pointer to any of the above.
;               o An array of pointers to 1D HEALPix or 2D image
;                 arrays.
;
; OPTIONAL INPUTS:
;       Range:  Sets the range of image intensities to map to the
;               displayed colour table. Options are:
;               o  Scalar: maximum image intensity to plot. The
;                   minimum to plot is the minimum in the data.
;               o  [minimum, maximum]
;               o  [minimum, maximum, beta] uses ASINH transfer function
;               o 'AUTO' or 'A' or '*': auto-scale image based on mode and
;                 robust estimate of standard deviation (separately
;                 for each channel).
;               o 'Full' or 'F': full range of image.
;               Default is Auto.
;
;       Proj:   'GRID', 'NPOLE', or 'SPOLE' if HEALPix array is
;               supplied. Default: 'GRID' = HPX projection.
;
;       Order:  'NESTED' or 'RING' if HEALPix array is
;               supplied. Default: header value if any, otherwise
;               'RING'
;
; KEYWORD PARAMETERS:
;       NPOLE:  Set for 'NPOLE' projection (alternative to Proj input)
;
;       SPOLE:  Set for 'SPOLE' projection (alternative to Proj input)
;
;       RING:   Set for RING order   (alternative to Order input)
;
;       NESTED: Set for NESTED order (alternative to Order input)
;
;       COLUMN: Single value or array of:
;               o Numbers or names of the binary table column to read,
;               o the plane index (coordinate on 3rd dimension) if
;                 the input is a 3-D array. (Dimensions higher than 3
;                 are "collapsed" into the 3rd dimension).
;               o Coordinate on 2nd dimension if the input is an array
;                 of (1-D) HEALPix arrays
;               If not specified, all the columns/planes/arrays are
;               read in.
;
;       EXTENSION:  FITS extension containing data to read (primary
;               HDU is considered extension 0). Default: extension 0 if
;               it contains data, otherwise extension 1.
;
;       WRAP:   = 0 (unset): intensities outside Range are displayed
;                            as min or max colour, as appropriate.
;               < 0: Intensities greater than max set by Range use
;                    "wrapped" colours, starting again from the bottom
;                    of the colour table. Intensities less than the
;                    min set by Range are displayed as the min colour.
;                >0: colour table is wrapped at both ends;
;
;       ROLL:   Enable rolling of image through +/-180 longitude
;               (enabled automatically if program determines the image
;               is a full-sphere HPX projection).
;
;       NAME:   Title for image.
;               Default: constructs one from filename and/or header.
;
;       LOG:    Set (/LOG or LOG=1) to give logfile a unique name of
;               the form "ximview_n.log" where n is an integer. (n=m+1
;               where "ximview_m.log" is the file in the default
;               directory with the largest index m).
;               OR: set = <string>, the name of the logfile.
;
;       OVERWRITE: if set, replaces existing image with new one,
;               allowing change of image size, astrometry, etc.
;
;       MACSTYLE: If set, zoom by ctrl-click [in]/ctrl-shift-click
;                 [out], and mark point by shift-click.
;
;       TEMPORARY: if set, allows the input array to be overwritten to 
;               save space.
;
;       VERBOSE: if set, prints diagnostics
;
;       EXPERT: sets the !ximview_expert system parameter to 1, which
;               suppresses some information and obvious warning
;               messages. Once set, this lasts for the rest of the IDL
;               session. 
;
; OUTPUTS:
;       Logfile: see LOG above.
;
; COMMON BLOCKS:
;       GR_GLOBAL:    Contains global graphics state parameters
;       GRID_GSCROLL: Used by low-level GSCROLL library.
;       XY2PIX:       HEALPix low-level library. Only used when
;                     displaying HEALPix images.
;
; SIDE EFFECTS:
;       If your plot device Visual Class is "DirectColor" a private
;       colour map is (usually) used when the cursor is in Ximview's
;       image display area, causing all elements on your screen to
;       change colour while the cursor remains there. To avoid such
;       "flashing", set DEVICE, TRUE_COLOR=24 at the start of your
;       session.
;
;       A similar effect may occur for Visual Class "PseudoColor"
;       except that colours will change only when the widget gains or
;       loses "keyboard focus" (i.e. it becomes the active window).
;
; RESTRICTIONS:
;       Only one instance of Ximview is allowed at a time.
;
;       XIMVIEW uses the HEALPix IDL library for displaying HEALPix
;       images. For best behaviour the XIMVIEW directories should
;       occur before the HEALPix ones.
;
; PROCEDURE:
;       Ximview prints some basic instructions when it starts up.
;       Detailed instructions are available from the HELP button.
;
;       For a description of HEALPix data, and to download the
;       software see:
;       o  http://healpix.jpl.nasa.gov/
;       o  Gorski et al. 2005, Astrophysical Journal, vol 622, p. 759
;
;       HEALPix datasets are converted to 2D arrays in one of the
;       projections described by Calabretta and Roukema (2007, Monthly
;       Notices of the Royal Astronomical Society, vol 381, p. 865.)
;
; EXAMPLE:
;       Display test-card image, scaled min (=0) to 100:
;
;               XIMVIEW, DIST(200,200), 100
;
;       To display the first three HEALPix arrays stored in the first
;       extension to file 'wmap_band_iqumap_r9_3yr_K_v2.fits'
;       (available from http://lambda.gsfc.nasa.gov/), with auto-scaling:
;
;               XIMVIEW, 'wmap_band_iqumap_r9_3yr_K_v2', '*', COL=[1,2,3]
;
; MODIFICATION HISTORY:
;       Written by:     J. P. Leahy, Jan-Feb 2008 (to v.0.3)
;       March 2008      v0.4: added RGB and HSV display, CUT4 files,
;                             better colour handling, numerous minor
;                             improvements (see Release Notes).
;       April 2008      v0.4.2: Bug fixes
;       July  2008      v0.5: Bug fixes
;       November 2008   v0.6: Bug fix (HP2HPX), added scale to PNG output.
;       August 2009     v0.6.1: Bug fix in grid2hp_index & related
;                               progs.
;       June 2013       v0.7: consolidated various bug fixes.
;       August 2013     v0.8: launch without input, scroll-wheel zoom,
;                             revised LUT handling, /HELP option, many
;                             bug fixes, deviant IDL astrolib routines
;                             merged into main distribution, files/header
;                             array introduced, changed tab labelling,
;                             sorted subroutines in this file.
;      September 2013  v0.8.1 Bug fixes; started to implement scroll
;                             bars on choice matrices in popup
;                             windows.
;      June 2014       v0.8.2 Bug fixes to numunit & set_print_fmt
;      July 2016       v0.9   Some updates to tab labelling, part
;                             implementation of draw_grid and VO,
;                             catalogue plotting, set_coord,
;                             programmable buttons, ability to delete
;                             last tab without closing program.
;      Sept-Oct 2020   v0.9.1 coordinate grids done, multiple overlay
;                             colours, non-fatal messages moved to popup
;                             windows, numerous tweaks and bug fixes.
;-
COMPILE_OPT IDL2

; Global parameters describing graphics state:
COMMON gr_global, windev, redraw_req, colmap, badcol, syscol, ibot, itop, $
   m2_col, m3_col, grid_col, lab_col, geo_col

start = SYSTIME(1)
version = '0.9.1'

IF KEYWORD_SET(help) THEN BEGIN
    PRINT, 'Ximview version '+version
    PRINT, ''
    PRINT, 'Syntax:'
    PRINT, '    XIMVIEW, [/Verbose]'
    PRINT, '      --- Launches GUI, from where data files can be loaded.'
    PRINT, '          Other keyword options (below) can be specified, but'
    PRINT, '          can also but set via the GUI.'
    PRINT, ''
    PRINT, '    XIMVIEW, input, [range, Column=, Extension=, NAme=, '
    PRINT, '         {proj | /NPole | /Spole}, {order | /RIng | /NEsted},'
    PRINT, '         /MAcstyle, /Wrap, /ROll, Log=, /Temporary, /Overwrite,'
    PRINT, '         /Verbose ]'
    PRINT, '      --- "input" can be the name of a FITS file, or an IDL array,'
    PRINT, '          structure, pointer, or pointer array.'
    PRINT, ''
    PRINT, '      See Ximview "Help" menu for full documentation.'
    RETURN
ENDIF

temporary = KEYWORD_SET(temporary)
restore_size = 0
ntab = 0

; Are we launching ximview widget now?
launching = XREGISTERED('ximview', /NOSHOW) EQ 0  
first = launching   ; updated below if widget launched but no data loaded.

verbose = KEYWORD_SET(verbose)
IF verbose THEN BEGIN
    ON_ERROR, 0                     ; For debug:
    error_status = 0
    HELP, /MEMORY
ENDIF ELSE BEGIN
    ON_ERROR, 2
    CATCH, error_status
ENDELSE

overwrite = KEYWORD_SET(overwrite)

BAIL:  ; in debug mode, tidy up by setting error_status manually, then GOTO, BAIL

IF error_status NE 0 THEN BEGIN
; Should arrive here if error occurs before XMANAGER is called and
; starts to handle errors
    CATCH, /CANCEL
    HELP, /LAST_MESSAGE

    IF temporary && N_ELEMENTS(data) GT 0 THEN FOR i=0,ntab-1 DO $
      IF PTR_VALID(data[i]) THEN PTR_FREE, data[i]

    IF restore_size THEN $ ; Restore any dummy axes in input
      input = REFORM(input, input_size[1:input_size[0]], /OVERWRITE)

    IF launching THEN BEGIN
        ximview_tidy, state

        IF N_ELEMENTS(tlb) GT 0 THEN WIDGET_CONTROL, tlb, /DESTROY

        MESSAGE, 'Problem initializing Ximview'
    ENDIF ELSE BEGIN            ; Clear up any pointers and pixmaps
                                ; related to new tabs
        tabarr = *state.TABARR
        ntab = N_ELEMENTS(tabarr)
        FOR i=noldtab,ntab-1 DO BEGIN
            id = i
            str = tabarr[i]
            deadtab = str.BASE
            index   = str.SCREEN
            WIDGET_CONTROL, deadtab, /DESTROY
            str.SCREEN = -1
            IF PTR_VALID(str.LUT) THEN PTR_FREE, str.LUT
            IF PTR_VALID(str.BYTE_PTR) THEN PTR_FREE, str.BYTE_PTR
            IF str.TEMPORARY && PTR_VALID(str.IM_PTR) THEN $
              PTR_FREE, str.IM_PTR

            current = WIDGET_INFO(state.TABS, /TAB_CURRENT)
            gscroll_newscreen, current, tabarr, mode.ZOOM_FACTOR, $
              mode.X_CENTRE, mode.Y_CENTRE, mode.XHALF, mode.YHALF, done, 1B
        ENDFOR
        *state.TABARR = tabarr
        MESSAGE, /INFORMATIONAL, 'Problem loading new data'
        WIDGET_CONTROL, state.TABS, /SENSITIVE
    ENDELSE

    RETURN
ENDIF                           ; Catch block ends here

IF ~launching THEN BEGIN
    IF ~!ximview_expert THEN MESSAGE, /INFORMATIONAL, $
        'Data will be added to existing Ximview window'
                                ; Find Ximview
    xim_catch, tlb, state, mode
    current = WIDGET_INFO(state.TABS, /TAB_CURRENT)
ENDIF

fmt = "('XIMVIEW: ',F7.3,' seconds: ',A)"

; Sort out input parameters
nin = N_ELEMENTS(input)
IF nin EQ 0 && ~launching THEN BEGIN
    MESSAGE, /INFORMATIONAL, 'No data to add.'
    RETURN
ENDIF

ring = KEYWORD_SET(ring)
nest = KEYWORD_SET(nest)
order_set = N_ELEMENTS(order) GT 0

CASE order_set + ring + nest OF
    0: ; Do nothing, may not be healpix
    1: IF ~order_set THEN order = nest ? 'NESTED' : 'RING'
    ELSE:  MESSAGE, 'Please specify ordering only once.'
ENDCASE

; Should we auto-scale?
auto_scale = 0B
IF N_ELEMENTS(range) EQ 0 THEN auto_scale = 1B
IF SIZE(range,/TYPE) EQ 7 THEN BEGIN
    rtext = STRMID(STRTRIM(STRUPCASE(range),2),0,1)
    SWITCH rtext OF
        '*':
        'A': BEGIN 
            auto_scale = 1B
            BREAK
        END
        'F': BEGIN
            auto_scale = 0B
            BREAK
        END
        ELSE: MESSAGE, /INFORMATIONAL, 'Option ' + rtext + $
          ' not yet available.'
    ENDSWITCH
ENDIF

; Initial setup: system info, logfile, widget, common:
IF launching THEN BEGIN 
; Get system & directory info:
    is_unix = STRCMP(!version.OS_FAMILY, 'UNIX', 4, /FOLD_CASE)
    windev = is_unix ? 'X' : 'WIN'
    null = is_unix ? '/dev/null' : 'NUL:'
    DEFSYSV, '!ximview_expert', EXIST=is
    IF ~is THEN DEFSYSV, '!ximview_expert', 0B

    IF KEYWORD_SET(expert) THEN !ximview_expert = 1B

    ; Find help directory & add to !help_path if necessary
    info = ROUTINE_INFO('ximview', /SOURCE)
    helpdir = STRSPLIT(info.PATH, 'ximview.pro', /REGEX, /EXTRACT) $
        + 'docs'
    pathsep = path_sep(/SEARCH_PATH)
    dirs = STRSPLIT(!help_path, pathsep, /EXTRACT)
    ndir = N_ELEMENTS(dirs)
    IF dirs[ndir-1] NE helpdir THEN !help_path = !help_path + pathsep + helpdir

; Logfile:
    IF ~KEYWORD_SET(log) THEN logfile = 'ximview.log' ELSE BEGIN
        ltype = SIZE(log,/TYPE) ; ltype = 7 for a string
        logfile = ltype EQ 7 ? log : get_log_name()
    ENDELSE
    OPENW, loglun, logfile, /GET_LUN, ERROR = err
    IF err NE 0 THEN BEGIN
        MESSAGE, /INFORMATIONAL, !ERROR_STATE.MSG
        PRINT, 'Ximview: log file not opened.'
        PRINT, "You probably don't have write access to the current directory."
        PRINT, 'Log information will be discarded.'
        OPENW, loglun, null, /GET_LUN
    ENDIF
    PRINTF, loglun,  version, SYSTIME(), FORMAT = $
      "('XIMVIEW Version ',A,' started at ',A)"

; Inputs involving input data:

    roll = KEYWORD_SET(roll)

    proj_set = N_ELEMENTS(proj) GT 0
    npole = KEYWORD_SET(npole) & spole = KEYWORD_SET(spole)
    IF proj_set + npole + spole GT 1 THEN $
        MESSAGE, 'Please specify projection only once!'
    IF npole THEN proj = 'NPOLE'
    IF spole THEN proj = 'SPOLE'
    IF N_ELEMENTS(proj) EQ 0 THEN proj = ''
    coltab = 0
    noldtab = 0

    wd = FILE_EXPAND_PATH('')
    title = {head: ' ', unit: ' ', tail: ' '} ; Dummy value
    null = PTR_NEW()
    state = {zoomcol: 0L, zoomin: 0L, zoomout: 0L, progcol: 0L, $ ; Widget IDs
             tabs: 0L, readout: 0L, label: 0L, readlab: 0L, $     ; Widget IDs
             pad1: 0L, pad2: 0L, OV_button: 0L, blink: 0L,  $     ; Widget IDs
             imstats: 0L, peakfit: 0L, profile: 0L, frames: 0L, $ ; Widget IDs
             disp: 0L, datasets: 0L, $ ; Widget IDs
             maxzoom: 0L, $            ; maxzoom variable 
             tab_template: null, $ ; -> template for tab description structure
             tabarr: null, $       ; -> structure array for each tab
             lasttab: -1L, $          ; last tab # referred to in log file
             title: title,       $    ; structure with elements of read label
             loglun: loglun,     $    ; Logical unit number for log file
             version: version,   $    ; Ximview version #
             first: 1B,          $    ; still loading first image(s)
             proj: proj, roll: roll, verbose: verbose, $ ; Input parameters
             focus: 1B, $             ; True if we have keyboard focus
             mark_opt: 0B, $          ; Option for how mark is set
             zoom_opt: 0B, $          ; Option for how zoom is done
             in_count: 0S, $          ; number of input datasets read
             files: PTR_NEW(/ALLOCATE_HEAP), $ ; -> array of dataset info
             maxwin: [512S, 512S], $  ; max size of draw window (updated later)
             global_colour: 0B, $     ; True if all tabs should use same LUT
             xim_graph: PTR_NEW(xim_graph), $ ; -> graphics state for ximview
             blink_seq: PTR_NEW(-1), $        ; -> blink sequence array
             xsize: 0S, ysize: 0S, mbar: 0S, $    ; tlb geometry
             newwidth: 0S, newheight: 0S, $       ; tlb geometry
             imsize: LONARR(5), nside: 0L, ns4: 0L, $ ; image geometry
             ring0: PTR_NEW(/ALLOCATE_HEAP), $    ; -> scratch array
             roi: 0B, statbox: [33S, 33S], $      ; Parameters for Imstats
             maxbox: 7S, peak: 1S, $              ; Parameters for  Peakfit
             is_astrom: 0B, $                ; generic astrometry available
             pole_pix: DBLARR(4), $          ; Pixel coords of poles, or nan
             ll_fmt: "(' ',2F9.3)", $        ; format for long/lat
             prec: 0S, $                     ; # decimal places for arcsec
             astrom: null, csystem: '', $    ; astrometry
             prog: '', $ ;                   ; callable application name
             control: PTR_NEW(/ALLOCATE_HEAP)} ; control structure for 
                                             ; applications

; Define structure containing data for each tab.
;
; Note: absrange, sdev, mode, range, zero are in the units of the input image
; However, displayed values in the readout and in widgets such as XIMVIEW_SCALE
; are multiplied by tab_str.mult to give a sensible range and the stored unit 
; string is appropriate when this re-scaling has been done.
;
    tab_str = {base: 0L, draw: 0L, scale: 0L, $     ; Widget IDs
           window: 0L, scale_index: 0L, $       ; Draw window indices
           im_ptr: null, $                      ; -> original image
           byte_ptr: null, $                    ; -> scaled images
           lut: null, $                         ; -> colour scale structure
           rgb: REPLICATE(null, 3), $           ; -> RGB screens, if any.
           absrange: [0.0, 0.0], $              ; Full intensity range
           sdev: 0.0, $                         ; estimated "noise" rms
           noise_type: 0, $                     ; noise correlation type
           mode: 0.0, $                         ; modal intensity.
           range: [0.0, 0.0], $                 ; intensity range displayed
           beta: 1.0, $                         ; parameter for Asinh scaling
           zero: 0.0, $                         ; zero for display
           unit: 'unknown', $               ; intensity unit
           beam: {BEAM, bmaj: 0.0, bmin: 0.0, bpa: 0.0, beam_area: 0.0}, $
           threshold: !values.F_NAN, $          ; Parameter for Imstats
           mult: 1.0, $                         ; Multiplier for flux
           format: "(1X,F9.3,'  ')", $          ; Format for pixel value
           screen: 0, $                         ; tab #
           column: -1, $                        ; Column number in input file
           polcode: 0, $                        ; Polarization code
           freqcode: '', $                      ; Frequency descriptor
           tablab: '', $                        ; End of label for tab
           wrap: 0, $       ; controls display of intensities outside "range"
           trfunc: 0, coltab: 0, $              ; display options
           izero: 0B, $                         ; LUT level of zero
           decomposed: 0B, $                    ; is device decomposed?
           temporary: 0B, $                     ; Delete data at end
           collab: ['mono','',''], $            ; Labels for colour channels
           map_seq: 0B}                         ; input dataset sequence #

    state.TAB_TEMPLATE = PTR_NEW(tab_str)
    tablab =  nin GT 0 ? 'Loading...' : ' '
    IF is_gdl() THEN tablab = [tablab,'Temporary']
; Save current graphics state in structure old_graph:
    swap_lut, xim_graph, dummy, old_graph

; make and realize the main widget
    ntab = is_gdl() ? 2 : 1
    make_ximview, 'Welcome to Ximview', ' ', state, tlb, ntab, tablab
    str = (*state.TABARR)[0]  ; NB at this stage the *only* tab!
    noldtab = 1

; Temporarily set LUT to old values
    lut = {R: old_graph.RED, G: old_graph.GREEN, B: old_graph.BLUE, $
           LINE: old_graph.P.COLOR, ABSENT: old_graph.P.BACKGROUND}
    *str.LUT = lut
    WSET, str.WINDOW
    
;set up gr_global common:
    sys_cols = WIDGET_INFO(tlb, /SYSTEM_COLORS)
    nsyscol = N_TAGS(sys_cols)
                                ; Find grey levels amoung system colours
    cols = [0]
    FOR i=0,nsyscol-1 DO BEGIN
        col = sys_cols.(i)      ; NB: integer triplet not byte triplet
        IF col[0] NE -1 THEN BEGIN
            IF col[0] EQ col[1] && col[0] EQ col[2] THEN cols = [cols, col[0]]
        ENDIF
    ENDFOR
    sys_cols = 0
    cols  = cols[ UNIQ( cols, BSORT(cols) ) ]
    IF verbose THEN PRINT, 'System grey-levels:', cols

    IF is_unix THEN greys = WHERE(cols NE 0B AND cols NE 255B, ngrey) $
               ELSE ngrey = 0
    IF ngrey GT 0 THEN greys = cols[greys]

; Set grey levels for non-data pixels. colour index = grey level so
; that the CI does not need to be changed for decomposed (RGB) images:
    IF ngrey LT 1 THEN greys = [106] ; Default off-sky level
    IF ngrey LT 2 THEN greys = [greys, 192] ; Default bad pixel level
    
; Find out about colour maps: Do we have one? Do we have to
; re-draw image if the map is changed?
    DEVICE, DECOMPOSED = 0  ; Affects results under TrueColor
    colmap     = COLORMAP_APPLICABLE(redraw_req)
    redraw_req = colmap && redraw_req

    IF redraw_req THEN greys = greys[0:1] ; only keep non-data levels.
    syscol = BYTE(greys)
    badcol = syscol[1]
    
    white = !D.table_size - 1
    ibot = redraw_req ? 0B : 1B
    itop = (white + ngrey - 4) < (white - 2)

; Placeholders for graphics colours:
    m2_col = 0B
    m3_col = m2_col
    grid_col = m2_col
    lab_col = m2_col
    geo_col = m2_col

    IF verbose THEN PRINT, SYSTIME(1)-start, "Graphics state set", FORMAT=fmt
    IF verbose THEN HELP, /MEMORY
                                ; End of "launching" set up.
ENDIF ELSE BEGIN                ; Loading new data into existing widget

    IF KEYWORD_SET(expert) THEN !ximview_expert = 1B

    first = state.FIRST ; true if we did not load any data originally.

    IF ~first && overwrite && (nin GT 0) THEN BEGIN
       IF N_ELEMENTS(*state.TABARR) GT 1 THEN BEGIN
          msg = 'Delete all tabs but the last before trying to overwrite'
          ok = DIALOG_MESSAGE(msg, /ERROR, DIALOG_PARENT=tlb)
       ENDIF
      ; One tab left. Clear it by sending event to ximview_deltab
       newevent = {ID: state.PAD2, TOP: tlb, HANDLER: 0L, TAB: 0L}
       ximview_deltab, newevent
      ; Retrieve state structure: reset by deleting last tab.
       WIDGET_CONTROL, tlb, GET_UVALUE = state
    ENDIF
    proj   = state.PROJ
    roll   = state.ROLL
    first = state.FIRST ; may have been reset by ximview_deltab

    str = (*state.TABARR)[0]
    
    prep_screen, state, mode, otabarr, old_graph
    noldtab = N_ELEMENTS(otabarr)
    otabarr = 0
                                ; grab graphics focus if we don't
                                ; already have it:
    IF ~state.FOCUS THEN BEGIN
        state.FOCUS = 1
        WIDGET_CONTROL, tlb, SET_UVALUE = state
        WIDGET_CONTROL, str.DRAW, /INPUT_FOCUS
    ENDIF
ENDELSE

IF nin GT 0 THEN BEGIN ; we have some data to load

;WIDGET_CONTROL, /HOURGLASS

; Suppress any dummy axes in input
    input_size = SIZE(input)
    restore_size = ~temporary && input_size[0] GT 0
    IF input_size[0] GT 0 THEN input = REFORM(input,/OVERWRITE)

                                ; Parse the input parameter
    parse_input, input, ORDER = order, PROJ = proj, TEMPORARY = temporary, $
        COLUMN = column, data, header, newname, howto, tablab, scale_pars, $
        namestr, polcodes, freqcodes, GET_SCALE = auto_scale, $
        EXTENSION = exten, VERBOSE = verbose

    IF verbose THEN BEGIN
        PRINT, SYSTIME(1)-start, "Input converted", FORMAT=fmt
        HELP, /MEMORY
    ENDIF

; create widget tab(s), and stash data in widget uservalues/heap
    make_tabs, input, proj, column, wrap, roll, name, temporary, state, $
        noldtab, data, header, newname, howto, tablab, namestr, $
        polcodes, freqcodes, ntab, imsize, line, extradims, mismatch

    IF mismatch THEN MESSAGE, 'New data does not match data already loaded'

; Scale the image
    scale_tabs, ntab, noldtab, column, state, auto_scale, scale_pars, howto, $
                extradims, input, range, start

    start_gscroll, ntab, state, imsize, noldtab, mode, start

    tabarr = state.TABARR

; Draw initial screens:
    fill_screens, ntab, noldtab, tabarr, mode, state, start

                                ; Update readout label
    title = state.TITLE
    title_string = title.HEAD + form_unit((*tabarr)[noldtab].UNIT) + title.TAIL
    WIDGET_CONTROL, state.READLAB, SET_VALUE = title_string

    state.FIRST = 0B

; Store global data
    WIDGET_CONTROL, state.TABS, SET_UVALUE = mode
ENDIF ELSE WIDGET_CONTROL, state.TABS, GET_UVALUE = mode

; Set Mac cursor options if requested:    
macstyle = KEYWORD_SET(macstyle)
IF macstyle THEN BEGIN
   state.zoom_opt = 1B
   state.mark_opt = 2B
ENDIF

; Store global data
WIDGET_CONTROL, tlb, SET_UVALUE = state

IF first && ~state.FIRST THEN BEGIN
   start_msg = ['Click cursor to define centre of zoom/scroll, then:',$
                '      Mouse button 1 to drag image', $
                '                   2 to mark point and record value', $
                '      zoom with scroll wheel or buttons at screen left', '', $
                'Maxfit and Imstats buttons work around last marked point.']
   IF macstyle THEN $
      start_msg = ['Click cursor to define centre of zoom/scroll, then:',$
                '      Hold down mouse button to drag image', $
                '      Shift-click to mark point and record value', $
                '      zoom with Ctrl-click [in]/Ctrl-shift click [out]', $
                '                or buttons at screen left', '', $
                'Maxfit and Imstats buttons work around last marked point.']
   IF ~!ximview_expert THEN $
      ok = DIALOG_MESSAGE(start_msg, /INFORMATION, $
                          TITLE = 'XIMVIEW', DIALOG_PARENT = tlb)
   line = [line, title_string]
   PRINTF, state.LOGLUN, line, FORMAT="(A)"
   PRINT, line, FORMAT="(A)"
ENDIF

; Enable buttons etc:
*state.BLINK_SEQ = INDGEN(ntab+noldtab)
WIDGET_CONTROL, tlb, /SENSITIVE
IF mode.OVERVIEW THEN WIDGET_CONTROL, state.ZOOMCOL, SENSITIVE = 0
sens =  ntab + noldtab GE 2
WIDGET_CONTROL, state.BLINK,  SENSITIVE = sens
WIDGET_CONTROL, state.FRAMES, SENSITIVE = 1

WIDGET_CONTROL, state.PROGCOL, SENSITIVE = state.prog NE ''

; Turn on events only if we're really ready
sens =  ~state.FIRST
WIDGET_CONTROL, state.DISP, SENSITIVE = sens
WIDGET_CONTROL, state.TABS, SENSITIVE = sens
WIDGET_CONTROL, str.DRAW, SENSITIVE = sens
WIDGET_CONTROL, state.OV_BUTTON, SENSITIVE = sens
WIDGET_CONTROL, state.IMSTATS, SENSITIVE = sens
WIDGET_CONTROL, state.PEAKFIT, SENSITIVE = sens
WIDGET_CONTROL, state.PROFILE, SENSITIVE = sens

CATCH, /CANCEL

IF launching THEN BEGIN
    IF verbose THEN BEGIN
        PRINT, SYSTIME(1)-start, "Starting event loop", FORMAT=fmt
        HELP, /MEMORY
    ENDIF
    XMANAGER, 'ximview', tlb, CLEANUP='ximview_cleanup', /NO_BLOCK
ENDIF ELSE WIDGET_CONTROL, state.TABS, /SENSITIVE

restore_lut, *state.XIM_GRAPH, old_graph

END
