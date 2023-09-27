;-----------------------------------------------------------------------------
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
;---------------------------------------------------------------------------
;
; Widget creation routines for main ximview widget and its tab screens
;
; Contents:
;   MATCH_TAB      Function to check if proposed new data matchs other
;                  tabs in size and astrometry
;   UPDATE_LABELS  Attempts to generate a set of meaningful and unique
;                  tab labels given header info.
;   CREATE_TAB     Does the actual widget creation for a tab
;   MAKE_TABS      Wrapper routine to create a set of tabs. Calls all
;                  the previous routines as well as PARSE_HEADER
;   MAKE_RGB_TAB   Specialist version of MAKE_TABS for 3-colour tab
;   SCALE_TABS     Applies scaling and transfer function to images on
;                  tabs, and also updates units
;   MAKE_XIMVIEW   Creates the main Ximview widget.
;
FUNCTION match_tab, state, imsize, proj, is_astrom, astrom, csystem
;
; Check that newly loaded data matches other tabs in size and astrometry
;
; Inputs
;   state:     State structure describing existing tabs
;   imsize:    SIZE array for new data
;   proj:      HEALPix projection (if any)
;   is_astrom: General astrometry is available
;   astrom:    Astrolib WCS astrometry structure
;   csystem:   String describing coordinate system
;
COMPILE_OPT IDL2, HIDDEN

bad1 = ~ARRAY_EQUAL(imsize[1:2], state.IMSIZE[1:2])
IF bad1 THEN MESSAGE, /INFORMATIONAL, $
  'New array size does not match data already loaded'

bad2 = proj NE state.PROJ
IF bad2 THEN  MESSAGE, /INFORMATIONAL, $
  'New array HP projection:' + proj + ' does not match data already loaded'

bad3 = is_astrom NE state.IS_ASTROM
IF bad3 THEN MESSAGE, /INFORMATIONAL, 'Warning: ' +  (is_astrom ? $
    'New array has astrometry but old does not' : $
    'New array has no astrometry but old does')

IF is_astrom && state.IS_ASTROM THEN BEGIN         ; check astrom structure
    good = 1
    ostrom = *state.ASTROM
    FOR i=0,4 DO good *= ARRAY_EQUAL(astrom.(i), ostrom.(i))
    good *= ARRAY_EQUAL(astrom.PV1, ostrom.PV1)
    good *= ARRAY_EQUAL(astrom.PV2, ostrom.PV2)
    IF FINITE(astrom.EQUINOX) && FINITE(ostrom.EQUINOX) THEN $
         good *= astrom.EQUINOX  EQ ostrom.EQUINOX
    IF astrom.CTYPE[0] NE ostrom.CTYPE[0] THEN BEGIN
        MESSAGE, /INFORMATIONAL, 'Warning: new data coordinate type: ' + $
          astrom.CTYPE[0]
        MESSAGE, /INFORMATIONAL, '         does not match previous:  ' + $
          ostrom.CTYPE[0]
    ENDIF
    bad4 = ~good
    IF bad4 THEN MESSAGE, /INFORMATIONAL, $
        'New array astrometry cannot match data already loaded'
ENDIF ELSE bad4 = 0B

IF is_astrom && csystem NE state.CSYSTEM THEN  BEGIN
    MESSAGE, /INFORMATIONAL, 'Warning: new array coordinate system: '+csystem

    MESSAGE, /INFORMATIONAL, '         does not match previous one: ' + $
      state.CSYSTEM
ENDIF

RETURN, ~(bad1 || bad2 || bad4)

END

PRO update_labels, state, newname, newtabs, tab_lab, topname
;
; Sets top title and tab labels when new data tabs are added
;
; Aim:
;   Labels should all be different & pattern of labels should be roughly
;   the same for all tabs. Need to maintain this as we load new datasets
;
;   format of label is <name:> <polarization> <freq> <other>
;   Each element can be blank.
; 
;   <name> is the user-supplied label for dataset on input.
;   If present for some, but not others, missing items have name = 'File n'
;   Name field is used if supplied for at least one input, or if the other
;   fields do not distinguish tabs. Note if the user specifies the same name
;   for more than one input there will probably be repeated labels anyway.
;
;   Polarization is used iff there is more than one polarization amongst the
;   tabs. Extended AIPS definition of pol includes spectral index, tau, etc.
;
;   Same goes for <frequency> (this may actually be wavelength, energy, band, 
;   velocity, etc).
;  
;   <other> is any other tab-distinguishing info extracted from the header
;   by parse_input. Default value 'Map n' is ignored unless all else fails.
;
; Inputs:
;    state:   structure with main global variables
;    newname: Name string built by parse_input
; Input/output
;    newtabs: On input, array of structures to describe new tabs.
;             on output, array describing all tabs.
;    tab_lab: On input, array of tab labels from parse_input
;             on output, updated
; Output:
;    topname: New top title
;
COMPILE_OPT IDL2, HIDDEN
ON_ERROR, 2

scodes = ['YX','XY','YY','XX', 'LR', 'RL', 'LL','RR', '', 'I','Q','U','V', $
          'Pol Int', 'Frac Pol', 'Pol Ang', 'alpha', 'tau', 'RM']
scodes[0:12] += ' Pol'

first = state.FIRST
ncol = N_ELEMENTS(newtabs)
IF first THEN tab_arr = newtabs ELSE BEGIN
    oldtabs = *state.TABARR
    oldlabs = get_tab_uvals(state.TABS)
    tab_arr = [oldtabs, newtabs]
ENDELSE
oldtabs = 0  &  newtabs = 0
n_tab = N_ELEMENTS(tab_arr)
tab_lab = STRARR(n_tab)

; Exclude tabs such as RGB comparison and polarization that borrow data
; from other tabs
real = WHERE(tab_arr.MAP_SEQ GE 1, ntab, COMPLEMENT=special)
IF ntab LT n_tab THEN tab_lab[special] = oldlabs[special]

tabarr = tab_arr[real]

file_str = *state.files
nfiles = N_ELEMENTS(file_str) ; should equal state.IN_COUNT
names  = file_str.NAME
freqs  = tabarr.FREQCODE
fid    = tabarr.MAP_SEQ - 1 ; index from files array to tab array

test = WHERE(~STRCMP(names, 'File', 4), nname)
test = WHERE(tabarr.POLCODE NE tabarr.POLCODE[0], npol)
test = WHERE(freqs NE freqs[0], nfreq)
test = WHERE(~STRCMP(tabarr.TABLAB, 'Map', 3), nother) 

IF nfreq GT 0 THEN FOR itab = 0,ntab-1 DO BEGIN
                               ; Eliminate repetition of 'band'
    freqs[itab] = STRJOIN(STRSPLIT(freqs[itab],'band', /REGEX, /EXTRACT))
ENDFOR

REDO:

tablab = STRARR(ntab)
IF npol GT 0   THEN tablab += scodes[tabarr.POLCODE-8]
IF nfreq GT 0  THEN tablab += ' ' + freqs
IF nother GT 0 THEN tablab += ' ' + tabarr.TABLAB
tablab = STRTRIM(STRCOMPRESS(tablab), 2)
IF nname GT 0 THEN BEGIN
    aa = WHERE(tablab NE '', na, COMPLEMENT = bb)
    IF na GT 0    THEN tablab[aa] = names[fid[aa]]+': '+tablab[aa]
    IF na LT ntab THEN tablab[bb] = names[fid[bb]]
ENDIF

test = UNIQ(tablab, SORT(tablab))
blank = WHERE(tablab EQ '',nblank)
bad = N_ELEMENTS(test) LT ntab || nblank GT 0

; If labels are not yet unique, use file number and then channel number:
IF bad && nname EQ 0 && nfiles GT 1 THEN BEGIN
    nname = 1
    GOTO, REDO
ENDIF
IF bad && nother EQ 0 THEN BEGIN 
    nother = 1
    GOTO, REDO
ENDIF  

tab_lab[real] = tablab

                                ; Work out new plot title and tab
                                ; prefix from name structure
files  = file_str.FILE            
IF nfiles GT 1 THEN BEGIN
    newname = STRJOIN(files,', ') + ': '
    len1 = STRLEN(newname)
    start = WHERE( STRCMP( TAG_NAMES(file_str), 'file', /FOLD_CASE) )
    nntail = ''
    FOR i = 1,6 DO BEGIN
        param = file_str.(i + start)
        par   = param[ UNIQ( param, BSORT(param) ) ]
        npar  = N_ELEMENTS(par)
        IF npar EQ 1 THEN nntail += par 
    ENDFOR
    nntail = STRTRIM(STRCOMPRESS(nntail),2)
    len2 = STRLEN(nntail)

    IF len1+len2 GT 80 THEN newname = files[0]+' and others: '
    newname += nntail
ENDIF             ; Else newname supplied by parse_input should be OK.
topname = newname

tab_arr[real] = tabarr
newtabs = tab_arr

END
;
PRO create_tab, base, tablab, uval, uname, tabarr, itab, xsize, ysize, lut
; Creates the widgets that populate a Ximview tab and fills in the mandatory
; fields of the tab descriptor structure
;
COMPILE_OPT IDL2, HIDDEN
str = tabarr[itab]

is_lut = N_ELEMENTS(lut) GT 0
IF str.LUT THEN PTR_FREE, str.LUT ; should never happen
str.LUT  = is_lut ? PTR_NEW(lut) : PTR_NEW(/ALLOCATE_HEAP)
       
str.BASE = WIDGET_BASE(base, TITLE = tablab, UNAME = uname, $
                       UVALUE = uval, /COLUMN, XPAD = 0, YPAD = 0, SPACE = 0)
str.DRAW = WIDGET_DRAW(str.BASE, XSIZE = xsize, YSIZE = ysize, RETAIN = 2, $
                       EVENT_FUNC = 'ximview_scroll', $
                       /BUTTON_EVENTS, /MOTION_EVENTS, /WHEEL_EVENTS)
str.SCALE = WIDGET_DRAW(str.BASE, XSIZE = xsize, YSIZE = 45, RETAIN = 2)
str.SCREEN = itab

; Get the window IDs for the draw widgets, stored as their widget values:
WIDGET_CONTROL, str.DRAW, GET_VALUE = index
str.WINDOW = index
WIDGET_CONTROL, str.SCALE, GET_VALUE = index
str.SCALE_INDEX = index

tabarr[itab] = str

END
;
PRO make_tabs, input, proj, column, wrap, roll, name, temporary, state, $
               noldtab, data, header, newname, howto, tablab, $
               namestr, polcodes, freqcodes, ntab, T, line, extradims, mismatch
; Wrapper routine:
;   Get astrometry & unit info (parse_header)
;   Update file list
;   Create tab structures to describe tabs & fill in from header
;   Create tabs on widget
;   Store data on the heap
;
; INPUTS
;  input, proj, column, wrap, roll, name, temporary:
;              Ximview inputs, defaults filled in; 
;              column converted to integer indices if needed.
;  state:      top-level control structure
;  data:       data or pointer to data, unless still in input.
;  header:     FITS header for data
;  newname:    name constructed by PARSE_INPUT from namestr
;  howto:      control parameter for accessing data
;  tablab:     (partial) labels for tabs to be created, from PARSE_INPUT
;  namestr:    structure with info from header
;  polcodess:   polarization codes for each tab to be created
;  freqcodes:  frequency codes (strings) for each tab
;
; INPUT/OUTPUT
;  noldtab:    Number of tabs already created (>= 1). 
;              Just a dummy if this is the first time data has been loaded;
;              if so, set to zero on output.
;
; OUTPUTS:
;  ntab:       Number of tabs made
;  T:          IDL SIZE array for data
;  line:       String describing data to print in log.
;  extradims:  true if data has more than 2 dimensions.
;  mismatch:   true if new data does not match existing tabs in size or
;              astrometry
;
COMPILE_OPT IDL2, HIDDEN
COMMON gr_global
ON_ERROR, 0
verbose = state.VERBOSE
first = state.FIRST

ntab = N_ELEMENTS(column)

IF N_ELEMENTS(name) EQ 0 THEN name = ''
log_header = newname + ' ' + name
ID_string = STRING(log_header, FORMAT = "('Dataset: ',A)")
PRINTF, state.LOGLUN, ID_string
PRINT, ID_string

IF howto LT 2 THEN BEGIN        ; swap to pointer usage if possible
    IF ntab EQ 1 THEN BEGIN
        temp = PTR_NEW(/ALLOCATE_HEAP)
        CASE howto OF
            0: *temp = input
            1: *temp = TEMPORARY(data)
        ENDCASE
        data  = temp
        howto = 2
        temporary = 1B ; we can delete heap variable
    ENDIF ELSE IF ~first THEN BEGIN  
                       ; should never get here as PARSE_INPUT always uses
                       ; pointer arrays for multi-column
        used = PTR_VALID((*state.TABARR).IM_PTR)
        IF MIN(used) EQ 0 THEN $
          MESSAGE, 'Internal error: bulk data slot already in use'
    ENDIF
ENDIF

IF howto LT 2 THEN void = $
   DIALOG_MESSAGE('Still using bulk data slot!',/CENTER)

CASE howto OF
    0: T = SIZE(input) ; input data was 2+D array
    1: T = SIZE(data)  ; input data converted to 2-D array
    2: T = SIZE(*data) ; 2-D array is stored via a pointer
    3: T = SIZE(*data[0]) ; array of 2-D images via pointer array
    ELSE: MESSAGE, 'Unknown howto value!'
ENDCASE

IF howto EQ 3 THEN BEGIN
    extras = N_ELEMENTS(data)
    IF extras EQ 1 THEN howto = 2 ELSE BEGIN
        IF T[0] NE 2 THEN MESSAGE, 'Pointer array should be to 2-D images'
    ENDELSE
ENDIF

extradims = T[0] GT 2
intype    = T[T[0]+1] ; IDL data type

IF extradims THEN BEGIN
    IF howto NE 0 THEN $
      MESSAGE, STRING(howto, FORMAT = $
                      "('Internal error: howto =',I2,' and >2 dimensions')")
    T = SIZE(input[*,*,0])
ENDIF

; Extract info from header
parse_header, T, header, column, roll, state.VERBOSE, $
              astrom, is_astrom, pole_pix, csystem, proj, unit, beam, title, $
              nside, ns4, ll_fmt, prec, boxsize 

IF first THEN BEGIN
                                ; Store info in state structure
    state.PROJ  = proj       &  state.ROLL = roll
    state.IMSIZE = T         &  state.NSIDE = nside  &  state.NS4  = ns4
    str = (*state.TABARR)[0]
    mismatch = 0B
    line = ['']
    IF nside GT 0 THEN line = [line, 'Seems to be HEALPix with N_side' + $
                              STRING(nside,FORMAT="(I5)"), '']
    state.TITLE = title
    state.STATBOX = boxsize > 33
ENDIF ELSE BEGIN             ; Check for agreement with existing maps:
    mismatch = ~match_tab(state, T, proj, is_astrom, astrom, csystem)
    IF mismatch THEN BEGIN
        IF howto GE 2 && temporary THEN  FOR i=0,ntab-1 DO $
          IF PTR_VALID(data[i]) THEN PTR_FREE, data[i]
        RETURN
    ENDIF
    str = *state.TAB_TEMPLATE
ENDELSE

; Update labelling if we have gained astrometry:
IF is_astrom  && ~state.IS_ASTROM THEN BEGIN ; update astrom stuff
    IF N_ELEMENTS(line) EQ 0 THEN line = ''
    line = [line, 'Coordinate system: '+csystem, '']                           
    state.CSYSTEM = csystem  &  state.IS_ASTROM = is_astrom
    IF SIZE(astrom, /TYPE) EQ 8 THEN BEGIN
        state.ASTROM = PTR_NEW(astrom)
        state.pole_pix = pole_pix
        state.LL_FMT = ll_fmt
        state.PREC = prec
    ENDIF
    state.TITLE = title
ENDIF
astrom = 0

; Update list of files read in
namestr.HEADER = PTR_NEW(header)
newfile = 1B
state.IN_COUNT += 1
IF first THEN *state.FILES = [namestr] ELSE BEGIN
          ; Check this really is a new dataset, not one already read.
    mswin = windev EQ 'WIN'
    old1 = WHERE(STRCMP(namestr.FILE, (*state.FILES).FILE, $
                 FOLD_CASE=mswin), hit1)
    old2 = WHERE(namestr.CHECKSUM EQ (*state.FILES).CHECKSUM, hit2)
    IF ~hit1 || ~hit2 || old1 NE old2 THEN BEGIN
        *state.FILES = [*state.FILES, namestr]
    ENDIF ELSE BEGIN
        state.IN_COUNT -= 1
        newfile = 0B
    ENDELSE
ENDELSE 
IF newfile THEN BEGIN
    fnum = state.IN_COUNT
    IF name EQ '' THEN name = 'File '+STRTRIM(STRING(fnum),2)
    (*state.FILES)[fnum-1].NAME = name
    cat_entry = name+': '+namestr.FILE
    ext = namestr.EXTENSION
    IF ext GT 0 THEN cat_entry += ' Ext '+STRTRIM(STRING(ext),2)
    void = WIDGET_BUTTON(state.DATASETS, VALUE = cat_entry, $
                         EVENT_PRO='ximview_header', UVALUE = fnum)
ENDIF

; Make tab array for new tabs
IF N_ELEMENTS(wrap) THEN str.WRAP = wrap
str.MAP_SEQ  = newfile ? state.IN_COUNT : old1+1 

IF howto LE 1 THEN temporary = 1B ; always safe to delete stored version.
str.TEMPORARY  = temporary
newtabs = REPLICATE(str, ntab)
newtabs.UNIT     = unit
newtabs.BEAM     = REPLICATE(beam,ntab)
newtabs.NOISE_TYPE = 0 ; unknown noise type; not stored in header
newtabs.POLCODE  = polcodes
newtabs.FREQCODE = freqcodes
newtabs.COLUMN   = column
newtabs.TABLAB   = tablab

; Store data:
CASE howto < 2 OF
    0: WIDGET_CONTROL, state.LABEL, SET_UVALUE = input
    1: WIDGET_CONTROL, state.LABEL, SET_UVALUE = data, /NO_COPY
    2: ; Data stored via pointer in state structure.
ENDCASE

IF intype EQ 1 THEN BEGIN ; input is already in bytes
    IF howto LE 1 THEN BEGIN ; this never happens (I hope!)
        WIDGET_CONTROL, state.LABEL, GET_UVALUE = data
        newtabs.BYTE_PTR = PTR_NEW(TEMPORARY(data))
    ENDIF ELSE newtabs.BYTE_PTR = data
ENDIF ELSE newtabs.BYTE_PTR = PTRARR(ntab, /ALLOCATE_HEAP)

newtabs.IM_PTR   = howto LT 2 ? PTR_NEW() : data

                        ; Update top title and tab labels, and
                        ; merge old & new arrays
update_labels, state, newname, newtabs, tablab, topname
tab_arr = TEMPORARY(newtabs)    ; Now contains old and new

WIDGET_CONTROL, state.LABEL, SET_VALUE = topname

                        ; Change labels of old tabs:
FOR itab=0,noldtab-1 DO BEGIN
    IF ~is_gdl() THEN $
       WIDGET_CONTROL, tab_arr[itab].BASE, BASE_SET_TITLE = tablab[itab]
    WIDGET_CONTROL, tab_arr[itab].BASE, SET_UVALUE = tablab[itab]
ENDFOR

nnewtab = first ? ntab-1 : ntab

xsize = !D.x_vsize  & ysize = !D.y_vsize
FOR icol = 0,nnewtab-1 DO BEGIN    ; Make the new tabs
    itab = icol + noldtab
    tab_arr[itab].LUT = PTR_NEW() ; otherwise create_tab will free old ptr
                                  ; which may be in use for first tab.
    create_tab, state.TABS, tablab[itab], tablab[itab], 'mono', $
      tab_arr, itab, xsize, ysize
ENDFOR
IF first THEN noldtab = 0         ; we want to start with tab 0 from now on

*state.TABARR = tab_arr

END
;
PRO make_rgb_tab, str, ntab, rgbptr, collab, tablab, state
;
; Creates a tab for an RGB image, fills in the tab structure
; Fills in the image and labels the plot. NB special tab as does not
; correspond to an input dataset.
;
; Inputs
;     rgbptr:       Array of 3 pointers to the rgb byte images
;     collab:       Array of 3 labels for the rgb channels
;     tablab:       Tab title label
;     state:        Usual system structure
;
; INPUT/OUTPUT
;     str:          Structure to describe tab to be created
; OUTPUT
;     ntab:         Number of existing tabs / index of new tab
;
COMPILE_OPT IDL2, HIDDEN

FIRST = state.FIRST

str.DECOMPOSED = 1
str.RANGE      = [0, 255]
str.ABSRANGE   = str.RANGE
str.RGB        = rgbptr
str.UNIT       = '   R   G   B' ; readout gives RGB byte values for now.
str.COLLAB     = collab
str.MAP_SEQ    = -1             ; indicates that data stored in other tabs.
; Cancel existing pointers
str.IM_PTR = PTR_NEW()  &  str.BYTE_PTR = PTR_NEW()
ramp = LINDGEN(255)
lutstr = {r:ramp, g: ramp, b: ramp, line: !P.color, absent: !P.background}
; For the sake of DirectColor visuals, not that they seem to work.

tabarr = first ? str : [TEMPORARY(*state.TABARR),str]
ntab   = N_ELEMENTS(tabarr) - 1 

IF first THEN BEGIN ; Initial tab already created by make_ximview
    *str.LUT = lutstr
ENDIF ELSE BEGIN
; Create a tab for the RGB image:
    tabarr[ntab].LUT = PTR_NEW()
    xsize = !D.x_vsize  & ysize = !D.y_vsize
    create_tab, state.TABS, tablab, tablab, tablab, tabarr, ntab, $
                xsize, ysize, lutstr
ENDELSE

*state.TABARR = tabarr
str = tabarr[ntab]

END
;
PRO scale_tabs, ntab, noldtab, column, state, auto_scale, scale_pars, howto, $
                extradims, input, range, start
;
; Scales a bunch of new images
;
COMPILE_OPT IDL2, HIDDEN
COMMON gr_global
ON_ERROR, 2
verbose = state.VERBOSE
fmt = "(F7.3,' seconds: ',A)"

nr = N_ELEMENTS(range)
no_range = nr EQ 0 && ~auto_scale
IF nr EQ 3 THEN BEGIN
   beta = range[2]
   trfunc = 1
   range = range[0:1]
ENDIF

tabarr = state.TABARR

first = state.FIRST
IF ~first THEN BEGIN
    laststr = (*tabarr)[noldtab-1]
    lutstr = *laststr.LUT
    decomp = laststr.DECOMPOSED
    oldcoltab = laststr.COLTAB
    oldzero = laststr.IZERO
ENDIF ELSE BEGIN
    oldcoltab = -1
    oldzero = !values.F_NAN
ENDELSE

FOR icol = 0,ntab-1 DO BEGIN
    itab = icol + noldtab
    incol = column[icol] - 1
    str = (*tabarr)[itab]
    bptr = str.BYTE_PTR  &  iptr = str.IM_PTR &  wrap = str.WRAP
    IF nr NE 3 THEN BEGIN
       trfunc = str.TRFUNC
       beta = str.BETA
    ENDIF ELSE BEGIN
       str.TRFUNC = trfunc
       str.BETA = beta
    ENDELSE
    zero = str.ZERO  & izero = str.IZERO

    IF auto_scale THEN BEGIN
        set_scale, scale_pars[*, icol], str
        range = str.RANGE
    ENDIF ELSE IF no_range THEN range = 'Full' 
                                  ; tells scale_image to use full data range

    CASE howto < 2 OF
        0: *bptr = extradims ? scale_image(input[*,*,incol], range, wrap, $
                                           trfunc, zero, beta, $
                                           ABS_RANGE = ar, VERBOSE=verbose) $
                             : scale_image(input, range, wrap, trfunc, zero, $
                                           beta, ABS_RANGE = ar,  $
                                           VERBOSE=verbose)
        1: BEGIN
            WIDGET_CONTROL, state.LABEL, GET_UVALUE = data, /NO_COPY
            *bptr = scale_image(data, range, wrap, trfunc,zero, beta, $
                                ABS_RANGE = ar, VERBOSE=verbose)
            WIDGET_CONTROL, state.LABEL, SET_UVALUE = data, /NO_COPY
        END
        2: *bptr = scale_image(*iptr, range, wrap, trfunc, zero, beta, $
                               ABS_RANGE = ar, VERBOSE=verbose)
    ENDCASE

    str.ABSRANGE = ar
    str.RANGE    = range

    absmax = 0.1d0*MAX(ABS(ar)) ; will put max in range [10-10000)
    test = NUMUNIT(absmax, str.UNIT, MULTIPLIER = mult, OUT_UNIT = ounit, $
                   /FORCE)
    IF mult EQ 0.0 THEN BEGIN
        PRINT, 'Absmax:', 10d0*absmax,' ',str.UNIT
        MESSAGE, 'Zero multiplier'
    ENDIF
    str.UNIT = ounit
    str.MULT = mult
    
    set_print_fmt, str
    
    str.IZERO = scale_image(zero, range, wrap, trfunc, zero, beta)

    (*tabarr)[itab] = str

; If necessary set colour table and remember it
    IF  str.COLTAB NE oldcoltab || $
       (str.COLTAB EQ 4 &&  oldzero NE str.IZERO) THEN BEGIN 
        ximview_lut, str.COLTAB, str.IZERO, decomp
        TVLCT, r, g, b, /GET
        lutstr = {r:r, g:g, b:b, line: !P.color, absent: !P.background}
    ENDIF
    
    *str.LUT = lutstr
    str.DECOMPOSED = decomp
    oldtab = str.coltab
    oldzero = str.izero

ENDFOR

IF verbose THEN MESSAGE, /INFORMATIONAL, $
  STRING(SYSTIME(1)-start, 'Images scaled', FORMAT = fmt)
IF verbose THEN HELP, /MEMORY

                               ; Set off-sky pixels to absent:
IF state.NSIDE NE 0 && state.IS_ASTROM THEN BEGIN
    fill_gores, state.NSIDE, state.IMSIZE, state.ASTROM, (*tabarr).BYTE_PTR

    IF verbose THEN MESSAGE, /INFORMATIONAL, $
      STRING(SYSTIME(1)-start, "absent pixels noted", FORMAT = fmt)
ENDIF

END
;
PRO make_ximview, name, title, state, tlb, ntab, tablab
;
; Sets up widget geometry and installs defaults
;
; INPUTS:
;    name    Initial title
;    title   Initial title string for cursor readout
;    state   top-level control structure for ximview
;    ntab    Number of tabs
;    tablab  array of labels for each tab  
;
; OUTPUTS:
;    tlb     Widget index for top-level base
;
COMPILE_OPT IDL2, HIDDEN

is_unix = STRCMP(!version.OS_FAMILY,'UNIX', 4,/ FOLD_CASE)

xsize = 512  & ysize = 512

tlb = WIDGET_BASE(TITLE = 'XIMVIEW', MBAR = menu, /COLUMN, $
                  UNAME = 'XIMVIEW', /TLB_SIZE_EVENTS, /KBRD_FOCUS_EVENTS)

; Funny business to find system default font name

dummy = WIDGET_BUTTON(tlb, VALUE = 'Redshirt')
devfont = WIDGET_INFO(dummy, /FONTNAME)
WIDGET_CONTROL, dummy, /DESTROY
fixedfont = is_unix ? devfont : 'lucida console*10'
;
; Set up menu buttons:
;
file     = WIDGET_BUTTON(menu, VALUE = 'File', /MENU)
options  = WIDGET_BUTTON(menu, VALUE = 'Options', /MENU)
display  = WIDGET_BUTTON(menu, VALUE = 'Display', /MENU)
frames   = WIDGET_BUTTON(menu, VALUE = 'Tabs', /MENU)
analysis = WIDGET_BUTTON(menu, VALUE = 'Analysis', /MENU)
datasets = WIDGET_BUTTON(menu, VALUE = 'Headers', /MENU)
help     = WIDGET_BUTTON(menu, VALUE = 'Help', /HELP, /MENU)

desc = ['0\Load FITS\ximview_fitsload', $
        '0\Load image file\ximview_rdimage', $
        '1\New logfile\ximview_newlog', $
           '0\Overwrite old file', '0\New sequence #', '2\Named...', $
        '1\Write PNG image\ximview_2png', $
           '0\Image and scale bar', '2\Image only', $
        '0\Reset\ximview_reset']
desc2 = ['0\State info\ximview_dump', $  ; Debug items
         '0\Tab info\ximview_dump', $
         '0\Mode info\ximview_dump', $ 
         '0\LUT info\ximview_dump', $
         '0\Astrometry\ximview_dump']
desc3 = ['0\Exit\ximview_exit']
IF state.VERBOSE THEN desc = [desc,desc2,desc3] ELSE desc = [desc,desc3]
file_pd = CW_PDMENU(file, desc, /MBAR, /RETURN_NAME)
desc = ['1\Mark point', '0\Middle mouse button','0\Right mouse button', $
        '2\Shift-click','1\Zoom in/out','0\Mouse wheel','2\Ctrl-/ctrl-shift- click']
opt_pd  = CW_PDMENU(options, desc, /MBAR, /RETURN_FULL_NAME)
cs   = colour_schemes(-1)
ncs  = N_ELEMENTS(cs)
cs   = ['0\'+cs[0:ncs-2], '2\'+cs[ncs-1]]
fun  = scale_funs(-1)
nfun = N_ELEMENTS(fun)
fun  = ['0\'+fun[0:nfun-2], '2\'+fun[nfun-1]]
desc = ['0\Adjust scaling\ximview_scale', $
        '0\Auto scale all tabs\ximview_autoscale_all', $
        '1\Colour table\ximview_newlut', cs, $
;        '1\Colour handling\ximview_colour', $
;        '0\Same for all', '2\Separate on each tab', $
        '1\Grid\ximview_grid', '0\On/Off', '2\Grid interval', $
        '0\Set view centre\ximview_goto', $
        '0\Restore default screen size\ximview_event']
disp_pd = CW_PDMENU(display, desc, /MBAR, /RETURN_FULL_NAME)
desc = ['0\Blink setup\ximview_blink', $
        '0\Red-Green-Blue\ximview_rgb', $
        '1\Polarization\ximview_pol', '2\Colour', $ ; '0\Vector', '2\LIC', $
        '1\Delete tab\ximview_deltab', '0\Current tab', '2\Specify']
fram_pd = CW_PDMENU(frames, desc, /MBAR, /RETURN_FULL_NAME)
desc = ['1\Imstats\ximview_imstat_opts', $
         '0\Box', '0\Region of Interest', '0\Set box size', '2\Set threshold', $
        '1\Peakfit\ximview_maxfit_opts', $
           '0\Find extremum', '0\Find maximum', '0\Find minimum',  $
           '2\Set box size', $
        '1\Clear marks\ximview_clear_mark', $
        '0\Last','2\All', $
        '0\Set image properties\ximview_setprop2', $
; '0\polarization\ximview_setpol', '0\frequency\ximview_setfreq',$
        '0\Set coordinate system\ximview_setcoord', $
        '1\Catalogue\ximview_vo','0\On/Off','0\Process catalogue', $
        '2\Load new catalogue', $
        '0\Profile options\ximview_profile_opts', $
        '1\Ruler\ximview_geom','0\Use marked points', '2\Interactive', $
        '0\Protractor\ximview_geom']
;       '0\Circles', 
;       '0\Mark direction', $
;       '1\Instrument FOVs','2\Planck']
ana_pd  = CW_PDMENU(analysis, desc, /MBAR, /RETURN_FULL_NAME)
desc = ['0\Help\ximview_help', '0\Release Notes\ximview_help', $
        '0\About\ximview_help']
help_pd = CW_PDMENU(help, desc, /MBAR, /HELP, /RETURN_FULL_NAME)

;
; Set up zoom buttons: 
;
base1 = WIDGET_BASE(tlb, /ROW)

leftcol = WIDGET_BASE(base1,/COLUMN, /ALIGN_CENTER)
zoomcol = WIDGET_BASE(leftcol,/COLUMN, EVENT_FUNC = 'ximview_zoom', $
                      /ALIGN_CENTER, UVALUE ='not button')
zoomt = WIDGET_LABEL(zoomcol, VALUE = 'Z')
zoomt = WIDGET_LABEL(zoomcol, VALUE = 'o')
zoomt = WIDGET_LABEL(zoomcol, VALUE = 'o')
zoomt = WIDGET_LABEL(zoomcol, VALUE = 'm')
zoomin    = WIDGET_BUTTON(zoomcol, VALUE = 'in', UVALUE = 'IN', UNAME = 'IN')
;,$ TOOLTIP = 'Zoom in x 2')
zoomreset = WIDGET_BUTTON(zoomcol, VALUE = '1:1', UVALUE = '1:1', $
                          TOOLTIP= 'Reset to 1 image pixel per screen pixel')
zoomout = WIDGET_BUTTON(zoomcol, VALUE = 'out', UVALUE = 'OUT', UNAME = 'OUT')
;,$ TOOLTIP = 'Zoom out x 2')

; Set up programmable buttons:
progcol = WIDGET_BASE(leftcol,/COLUMN, EVENT_PRO = 'ximview_prog', $
  /ALIGN_CENTER, UVALUE ='not button')
zoomt = WIDGET_LABEL(progcol, VALUE = ' ') ; spacer
buttonA = WIDGET_BUTTON(progcol, VALUE = 'A', UVALUE = 'A')
buttonB = WIDGET_BUTTON(progcol, VALUE = 'B', UVALUE = 'B')
buttonC = WIDGET_BUTTON(progcol, VALUE = 'C', UVALUE = 'C')

; Set up main image display area: top title and tabs
rightcol = WIDGET_BASE(base1,/COLUMN)
label_base = WIDGET_BASE(rightcol,/COLUMN)
label    = WIDGET_LABEL(label_base, value=name, /DYNAMIC_RESIZE)
                                ; Set up a separate tab with draw &
                                ; scale windows for each image to
                                ; display
tabs     = is_unix ? WIDGET_TAB(rightcol, EVENT_PRO = 'ximview_tab', $
                                MULTILINE = 20) $
                   : WIDGET_TAB(rightcol, EVENT_PRO = 'ximview_tab')

tab_arr = REPLICATE(*state.TAB_TEMPLATE, ntab)
FOR itab = 0,ntab-1 DO create_tab, tabs, tablab[itab], tablab[itab], 'mono', $
  tab_arr, itab, xsize, ysize
;
; Set up cursor readout area below image tabs:
readlab = WIDGET_LABEL(rightcol, value = title, /ALIGN_LEFT, $
                       FONT = fixedfont, /DYNAMIC_RESIZE)
readout = WIDGET_LABEL(rightcol, value = 'No pixel assigned', $
                       FONT = fixedfont, /ALIGN_LEFT, /DYNAMIC_RESIZE)
;
; Set up main action buttons at bottom of widget
button_row = WIDGET_BASE(tlb, /ROW)
overview   = WIDGET_BUTTON(button_row, VALUE = 'Overview', FONT = bfont, $
                           EVENT_PRO = 'ximview_overview', $
             TOOLTIP = 'Return to initial view of whole image in one panel')

; Use pad1 as a convenient widget when we want explicitly to request scaling:
pad1       = WIDGET_BASE(button_row, EVENT_PRO = 'ximview_scale')
blink      = WIDGET_BUTTON(button_row, VALUE='Blink on/off', FONT = bfont, $
                           EVENT_PRO = 'ximview_blink')
imstats    = WIDGET_BUTTON(button_row, VALUE = 'Imstats', FONT = bfont, $
                           EVENT_FUNC = 'ximview_imstats', $
                    TOOLTIP = 'Prints statistics for box around marked point')
peakfit    = WIDGET_BUTTON(button_row, VALUE = 'Peakfit', FONT = bfont, $
                           EVENT_PRO  ='ximview_maxfit',  $
                       TOOLTIP = 'Fits peak near marked point')
profile    = WIDGET_BUTTON(button_row, VALUE = 'Profile', FONT = bfont, $
                           EVENT_PRO = 'ximview_profile', $
                           TOOLTIP = 'Plot 1D profile along slice')
; pad2 is another convenient peg:
pad2  = WIDGET_BASE(button_row, EVENT_PRO = 'ximview_deltab')
exit  = WIDGET_BUTTON(button_row, VALUE = 'Exit', FONT=bfont, $
                      EVENT_PRO = 'ximview_exit')
;
; basic widget build finished.
;
; Adjust button sizes to look neat:
;
ogeom = WIDGET_INFO(overview,/GEOMETRY)
WIDGET_CONTROL, exit, XSIZE=ogeom.XSIZE

bgeom = WIDGET_INFO(base1, /GEOMETRY)
brgeom = WIDGET_INFO(button_row, /GEOMETRY)
xgeom = WIDGET_INFO(exit,/GEOMETRY)
length = xgeom.XOFFSET + xgeom.XSIZE + brgeom.XPAD
dx  = FIX(bgeom.XSIZE - bgeom.XPAD - length)
WIDGET_CONTROL, pad1, XSIZE = dx/2
WIDGET_CONTROL, pad2, XSIZE = dx - (dx/2)
;
; Record geometry details in mode structure
;
mode = {drag: 0B, $             ; true if currently dragging
        pan: 1B,  $             ; true if panning enabled
        overview: 1B, $         ; true if in overview mode
        roll: state.ROLL, $     ; true if rolling is currently enabled
        new_view: 1B, $         ; true if re-load of pixmaps is needed
        done: 1B, $             ; true if all pixmaps up to date
        blink: 0B, $            ; true if blinking is going on
        tab: 0B, $              ; currently visible tab
        bbase: -1L, $           ; base widget of blink tab
        bwin:  -1L, $           ; Window index of main blink window
        bswin: -1L, $           ; Window index of blink scalebar window
        period: 0.5, $          ; current blink time/tab (seconds)
        do_scale: 1, $          ; blink scale as well as image
        zoom_factor: 1.0, $     ; current zoom factor
        zoom: 0, $              ; log_2 of zoom factor
        zfac: 1, $              ; zoom factor or 1/zoom factor
        resamp: [1.0, 1.0], $   ; resampling factors for overview mode.
        corner: [0, 0], $       ; screen coord of image BLC in overview mode
        x_centre: 0, y_centre: 0, $ ; image coord of pixel at screen centre
        xpix: 0, ypix: 0, $     ; image coords of pixel under cursor
        xpt: -1., ypt: -1., $   ; image coords of marked pixel
        xpt1: -1., ypt1: -1., $ ; previous marked pixel
        xpt2: -1., ypt2: -1., $ ; 2-previous marked pixel
        xhalf: 0, yhalf: 0, $   ; screen coord of screen centre
        oxtv: xsize/2, oytv: ysize/2, $ ; screen coords of last cursor pos.
        grid_plot: 0B, $        ; If true, plot coordinate grid
        grid_calc: 1B, $        ; If true, auto-set grid increment
        grid_step: [-1d0, -1d0], $  ; Intervals between grid lines [degrees]
        catplot: 1B, $          ; If true, plot catalogue(s).
        catalog: PTR_NEW(/ALLOCATE_HEAP) $ ; Pointers to pointers (!) 
                                           ; to catalogues
       }

tlb_geom = WIDGET_INFO(tlb, /GEOMETRY)
DEVICE, GET_SCREEN_SIZE = screen

xmargin = (screen[0]- tlb_geom.SCR_XSIZE)/2
ymargin = (screen[1]- tlb_geom.SCR_YSIZE)/2
xoff = 100 < xmargin

WIDGET_CONTROL, tlb, TLB_SET_XOFFSET = xoff
WIDGET_CONTROL, tlb, /REALIZE
WIDGET_CONTROL, tlb, SENSITIVE = 0   ; Wait until finished setting up
WIDGET_CONTROL, tabs, SET_UVALUE = mode

FOR itab = 0, ntab-1 DO BEGIN
    WIDGET_CONTROL, tab_arr[itab].draw, GET_VALUE = index
    tab_arr[itab].window = index
    WIDGET_CONTROL, tab_arr[itab].scale, GET_VALUE = index
    tab_arr[itab].scale_index = index
ENDFOR

*state.BLINK_SEQ = LINDGEN(ntab)

tlb_geom = WIDGET_INFO(tlb, /GEOMETRY)
;
; Fill out state structure:
;
IF is_unix THEN BEGIN
    state.XSIZE     = tlb_geom.SCR_XSIZE ; Use SCR values as they include menu
    state.YSIZE     = tlb_geom.SCR_YSIZE ; bar and also the "size" returned by
    state.NEWWIDTH  = tlb_geom.SCR_XSIZE ; re-size events reflects SCR
    state.NEWHEIGHT = tlb_geom.SCR_YSIZE ; values under Linux.
    state.MBAR      = tlb_geom.SCR_YSIZE - tlb_geom.YSIZE
ENDIF ELSE BEGIN
    state.XSIZE     = tlb_geom.XSIZE ; But under MS WIndows they don't
    state.YSIZE     = tlb_geom.YSIZE
    state.NEWWIDTH  = tlb_geom.XSIZE
    state.NEWHEIGHT = tlb_geom.YSIZE
ENDELSE
state.TABS    = tabs       &  state.TABARR  = PTR_NEW(tab_arr)
state.READOUT = readout    &  state.LABEL   = label
state.READLAB = readlab    &  state.ZOOMCOL = zoomcol
state.ZOOMIN  = zoomin     &  state.ZOOMOUT = zoomout        
state.PAD1    = pad1       &  state.PAD2    = pad2
state.OV_BUTTON = overview &  state.BLINK   = blink
state.IMSTATS = imstats    &  state.PEAKFIT = peakfit & state.PROFILE = profile
state.DISP = display       &  state.progcol = progcol
state.FRAMES  = frames     &  state.DATASETS = datasets

IF tlb_geom.SCR_YSIZE LT tlb_geom.YSIZE THEN BEGIN ; widget does not fit!
   state.NEWHEIGHT = state.YSIZE - 128
   ximview_resize, tlb, state
                                ; try again
   tlb_geom = WIDGET_INFO(tlb, /GEOMETRY)
   IF is_unix THEN state.MBAR = tlb_geom.SCR_YSIZE - tlb_geom.YSIZE
   IF state.MBAR LT 0 THEN MESSAGE, 'Cannot fit Ximview on the screen'
ENDIF

END
