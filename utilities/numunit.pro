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
;  Module NUMUNIT
;  
;  J. P. Leahy 2008
;
;  Contents:
;  get_power:  utility routine: interprets "raised to the power"
;  lcunit:     utility routine: converts unit names to standard form
;  numunit:    main function.
; 
;  A standard IDL documentation block follows the declaration of
;  numunit.
;
PRO get_power, string, power, irest, valid
;
; Gets the exponent of a valid FITS unit "raised to the power of" expression
;
; Input: string: contains unit from the point the "power of" part starts
; Output: irest: first character position in string after the
;                expression
;         power: the power
;         valid: 2 if there is a valid integer power expression
;                1 if there is valid fractional power expression
;                0 if there is not a valid expression
;
COMPILE_OPT IDL2, HIDDEN

irest = 0  & power = 1 & valid = 0
IF STRLEN(string) EQ 0 THEN RETURN

IF STRCMP(string,'**',2) THEN ksta = 2 ELSE $
  IF STRCMP(string,'^',1) THEN ksta = 1 ELSE ksta = 0

brace = STRCMP(STRMID(string,ksta),'(',1)
IF brace THEN BEGIN
    ksta = ksta + 1
    kend = STREGEX(string,')') - 1
    IF kend LT ksta THEN MESSAGE, /INFORMATIONAL, $
        'Mismatched parentheses in unit string'

    ktext = STRMID(string,ksta,kend-ksta+1)
    dot = STREGEX(ktext,'\.',/BOOLEAN)
    solidus = STREGEX(ktext,'/')
    ratio = solidus GT -1
    irest = kend + 2
ENDIF ELSE BEGIN
    dot = 0  & ratio = 0
; Find first char not 0 to 9 or + or -
    klen = STREGEX(STRMID(string,ksta),'[^0-9+-]')
    IF klen GE 0 THEN BEGIN 
      ktext = STRMID(string,ksta,klen)
      irest = ksta + klen
    ENDIF ELSE BEGIN ; nothing left after exponent.
      ktext = STRMID(string,ksta)
      irest = STRLEN(string)
    ENDELSE
ENDELSE
valid = 2 - dot - ratio

IF irest EQ 0 || valid EQ 0 THEN BEGIN
    valid = 0
    irest = 0
    RETURN
ENDIF

CATCH, error_code
IF error_code NE 0 THEN BEGIN ; Problem converting string.
    CATCH, /CANCEL
    power = 1
    valid = 0
    irest = 0
    PRINT, 'unconvertible ktext: "'+ktext+'"'
ENDIF

IF valid EQ 2 THEN power = FIX(ktext) ELSE IF dot THEN $
  power = DOUBLE(ktext) ELSE BEGIN
    p1 = DOUBLE(STRMID(ktext,0,solidus))
    p2 = DOUBLE(STRMID(ktext,solidus+1))
    power = p1/p2
ENDELSE

RETURN

END

FUNCTION lcunit, unit_in, uc, prefix_allowed, matched
;
; Input purports to be a prefixless unit. Returns standard abbreviation
; of unit if recognised, otherwise returns the unit unchanged.
;
; INPUTS:
;         unit_in   String containing possible unit (no whitespace)
;         uc        If true, unit_in may have been upper-cased
; OUTPUTS:
;        prefix_allowed  true if the unit can take multiplier prefixes
;                        If unit is not matched, set to false
;        matched         0: no match; 1: recognised; 2: ambiguous
;                        3: is a standard function
;
COMPILE_OPT IDL2, HIDDEN

; Units allowed by FITS standard that accept multiplier prefixes:
standard = ['m','g','s','rad','sr','K','A','mol','cd', $
         'Hz','J','W','V','N','Pa','C','Ohm','S','F','Wb','T','H','lm','lx', $
         'a','yr','eV','pc','Jy','mag','R','G','barn','bit','byte']
; Full names of units in standard list accepting prefixes:
name  = ['metre', 'gram', 'second', 'sec', 'radian', 'steradian', 'str', $
         'kelvin', 'ampere', 'mole', 'candela', 'hertz', $
         'joule', 'watt', 'volt', 'newton', 'Pascal','coulomb', 'ohm', $
         'siemens', 'farad', 'weber', 'tesla', 'henry', 'lumen', 'lux', $
         'year', 'electronvolt','parsec', 'jansky', 'magnitude', 'rayleigh', $
         'gauss']
; Corresponding symbols:
name_sym = ['m','g','s','s','rad','sr','sr','K','A','mol','cd', $
         'Hz','J','W','V','N','Pa','C','Ohm','S','F','W','T','H','lm','lx', $
         'yr','eV','pc','Jy','mag','R','G']

; Sexagesimal units:
std_sex  = ['deg','arcmin','arcsec','mas','min','h']
name_sex = ['degree','arcminute','arcsecond','minute','hour']
nns_sym  = ['deg','arcmin','arcsec','min','h','d']

; Units allowed by FITS standard that do not accept prefixes:
std_nopref = [std_sex,'d','erg','Ry', $
              'solMass','u','solLum','Angstrom','solRad','AU','lyr','count', $
              'ct','photon','ph','pixel','pix','D','Sun','chan','bin', $
              'voxel','adu','beam','unknown']
name_nopref = [name_sex,'day','rydberg','Msun','Lsun','Rsun','lightyear', $
               'debye', 'channel']
nnp_sym     = [nns_sym,'d','Ry','solMass','SolLum','solRad','lyr','D','chan']

functions = ['log','ln','exp','sqrt']

first = 1B

; split off comments/subscripts to units
subs_idx = STREGEX(unit_in+',','[,|_]')
unit = STRMID(unit_in,0,subs_idx)
subscript = STRMID(unit_in,subs_idx)

REDO: 

matched = 1
length = STRLEN(unit)

IF length EQ 1 THEN BEGIN
; Unit names that are ambiguous when uppercased:
    ambig = ['G','S','A','H','D']
    test = WHERE(STRUPCASE(unit) EQ ambig)
    IF test NE -1 THEN BEGIN
        IF uc THEN matched = 2
        prefix_allowed = uc ? test NE 4 : unit NE 'h' && unit NE 'd'
        RETURN, unit_in
    ENDIF
ENDIF

IF uc THEN BEGIN
   unituc = STRUPCASE(unit)

   prefix_allowed = 1B
   test = WHERE(unituc EQ STRUPCASE(standard)) 
   IF test NE -1 THEN RETURN, standard[test]+subscript
   test = WHERE(unituc EQ STRUPCASE(name))
   IF test NE -1 THEN RETURN, name_sym[test]+subscript

   prefix_allowed = 0B
   test = WHERE(unituc EQ STRUPCASE(std_nopref))
   IF test NE -1 THEN RETURN, std_nopref[test]+subscript
   test = WHERE(unituc EQ STRUPCASE(name_nopref))
   IF test NE -1 THEN RETURN, nnp_sym[test]+subscript

   matched = 3
   test = WHERE(unituc EQ STRUPCASE(functions))
   IF test NE -1 THEN RETURN, functions[test]+subscript
ENDIF ELSE BEGIN
   prefix_allowed = 1B
   test = WHERE(unit EQ standard) 
   IF test NE -1 THEN RETURN, standard[test]+subscript
   test = WHERE(unit EQ name)
   IF test NE -1 THEN RETURN, name_sym[test]+subscript

   prefix_allowed = 0B
   test = WHERE(unit EQ std_nopref)
   IF test NE -1 THEN RETURN, std_nopref[test]+subscript
   test = WHERE(unit EQ name_nopref)
   IF test NE -1 THEN RETURN, nnp_sym[test]+subscript

   matched = 3
   test = WHERE(unit EQ functions)
   IF test NE -1 THEN RETURN, functions[test]+subscript
ENDELSE

IF first THEN BEGIN
   last = STRMID(unit_in,length-1,1)
   plural = uc ? last EQ 'S' : last EQ 's'
   IF plural THEN BEGIN
                                ; may be plural; strip S & try again
      unit = STRMID(unit_in,0,length-1)
      first = 0B
      GOTO, REDO
   ENDIF
ENDIF

matched = 0
RETURN, unit_in

END
;
PRO analyse_unit, un, oldpref, un0, set_pref
;
; Analyses text part of unit into elements consisting of multiplier
; prefices and remaining unit. Converts prefices and units into standard form.
;
;
; INPUT/OUTPUT
;     un       On input, text string purporting to be a FITS-compliant unit. 
;              on output, converted to FITS standard if needed.
; OUTPUT
;     oldpref  scaling prefix in FITS standard form
;     un0      First element of unit in fits standard form, without prefix
;     set_pref If true (1B) try to set a scaling prefix on the output unit

; Constants and defaults
oldpref = ''
breaks = ' *./^0123456789+-()'
functions = ['log','ln','exp','sqrt']

; SI prefixes (powers of 1000):
prefix = ['y','z','a','f','p','n','u','m','','k','M','G','T','P','E','Z','Y']
; Non-standard multipler prefixes (powers of 10):
altpref = ['c','d','','da','h']
pf1 = [prefix,altpref]

uc = un EQ STRUPCASE(un) ; If true, unit is in CAPS (very confusing!)
sure_uc = 0B  ; set to 1 if we are sure units have been uppercased.

; Break unit string into words, each of which should have the structure 
; <prefix>unit, unless they are a function specifier:

word = STRSPLIT(un, breaks, /EXTRACT, COUNT = count)
IF count EQ 0 THEN BEGIN   ; null unit
    set_pref = 0B
    un0 = un
    RETURN
ENDIF

IF uc THEN BEGIN
; Spelled-out prefixes:
    prefuc = ['YOCTO','ZEPTO','ATTO','FEMTO','PICO','NANO','MICRO','MILLI', $
              '', 'KILO','MEGA','GIGA','TERA','PETA','EXA','ZETTA','YOTTA']
    apfuc  = ['CENTI','DECI','','DECA','HECTO']
    prefs = [prefuc, apfuc]
    preflen = STRLEN(prefs)
    pref3 = STRMID(prefs,0,3) ; first 3 letters of prefix
                              ; unambiguous except for DEC
ENDIF

; Loop over words in unit string, separating out prefix and unit and 
; converting each to standard form if recognised. 
; Remember result for first word.
FOR j = 0,count-1 DO BEGIN 
    matched = 0
    pref0 = ''
    gotpref = -1
    current = word[j]
    unlen = STRLEN(current)
    IF current EQ 'Pa' THEN BEGIN   ; awkward case: not Peta annum !
       gotpref = 1
       matched = 1
       prefix_allowed = 1B
    ENDIF

    IF uc && current EQ 'EV' THEN BEGIN ; eV or Exa volt ?
        gotpref = 1
        matched = 1
        prefix_allowed = 1B
        IF N_ELEMENTS(evs) EQ 0 THEN evs = j ELSE evs = [evs,j]
    ENDIF
    
    IF uc && unlen GE 4 THEN BEGIN       ; look for spelled-out prefix
        gotpref = WHERE(STRMID(current,0,3) EQ pref3, nmatch)
        IF nmatch EQ 2 THEN CASE STRMID(current,0,4) OF
            'DECI': gotpref = gotpref[0]
            'DECA': gotpref = gotpref[1]
            ELSE: gotpref = -1
        ENDCASE
        IF (gotpref GT -1) && $
            (STRMID(current,0,preflen[gotpref]) EQ prefs[gotpref]) THEN BEGIN
            pref0 = pf1[gotpref]
            current = STRMID(current,preflen[gotpref]) 
        ENDIF ELSE gotpref = -1
    ENDIF
    ; If gotpref > -1, current is now a pure unit, otherwise 
    ; we have at most a two-character prefix.
    IF gotpref EQ -1 && unlen GE 3 THEN BEGIN
        p2 = STRMID(current,0,2)
        deca = uc ? p2 EQ 'DA' : p2 EQ 'da'
        IF deca THEN BEGIN
            trial = lcunit(STRMID(current,2), uc, prefix_allowed, matched)
            IF matched GT 0 && prefix_allowed THEN BEGIN
                pref0 = 'da'
                gotpref = 1
                current = trial[0]
            ENDIF ELSE matched = 0
        ENDIF
    ENDIF
    IF gotpref EQ -1 && unlen GE 2 THEN BEGIN ; one-char prefix or pure unit
        p1 = STRMID(current,0,1)
        IF uc THEN BEGIN
            gotpref = WHERE(p1 EQ STRUPCASE(pf1),nmatch)
                         ; choose uppercase prefix if ambiguous:
            IF nmatch EQ 2 THEN gotpref = gotpref[1]
        ENDIF ELSE gotpref = WHERE(p1 EQ pf1, nmatch)
        IF nmatch GT 0 THEN BEGIN
            trial = lcunit(STRMID(current,1), uc, prefix_allowed, matched)
            IF matched GT 0 && prefix_allowed THEN BEGIN
                pref0 = pf1[gotpref]
                current = trial[0]
            ENDIF ELSE matched = 0
        ENDIF
    ENDIF  
    IF matched EQ 0 THEN BEGIN
        trial = lcunit(current, uc, prefix_allowed, matched)
        IF matched GT 0 && (pref0 EQ '' || prefix_allowed) THEN $
            current = trial[0] ELSE matched = 0
    ENDIF
    IF matched EQ 0 THEN BEGIN
       current = word[j]
       pref0 = ''
       prefix_allowed = 0B
    ENDIF
    IF j EQ 0 THEN BEGIN
        oldpref = pref0[0]
        un0 = current
        set_pref = prefix_allowed
        is_func = matched EQ 3
    ENDIF
    sure_uc = sure_uc || (uc && pref0+current NE word[j])
    word[j] = pref0 + current
ENDFOR

; Decide eV/EV ambiguity: assume eV if we're sure of uppercasing 
IF N_ELEMENTS(evs) GT 0 THEN BEGIN
    IF sure_uc THEN BEGIN
        word[evs] = 'eV' 
        IF evs[0] EQ 0 THEN un0 = 'eV'
    ENDIF ELSE IF evs[0] EQ 0 THEN BEGIN
        oldpref = 'E'
        un0 = 'V'
    ENDIF
ENDIF
                                ; Now put string back together
start = STRSPLIT(un, breaks, LENGTH = len)
IF start[0] GT 0 THEN new = STRMID(un,0,start[0]) ELSE new = ''
FOR i = 0,count-2 DO BEGIN
    gap0 = start[i] + len[i]
    fill = STRMID(un, gap0, start[i+1] - gap0)
    new += word[i] + fill
ENDFOR
new = new + word[count-1]
gap0 = start[count-1]+len[count-1]
lengap = STRLEN(un) - gap0
IF lengap GT 0 THEN new += STRMID(un,gap0,lengap)
un = new

END
;
FUNCTION numunit, numin, unit, errin, PRECISION = prec_in, DECIMALS = dec_in, $
    FORCE = force, MULTIPLIER = multiplier, OUT_UNIT = out_unit
;+
; NAME:
;       NUMUNIT
;
; PURPOSE: 
;       Parses (number, unit) pairs to a string suitable for
;       printing. A multiplier prefix is chosen for the unit to put
;       the number in the range [1,1000). Alternatively, parses
;       (number, error, unit) sets. Assumes unit is encoded in
;       the FITS convention (See the FITS standard V3.0 at
;       fits.gsfc.nasa.gov), generalised in that

;       o  "<unit>s" is acceptable as well as "<unit>" 
;          (unless unit is a single character),
;       o  Spelled-out names of some units are allowed.
;       o  unit names which have been upper-cased are allowed,
;          provided they are not too ambiguous. For upper-case unit
;          descriptors, we also search for spelled-out multiplier
;          prefixes, eg. "MILLI".
;       o  If a unit name contains a comma or an underscore (e.g. "K_CMB"), 
;          this is assumed to denote a comment or subscript and only the text
;          preceding it is considered as the unit proper. 
;       o  The output unit name is in the FITS standard form
;          (with any comment/subscript preserved unchanged).
;
;       If the leading component of the unit (E.g. "count" in
;       "count/s") is one listed in the FITS standard as being
;       inappropriate for multiplier prefixes, or if it is not
;       recognised, the unit is returned unchanged unless /FORCE is
;       specified, in which case a numerical prefix is supplied.
;
;       Optional outputs give the conversion factor between input and
;       output units (a power of ten) and the output unit. These can
;       be used for formatting tabular data etc.
;
; CATEGORY:
;       Input/Output, String Processing
;
;
; CALLING SEQUENCE:
;
;       Result = NUMUNIT(Number, Unit)
;
; INPUTS:
;       Number:  Number to be formatted (assumed float or double)
;       Unit:    String containing units of number, may already
;                contain a prefix
;
; OPTIONAL INPUT
;
;       error:   Error in number.
;
; KEYWORD PARAMETERS:
;       PRECISION:  Number of significant figures (integer, < 20).
;                   Default: 3. Ignored if ERROR is specified
;       DECIMALS:   Number of decimal places (integer, < 19)
;                   (alternative to PRECISION).
;       FORCE:      Re-scale even units which should not have
;                   prefixes, by including "10^n" in the unit string.
;
; OUTPUTS:
;      Result is a string containing rescaled number and unit, with
;      error if specified.
;
; OPTIONAL OUTPUTS:
;       MULTIPLIER: Scale factor between input and output number:
;                   output = input * multiplier
;       OUT_UNIT:   String containing just the output unit (including
;                   prefix if any). 
;
; EXAMPLE:
;       Simple use:
;
;       IDL> PRINT, NUMUNIT(0.00001281,'  Ohms')
;       12.8 uOhm
;   
;       The illegitimate terminal "s" has been removed. (ohm-seconds
;       should be coded as "Ohm s", "Ohm.s" or "Ohm*s"). 
;
;       IDL> PRINT, NUMUNIT(1.7377e7,'JY/PIXEL', DEC=2)
;         17.38 MJy/pixel
;
;       (Such upper-case units may occur in old FITS headers)
;       Note that the leading blanks are not stripped when DECIMALS is
;       set, to facilitate lining results up in columns.
;
;       IDL> PRINT, NUMUNIT(!pi,'Hz', 0.345, DEC = 5)
;          3.14 +/- 0.34 Hz
;      
;       (DEC and PREC c parameters are ignored if error is
;       specified. Errors are reported to 2 significant figures.)
;
;       To choose a format and unit for a table, run NUMUNIT on the
;       maximum value to ensure there are no overflows, and specify
;       plenty of precision to allow for smaller values:
;
;       IDL> data = 100.0*EXP(RANDOMN(seed, 100))
;       IDL> PRINT, NUMUNIT(MAX(data),'Hz', PREC=5, MULT=mult, OUT=unit)
;       1.9922 kHz
;       IDL> PRINT, unit, mult
;       kHz   0.00100000
;       IDL> PRINT, 'Frequency', '('+unit+')', FORMAT = "(A)"
;       Frequency
;       (kHz)
;       IDL> PRINT, mult*data, FORMAT = "(F7.4)"
;        0.0763
;        0.0816
;       ...etc
;
; MODIFICATION HISTORY:
;       Written by:     J. P. Leahy, Jan-Feb 2008
;       Version 2 July 2013: simplified logic & more thorough conversion
;                            to standard forms.
;       Version 2.1 Sept 2020: added error input.
;       Version 2.2 Aug  2022: fixed bug in format statement
;
;-
COMPILE_OPT IDL2, HIDDEN
ON_ERROR, 0

; Constants & defaults:

; SI prefixes (powers of 1000):
prefix = ['y','z','a','f','p','n','u','m','','k','M','G','T','P','E','Z','Y']

prefexp = ''
prefout = ''
untail = ''
power = 0 ; Default prefix power of ten

; Check inputs:
nerr = N_ELEMENTS(errin)
IF N_ELEMENTS(numin) NE 1 || N_ELEMENTS(unit) NE 1 || nerr GT 1 THEN $
   MESSAGE, 'Inputs must be scalar' 

; Get power of ten: note that we use original-precision input because
; rounding errors in converting floats to doubles can give the wrong
; answer if we use the converted value "number" 
inpower = numin EQ 0.0 ? 0 : FLOOR(ALOG10(ABS(numin)))

do_err = nerr EQ 1

; Set required precision
IF do_err THEN BEGIN
   number = DOUBLE([numin,errin])
   precset = 1B
   eprec =  2  ; 2 sig. figs for error
   epower = FLOOR(ALOG10(ABS(errin)))
   precision = inpower - epower + eprec
ENDIF ELSE BEGIN
   number = DOUBLE(numin)
   decset  = N_ELEMENTS(dec_in) GT 0
   IF decset THEN decimals = CEIL(dec_in)
   IF decset  THEN BEGIN
      IF decimals GT 19 THEN MESSAGE, 'Max 19 decimal places'
      IF decimals LT 0 THEN MESSAGE, '# decimal places should be >= 0'
   ENDIF
   precset = N_ELEMENTS(prec_in) GT 0
   IF precset THEN precision = CEIL(prec_in)
   IF precset THEN BEGIN
      IF precision GT 20 THEN MESSAGE, 'Max 20 significant figures'
      IF precision LT 0 THEN MESSAGE, '# significant figures should be >= 0'
   ENDIF
   IF decset AND precset THEN MESSAGE, 'Set only one of DECIMALS or PRECISION'

   IF decset + precset EQ 0 THEN BEGIN
      precset = 1B
      precision = 3
   ENDIF
ENDELSE

; Set default outputs
premult = 0B
numout = number

force = KEYWORD_SET(force)

un = STRTRIM(unit, 2)

; Now do the interpretation:

; Check for power of ten in unit string. If present, adjust number and
; remove from unit string
tenpos = STREGEX(un,'10')
IF tenpos NE -1 THEN BEGIN
   get_power, STRMID(un,2+tenpos), kpow, irest, valid
   IF valid EQ 2 THEN BEGIN
      number = number * (10d0)^kpow
      inpower += kpow
      un = STRTRIM(STRMID(un,tenpos+irest+2),1)
   ENDIF
ENDIF

; Sort out text part of unit
analyse_unit, un, oldpref, un0, set_pref

word0 = oldpref+un0
wlen = STRLEN(word0)
unlen = STRLEN(un)
IF un NE '' THEN start = STREGEX(un, word0) ELSE start = 0
unnose = STRMID(un,0,start)
untail = STRMID(un,start+wlen)

IF set_pref THEN BEGIN ; Do some checks to see if we can really re-scale:

; See if first unit is included in a (nested) group & if any are raised
; to powers:
   powered = 0B
    
   b1 = STRPOS(unnose,'(')
   group = un
   glen = unlen
   WHILE b1 NE -1 DO BEGIN
      b2 = STRPOS(group,')',/REVERSE_SEARCH) ; matching closing parenthesis
      IF b2 LT b1 THEN MESSAGE, /INFORMATIONAL, $ 
                                'Mismatched parentheses in unit string'
        
      IF b2 LT glen-1 THEN BEGIN
         tail = STRMID(group,b2+1)
         get_power, tail, unpower, irest, valid
         IF valid && unpower NE 1.0 THEN powered = 1B
      ENDIF
      group = STRMID(group,b1+1,b2-b1-1)
      start = STREGEX(group, word0)
      gnose = STRMID(group,0,start)
      b1 = STRPOS(gnose,'(')
   ENDWHILE

; Check to see if leading unit is raised to a funny power
   get_power, untail, unpower, irest, valid
   IF valid && unpower NE 1.0 THEN powered = 1B
   
   IF powered THEN BEGIN
      un0 = word0
      GOTO, NOSCALE
   ENDIF
     
   IF oldpref NE '' THEN BEGIN  ; interpret scaling prefix as power of 10
      pref = WHERE(oldpref EQ prefix)
      pref = pref[0]
      IF pref NE -1 THEN power = 3*(pref - 8) ELSE BEGIN
         pref = WHERE(oldpref EQ altpref)
         IF pref NE -1 THEN power = pref - 2
      ENDELSE
      IF pref EQ -1 THEN MESSAGE, $
         'Internal error: scaling prefix set but not recognised'
   ENDIF
   
   GOTO, PREFSET
ENDIF ELSE BEGIN
   ; might be sexagesimal:
   std_sex  = ['deg','arcmin','arcsec','mas','min','h']
   null = WHERE(un0 EQ std_sex,nsex)
   sexagesimal = nsex EQ 1
   ; nothing is done with this yet.
ENDELSE
; End of setting text part of unit
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

NOSCALE:
    IF force || do_err THEN premult = 1B ELSE BEGIN
       prefout = ''
       numout = number
       IF inpower LT 0 OR inpower GE 3 THEN BEGIN
          IF precset THEN decimals = precision - 1
          fmt = STRING(8+decimals, decimals, FORMAT = "('(E',I2,'.',I1,')')")
          GOTO, FINISH
       ENDIF
    ENDELSE

PREFSET:
    number = number * 10d0^(power)
    pon3 = FLOOR((inpower+power)/3d0)
    outpower = 3*pon3
    IF premult THEN BEGIN                          ; Put x 10^x in front of unit
        pref_index = 8
        IF outpower EQ 0 THEN prefexp = '' ELSE BEGIN
            prefexp = 'x 10^' + STRTRIM(STRING(outpower),1) + ' '
        ENDELSE
    ENDIF ELSE BEGIN                                ; Use SI prefixes
        pref_index = pon3 + 8
        IF pref_index LT 0 OR pref_index GT 16 THEN BEGIN
            number = number * 10d0^(-power)
            GOTO, NOSCALE                           ; out of prefixes
        ENDIF
        prefout = prefix[pref_index]
    ENDELSE
    numout = number * 10d0^(-outpower)

    outlog = numout[0] EQ 0.0 ? 0 : FLOOR(ALOG10(ABS(numout[0])))
    IF precset THEN decimals = precision - outlog - 1

                                ; Watch out for rounding up:
    frac = numout * 10d0^(decimals)
    frac = frac - FIX(frac)
    IF ABS(frac[0]) GT 0.5 && FIX(ABS(numout[0]))+1 GE 1000 THEN BEGIN
        IF premult THEN BEGIN
           numout /= 1d3
           outpower += 3
           prefexp = 'x 10^' + STRTRIM(STRING(outpower),1) + ' '
        ENDIF ELSE BEGIN
           IF pref_index + 1 GT 16 THEN BEGIN
              number = number * 10d0^(-power)
              GOTO, NOSCALE
           ENDIF
           numout  /= 1d3
           prefout = prefix[pref_index + 1]
        ENDELSE
        IF precset THEN decimals = $
           precision - FLOOR(ALOG10(ABS(numout[0])+1)) - 1
    ENDIF

; Create suitable output format for number, and error if present.
    IF decimals LE 0 THEN BEGIN
       fmt = '(I4'
       IF do_err THEN fmt += epower EQ 2 ? ",' +/-',I4" : ",' +/-',I3" 
    ENDIF ELSE BEGIN
       lead = STRTRIM(STRING(5+decimals),2)
       elead = STRTRIM(STRING(3+decimals),2)
       ndp = STRTRIM(STRING(decimals),2)
       fmt = '(F'+lead+'.'+ndp
       IF do_err THEN fmt += ",' +/-',F"+elead+'.'+ndp
    ENDELSE
    fmt += ')'
    IF decimals LE 0 THEN numout = DOUBLE(ROUND(numout)) 
    
FINISH:

out_unit = prefexp + unnose + prefout + un0 + untail

numstring = STRING(numout, FORMAT = fmt) 
IF precset THEN numstring = STRTRIM(numstring, 2)
IF do_err THEN numstring = '('+numstring+')'
numstring += ' ' + out_unit

multexp = numin NE 0.0 ? ROUND(ALOG10(numout[0]/numin)) : 0
multiplier = 10.0^multexp

RETURN, numstring

END
