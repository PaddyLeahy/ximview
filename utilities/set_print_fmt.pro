; SET_PRINT_FMT
;
; Copyright (C) 2014 J. P. Leahy 
;
;    Thes file is part of XIMVIEW
;
;    XIMVIEW is free software: you can redistribute it and/or modify
;    it under the terms of the GNU General Public License as published by
;    the Free Software Foundation, either version 3 of the License, or
;    (at your option) any later version.
;
;    This program is distributed in the hope that it will be useful,
;    but WITHOUT ANY WARRANTY; without even the implied warranty of
;    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
;    GNU General Public License for more details.
;
;    You should have received a copy of the GNU General Public License
;    along with this program.  If not, see <https://www.gnu.org/licenses/>.
;
PRO set_print_fmt, tabstr
;
; Creates the format string for printing out pixel values in Ximview, based
; on information stored in the standard tab structure
; 
is_unix = STRCMP(!version.OS_FAMILY,'UNIX', 4,/ FOLD_CASE)

absrange = tabstr.ABSRANGE
sdev     = tabstr.SDEV
mult     = tabstr.MULT

; Find (max) number of figures before decimal point:
nf1 = CEIL(ALOG10(mult*absrange[1]))
nf2 = 0
IF absrange[0] LT 0.0 THEN $
  nf2 = 1 + CEIL(ALOG10(-mult*absrange[0]))
nf = MAX([nf1,nf2])

; dynamic range
; TODO: if error is available, replace sdev with error
drr = sdev NE 0.0 ? ALOG10(MAX(ABS(absrange))/sdev) : 4
dr = FINITE(drr) ? CEIL(drr) : 4

; number of decimal places needed to display numbers with values 
; of order the standard deviation (and/or error):
ndp = CEIL(dr) - nf + 2
IF ndp LE 0 THEN fmt = "F6.0,5(' ')" ELSE BEGIN
      ; number of characters required in display
    nn = nf + 1 + ndp 
    IF nn LE 11 THEN BEGIN ; F format is possible
                     ; If possible, put decimal point at character # 6
        ns = 5 - ndp
        tail = ns GT 0 ? ",'"+STRJOIN(REPLICATE(' ',ns),'')+"'" : ''
        ns = ns GE 0 ? 5 - nf : 11 - nn
        fmt = 'F'+STRTRIM(STRING(nn+ns),2)+'.'+STRTRIM(STRING(ndp),2)+tail
    ENDIF ELSE fmt = is_unix ? 'E11.4' : 'E11.3'
ENDELSE
tabstr.FORMAT = '(1X,'+fmt+')'

END
