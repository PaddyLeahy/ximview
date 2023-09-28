; PADSTRING
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
FUNCTION padstring, instring, len, just, FORMAT=format
;
; Pads a string with blanks
; 
; INPUTS
;      instring string to pad
;      len      (integer) required length in characters
;      just     (string): justification code:
;               'L':  justify left, add spaces on right
;               'R':  justify right, add spaces on left
;               'C': centre, if odd number of pads, put
;                     extra space on right
;               'CR': centre, any extra space on left.
;
; KEYWORD INPUT
; TODO: implement format option
;      FORMAT   if set, instring is the middle of a format statement and 
;               the output is also a format statement that will print the 
;               required number of blanks on either side. 
;
len_in = STRLEN(instring)
IF len_in GT len THEN MESSAGE, 'not enough space for input string'

CASE STRMID(just,0,1) OF
    'L': fmt ='(A-'+STRTRIM(STRING(len), 2)+')'
    'R': fmt ='(A'+STRTRIM(STRING(len), 2)+')'
    'C': BEGIN
       pad = INTARR(2)
       extra = len - len_in
       pad[0] = extra / 2
       pad[1] = extra - pad[0]
       IF STRMID(just,1,1) EQ 'R' THEN pad = REVERSE(pad)
       fmt  = pad[0] GT 0 ? '(' + STRTRIM(STRING(pad[0]), 2) + "(' ')," : '('
       fmt += pad[1] GT 0 ? 'A,' + STRTRIM(STRING(pad[1]), 2) + "(' '))" : 'A)'
    END
    ELSE: MESSAGE, 'Justification code "'+just+'" not recognised.'
ENDCASE

RETURN, STRING(instring, FORMAT=fmt)

END
