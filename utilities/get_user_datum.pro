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
FUNCTION get_user_datum, event, question, field_size, default
;+
; NAME:
;       GET_USER_DATUM
;
; PURPOSE:
;       Launches a modal dialog to get a piece of info from the user
;
; CATEGORY:
;       Widget Routines, Compound
;
; CALLING SEQUENCE:
;
;       Result = GET_USER_DATUM(Event, Question, Field_size, Default)
;
; INPUTS:
;       Event:      Unused. Present for compatibility with old version.
;       Question:   Text for label may be array
;       Field_size: Size of text box(es) for answer, in characters.
;
; OPTIONAL INPUTS:
;       Default:    Default value placed in box(es).
;
; OUTPUTS:
;       Result is either the contents of the text box or 'Cancel' if
;       Cancel button was pressed.
;
; EXAMPLE:
;
;       name = GET_USER_DATUM('Please enter your name',30)
;
; MODIFICATION HISTORY:
;       Written by:      J P Leahy March 2008 
;                        (complete re-write of earlier version).
;-
COMPILE_OPT IDL2, HIDDEN

nq = N_ELEMENTS(question)
IF nq EQ 0 THEN MESSAGE, 'no question to ask'
ns = N_ELEMENTS(field_size)
IF ns EQ 1 && nq GT 1 THEN field_size = REPLICATE(field_size,nq) ELSE $
   IF ns NE nq THEN MESSAGE, 'should be same number of questions & field sizes'

nd = N_ELEMENTS(default)
IF nd EQ 0 THEN BEGIN 
   default = ''
   nd = 1
ENDIF
IF nd EQ 1 && nq GT 1 THEN default = REPLICATE(default,nq) ELSE $
   IF nd NE nq THEN MESSAGE, 'Should be same number of questions & defaults'

; remove illegitimate character to make tagnames
reserved_char = '!"$%^&*(){}[]<>-+=,;:./\|?@''~# '
tagname = STRJOIN(STRSPLIT(question,reserved_char,/EXTRACT),'_')
FOR iq=0,nq-1 DO BEGIN
; Escape backslashes and commas in the input strings:
   qq = STRJOIN(STRSPLIT(question[iq],'\',/EXTRACT,/PRESERVE_NULL),'\\')
   question[iq] = STRJOIN(STRSPLIT(qq,',',/EXTRACT,/PRESERVE_NULL),'\,')
   dd = STRJOIN(STRSPLIT(default[iq],'\',/EXTRACT,/PRESERVE_NULL),'\\')
   default[iq] = STRJOIN(STRSPLIT(dd,',',/EXTRACT,/PRESERVE_NULL),'\,')
ENDFOR
qlen = STRLEN(question)
fmt = "(A-"+STRTRIM(STRING(MAX(qlen)),2)+",':')"
qq = STRING(question, FORMAT=fmt)
level = nq GT 1 ? [REPLICATE('0',nq-1),'2'] : ['2']
fs = STRTRIM(STRING(field_size),2)
desc = ['1,BASE,,ROW']
FOR iq=0,nq-1 DO desc = [desc, $
       level[iq]+', TEXT, '+default[iq]+', TAG='+tagname[iq]+', WIDTH=' $
                         +fs[iq] + ',LABEL_LEFT='+qq[iq] ] 
desc = [desc,'1, BASE,,ROW', '0, BUTTON, Accept, QUIT, TAG=OK', $
        '2, BUTTON, Cancel, QUIT']
;nline = N_ELEMENTS(desc)
;FOR i=0,nline-1 do print, desc[i]+'"'

result = CW_FORM(desc, TITLE = 'Enter data',/COLUMN)
IF nq EQ 1 THEN RETURN, result.OK ? result.(0) : 'Cancel'  ELSE BEGIN
   RETURN, result.OK ? result : 'Cancel'  
ENDELSE

END

