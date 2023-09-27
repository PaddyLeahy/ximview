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
FUNCTION divup, I, J
;+
; NAME:
;       DIVUP
;
; PURPOSE:
;       Returns integers rounded away from zero, instead of towards
;       zero as in normal integer division. Arguments can be two
;       scalars, two equal-length arrays, or a scalar and an array.
;
; CATEGORY:
;       Mathematics
;
;
; CALLING SEQUENCE:
;
;       Result = DIVUP(I, J)
;
; INPUTS:
;       I:  Numerator (may be an array)
;       J:  Denominator (may be an array)
;
; OUTPUTS:
;       Rounded-up ratio.
;
; EXAMPLE:
;
;          PRINT, DIVUP(1,2), DIVUP(81,9), DIVUP([-1,2,-4],[-3,-4,3])
;       
;       IDL prints:   1           9           1          -1          -2
;
; MODIFICATION HISTORY:
;       Written by:     J. P. Leahy November 2007
;       March 2008:     Modified to accept arrays.
;-
COMPILE_OPT IDL2, HIDDEN

pos = I GE 0
negs = WHERE(~pos)
sign = FIX(pos)
IF negs[0] NE -1 THEN sign[negs] = -1
RETURN, sign * (ABS(I) + ABS(J) - 1) / J
END
