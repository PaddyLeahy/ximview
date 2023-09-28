; GET_JUMPS
;
; Copyright (C) 2020, 2023 J. P. Leahy 
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
PRO get_jumps, xx, yy, jump, njump
;
; Finds jumps in grid lines caused by going over the 
; edge of a map or over a gore.
;
; INPUTS
;    xx, yy   Set of coordinates for points on a grid line
;
; OUTPUTS
;    jump     List of indices of the first point in a contiguous segment
;    njump    number of jumps (number of segments is njump + 1).
;
dx = xx - SHIFT(xx,1)
dy = yy - SHIFT(yy,1)
dxy = SQRT(dx^2 + dy^2)
medstep = MEDIAN(dxy)

; For debugging
;PRINT, 'Median and maximum step:', medstep, MAX(dxy)

jump = WHERE(dxy GT (5*medstep > 1), njump)
IF njump EQ 0 THEN jump = [0, N_ELEMENTS(dxy)] ELSE $
   jump = [0, jump, N_ELEMENTS(dxy)]

END
