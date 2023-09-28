; ROTATE_TRANSFORM
;
; Copyright (C) 2020 J. P. Leahy 
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
FUNCTION rotate_transform, vec, angle_in
; constructs rotation matrix for active rotation around axis vec
; by angle angle_in. (Right-hand sense).
;
; vec MUST be normalised on input (not checked, for speed).
;
angle = DOUBLE(angle_in[0])
sa = SIN(angle)  &  ca = COS(angle)
ica = 1d - ca
      
diag = (vec*vec)*ica + ca
sinpart = REVERSE(REFORM(vec))*sa

cospart = [vec[0]*vec[1:2],vec[1]*vec[2]]*ica

cps = cospart + sinpart
cms = cospart - sinpart

entry = [[diag[0],  cms[0],  cps[1]], $
         [ cps[0], diag[1],  cms[2]], $
         [ cms[1],  cps[2], diag[2]]]

RETURN, entry

END

