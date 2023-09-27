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
FUNCTION get_pix_area, pixxy, astrom
;
; Calculates (approximate) pixel area. Includes several methods, not
; all actually used at the moment.
;
COMPILE_OPT IDL2, HIDDEN

s2r = !dpi / (3600d0*180d0)
r2d = 180d0 / !dpi
pi2 = 0.5d0*!dpi

xoff = 0.5d0*[-1, -1, 1, 1]  &  yoff = 0.5d0*[-1, 1, 1, -1]
XY2AD, pixxy[0] + xoff, pixxy[1] + yoff, astrom, longs, lats

                                ; Check for wrap in longitude:
range = MINMAX(longs)
IF range[1]-range[0] GT 180d0 THEN longs[WHERE(longs GT 180d0)] -= 360d0

                                ; Jacobian approach: accurate if pixel
                                ; boundaries follow coordinate lines,
                                ; and if pixel doesn't straddle +/-180 deg.
theta = (90d0 - MEAN(lats))/r2d
dpdx =  0.5d0 * ( (longs[2]-longs[1]) + (longs[3]-longs[0]) )
dpdy =  0.5d0 * ( (longs[1]-longs[0]) + (longs[2]-longs[3]) )
dtdx = -0.5d0 * ( (lats[2] - lats[1]) + (lats[3] - lats[0]) )
dtdy = -0.5d0 * ( (lats[1] - lats[0]) + (lats[2] - lats[3]) )
omega = SIN(theta) * ABS(dpdx*dtdy - dpdy*dtdx) / r2d^2

; PRINT, 'Omega:', omega

; Get great-circle distances between all corners
GCIRC, 2, longs[0], lats[0], longs[1:3], lats[1:3], distance1
GCIRC, 2, longs[1], lats[1], longs[2:3], lats[2:3], distance2
GCIRC, 2, longs[2], lats[2], longs[3], lats[3], distance3
distance = s2r * [distance1, distance2, distance3]

a = distance[0] ; corner 0 -> corner 1 
b = distance[3] ; corner 1 -> corner 2
c = distance[1] ; corner 2 -> corner 0
d = distance[2] ; corner 0 -> corner 3
e = distance[5] ; corner 2 -> corner 3

IF (a-e) LT 0.001*a AND (b-d) LT 0.001*b THEN BEGIN 
                                ; pretty regular, probably safe to use
                                ; small-angle approximation

    cc1 = ACOS( (a*a + b*b - c*c) / (2d0*a*b))
    sincc1 = SIN(cc1)
    aa = ASIN(sincc1 * (a/c))
    bb = ASIN(sincc1 * (b/c))
    area1 = 0.5d0 * a*b*sincc1

    cc2 = ACOS( (d*d + e*e - c*c) / (2d0*e*d))
    sincc2 = SIN(cc2)
    dd = ASIN(sincc2 * (d/c))
    ee = ASIN(sincc2 * (e/c))
    area2 = 0.5d0 * d*e*sincc2

;    PRINT, 'Angles: ', (bb+ee)*r2d, cc1*r2d, (aa+dd)*r2d, cc2*r2d
;    PRINT, 'areas:', area1, area2, area1+area2
ENDIF ELSE BEGIN 
                                ; Exact result for infinite precision
                                ; and great-circle boundaries,
                                ; but subject to round-off errors  

    s1 = 0.5d0 * (a + b + c)    ; semi-perimeter
    s2 = 0.5d0 * (d + e + c)

    sina = SIN(a)  &  sinb = SIN(b)  & sinc = SIN(c)
    sind = SIN(d)  &  sine = SIN(e)

; Get corner angles avoiding cosines for sides (probably small angles)
    cc1 = 2d0*ASIN( SQRT( SIN(s1-a) * SIN(s1-b) / (sina * sinb) ) )
    sincc1 = SIN(cc1)
    aacomp = ACOS(sincc1 * (sina/sinc) )
    bbcomp = ACOS(sincc1 * (sinb/sinc) )

    cc2 = 2d0*ASIN( SQRT( SIN(s2-d) * SIN(s2-e) / (sind * sine) ) )
    sincc2 = SIN(cc2)
    ddcomp = ACOS( sincc2 * (sind/sinc) )
    eecomp = ACOS( sincc2 * (sine/sinc) )

    area1 = cc1 - aacomp - bbcomp ; Spherical triangle area rule
    area2 = cc2 - ddcomp - eecomp  

;    aa = pi2-aacomp  &  bb = pi2-bbcomp  
;    dd = pi2-ddcomp  &  ee = pi2-eecomp
;    PRINT, 'Angles: ', (bb+ee)*r2d, cc1*r2d, (aa+dd)*r2d, cc2*r2d
;    PRINT, 'areas:', area1, area2, area1+area2
ENDELSE

RETURN, area1 + area2
END

