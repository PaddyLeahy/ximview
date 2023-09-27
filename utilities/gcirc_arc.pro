FUNCTION gcirc_arc, ra1, dec1, ra2, dec2, interval, dis
;
; Returns the coordinates of a set of points along the great circle
; arc connecting two sky positions
;
; INPUTS
;  ra1, dec1, ra2, dec2  Coordinates of endpoints (degrees)
;  interval              Interval between circle points (degrees)
;
; OPTIONAL OUTPUT
;  dis                   Great-circle distance in arcsec
;
; RETURNS
;    [[ra],[dec]]        Array of points along circle
;  
gcirc, 2, ra1, dec1, ra2, dec2, dis
dis_deg = dis / 3600d0               ; convert to degrees
npt = CEIL(dis_deg / interval)
dis_rad = dis_deg * !dpi / 180d0 ; convert to radians
ang2vec, [dec1,dec2], [ra1, ra2], vecs, /ASTRO
axis = CROSSP(vecs[0,*],vecs[1,*])
axis = norm_vec(axis,/OVERWRITE)
rots = DINDGEN(npt) * dis_rad/(npt-1) 
vec_circ = DBLARR(npt,3)

FOR i=0,npt-1 DO BEGIN
   mat = ROTATE_TRANSFORM(axis,rots[i])
   vec_circ[i,*] = mat ## vecs[0,*]
ENDFOR

vec2ang, vec_circ, dec, ra, /ASTRO

RETURN, [[ra],[dec]]
END
