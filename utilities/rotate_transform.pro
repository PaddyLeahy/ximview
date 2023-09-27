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

