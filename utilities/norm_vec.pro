FUNCTION norm_vec, vec, OVERWRITE=overwrite
;
; Normalises an array of vectors with shape [nvec,ndim]
;
svec = SIZE(vec)
IF svec[0] GT 0 THEN BEGIN
   ndim = svec[svec[0]] ; number of elements in each vector
   dims = svec[1:svec[0]]
   length = svec[svec[0]+2]
ENDIF ELSE MESSAGE, 'Input is a scalar'

overwrite = KEYWORD_SET(overwrite)

norms = SQRT(TOTAL(vec^2,svec[0]))
IF overwrite THEN BEGIN
   vec = REFORM(vec,length/ndim,ndim,/OVERWRITE)
   FOR i=0,ndim-1 DO vec[*,i] /= norms
   vec = REFORM(vec,dims,/OVERWRITE)
   RETURN, vec
ENDIF ELSE BEGIN
   normvec = REFORM(vec,length/ndim,ndim)
   FOR i=0,ndim-1 DO normvec[*,i] /= norms
   normvec = REFORM(normvec,dims,/OVERWRITE)
   RETURN, normvec
ENDELSE

END
