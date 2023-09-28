; CUBIC
;
; Copyright (C) 2020, 2023 J. P. Leahy 
;
;    This program is free software: you can redistribute it and/or modify
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
FUNCTION cubic, a, b, c, d, CHECK = check
;
; Finds the roots of cubic equations with real coefficients:
; 
;   a*x^3 + b*x^2 + c*x + d = 0
;   
; Caveat emptor: based on Wikipedia algorithm
;
; Input parameters can be arrays (all same size), in which case array
; of equations is solved
;
; Also, b, c, d may be scalar, but a an array.
;
; Returns a real array if all roots are real, otherwise a complex array with 
; the real root as the first element (i.e. imaginary part = 0).
; 
check = KEYWORD_SET(check) ; currently has no effect

na = N_ELEMENTS(a)
nb = N_ELEMENTS(b)
nc = N_ELEMENTS(c)
nd = N_ELEMENTS(d)
ok = na EQ nb AND na EQ nc AND na EQ nd
semiscalar = nb EQ 1 and nc EQ 1 and nd EQ 1
IF ~(ok || semiscalar) THEN $
   MESSAGE, 'Coefficient arrays should have equal sizes'

not_cubic = WHERE(a EQ 0,nnot)
IF nnot GT 0 THEN BEGIN
   IF na EQ 1 THEN MESSAGE, 'Not a cubic: first coefficient zero'
   MESSAGE, STRING(nnot, FORMAT =$ 
        "(I,' coefficient sets are not cubics: first coefficents zero')")
ENDIF

; Trap simple case with zero root:
IF nd EQ 1 && d EQ 0 THEN BEGIN ; cubic = x*quadratic
   delta = b^2 - 4*a*c
   imag = WHERE(delta LT 0d0,nimag,/NULL)
   IF nimag EQ 0 THEN BEGIN
      roots = DBLARR(3,na) 
      rtdelta = SQRT(delta)
   ENDIF ELSE BEGIN
      roots = DCOMPLEXARR(3,na)
      rtdelta = SQRT(DCOMPLEX(delta))
   ENDELSE
   roots[0,*] = 0d0 
   roots[1,*] = (-b+rtdelta) / (2*a)
   roots[2,*] = (-b-rtdelta) / (2*a)
   RETURN, roots
ENDIF

a2 = a*a
b2 = b*b
c2 = c*c
discriminant = 18d0*a*b*c*d - 4d0*b2*b*d + b2*c2 - 4d0*a*c2*c - 27d0*a2*d*d

del0 = b2 - 3d0*a*c 
del1 = 2d0*b^3 - 9d0*a*b*c + 27d0*a2*d

; Generic case: one real root, two complex:
gen = WHERE(discriminant LT 0d0, ngen, COMPLEMENT=special)
cru = [DCOMPLEX(1d0,0d0),DCOMPLEX(-0.5d0, SQRT(0.75d0)), $
       DCOMPLEX(-0.5d0, -SQRT(0.75d0))] ; cube roots of unity
IF ngen GT 0 THEN BEGIN 
   roots = DCOMPLEXARR(3,na)
   C3 = DBLARR(ngen)
   simple = WHERE(del0[gen] EQ 0d0, nsimp, COMPLEMENT=full)
   IF nsimp GT 0 THEN BEGIN
      C3[simple] = del1[gen[simple]]
   ENDIF 
   nf = ngen - nsimp
   IF nf GT 0 THEN BEGIN
      gf = gen[full]
      rt = SQRT(-27d0*a2[gf]*discriminant[gf])
      pos1 = WHERE(del1[gf] GT 0d0, np, COMPLEMENT=neg1)
      IF np GT  0 THEN C3[full[pos1]] = 0.5d0*(del1[gf[pos1]] - rt[pos1])
      IF np LT nf THEN C3[full[neg1]] = 0.5d0*(del1[gf[neg1]] + rt[neg1])
   ENDIF
   negs = WHERE(C3 LT 0d0, nneg)
   C1 = ABS(C3)^(1d0/3d0)
   IF nneg GT 0 THEN C1[negs] *= -1
   FOR ii = 0,2 DO BEGIN
      CC = C1*cru[ii]
      IF semiscalar THEN $
         roots[ii,gen] = - (b + CC + del0[gen]/CC)/(3d0*a[gen]) ELSE $
         roots[ii,gen] = - (b[gen] + CC + del0[gen]/CC)/(3d0*a[gen])
   ENDFOR
ENDIF ELSE roots = DBLARR(3,na)

IF ngen EQ na THEN GOTO, do_check
; Remaining cases are 3 real roots

; check for degenerate case with d = 0 (implies x*quadratic = 0)
IF nd GT 1 THEN BEGIN ; scalar & semiscalar cases already covered
   quad = WHERE(d[special] EQ 0, nquad, complement=true_cubic)
   quad = special[quad]
   true_cubic = special[true_cubic]
   IF nquad GT 0 THEN BEGIN
      roots[0,quad] = 0d0
      rtdelta = SQRT(b[quad]^2 - 4*a[quad]*c[quad])
      roots[1,quad] = (-b[quad]+rtdelta) / (2*a[quad])
      roots[2,quad] = (-b[quad]-rtdelta) / (2*a[quad])
   ENDIF
   IF nquad+ngen EQ na THEN GOTO, do_check
ENDIF ELSE true_cubic = special

; use trig solution for non-zero discriminant
trig = WHERE(discriminant[true_cubic] GT 0d0, ntrig)
IF ntrig GT 0 THEN BEGIN 
   trig = true_cubic[trig]
   p  = - del0[trig] / (3d0*a2[trig])
   q  = del1[trig] / (27d0*a[trig]^3)
   ff = 2*SQRT(-p/3d0)
   acos1 = ACOS(3d0*q/(p*ff))
   tk = DBLARR(3,ntrig)
   FOR j=0,ntrig - 1 DO BEGIN ; Check!!
;      i = trig[j]
      tk[*,j] = ff[j]*COS((acos1[j]-[0, 2*!dpi, 4*!dpi])/3d0)
   ENDFOR
   bb = semiscalar ? b : b[trig]
   tk2x = bb / (3d0*a[trig])
   roots[0,trig] = tk[0,*] - tk2x
   roots[1,trig] = tk[1,*] - tk2x
   roots[2,trig] = tk[2,*] - tk2x
ENDIF 

; Special cases:
spec = WHERE(discriminant[true_cubic] EQ 0d0, nspec)
IF nspec GT 0 THEN BEGIN
   spec = true_cubic[spec]
   one = WHERE(del0[spec] EQ 0d0, n1, COMPLEMENT = two)
   s1 = spec[one]
   IF n1 GT 0 THEN BEGIN
      bb = semiscalar ? b : b[s1]
      x1 = -bb / (3d0*a[s1])
      roots[0,s1] = x1
      roots[1,s1] = x1
      roots[2,s1] = x1
   ENDIF
   IF n1 LT nspec THEN BEGIN
      s2 = spec[two]
      IF semiscalar THEN BEGIN
         x1 = (9*a[s2]*d - b*c) / (2d0 * del0[s2])
         x2 = (4*a[s2]*b*c - 9*a2[s2]*d - b^3)/(a[s2]*del0[s2])
      ENDIF ELSE BEGIN
         x1 = (9*a[s2]*d[s2] - b[s2]*c[s2]) / (2d0 * del0[s2])
         x2 = (4*a[s2]*b[s2]*c[s2] - 9*a2[s2]*d[s2] - b[s2]^3)/(a[s2]*del0[s2])
      ENDELSE
      roots[0,s2] = x2
      roots[1,s2] = x1
      roots[2,s2] = x1
   ENDIF
ENDIF

do_check:
; check alleged roots really are
test = DCOMPLEXARR(3,na)
grad = DCOMPLEXARR(3,na)
frac_err = DCOMPLEXARR(3,na)
FOR i=0,2 DO BEGIN
   test[i,*] = ((a*roots[i,*] + b)*roots[i,*] + c)*roots[i,*] + d
   grad[i,*] = (3*a*roots[i,*] + 2*b)*roots[i,*] + c
ENDFOR
non_zero = WHERE(test NE COMPLEX(0d0), nnz)
frac_err[non_zero] = test[non_zero] / (grad[non_zero] * roots[non_zero])
IF MAX(frac_err[0,*]) GT 1d-2 THEN STOP
; This REALLY shouldn't ever happen!

nbc = 0
IF ngen EQ 0 THEN badr = WHERE(ABS(frac_err) GT 1e-3, nbad) ELSE BEGIN
   badr = WHERE(ABS(frac_err[0,*]) GT 1e-3, nbad)
   badc = WHERE(ABS(frac_err[1:2,gen]) GT 1e-3,nbc)
   IF nbc GT 0 THEN check = 1B
ENDELSE
nbad2 = 0
IF ngen LT na THEN badr2 = WHERE(ABS(frac_err[1:2,special]) GT 1e-3, nbad2)

IF nbad+nbad2 GT 0 THEN check = 2B

; TODO: finish working through this case. What's badr2 for ?
IF check GE 2B THEN BEGIN
   IF nbad LT 10 THEN BEGIN
      badr = gen[badr]
      badr2 = 0
      FOR jj=0,nbad-1 DO BEGIN
         kk = badr[jj]
         PRINT, 'Discriminant:', discriminant[kk/3], $
                ', cubic coefficients:', a[kk/3], b[kk/3], c[kk/3], d[kk/3], $
                ', cubic fractional error:', frac_err[0,kk/3]
      ENDFOR
   ENDIF ELSE PRINT, 'CUBIC:', nbad, ' Bad real solutions and', nbc, $
                     ' Bad complex solutions'
ENDIF

IF na EQ 1 THEN roots = REFORM(roots,3)
RETURN, roots

END
