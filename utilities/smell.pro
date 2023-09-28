; SmEll
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
;+
; NAME:
;
;         SmEll
;
; PURPOSE:
;
;         Finds the smallest-area ellipse enclosing a given set of
;         points in 2D
;
; CATEGORY:
;
;         Geometry
;
; CALLING SEQUENCE:
;
;         ellipse = SmEll(points, [SEED=seed, /PLOT])
;
; INPUTS:
;        
;        points:   n by 2 array of x,y coordinates for points
;
; OPTIONAL INPUTS:
;
;        None
;
; KEYWORD PARAMETERS:
;
;        seed  Starting seed for random number generator. If omitted
;              RANDOMU initialises itself using the machine time. Just
;              used for shuffling hull points at the start.
;
;        plot  If set, plots input points, hull points, and fitted circle
;
; OUTPUTS:
;
;        returns a structure containing the polynomial coefficients of
;        the ellipse in centred form
;            .sup    Number of support points (-1 if unknown)
;            .c      Coordinates of ellipse centre
;            .m      2x2 matrix encoding axis ratio and orientation
;            .z      size parameter
;
;         Some other fields used internally are present for sup = 4
;
;         If p = (x,y)^T is the column vector representation of a point on
;         the ellipse, then the equation of the ellipse is
;                    (p - c)^T m (p-c) - z = 0
;
;         In IDL the vectors are stored as rows so the LHS is represented as
;               (p-c) ## m ## TRANSPOSE(p-c) - z
;         or equivalently but slightly more efficiently:
;               (p-c) ## MATRIX_MULTIPLY(p-c,m, /ATRANSPOSE) - z
;
;         This quantity is positive for points outside the ellipse,
;         and negative for points strictly inside, so this is a fast
;         test for inclusion of any point within the ellipse.
;
; OPTIONAL OUTPUTS:
;
;        None
;
; COMMON BLOCKS:
;
;        SmEll_common  contains last-calculated ellipse structure
;
; SIDE EFFECTS:
;
;        If /PLOT is specified, plots input and fitted ellipse on the
;        current plot device. 
;        
; RESTRICTIONS:
;
;        None
;
; PROCEDURE:
;
;        First finds the convex hull of the input list using
;        QHULL. Then finds the minimum enclosing ellipse using
;        recursive algorithm of Welzl as described by 
;        Gaertner & Schoenherr (1997, 1998). This code does not use
;        the rational arithmetic of G&S 1998 and so is subject to
;        rounding errors in marginal cases. In the "difficult" case of
;        4 support points the G&S algorithm does not explicitly solve
;        for the ellipse, and the final minimisation of the ellipse area
;        is found by solving a cubic equation whose zeros correspond
;        to the extrema of det(m/z), where area = pi/SQRT(det(m/z))
;
; EXAMPLE:
;
;        The following example uses ELLIPSE_BOX and OPLOT_ELLIPSE which are
;        included in this file and so available after SmEll has run:
;
;        IDL> points = RANDOMN(seed,200,2)
;        IDL> PLOT, points[*,0], points[*,1], PSYM = 1, /ISOTROPIC
;        IDL> ellipse = SmEll(points)
;        IDL> box = ellipse_box(ellipse)
;        IDL> oplot_ellipse, ellipse, COLOR = 1
;        IDL> id = [INDGEN(4),0] ; repeat first corner to link up box
;        IDL> OPLOT, box[0,id], box[1,id], COLOR = 2
;
;        The same results can be obtained with
;
;        IDL> ellipse = SmEll(RANDOMN(seed,200,2), /PLOT)
;
; MODIFICATION HISTORY:
;        
;        Written by  J. P. Leahy, University of Manchester 
;        Version 0 25/08/2020
;        Version 1  5/09/2020
;        Bug fix 27/09/2023 - Removed innappropriate use of SIGNUM
;-
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Starts with some ancillary functions
;
FUNCTION centrep2param, centrep
; 
; converts from centred conic to parametric representation.
;
; INPUT
;    centrep   Standard ellipse structure giving centred form of
;              polynomial
;
; RETURNS
;
;     [x_centre, y_centre, semi-major axis, semi-minor axis, cos(tau), sin(tau)]
;
;     where tau is the angle of the major axis relative to the x axis.
;
; First extract general conic coefficients:
centre = centrep.c
z = centrep.z
r = centrep.M[0,0]
s = centrep.M[1,1]
t = centrep.M[1,0] ; M should be symmetric so also = M[0,1]
w = centre ## MATRIX_MULTIPLY(centre, centrep.m, /ATRANSPOSE) - z
w = w[0] ; convert to scalar
uv = - centrep.m # centre
u = uv[0]
v = uv[1]

; Find semi-major and semi-minor axes, (a,b):
detM = r*s - t^2
det0 = w*detM + 2*t*u*v - s*u^2 - r*v^2

; eigenvalues of M:
lambda = (r+s) + [1,-1]*SQRT((r+s)^2 - 4*detM)
lambda /= 2
; PRINT, 'Eigenvalues of M:', lambda

IF z NE 0d0 THEN BEGIN
   ab = SQRT(-det0/(detM*lambda))
   a = MAX(ab, MIN=b)
ENDIF ELSE BEGIN ; degenerate case, at least one axis has zero length
   a = 0.5d0*SQRT(r + s)
   b = 0d0
ENDELSE
; PRINT, 'Semi-major and minor axes:', a, b

; Find rotation angle tau:
IF z EQ 0d0 THEN BEGIN
   cos2tau = (s-r) / (4*a*a)
ENDIF ELSE IF b NE a THEN BEGIN
   factor = (a*b)^2/((b-a)*(a+b)) / z ; NB: negative since b < a & z >= 0
   cos2tau = (r-s)*factor
ENDIF ELSE BEGIN ; nominal tau = 0 for circle
   cos2tau = 1
ENDELSE
IF cos2tau GT 1d0 THEN MESSAGE, 'Wrong cos(2tau)'
; Choose tau in range +/-(pi/2):
costau = SQRT((1+cos2tau)/2)
;sintau = -SIGNUM(t)*SQRT((1-cos2tau)/2) ; WRONG: don't want to
;set to 0 for exact zero t
sintau = SQRT((1-cos2tau)/2)
IF t GE 0d0 THEN sintau *= -1

RETURN, [centre, a, b, costau, sintau]
END
;
PRO parcheck, centre, a, b, ct, st
;
; Unpacks parameters of parametric representation of ellipse, if needed
;
nec = N_ELEMENTS(centre)
CASE nec OF
   2: 
   6: BEGIN
      a = centre[2]
      b = centre[3]
      ct = centre[4]
      st = centre[5]
      centre = centre[0:1]
   END
   ELSE: MESSAGE, 'wrong number of parameters'
ENDCASE
END
;
FUNCTION param2points, ci, a, b, ct, st, RANDOM = random, SEED = seed
;
; returns point set of ellipse given its parametric representation.
; By default, points are every degree in paramter t, but may be a
; random list.
;
do_random = N_ELEMENTS(random) GT 0 

centre = ci
parcheck, centre, a, b, ct, st

Rt = [[ct, -st],[st,ct]]
ab = [[a, 0],[0,b]]
m1 = Rt ## ab

; curve parametric angle:
t = do_random ? RANDOMU(seed,random,/DOUBLE) : DINDGEN(361)/360d0
t *= 2d0*!dpi

vec = [[COS(t)],[SIN(t)]] # m1

FOR i=0,1 DO vec[*,i] += centre[i]

RETURN, vec
END
;
FUNCTION ellipse_box, ellipse
;
; Finds the range of x and y surrounding an ellipse
;
z = ellipse.z
r = ellipse.m[0,0]/z
s = ellipse.m[1,1]/z
t = ellipse.m[0,1]/z
x01 = ellipse.c[0] + [-1,1]/SQRT(r - t^2/s)
y01 = ellipse.c[1] + [-1,1]/SQRT(s - t^2/r)
coords = [[x01[0], y01[0]], [x01[0], y01[1]], $
              [x01[1], y01[1]], [x01[1], y01[0]]]
RETURN, coords
END
;
PRO oplot_ellipse, ellipse_in, COLOR = color, PSYM = PSYM, LINE=line
;
; Plots ellipse specified by structure containing centred version of
; conic polynomial
;
IF N_ELEMENTS(color) EQ 0 THEN color = 13
IF N_ELEMENTS(psym) EQ 0 THEN PSYM = 0
IF N_ELEMENTS(line) EQ 0 THEN LINE = 0
nin = N_ELEMENTS(ellipse_in)
ellipse = nin EQ 1 ? ellipse_in :  conic2centrep(ellipse_in) 
param = centrep2param(ellipse) 
xy = param2points(param)
OPLOT, xy[*,0], xy[*,1], COLOR = color, PSYM = psym, LINE = line
END
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; The following functions are used in the main code
;
FUNCTION shuffle, n, seed
;
; Produces a randomly shuffled list of the integers 0 to n-1
;
id_random = LONG(n*RANDOMU(seed,n))
id = LINDGEN(n)
FOR i=n-1,0,-1 DO BEGIN
   k = id[i]
   j = id_random[i]
   id[i] = id[j]
   id[j] = k
ENDFOR
RETURN, id
END
;
FUNCTION conic2centrep, conic
;
; Convert from conic polynomial to centred form
;
r = conic[0]
s = conic[1]
t = conic[2]
u = conic[3]
v = conic[4]
w = conic[5]
MM = [[r, t],[t, s]]
detM = r*s - t^2
IF detM EQ 0d0 THEN MESSAGE, 'zero determinant: cannot convert to centred form'
invM = [[s, -t],[-t, r]]/detM
uv = [u,v]
centre = - invM # uv
z = centre ## MATRIX_MULTIPLY(centre, MM, /ATRANSPOSE) - w
RETURN, {sup: -1, c: centre, m: MM, z: z[0]}
END
;
FUNCTION line, pt
;
; Returns coefficients of the line defined by two points p1, p2 in the form
; c1*x + c2*y + c3 = [p1,p2,p] = 0  where p is point (x,y)
;
RETURN, [pt[0,1]-pt[1,1],pt[1,0]-pt[0,0],pt[0,0]*pt[1,1] - pt[1,0]*pt[0,1]]
END
;
FUNCTION hyperbola, p12,p34
;
; Calculates the degenerate hyperbola through three points from the
; line coefficients of each pair
;
r = 2*p12[0]*p34[0]
s = 2*p12[1]*p34[1]
t = (p12[0]*p34[1] + p34[0]*p12[1])
u = (p12[0]*p34[2] + p12[2]*p34[0])
v = (p12[1]*p34[2] + p12[2]*p34[1])
w = 2*p12[2]*p34[2]

RETURN, [r,s,t,u,v,w]
         
END
;
FUNCTION orient, p1, p2, p3
;
; finds the relative orientation of 3 points, which should not be colinear
;
orient = SIGNUM(-p3[0]*p2[1] + p1[0]*p2[1] + p2[0]*p3[1] $
                -p1[0]*p3[1] + p3[0]*p1[1] - p2[0]*p1[1])
IF orient EQ 0 THEN MESSAGE, 'Colinear points' 
RETURN, orient
END
;
PRO swap_pt, pt, i, j
;
; Swaps two points in a list
;
temp = pt[i,*]
pt[i,*] = pt[j,*]
pt[j,*] = temp
END
;
PRO hypgen, pt_in, hyp1, hyp2
;
; Generate a pair of degenerate hyperbolae passing through 4 points
;
pt = pt_in
; Place points in counter-clockwise order, so get orientation of
; triplets:
o1_24 = orient(pt[1,*],pt[3,*],pt[0,*])
o3_24 = orient(pt[1,*],pt[3,*],pt[2,*])
IF o1_24 EQ o3_24 THEN IF orient(pt[2,*],pt[1,*],pt[0,*]) NE o3_24 THEN $
   swap_pt, pt, 0, 1 ELSE swap_pt, pt, 1, 2

p12 = line(pt[0:1,*])
p34 = line(pt[2:3,*])
hyp1 = hyperbola(p12,p34)

p23 = line(pt[1:2,*])
p41 = line([[pt[3,0],pt[0,0]],[pt[3,1],pt[0,1]]])
hyp2 = hyperbola(p23,p41)

END
;
FUNCTION conic, pt, cc
;
; Evaluates conic polynomial for a particular set of points
;
x = REFORM(pt[*,0])
y = REFORM(pt[*,1])
RETURN, cc[0]*x^2 + cc[1]*y^2 + 2*(cc[2]*x*y + cc[3]*x + cc[4]*y) + cc[5]
END
;
FUNCTION detA, tau, cc, z
;
; Evaluates determinant of A = M/z, equal to pi/(area)^2 for ellipse 
; Function is only needed and only works for four-point case. NB, det A
; is negative for hyperbola, so best not to work out area directly
;
COMMON SmEll_common, ellipse
nt = N_ELEMENTS(tau)
detA = DBLARR(nt)

hyp1 = ellipse.hyp1 
hyp2 = ellipse.hyp2 

lambda = ellipse.lambda0 - tau*hyp2[0]
mu     = ellipse.mu0     + tau*hyp1[0]

FOR it=0,nt-1 DO BEGIN
   cc = lambda[it]*hyp1 + mu[it]*hyp2
   ; Ensure cc is positive definite
   IF cc[0] LT 0 THEN cc *= -1
   detM = cc[0]*cc[1] - cc[2]^2

   Z1 = cc[3]*cc[1] - cc[4]*cc[2]
   Z2 = cc[4]*cc[0] - cc[3]*cc[2]
   ZZ = cc[3]*Z1 + cc[4]*Z2
   z = ZZ/detM - cc[5]
   detA[it] = detM/z^2
ENDFOR

RETURN, detA
END
;
FUNCTION delta, tau, CUBIC = cubic
;
; Evaluates parameter delta which is a cubic in tau proportional to
; the derivative of the ellipse area with respect to tau.
; Only needed for four support points.
; Optionally returns the cubic coefficients
;
COMMON SmEll_common, ellipse

cubic = KEYWORD_SET(cubic)
IF cubic THEN  tau = INDGEN(4)

nt = N_ELEMENTS(tau)
delta = DBLARR(nt)

hyp1 = ellipse.hyp1
hyp2 = ellipse.hyp2

lambda = ellipse.lambda0 - tau*hyp2[0]
mu     = ellipse.mu0     + tau*hyp1[0]

dc = hyp1[0]*hyp2 - hyp2[0]*hyp1 ; d(cc)/d(tau)

FOR it=0,nt-1 DO BEGIN
   cc = lambda[it]*hyp1 + mu[it]*hyp2

   detM = cc[0]*cc[1] - cc[2]^2
   dd = dc[0]*cc[1] + dc[1]*cc[0] -2*dc[2]*cc[2] ; d(detM)/d(tau)

   Z1 = cc[3]*cc[1] - cc[4]*cc[2]
   Z2 = cc[4]*cc[0] - cc[3]*cc[2]
   ZZ = cc[3]*Z1 + cc[4]*Z2

; derivatives wrt tau parameter of G&S97

   dZ1 = dc[3]*cc[1] + dc[1]*cc[3] - dc[4]*cc[2] - dc[2]*cc[4]
   dZ2 = dc[4]*cc[0] + dc[0]*cc[4] - dc[3]*cc[2] - dc[2]*cc[3]
   dZZ = dc[3]*Z1 + cc[3]*dZ1 + dc[4]*Z2 + cc[4]*dZ2 ; d(ZZ)/d(tau)

   delta[it] = 3*dd*ZZ + detM*(2*detM*dc[5] - dd*cc[5] - 2*dZZ)

ENDFOR

IF cubic THEN BEGIN
; If we were asked for the cubic coefficients, work them out 
; delta contains values for tau = 0, 1, 2, 3 so polynomial 
; delta = (a t^3 + b t^2 + c t + d)
;can be  written
;         { 1, 1, 1}{a}   {delta(1) - d}
;         { 8, 4, 2}{b} = {delta(2) - d}
;         {27, 9, 3}{c}   {delta(3) - d}
; Recover by inverting the matrix. Matrix inverse is integer if scaled
; by 6, so we do that:
   d = delta[0]
   abc = MATRIX_MULTIPLY( delta[1:3] - d, [[3,-3,1],[-15,12,-3],[18,-9,2]], $
                          /ATRANSPOSE)
   delta = [REFORM(abc), 6*d]            ; actually 6 delta
ENDIF

RETURN, delta
END
;
FUNCTION in_E, pt, ellipse_in
;
; Returns true (1B) if point (xp,yp) is in the ellipse, otherwise 0B. 
;
COMMON SmEll_common, ellipse
junk = ''


tol = 0d0 ; No need for finite tolerance

; Check inputs
IF N_ELEMENTS(pt) NE 2 THEN MESSAGE, 'One point at a time!'

IF N_ELEMENTS(ellipse_in) EQ 1 THEN ellipse = ellipse_in
IF N_ELEMENTS(ellipse) NE 1 THEN MESSAGE, 'No ellipse to test!'

CASE ellipse.sup OF
   0: in_ellipse = 0B           ; nothing is in a null ellipse!
   1: in_ellipse = ARRAY_EQUAL(ellipse.c,pt)
   4: BEGIN     
      ; Find coefficients of conic passing through test point
      hyp1 = ellipse.hyp1
      hyp2 = ellipse.hyp2
      mu     = -conic(pt,hyp1)
      lambda =  conic(pt,hyp2)
      cc  = lambda[0]*hyp1 + mu[0]*hyp2
; Check to see if conic is an ellipse
      detM = cc[0]* cc[1] - cc[2]^2
      IF detM GT 0 THEN BEGIN ; 5-point conic is ellipse
; we still need to know if the point is inside or outside the
; as-yet-undetermined minimum-area ellipse through the four support
; points.
         dc = hyp1[0]*hyp2 - hyp2[0]*hyp1
         ellipse.lambda0 = lambda
         ellipse.mu0 = mu
         rho = conic(pt, dc)
         in_ellipse = rho*delta(0d0) LE tol
      ENDIF ELSE BEGIN ; 5-point conic is parabola or hyperbola

; per G&S97, our point is in or out of minimum-area ellipse if it is in or
; out of any ellipse through the four points, for instance the one we
; already computed and stored in the ellipse variable

         ellipse.sup = 3 ; turn off testing for 4-point case
         in_ellipse = in_E(pt, et)
         ellipse.sup = 4 ; restore correct value
      ENDELSE

   END
   ELSE: BEGIN
      dp = REFORM(pt) - ellipse.c ; Eliminate dummy first dimension
      conic_p = dp ## MATRIX_MULTIPLY(dp, ellipse.m, /ATRANSPOSE) - ellipse.z
      in_ellipse = conic_p LE tol
      IF ellipse.sup EQ 2 THEN BEGIN 
                                ; Check point is within semi-major
                                ; axis of ellipse centre:
         d2 = ellipse.m[0,0] + ellipse.m[1,1]
         in_ellipse *= 4*TOTAL(dp^2) -d2 LE tol
      ENDIF 
   END
ENDCASE

RETURN, in_ellipse
END
;
FUNCTION e3_comp, pt
;
; Computes minimum-area ellipse passing through three points
;
centre = MEAN(pt, DIMENSION = 1)

dv = pt
FOR i=0,1 DO dv[*,i] -= centre[i]

invM = MATRIX_MULTIPLY(dv, dv, /ATRANSPOSE)
detM = (invM[0,0]*invM[1,1]-invM[0,1]*invM[1,0])

; Avoid division by scaling M and z by det(M):
MM = 3*[ [invM[1,1],-invM[1,0]], [-invM[0,1],invM[0,0]] ]
z = 2d0 * detM

RETURN, {sup: 3, c: centre, m: MM, z: z}
END
;
FUNCTION e4_comp, pt
;
; Pretends to compute minimum-area ellipse passing through four
; points, but actually just returns the corresponding degenerate
; hyperbolae and a sample ellipse, not necessarily minimal.
;
; get degenerate hyperbolae
hypgen, pt, hyp1, hyp2
r1 = hyp1[0]
s1 = hyp1[1]
t1 = hyp1[2]
alpha = r1*s1 - t1^2            ; det(m1)

r2 = hyp2[0]
s2 = hyp2[1]
t2 = hyp2[2]
gamma = r2*s2 - t2^2            ; det(m2)

beta = r1*s2 + r2*s1 - 2*t1*t2

; Some ellipse (not minumum area) through the four points
lambda = 2*gamma - beta
mu     = 2*alpha - beta

conic = lambda*hyp1 + mu*hyp2
; make sure ellipse is positive definite
IF conic[0] LT 0d0 THEN conic *= -1
cr = conic2centrep(conic)

RETURN, {sup: 4, c: cr.c, m: cr.m, z: cr.z, $
         lambda0: lambda, mu0: mu, hyp1: hyp1, hyp2: hyp2}
END
;
FUNCTION e5_comp, pt
;
; Calculates two degenerate hyperbolae through the first four points
; per Fig. 2 of GS97, and then the unique conic through those points
; and a fifth point. Usually this was already done by a previous call
; to e4_comp / in_E which found that the 5th point was outside the
; 4-point ellipse.
;
COMMON SmEll_common, ellipse

                                ; process first four points if not
                                ; done before:
IF N_TAGS(ellipse) NE 8 THEN ellipse = e4_comp(pt[0:3,*])
q = pt[4,*]

lambda =  conic(q,ellipse.hyp2)
mu     = -conic(q,ellipse.hyp1)
cc  = lambda[0]*ellipse.hyp1 + mu[0]*ellipse.hyp2
IF cc[0] LT 0 THEN cc *= -1     ; Ensure positive definite

ellipse = conic2centrep(cc)

ellipse.sup = 5

RETURN, ellipse
END
;
FUNCTION SmEll_core, qpt, rpt
;
; Recursive part of SmEll algorithm as given in G&S98
;
; Note: we don't check for degenerate cases where 3 or more
; points are in a line, since in this implementation we are only
; looking at points in the convex hull, which are never colinear.
; Thus the only degenerate ellipses are those with 0,1, or 2 support points.
;
COMMON SmEll_common, ellipse

size_q = SIZE(qpt)
nq = size_q[1]
size_r = SIZE(rpt)
nr = size_r[1]

IF nq EQ 0 THEN BEGIN ; Actual computation if set Q is empty
   CASE nr OF
      0: ellipse = {sup: 0, z: 0d0}
                                ; Single-point degenerate ellipse
      1: ellipse = {sup: 1, c: rpt, m: [[0d0,1d0],[1d0,0d0]], z: 0d0}
      2: BEGIN ; "ellipse" is line segment between points
         centre = MEAN(rpt, DIMENSION=1, /DOUBLE)
         dp = rpt[1,*] - rpt[0,*]
         t = -dp[0]*dp[1]
         m = [[dp[1]^2, t],[t,dp[0]^2]] ; trace = (2 * semi-major axis)^2
         ellipse = {sup: 2, c: centre, m: m, z: 0d0}
      END
      3: ellipse = e3_comp(rpt)
      4: ellipse = e4_comp(rpt)
      5: ellipse = e5_comp(rpt)
      ELSE: MESSAGE, 'Too many points (>5) in R'
   ENDCASE

   RETURN, ellipse
ENDIF 
ellipse = SmEll_core([], rpt)
IF nr EQ 5 THEN RETURN, ellipse

nq1 = nq - 1
FOR ip=0,nq1 DO BEGIN
   in = in_E(qpt[ip,*], ellipse) 
   IF ~in THEN BEGIN
      id1 = ip EQ 0 ? !null : [0:ip-1] ; list of points up to current one
      qpt2 = qpt[id1,*]
                                ; insert current point at end of set R
                                ; (this is where a list would be
                                ; better!) 
      rpt2 = rpt EQ !null ? qpt[ip,*] $
                           : [[rpt[*,0],qpt[ip,0]],[rpt[*,1],qpt[ip,1]]]
      ellipse = SmEll_core(qpt2, rpt2)
      
      ; Move current point to front of P
      id2 = ip EQ nq1 ? !null : [ip+1:nq1]
      idx = [ip,id1,id2]
      qpt = qpt[idx,*]
   ENDIF
ENDFOR
RETURN, ellipse
END
;
FUNCTION SmEll, pt_in, SEED = seed, PLOT = do_plot
;
; Main routine. See earlier dochead
;
COMMON SmEll_common, ellipse
junk = ''

; process inputs
debug = KEYWORD_SET(do_plot)

; Convert internally to double precision
IF N_ELEMENTS(pt_in) EQ 0 THEN BEGIN
   MESSAGE, /INFORMATIONAL, 'No points in input list'
   RETURN, {sup: 0, z: 0}
ENDIF
pts = DOUBLE(pt_in)

; Check inputs
size_pts = SIZE(pts)
IF size_pts[0] NE 2 || size_pts[2] NE 2 THEN $
   MESSAGE, 'Points array should be n x 2'
np = size_pts[1]

IF debug THEN BEGIN
   PLOT, pts[*,0], pts[*,1], /ISOTROPIC, PSYM = 1, /YNOZERO   
   t0 = SYSTIME(/SECONDS)
   t1 = t0
ENDIF

IF np GT 3 THEN BEGIN
  ; Get convex hull
   QHULL, pts[*,0], pts[*,1], tr
;  tr is actually a list of the corners of facets, but the facets are
;  in no particular order so we just take the first point in each
;  facet as our list of points on the convex hull
   idx = REFORM(tr[0,*])
   pth = pts[idx,*] 
   nh = N_ELEMENTS(idx)
   t1 = SYSTIME(/SECONDS)
ENDIF ELSE BEGIN ; 1 to 3 points (!)
   pth = pts
   nh = np
ENDELSE

; centre points to avoid loss of accuracy due to large coefficients
centre = MEAN(pts,DIMENSION = 1)
FOR i=0,1 DO pth[*,i] -= centre[i]

IF debug THEN BEGIN
   PRINT, nh, FORMAT="(I6,' points in the convex hull.')"
                                ; Sort hull points into anticlockwise
                                ; order for plotting
   theta = ATAN(pth[*,1],pth[*,0])
   pth = pth[SORT(theta), *]

   OPLOT, [pth[*,0],pth[0,0]]+centre[0], [pth[*,1],pth[0,1]]+centre[1], COL = 2
   PLOTS, pth[0,0]+centre[0], pth[0,1]+centre[1], COL = 1, PSYM = 1
   READ, PROMPT = 'continue:', junk
ENDIF

idx = SHUFFLE(nh, seed)
pth = pth[idx,*]

t2b = SYSTIME(/SECONDS)
ellipse = nh LE 2 ? SmEll_core([],pth) : SmEll_core(pth,[])

IF debug THEN BEGIN
   t2 = SYSTIME(/SECONDS)
   PRINT, (t1-t0)*1e3,  FORMAT="('  Time for QHULL:     ',F7.3,' msec')"
   PRINT, (t2-t2b)*1e3, FORMAT="('  Time for SmEll_core:',F7.3,' msec')"  
   PRINT, 'Number of support points for ellipse:', ellipse.sup
ENDIF

; check for implicit 4-point ellipse
IF ellipse.sup EQ 4 THEN BEGIN ; find minimum-area ellipse
                                
; Equivalently, find maximum of detA(tau), or zero of d(detA)/d(tau)
; (in fact quantity delta which just has same sign as gradient of detA).
; But delta is just a cubic in tau, as fourth-order term cancels. 

; Extract coefficients of delta(tau):
   cd = delta(0d0,/CUBIC) 

; Solve cubic equation
   roots = cubic(cd[0],cd[1],cd[2],cd[3], /CHECK)

   root_type = SIZE(roots,/TYPE)
   IF root_type EQ 6 OR root_type EQ 9 THEN BEGIN ; complex roots, 1st is real
      tau_min = REAL_PART(roots[0])
      detA_max = detA(tau_min,cc)
   ENDIF ELSE BEGIN             ; 3 real roots
      dets = detA(roots)
      detA_max = MAX(dets,imax)
      detA_max = detA(roots[imax],cc) ; get conic coefficients
   ENDELSE
   IF detA_max LT 0d0 THEN MESSAGE, 'Failed to find ellipse'

   ellipse = conic2centrep(cc)
   ellipse.sup = 4
ENDIF

; restore original offset
ellipse.c += centre
FOR i=0,1 DO pth[*,i] += centre[i]

IF debug THEN BEGIN
   READ, PROMPT = 'continue:', junk
   box = ellipse_box(ellipse)
   xr = [box[0,0],box[0,2]]
   yr = [box[1,0],box[1,1]]
   PLOT, pts[*,0], pts[*,1], /ISO, PSYM = 1, /YNOZERO
   OPLOT, pth[*,0], pth[*,1], PSYM = 1, COL = 2
   oplot_ellipse, ellipse, COL=1
   id = [indgen(4),0]
   OPLOT, box[0,id], box[1,id], COL = 3
ENDIF

RETURN, ellipse

END
