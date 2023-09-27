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
