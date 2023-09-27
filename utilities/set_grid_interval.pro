PRO set_grid_interval, range, lat_lon, interval, start, nline, SEXAGESIMAL=sexa
;
; Tries to choose a sensible gap between grid lines
;
; INPUTS
;     range   min-max coordinates to cover
;     lat_lon = 1 for longitude, 2 for latitude, 0 for neither.
;    
; KEYWORD INPUTS
;    sexa      true if we want grid lines at sexagesimal intervals,
;
; OUTPUTS
;    interval  spacing between grid lines in degrees
;    start     lowest value of coordinate at which to draw line
;    nline     number of lines to draw.
;   
sex = KEYWORD_SET(sexa)

ntarget = 5 ; aim for about ntarget grid lines in each direction
dc = 0
dr = range[1] - range[0]
IF dr LE 0 THEN MESSAGE, 'Range should be finite and positive'

                                ; Deal with simple cases
CASE lat_lon OF 
   1: IF dr GE 360d0 THEN BEGIN
      start = 0d0
      interval = 60d0
      nline = 6
      RETURN
   ENDIF
   2: IF dr GE 180d0 THEN BEGIN
      start = -90d0
      interval = 30d0
      nline = 7
      RETURN
   ENDIF
   ELSE:
ENDCASE

IF sex && lat_lon EQ 1 THEN factor = 15d0 ELSE factor = 1d0
ifactor = 1/factor
lower = range[0] * ifactor
upper = range[1] * ifactor
dr *= ifactor

log_int0 = ALOG10(dr/ntarget)
IF sex && log_int0 LT 0 THEN BEGIN ; fraction of degree or hour
   scale0 = 1d0
   scale1 = 60d0  ; go to minutes
   IF dr*scale1 LT ntarget THEN BEGIN
      scale0 = scale1
      scale1 = 3600d0            ; go to seconds
   ENDIF
   int1 = 1d0/scale1
   IF dr*scale1 LT ntarget THEN BEGIN ; revert to decimal for subarcsec
      scale0 = scale1
      int1 = 10d0^FLOOR(ALOG10(dr*scale0/ntarget)) / scale0
   ENDIF 
   multiplier = [1,2,3,4,6,10,20,30,40,60]
ENDIF ELSE BEGIN
   int1  = 10^FLOOR(log_int0)
   scale0 = 0.1d0/int1
   multiplier = [1,2,5,10]
   IF lat_lon GT 0 && int1 GE 1 THEN multiplier = [1,2,3,4,6,10]
ENDELSE
; order of magnitude for grid interval:

; To align grid lines with round numbers, initially put the start at
; the first round number below the window. Will be refined later 
start = FLOOR(lower*scale0)/scale0

REPEAT BEGIN
   interval = int1*multiplier[dc]
   dc += 1
   nl1 = FIX((upper - start)/interval) + 1
   IF nl1 LE 0 THEN PRINT, range, start, int1, interval
   lines = start + INDGEN(nl1)*interval
   in = WHERE(lines GE lower,nline)
ENDREP UNTIL nline LE ntarget
start = lines[in[0]]

IF sex THEN BEGIN
   start *= factor
   interval *= factor
ENDIF

IF lat_lon EQ 1 THEN start >= -180d0 ELSE IF lat_lon EQ 2 THEN start >= -90d0

IF nline EQ 0 THEN MESSAGE, 'Failed to find suitable grid interval'
END
