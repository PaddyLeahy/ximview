FUNCTION string2coord, coord, hour_default, error
;
; Converts a string in many formats into a coordinate angle in degrees
;
; INPUT
;  coord         scalar string variable
;  hour_default  true (1) if hexadecimal input interpreted as hour by
;                default. If false (0), explicit hour notation is rejected
;                as an error.
;
; OUTPUT
;  error         0 if string converted correctly,
;                1 string specifies hours, but hours not allowed
;                2 uninterpretable string
;
; Possible strings and their interpretation:
;  Pure number     angle in degrees
;  xx:xx:xx        hexadecimal angle in deg/hr, min, sec.
;  xxhxxmxxs       hexadecimal angle in hours/min/sec
;  xxdxx'xx"       hexadecimal angle in degrees/min/sec
; Hexadecimal numbers can be truncated at any point; unspecified
; minutes are seconds are set to zero. m/' and s/" can be used interchangably
;
ON_ERROR, 2
CATCH, error_status

IF error_status NE 0 THEN BEGIN
   CATCH, /CANCEL
   HELP, /LAST_MESSAGE
   error = 2
   RETURN, 0d0
ENDIF

error = 0

hours = hour_default
number = STREGEX(coord,'^[+-]?[0-9.]+$',/BOOLEAN) 
IF number THEN angle = DOUBLE(coord) ELSE BEGIN
   minute = 0d0 & sec = 0d0
   IF STREGEX(coord,':',/BOOLEAN) THEN BEGIN
      ra = STRSPLIT(coord,':',/EXTRACT, COUNT = ct)
      IF ct GE 2 THEN minute =  ra[1]
      IF ct GE 3 THEN sec = ra[2]
   ENDIF ELSE BEGIN             ; must use hms or d'" as separators
      hours = STREGEX(coord,'h',/BOOLEAN)
      degrees = STREGEX(coord,'d',/BOOLEAN)
      IF (~hours && ~degrees) || (hours && degrees) THEN error = 2
      IF (hours && ~hour_default) THEN error = 1
      IF error NE 0 THEN RETURN, 0d0

      ra = STRSPLIT(coord,'hd',/EXTRACT, COUNT = ct)
      IF ct GT 1 THEN  BEGIN
         minute = STRSPLIT(ra[1],'m''',/EXTRACT, COUNT = ct)
         IF ct GT 1 THEN sec = STRSPLIT(minute[1],'s"',/EXTRACT)
      ENDIF
   ENDELSE
   angle = TEN(DOUBLE(ra[0]),DOUBLE(minute[0]),DOUBLE(sec[0]))
   IF hours THEN angle *= 15d0
ENDELSE

RETURN, angle

END
