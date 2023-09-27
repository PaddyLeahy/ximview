PRO grid_labels, lines, inc, lat_lon, sex, label, lablen
;
; Makes label strings for a set of coordinate lines
; set labels to be as compact as possible
;
; INPUTS
;   lines   Coordinate values (degrees for lat, lon, dec,; hours for
;           RA)
;   inc     increment between lines
;   lat_lon = 1 for longitude axis, = 2 for latitude, = 0 for neither.
;   sex     if true, interpret as sexagesimal coordinate
;
; OUTPUTS
;   label   label string for each line
;   lablen  approximate length of label in device coordinates
;
ON_ERROR, 2
ncom = 0 ; keep a tally of embedded control characters
IF sex THEN BEGIN
   hfmt = lat_lon EQ 1 ? "(I02,'!Uh!N')" : "(I+03,'!M%')"
   ncom += lat_lon ? 4 : 2
                                ; The tricky thing is to avoid
                                ; rounding errors
   ncom += lat_lon ? 4 : 2
   IF inc GE 1d0 THEN BEGIN
      hours = ROUND(lines)
      label = STRING(ROUND(lines),FORMAT=hfmt) 
   ENDIF ELSE BEGIN             ; Include minutes in label
      mm = lat_lon EQ 1 ? '!Um!N' : "!M'"
      ncom += lat_lon EQ 1 ? 4 : 2
      abslines = ABS(lines)
      abshours = FIX(abslines)
      sgn = 2* (lines GE 0) - 1
      
      minutes = 60*(abslines - abshours)
      IF inc GE 1d0/60d0 THEN BEGIN
         minutes = ROUND(minutes)
         rnd = WHERE(minutes EQ 60, nr)
         IF nr GT 0 THEN BEGIN
            minutes[rnd] = 0
            abshours[rnd] += 1
         ENDIF
                                ; floating zero has sign but integer does not.
         hours = sgn*abshours
         label = STRING(hours,FORMAT=hfmt) + $
                 STRING(minutes,FORMAT="(I02)") + mm
      ENDIF ELSE BEGIN          ; Include seconds in label
         ss = lat_lon EQ 1 ? '!Us!N' : '!M"'
         ncom += lat_lon EQ 1 ? 4 : 2
         iminutes = FIX(minutes)
         seconds = 60*(minutes - iminutes)
         IF inc GE 1d0/3600d0 THEN BEGIN
            seconds = ROUND(seconds)
            rnd = WHERE(seconds EQ 60, nr)
            IF nr GT 0 THEN BEGIN
               seconds[rnd] = 0
               iminutes[rnd] += 1
               rnd2 = WHERE(iminutes EQ 60, nr2)
               IF nr2 GT 0 THEN BEGIN
                  iminutes[rnd2] = 0
                  abshours[rnd2] += 1
               ENDIF
            ENDIF
            hours = sgn*abshours
            label = STRING(hours,FORMAT=hfmt) + $
                           STRING(iminutes,FORMAT="(I02)") + mm + $
                           STRING(seconds,FORMAT="(I02)") + ss
         ENDIF ELSE BEGIN       ; Include fractional seconds in label
            ndec = (-FLOOR(ALOG10(3600*inc)))
            ndigit = 3+ndec
            sfmt = STRCOMPRESS(STRING("(F",ndigit,'.',ndec,",A)"),/REMOVE_ALL)
            hours = sgn*abshours
            label = STRING(hours,FORMAT=hfmt) + $
                           STRING(iminutes,mm,FORMAT="(I02,A)") + $
                           STRING(seconds,ss,FORMAT=sfmt)
         ENDELSE
      ENDELSE
   ENDELSE
ENDIF ELSE BEGIN
   ndec = (-FLOOR(ALOG10(inc))) > 0
   sdec = '.' + STRTRIM(STRING(ndec),2)
   CASE lat_lon OF ; set format
      0: BEGIN 
         ndigit = CEIL(ALOG10(MAX(ABS(lines))))
         IF MIN(lines) LT 0 THEN ndigit += 1  ; allow for sign
         IF ndec GT 0 THEN ndigit += 1 + ndec ; allow for decimal point
         sdigit = STRTRIM(STRING(ndigit),2)
         fmt = ndec GT 0 ? "(F"+sdigit+sdec : "(I" + sdigit
      END
      1: BEGIN
         ndigit = ndec GT 0 ? 4 + ndec : 3
         sdigit = STRTRIM(STRING(ndigit),2)
         fmt = "('!8l!3 = '," 
         fmt += ndec GT 0 ? "F0" + sdigit + sdec : "I0" + sdigit
         ncom += 4
      END
      2: BEGIN
         ndigit = ndec GT 0 ? 4 + ndec : 3
         sdigit = STRTRIM(STRING(ndigit),2)
         fmt = "('!8b!3 = '," 
         fmt += ndec GT 0 ? "F+0" +sdigit + sdec : "I+0" + sdigit
         ncom += 4
      END
   ENDCASE
   fmt += ")"               
   label = STRING(lines, FORMAT=fmt)
ENDELSE

lablen = (STRLEN(label[0]) - ncom)*!D.X_CH_SIZE ; correct for embedded commands

END
