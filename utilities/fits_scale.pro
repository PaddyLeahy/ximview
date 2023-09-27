PRO fits_scale, data, header
;
; Apply FITS header scaling parameters, if present, and convert blanks
; to NaNs
;
bitpix = SXPAR(header,'BITPIX')
double = ABS(bitpix) EQ 64
n_bad = 0
IF bitpix GT 0 THEN BEGIN
   blank = SXPAR(header,'BLANK', COUNT = N_blank)
   IF n_blank THEN BEGIN
      SXDELPAR, header, 'BLANK'
      bad = WHERE(data EQ blank, n_bad)
      NaN = double ? !values.D_NAN : !values.F_NAN
   ENDIF 
ENDIF 

Bscale = sxpar( header, 'BSCALE' , Count = N_bscale)
IF n_bscale THEN SXDELPAR, header, 'BSCALE' ELSE Bscale = 1d0 
Bzero = sxpar(header, 'BZERO', Count = N_Bzero )
IF n_Bzero  THEN SXDELPAR, header, 'BZERO' ELSE Bzero  = 0d0

IF n_bscale || n_Bzero THEN BEGIN
   data = bscale*data + bzero
   IF ~double THEN data = FLOAT(data)
ENDIF
IF n_bad GT 0 THEN data[bad] = NaN

END
