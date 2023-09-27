PRO testhp
;
; Routine to test hpx to hp and back pixel conversion
;
; Coordinates of facet BLC in units of nside for facets 0 through 12,
; where facet 12 is the repeat of facet 6 at the top right of grid.
fx = [1,0,3,2, 2,1,0,3, 2,1,4,3, 4]
fy = [2,1,4,3, 2,1,0,3, 1,0,3,2, 4]
order = 'RING'
nside = 512L
npf  = nside^2
npix = NSIDE2NPIX(nside)
ip = LINDGEN(npix)
ipr = REORDER(ip,/R2N) ; contains ring pixel index in nest order. 
FOR facet=0L,11 DO BEGIN
   i1 = facet*npf
   i2 = i1 + npf - 1L
   ipn =ip[i1:i2]
   hp2hpx, ipn, 'NESTED', nside, jface, ix, iy
   ipx = hpx2hp(nside, ix+nside*fx[jface], iy+nside*fy[jface], order, $ 
                R_OFF=tring) 
   IF ipx NE ipr[i1:i2] THEN $
      PRINT, jface, ix, iy, ip, ipx, FORMAT="(I2,2I5,2I8)"
ENDFOR
END
