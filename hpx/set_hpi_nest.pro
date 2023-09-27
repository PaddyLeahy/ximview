; -----------------------------------------------------------------------------
;
;  Copyright (C) 2007-2008   J. P. Leahy
;
;
;  This file is part of Ximview and of HEALPix
;
;  Ximview and HEALPix are free software; you can redistribute them and/or modify
;  them under the terms of the GNU General Public License as published by
;  the Free Software Foundation; either version 2 of the License, or
;  (at your option) any later version.
;
;  Ximview and HEALPix are distributed in the hope that they will be useful,
;  but WITHOUT ANY WARRANTY; without even the implied warranty of
;  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
;  GNU General Public License for more details.
;
;  You should have received a copy of the GNU General Public License
;  along with HEALPix; if not, write to the Free Software
;  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
;
;
; -----------------------------------------------------------------------------
FUNCTION set_hpi_nest, nside, x, iy
;
; J. P. Leahy 2008
;
; Sets up look-up table from facet pixel index to healpix index. 
; Bit manipulation code is from HEALPix routine ring2nest.pro 
;
; Arguments: nside: usual HEALPix parameter
;            x, iy: Arrays of x and y indices for pixels in facet 
;                   NB: May be short integers.
; Returns: 1-D long int index for square facet array.
;
COMPILE_OPT HIDDEN

;  arrays for (x,y) -> pixel number mapping
COMMON xy2pix, x2pix, y2pix
sz = SIZE(x2pix)
IF (sz[sz[0]+1] EQ 0) THEN init_xy2pix ; initiate x2pix and y2pix

ix = (nside-1) - x  
ix_low = ix MOD 128
ix_hi  = TEMPORARY(ix)/128
iy_low = iy MOD 128
iy_hi  = iy/128

ipf =  (x2pix(ix_hi)  + y2pix(iy_hi)) * 16384L + $
  (x2pix(ix_low) + y2pix(iy_low)) ; in {0, nside**2 - 1}

RETURN, ipf

END

