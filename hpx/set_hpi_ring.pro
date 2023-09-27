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
PRO set_hpi_ring, nside, ix, iy, row, ringnum, ringpos, offset
;
; J. P. Leahy 2008
;
; Set up look-up table for rings:
;
; Inputs:
;   nside: Usual HEALPix parameter
;   ix, iy: arrays of pixel indices in facet (may be short ints)
;   row:   whether desired facet is in top, middle or bottom row
;
; Outputs: (all nside x nside arrays)
;   ringnum: index for the ring  (short int)
;   ringpos: position of pixel along ring for first facet in row
;            (short int except in middle row)
;   offset:  offset between adjacent facets (modulo 4*nside)
;
COMPILE_OPT IDL2, HIDDEN

nss = FIX(nside) ; Convert to short integer
nsmin1 = nss - 1S
nx = N_ELEMENTS(ix)  & ny = N_ELEMENTS(iy)
IF nx NE ny THEN MESSAGE, 'Inconsistent x and y pixel arrays'
xminy = ix - iy

CASE row OF
    0: BEGIN                    ; top facets:
        ringpos = nsmin1 - iy - (xminy > 0S)/2S
        ringnum = nsmin1 + TEMPORARY(xminy)
        offset  = (ringnum + 1S) < nss
    END
    1: BEGIN                    ; Middle facets: 
        ringpos = ((9*nside - 1) - xminy)/2L - iy
        ringpos MOD= 4*nside
        ringnum = (nss + nsmin1) + TEMPORARY(xminy)
        offset = nx GT 1 ? REPLICATE(nside,nx) : nss
    END

    2: BEGIN                    ; bottom facets:
        ringpos = nsmin1 - ix + (xminy < 0S)/2S
        ringnum = (3S*nss - 1S) + TEMPORARY(xminy)
        offset = ((4S*nss - 1S) - ringnum) < nss
    END
ENDCASE

END
