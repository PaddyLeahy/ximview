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
PRO hp2hpx, pixel, order, nside, jface, ix, iy
;
; Converts HEALPix pixel number to (x,y) coordinates and facet number
;
nface = nside*nside
npix = N_ELEMENTS(pixel)

CASE order OF
    'NESTED': BEGIN
        jface = pixel / nface
        iface = pixel MOD nface
        xyind = set_xy_nest(nside, iface)
    END
    'RING': BEGIN
        tring = set_ring_offset(nside)
        tring = [tring, NSIDE2NPIX(nside)]
        row   = INTARR(npix, /NOZERO)
        minpix = MIN(pixel, MAX = maxpix)
        minrow = (WHERE(minpix LT tring))[0] - 1
        maxrow = (WHERE(maxpix LT tring))[0] - 1
        FOR i=minrow,maxrow-1 DO BEGIN
            in = WHERE(pixel GE tring[i] AND pixel LT tring[i+1])
            IF in[0] NE -1 THEN row[in] = i
        ENDFOR
        in = 0
        row[WHERE(pixel GE tring[maxrow])] = maxrow

        nrow = tring[row+1]
        nrow -= tring[row]
        nrow /= 4
        irow = pixel
        irow -= tring[row]
        jface = irow
        jface /= nrow  ; Now contains longitude quadrant

        iq   = irow MOD nrow
        nrow = 0

        rows = WHERE(row GE 3*nside-1)
        IF rows[0] NE -1 THEN jface[rows] += 8

        mid = WHERE(row GE nside AND row LT 3*nside-1)
        IF mid[0] NE -1 THEN BEGIN
            midrow  = row[mid]
            iq   = iq[mid]
            test1 = iq LE (midrow - nside) / 2
            test2 = iq LT (3*nside - midrow) / 2
            iq = 0

            wanted = midrow GE nside AND midrow LT 2*nside
            rows = WHERE(wanted AND test1)
            IF rows[0] NE -1 THEN jface[mid[rows]] += 4

            rows = WHERE(~test1 AND ~test2)
            IF rows[0] NE -1 THEN BEGIN
                jface[mid[rows]] += 1
                jface[mid[rows]] MOD= 4
                jface[mid[rows]] += 4
            ENDIF

            wanted = midrow GE 2*nside AND midrow LT 3*nside-1
            rows = WHERE(wanted AND test2)
            IF rows[0] NE -1 THEN jface[mid[rows]] += 4
            rows = WHERE(wanted AND test1 AND ~test2)
            IF rows[0] NE -1 THEN jface[mid[rows]] += 8
        ENDIF
        rows = 0  &  mid = 0  &  quad = 0

;        PIX2VEC_RING, nside, pixel, vec
;        VEC2PIX_NEST, nside, vec, rpix
;        rface = rpix / nface
;        IF ~ARRAY_EQUAL(jface, rface) THEN PRINT, 'Oops'
; pixel, row, irow, jface, rface

        frow = jface / 4        ; facet row
        nss = FIX(nside)        ; Convert to short integer
        nsmin1 = nss - 1S
        void = HISTOGRAM(frow, MIN = 0, BINSIZE = 1, REVERSE_INDICES=ri)
        frow = 0
        fcol = jface MOD 4      ; facet sequence in row
        ix = INTARR(npix)
        iy = INTARR(npix)
        nbin = N_ELEMENTS(void)
        IF ri[0] NE ri[1] THEN BEGIN
            idx = ri[ri[0]:ri[1]-1]
            offset = row[idx]
            offset += 1S
            offset <= nss
            offset *= fcol[idx]
            ringpos = - TEMPORARY(offset)
            ringpos += irow[idx]
            xminy = row[idx]
            xminy -= nsmin1
            iy[idx] = nsmin1 - (xminy > 0S)/2S
            iy[idx] -= ringpos
            ix[idx] = xminy
            ix[idx] += iy[idx]
            ringpos = 0  &  xminy = 0
        ENDIF
        IF nbin eq 1 THEN GOTO, done
        IF ri[1] NE ri[2] THEN BEGIN
            idx = ri[ri[1]:ri[2]-1]
            offset = fcol[idx]
            offset *= nside
            ringpos = - TEMPORARY(offset)
            ringpos += irow[idx]
            negs = WHERE(ringpos GT nside)
            IF negs[0] NE -1 THEN ringpos[negs] -= 4*nside
            negs = 0
            xminy = row[idx]
            xminy -= (nss + nsmin1)
            iy[idx] = ((9L*nside - 1L) - xminy)/2L - ringpos
            iy[idx] MOD= nside
            negs = WHERE(iy[idx] LT 0)
            IF negs[0] NE -1 THEN iy[idx[negs]] += nside
            ix[idx] = xminy
            ix[idx] += iy[idx]
            ringpos = 0  &  xminy = 0
        ENDIF
        IF nbin eq 2 THEN GOTO, done
        IF ri[2] NE ri[3] THEN BEGIN
            idx = ri[ri[2]:ri[3]-1]
            offset  = - row[idx]
            offset += (4*nss - 1)
            offset <= nss
            offset *= fcol[idx]
            ringpos = - TEMPORARY(offset)
            ringpos += irow[idx]
            xminy = row[idx]
            xminy -= (3S*nss - 1S)

            ix[idx] = nsmin1 + (xminy < 0S)/ 2S - ringpos
            iy[idx] = ix[idx]
            iy[idx] -= xminy
            ringpos = 0  &  xminy = 0
        ENDIF
done:
        idx = 0
    END
ENDCASE

END
