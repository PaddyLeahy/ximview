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
PRO cursor_grip
;
; Sets the device cursor to a gripping hand
;
COMPILE_OPT IDL2, HIDDEN

curdat = [32769, 28686, 18482, 18514, 20114,  2448, 2432,  2433,  $
            384,   640,  1088,  2112,  2112,  4128,  4128, 61503]
hotspot = [8, 8]
mask   = [32769, 61455, 63551, 63615, 65279, 65535, 65535, 65535, $
          65535, 65279, 64639, 63615, 63615, 61503, 61503, 61503]

DEVICE, CURSOR_IMAGE = curdat, CURSOR_XY = hotspot, $
	CURSOR_MASK = mask

END
