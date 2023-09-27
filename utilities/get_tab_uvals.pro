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
FUNCTION get_tab_uvals, parent, tabID
;+
; NAME:
;       GET_TAB_UVALS
;
; PURPOSE:
;       Returns list of UVALUES for all the widgets that are (direct)
;       children of the widget with ID parent. Also returns the
;       widget-IDs via the tabID argument. Main application is where
;       tab names are stored in the UVALUES.
;
; CATEGORY:
;       Widgets
;
; CALLING SEQUENCE:
;
;       Result = GET_TAB_UVALS(Parent, TabID)
;
; INPUTS:
;       Parent:  Widget ID of parent (usually a tab)
;
; OUTPUTS:
;       Result is an array of UVALUES
;
; OPTIONAL OUTPUTS:
;       TabID:   An array containing the widget IDs of the children.
;
; RESTRICTIONS:
;       The UVALUES must all be the same type!
;
; EXAMPLE:
;
;        label = ['Tom', 'Dick', 'Harry']
;        ntabs = N_ELEMENTS(tablab)
;        tabs = WIDGET_TAB(tlb)
;        FOR i=0,ntab-1 DO $
;          junk = WIDGET_BASE(tabs, TITLE=label[i], UVALUE=label[i])
;
;   ...later...
;        PRINT, GET_TAB_UVALS(tabs)
;
;   IDL prints: Tom Dick Harry
; 
; MODIFICATION HISTORY:
;       Written by:     J. P. Leahy, February 2008
;-
COMPILE_OPT IDL2, HIDDEN

ntab = WIDGET_INFO(parent, /TAB_NUMBER)
tab0 = WIDGET_INFO(parent, /CHILD)

WIDGET_CONTROL, tab0, GET_UVALUE = tabname
tabid = [tab0]
tabs = [tabname]

FOR itab = 1,ntab-1 DO BEGIN
    tab0 = WIDGET_INFO(tab0, /SIBLING)
    tabid = [tabid, tab0]
    WIDGET_CONTROL, tab0, GET_UVALUE = tabname
    tabs = [tabs,tabname]
ENDFOR

RETURN, tabs
END
