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
PRO item_capture, event
; Extracts item index from a query widget and deletes it.
;
WIDGET_CONTROL, event.ID, GET_UVALUE=container
WIDGET_CONTROL, container, SET_UVALUE=event.INDEX
WIDGET_CONTROL, event.TOP, /DESTROY

END

FUNCTION get_user_item, event, question, items
;+
; NAME:
;       GET_USER_ITEM
;
; PURPOSE:
;       Launches a modal dialog to get a piece of info from the user
;
; CATEGORY:
;      Widget Routines
;
; CALLING SEQUENCE:
;
;      Result = GET_USER_ITEM(Event, Question, Items)
;
; INPUTS:
;      Event:    Event structure. event.TOP is used for group leader 
;                and the uservalue slot of event.ID is borrowed to
;                store the datum. 
;      Question: Text for label
;      Items:    String array (list of items)
;
; OUTPUTS:
;      Result is the item selected from the list.
;
; EXAMPLE:
;
;      names = ['Tom', 'Dick', 'Harry']
;      who = GET_USER_ITEM(event, 'Who is your best friend?', names) 
;
; MODIFICATION HISTORY:
;       Written by:     J. P. Leahy, February 2008
;-
COMPILE_OPT IDL2, HIDDEN

WIDGET_CONTROL, event.ID, GET_UVALUE = save
query = WIDGET_BASE(GROUP_LEADER = event.TOP, TITLE='Select Item', $
                    /MODAL, /COLUMN)
label = WIDGET_LABEL(query, VALUE = question)
drops = WIDGET_DROPLIST(query, VALUE = items, $
                        EVENT_PRO = "item_capture", UVALUE = event.ID)
WIDGET_CONTROL, query, /REALIZE
XMANAGER, 'get_user_item', query
WIDGET_CONTROL, event.ID, GET_UVALUE=item
WIDGET_CONTROL, event.ID, SET_UVALUE = save

RETURN, item

END
