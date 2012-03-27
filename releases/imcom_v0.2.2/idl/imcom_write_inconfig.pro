PRO imcom_write_inconfig, INCONFIGFILE=inconfigfile, $
                          USERXY=userxy,             $
                          PSFFILE=psffile,           $
                          GIMFILE=gimfile,           $
                          ROTANGDEG=rotangdeg,       $
                          NOISE=noise,               $
                          GIMXFILE=gimxfile,         $
                          GIMYFILE=gimyfile,         $
                          GIMXSCALE=gimxscale,       $
                          GIMYSCALE=gimyscale,       $
                          DITHER=dither

;
;    IMCOM_WRITE_INCONFIG.PRO - writes all input image config files
;                               for the prototype IMCOM package
;    Copyright (C) 2011 Barnaby Rowe
;
;    This program is free software: you can redistribute it and/or modify
;    it under the terms of the GNU General Public License as published by
;    the Free Software Foundation, either version 3 of the License, or
;    (at your option) any later version.
;
;    This program is distributed in the hope that it will be useful,
;    but WITHOUT ANY WARRANTY; without even the implied warranty of
;    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
;    GNU General Public License for more details.
;
;    You should have received a copy of the GNU General Public License
;    along with this program.  If not, see <http://www.gnu.org/licenses/>.
;
;    DESCRIPTION:
;    Writes all input image configuration info (supplied by keywords)
;    to ASCII configuration files in the format in which IMCOM is 
;    guaranteed to function correctly.  For a description of these
;    keywords and the values they should take, please refer to the 
;    supplied user manual IMCOM_DOCUMENTATION.PDF
;
;    CALLED BY: IMCOM_WRITE_CONFIG.PRO
;

if not keyword_set(inconfigfile) then message, "ERROR: INCONFIGFILE keyword must be set"
message, "Writing configuration file to "+inconfigfile, /INFO
if not keyword_set(psffile) then message, "ERROR: PSFFILE keyword must be set"
if not keyword_set(gimfile) then message, "ERROR: GIMFILE keyword must be set"

if n_elements(rotangdeg) eq 0 and not keyword_set(rotangdeg) then message, "ERROR: ROTANGDEG keyword must be set"
if not keyword_set(noise) then message, "ERROR: NOISE keyword must be set"

if keyword_set(userxy) then begin
  if not keyword_set(gimxfile) then message, "ERROR: GIMXFILE keyword must be set if USERXY=1"
  if not keyword_set(gimyfile) then message, "ERROR: GIMYFILE keyword must be set if USERXY=1"
endif else begin
  if not keyword_set(gimxscale) then message, "ERROR: GIMXSCALE keyword must be set if USERXY=0"
  if not keyword_set(gimyscale) then message, "ERROR: GIMYSCALE keyword must be set if USERXY=0"
  if not keyword_set(dither) then message, "ERROR: DITHER keyword must be set if USERXY=0"
endelse
openw, lun, inconfigfile, /GET_LUN
printf, lun, 'PSFFILE     '+strtrim(psffile, 2)
printf, lun, 'GIMFILE     '+strtrim(gimfile, 2)
printf, lun, 'ROTANGDEG   '+strtrim(rotangdeg, 2)
printf, lun, 'NOISE       '+strtrim(noise, 2)
if keyword_set(userxy) then begin
  printf, lun, 'GIMXFILE    '+strtrim(gimxfile, 2)
  printf, lun, 'GIMYFILE    '+strtrim(gimyfile, 2)
endif else begin
  printf, lun, 'GIMXSCALE   '+strtrim(gimxscale, 2)
  printf, lun, 'GIMYSCALE   '+strtrim(gimyscale, 2)
  printf, lun, 'DITHER      '+strtrim(dither[0], 2)+' '+strtrim(dither[1], 2)
endelse
free_lun, lun
END
