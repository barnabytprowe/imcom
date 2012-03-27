PRO imcom_write_outconfig, OUTCONFIGFILE=outconfigfile, $
                           USERXY=userxy,               $
                           GAMFILE=gamfile,             $
                           OUTFILE=outfile,             $
                           KFILE=kfile,                 $
                           TFILE=tfile,                 $
                           SFILE=sfile,                 $
                           UFILE=ufile,                 $
                           OUTXFILE=outxfile,           $
                           OUTYFILE=outyfile,           $
                           OUTXSCALE=outxscale,         $
                           OUTYSCALE=outyscale,         $
                           OUTPOS=outpos,               $
                           NOUT=nout

;
;    IMCOM_WRITE_OUTCONFIG.PRO - writes all output image config files
;                                for the prototype IMCOM package
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
;    Writes all output image configuration info (supplied by keywords)
;    to ASCII configuration files in the format in which IMCOM is 
;    guaranteed to function correctly.  For a description of these
;    keywords and the values they should take, please refer to the 
;    supplied user manual IMCOM_DOCUMENTATION.PDF
;
;    CALLED BY: IMCOM_WRITE_CONFIG.PRO
;

if not keyword_set(outconfigfile) then message, "ERROR: OUTCONFIGFILE keyword must be set"
message, "Writing configuration file to "+outconfigfile, /INFO
if not keyword_set(gamfile) then message, "ERROR: GAMFILE keyword must be set"
if not keyword_set(outfile) then message, "ERROR: OUTFILE keyword must be set"
if not keyword_set(sfile) then message, "ERROR: SFILE keyword must be set"
if not keyword_set(ufile) then message, "ERROR: UFILE keyword must be set"
if not keyword_set(kfile) then message, "ERROR: KFILE keyword must be set"
if not keyword_set(tfile) then message, "ERROR: TFILE keyword must be set"
if keyword_set(userxy) then begin
  if not keyword_set(outxfile) then message, "ERROR: OUTXFILE keyword must be set if USERXY=1"
  if not keyword_set(outyfile) then message, "ERROR: OUTYFILE keyword must be set if USERXY=1"
endif else begin
  if not keyword_set(outxscale) then message, "ERROR: OUTXSCALE keyword must be set if USERXY=0"
  if not keyword_set(outyscale) then message, "ERROR: OUTYSCALE keyword must be set if USERXY=0"
  if not keyword_set(outpos) then message, "ERROR: OUTPOS keyword must be set if USERXY=0"
endelse

openw, lun, outconfigfile, /GET_LUN
printf, lun, 'GAMFILE     '+strtrim(gamfile, 2)
printf, lun, 'OUTFILE     '+strtrim(outfile, 2)
printf, lun, 'KFILE       '+strtrim(kfile, 2)
printf, lun, 'TFILE       '+strtrim(tfile, 2)
printf, lun, 'SFILE       '+strtrim(sfile, 2)
printf, lun, 'UFILE       '+strtrim(ufile, 2)
if keyword_set(userxy) then begin
  printf, lun, 'OUTXFILE    '+strtrim(outxfile, 2)
  printf, lun, 'OUTYFILE    '+strtrim(outyfile, 2)
endif else begin
  printf, lun, 'OUTXSCALE   '+strtrim(outxscale, 2)
  printf, lun, 'OUTYSCALE   '+strtrim(outyscale, 2)
  printf, lun, 'OUTPOS      '+strtrim(outpos[0], 2)+' '+strtrim(outpos[1], 2)
  printf, lun, 'NOUT        '+strtrim(nout[0], 2)+' '+strtrim(nout[1], 2)
endelse
free_lun, lun

END
