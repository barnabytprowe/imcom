PRO imcom_write_config3v0, CONFIGFILE=configfile, $
                           NEXP=nexp,             $
                           USERXY=userxy,         $
                           GIMFILE=gimfile,       $
                           GIMXFILE=gimxfile,     $
                           GIMYFILE=gimyfile,     $
                           GIMXSCALE=gimxscale,   $
                           GIMYSCALE=gimyscale,   $
                           PSFFILE=psffile,       $
                           PSFXSCALE=psfxscale,   $
                           PSFYSCALE=psfyscale,   $
                           ROTANGDEG=rotangdeg,   $
                           NOISE=noise,           $
                           DITHERS=dithers,       $
                           OUTFILE=outfile,       $
                           AFILE=afile,           $
                           BFILE=bfile,           $
                           QFILE=qfile,           $
                           LFILE=lfile,           $
                           PFILE=pfile,           $
                           KFILE=kfile,           $
                           TFILE=tfile,           $
                           SFILE=sfile,           $
                           UFILE=ufile,           $
                           GAMFILE=gamfile,       $
                           OUTXFILE=outxfile,     $
                           OUTYFILE=outyfile,     $
                           OUTXSCALE=outxscale,   $
                           OUTYSCALE=outyscale,   $
                           OUTPOS=outpos,         $
                           NOUT=nout

ON_ERROR, 2

if not keyword_set(configfile) then configfile = "imcom.config"
message, "Writing configuration file to "+configfile, /INFO

if not keyword_set(nexp) then message, "ERROR: NEXP keyword must be set"
if not keyword_set(psffile) then message, "ERROR: PSFFILE keyword must be set to specify PSF image filenames"
if not keyword_set(gimfile) then message, "ERROR: GIMFILE keyword must be set to specify image filenames"
if not keyword_set(psfxscale) then message, "ERROR: PSFXSCALE keyword must be set to specify PSF image pixel scales"
if not keyword_set(psfyscale) then message, "ERROR: PSFYSCALE keyword must be set to specify PSF image pixel scales"
if not keyword_set(rotangdeg) then begin
  message, "WARNING: ROTANGDEG keyword should be set to specify PSF rotation angle, assuming zero rotation for all images", /INFO 
  rotangdeg = dblarr(nexp)
endif else if size(rotangdeg, /DIM) ne nexp then message, "ERROR: float array ROTANGDEG must have NEXP elements"
if not keyword_set(noise) then message, "ERROR: NOISE keyword must be set to specify RMS noise in input images"
if size(noise, /DIM) ne nexp then message, "ERROR: float array NOISE must have NEXP elements"

if keyword_set(userxy) then begin
  if not keyword_set(gimxfile) then message, "ERROR: GIMXFILE keyword must be set to specify image pixel scales"
  if not keyword_set(gimyfile) then message, "ERROR: GIMYFILE keyword must be set to specify image pixel scales"
  if size(gimxfile, /DIM) ne nexp then message, "ERROR: string array GIMXFILE must have NEXP elements"
  if size(gimyfile, /DIM) ne nexp then message, "ERROR: string array GIMYFILE must have NEXP elements"
endif else begin
  if not keyword_set(gimxscale) then message, "ERROR: GIMXSCALE keyword must be set to specify image pixel scales"
  if size(gimxscale, /DIM) ne nexp then message, "ERROR: float array GIMXSCALE must have NEXP elements"
  if not keyword_set(gimyscale) then message, "ERROR: GIMYSCALE keyword must be set to specify image pixel scales"
  if size(gimyscale, /DIM) ne nexp then message, "ERROR: float array GIMYSCALE must have NEXP elements"
  if not keyword_set(dithers) then begin
    dithers = fltarr(nexp, 2)
    message, "WARNING: DITHERS keyword not set, assuming no offset between input images!", /INFO
  endif
  dsize = size(dithers, /DIM)
  if dsize[0] ne nexp then message, "ERROR: float array DITHERS must have 2 x NEXP elements"
endelse

instring = '.input'+strtrim(lindgen(nexp)+1L, 2)
inconfigfile = configfile+instring
outstring = '.output'
outconfigfile = configfile+outstring
if not keyword_set(gamfile) then message, "ERROR: GAMFILE keyword must be set to specify Gamma image filename"
if not keyword_set(outfile) then message, "ERROR: OUTFILE keyword must be set"

if keyword_set(userxy) then begin
  if not keyword_set(outxfile) then message, "ERROR: OUTXFILE keyword must be set if USERXY=1"
  if not keyword_set(outyfile) then message, "ERROR: OUTYFILE keyword must be set if USERXY=1"
endif else begin
  if not keyword_set(outxscale) then message, "ERROR: OUTXSCALE keyword must be set if USERXY=0"
  if not keyword_set(outyscale) then message, "ERROR: OUTYSCALE keyword must be set if USERXY=0"
  if not keyword_set(outpos) then message, "ERROR: OUTPOS keyword must be set if USERXY=0"
  if n_elements(outpos) ne 2 then message, "ERROR: float array OUTPOS must have 2 elements"
  if not keyword_set(nout) then message, "ERROR: NOUT keyword must be set to specify output image dimensions [n1out, n2out] if USERXY=0"
  if n_elements(nout) ne 2 then message, "ERROR: NOUT keyword must have 2 elements to specify output image dimensions [n1out, n2out]"
endelse

for iexp=0L, nexp-1L do begin

  if keyword_set(userxy) then begin
    imcom_write_inconfig, INCONFIGFILE=inconfigfile[iexp], $
                          USERXY=userxy,                   $
                          PSFFILE=psffile[iexp],           $
                          GIMFILE=gimfile[iexp],           $
                          ROTANGDEG=rotangdeg[iexp],       $
                          NOISE=noise[iexp],               $
                          GIMXFILE=gimxfile[iexp],         $
                          GIMYFILE=gimyfile[iexp]
  endif else begin
    imcom_write_inconfig, INCONFIGFILE=inconfigfile[iexp], $
                          USERXY=userxy,                   $
                          PSFFILE=psffile[iexp],           $
                          GIMFILE=gimfile[iexp],           $
                          ROTANGDEG=rotangdeg[iexp],       $
                          NOISE=noise[iexp],               $
                          GIMXSCALE=gimxscale[iexp],       $
                          GIMYSCALE=gimyscale[iexp],       $
                          DITHER=reform(dithers[iexp, *])
  endelse

endfor

if keyword_set(userxy) then begin
  imcom_write_outconfig, OUTCONFIGFILE=outconfigfile, $
                         USERXY=userxy,               $
                         GAMFILE=gamfile,             $
                         OUTFILE=outfile,             $
                         SFILE=sfile,                 $
                         UFILE=ufile,                 $
                         KFILE=kfile,                 $
                         TFILE=tfile,                 $
                         OUTXFILE=outxfile,           $
                         OUTYFILE=outyfile
endif else begin
  imcom_write_outconfig, OUTCONFIGFILE=outconfigfile, $
                         USERXY=userxy,               $
                         GAMFILE=gamfile,             $
                         OUTFILE=outfile,             $
                         SFILE=sfile,                 $
                         UFILE=ufile,                 $
                         KFILE=kfile,                 $
                         TFILE=tfile,                 $
                         OUTXSCALE=outxscale,         $
                         OUTYSCALE=outyscale,         $
                         OUTPOS=outpos,               $
                         NOUT=nout
endelse
if keyword_set(userxy) then ustring = '1' else ustring = '0'

openw, lun, configfile, /GET_LUN
printf, lun, 'PSFXSCALE   '+strtrim(psfxscale, 2)
printf, lun, 'PSFYSCALE   '+strtrim(psfyscale, 2)
printf, lun, 'NEXP        '+strtrim(nexp, 2)
printf, lun, 'USERXY      '+ustring
for i=0L, nexp-1L do printf, lun, 'INCONFIG'+strtrim(i+1L, 2)+'   '+inconfigfile[i]
printf, lun, 'OUTCONFIG   '+outconfigfile
if keyword_set(afile) then printf, lun, 'AFILE       '+strtrim(afile, 2)
if keyword_set(bfile) then printf, lun, 'BFILE       '+strtrim(bfile, 2)
if keyword_set(qfile) then printf, lun, 'QFILE       '+strtrim(qfile, 2)
if keyword_set(lfile) then printf, lun, 'LFILE       '+strtrim(lfile, 2)
if keyword_set(pfile) then printf, lun, 'PFILE       '+strtrim(pfile, 2)
free_lun, lun

END
