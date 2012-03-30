PRO make_testimages6v2_rt5, DITHERS=dithers,       $
                            PSFFILE=psffile,       $
                            OUTFILE=outfile,       $
                            AFILE=afile,           $
                            BFILE=bfile,           $
                            QFILE=qfile,           $
                            LFILE=lfile,           $
                            PFILE=pfile,           $
                            TFILE=tfile,           $
                            KFILE=kfile,           $
                            SFILE=sfile,           $
                            UFILE=ufile,           $   
                            CONFIGFILE=configfile, $
                            FSTRING=fstring,       $
                            PLOT=plot

if not file_test("dirs.sav") then begin
  psfdir = "./psfs"+path_sep()
  configdir = "./configs"+path_sep()
  fitsdir = "./fits"+path_sep()
  save, psfdir, configdir, fitsdir, FILE="dirs.sav"
endif else restore, FILE="dirs.sav"
if not keyword_set(psffile) then begin
  psffile = psfdir+'PSF_WFIRST_1.3_1.0UM.sav'
  restore, FILE=psffile
endif else restore, FILE=psffile
psf = psf / total(psf)

jdemscale = 0.18d0

noise = 0.d0
nscale1 = jdemscale / psfscale
nscale2 = jdemscale / psfscale
nscale3 = jdemscale / psfscale
nscale4 = jdemscale / psfscale
nscale5 = jdemscale / psfscale

cyclesparsec = (1.3d0 / 1.d-6 / (180.d0 / !dpi) / 3600.d0)
; this is the number of cycles per arcsec at the band limit of the
; telescope = D / \Lambda [cycles per radian]
outxscale = 0.5d0 / cyclesparsec   ; for output = Nyquist sampling
outyscale = 0.5d0 / cyclesparsec   ; for output = Nyquist sampling
nscaleout = outxscale / psfscale

npix1 = 32L
npix2 = 32L 
npix3 = 32L 
npix4 = 32L
npix5 = 32L 
npixout = fix(npix1 * jdemscale / outxscale)

nover = 5L
npoly = 13L
if not keyword_set(fstring) then fstring = 'rt5xrt5'

gimxscale = double([nscale1, nscale2, nscale3, nscale4, nscale5]) * psfscale
gimyscale = double([nscale1, nscale2, nscale3, nscale4, nscale5]) * psfscale
xc = 0.5d0 * [outxscale, outyscale] * double(npixout)
print, outxscale, nscaleout, npixout, xc[0], xc[1]

if not keyword_set(dithers) then begin
  dithers = [[([0.d0, 2.d0, 4.d0, 1.d0, 3.d0] / 5.d0) * gimxscale[0]], $
             [([0.d0, 1.d0, 2.d0, 3.d0, 4.d0] / 5.d0) * gimyscale[0]]]
endif
fits_write, fitsdir+'psf.'+fstring+'.fits', psf

entry_device = !d.name
set_plot, 'PS'
device, FILE="./plots/xy_locations."+fstring+".ps", $
        /COLOR, XSIZE=4, YSIZE=4, /ENCAP, /INCHES
make_testxy, [npixout, npixout], xout, yout, XSCALE=outxscale, $
                                             YSCALE=outyscale, $
                                             DITHER=[0.d0, 0.d0]
nnout = npixout * npixout
plot, xout, yout, PSYM=3,                                 $
                  XRANGE=[-1, 15] * outxscale, $
                  YRANGE=[-1, 15] * outyscale, $
                  XTITLE="x [arcsec]",                    $
                  YTITLE="y [arcsec]",                    $
               ;   TITLE=fstring,                        $
                  SYMSIZE=1.5, /NODATA, $
                  XSTYLE=1, YSTYLE=1,   $
                  POSITION=[0.15, 0.12, 0.97, 0.97]
tvlct, 255, 0, 0, 1
tvlct, 0, 255, 0, 2
tvlct, 0, 0, 255, 3
tvlct, 255, 255, 0, 4
tvlct, 255, 0, 255, 5
make_testxy, [npix1, npix1], x1, y1, XSCALE=gimxscale[0], $
                                     YSCALE=gimyscale[0], $
                                     DITHER=reform(dithers[0, *])
make_testxy, [npix2, npix2], x2, y2, XSCALE=gimxscale[1], $
                                     YSCALE=gimyscale[1], $
                                     DITHER=reform(dithers[1, *])
make_testxy, [npix3, npix3], x3, y3, XSCALE=gimxscale[2], $
                                     YSCALE=gimyscale[2], $
                                     DITHER=reform(dithers[2, *])
make_testxy, [npix4, npix4], x4, y4, XSCALE=gimxscale[3], $
                                     YSCALE=gimyscale[3], $
                                     DITHER=reform(dithers[3, *])
make_testxy, [npix5, npix5], x5, y5, XSCALE=gimxscale[4], $
                                     YSCALE=gimyscale[4], $
                                     DITHER=reform(dithers[4, *])
oplot, x1, y1, PSYM=1, COL=1, SYMSIZE=2., THICK=5
oplot, x2, y2, PSYM=1, COL=2, SYMSIZE=2., THICK=5
oplot, x3, y3, PSYM=1, COL=3, SYMSIZE=2., THICK=5
oplot, x4, y4, PSYM=1, COL=4, SYMSIZE=2., THICK=5
oplot, x5, y5, PSYM=1, COL=5, SYMSIZE=2., THICK=5
oplot, xout, yout, PSYM=7, SYMSIZE=0.7, THICK=2
device, /CLOSE_FILE
set_plot, entry_device


if not keyword_set(outfile) then outfile = fitsdir+"H." $
                                                  +strtrim(fstring, 2)+".fits"
if not keyword_set(afile) then afile = fitsdir+"A."+strtrim(fstring, 2)+".fits"
if not keyword_set(bfile) then bfile = fitsdir+"B."+strtrim(fstring, 2)+".fits"
if not keyword_set(qfile) then qfile = fitsdir+"Q."+strtrim(fstring, 2)+".fits"
if not keyword_set(lfile) then lfile = fitsdir+"L."+strtrim(fstring, 2)+".fits"
if not keyword_set(pfile) then pfile = fitsdir+"P."+strtrim(fstring, 2)+".fits"
if not keyword_set(tfile) then tfile = fitsdir+"T."+strtrim(fstring, 2)+".fits"
if not keyword_set(kfile) then kfile = fitsdir+"K."+strtrim(fstring, 2)+".fits"
if not keyword_set(sfile) then sfile = fitsdir+"S."+strtrim(fstring, 2)+".fits"
if not keyword_set(ufile) then ufile = fitsdir+"U."+strtrim(fstring, 2)+".fits"
if not keyword_set(configfile) then configfile = configdir+strtrim(fstring, 2)+".config"

; do two test setups, one with psf different in each case, one with const

diff_frac = !dpi * 1.87d0 / 2.d0 / 18.d0  ; (from Roger Smith & Hardy et al 2007)
                            ; diffusion length as fraction of pixel size

charge_diff2v0, psf, psfcd1, DIFFSCALE=double(nscale1) * diff_frac,        $
                             FILE=fitsdir+"psf.cd."+fstring+".1.fits"
pix_response2v0, psf, psfpix1, NSCALE=nscale1,                             $
                               FILE=fitsdir+"psf.pix."+fstring+".1.fits"
charge_diff2v0, psfpix1, psfcdpix1, DIFFSCALE=double(nscale1) * diff_frac, $
                FILE=fitsdir+"psf.cd.pix."+fstring+".1.fits"
;
charge_diff2v0, psf, psfcd2, DIFFSCALE=double(nscale2) * diff_frac,        $
                             FILE=fitsdir+"psf.cd."+fstring+".2.fits"
pix_response2v0, psf, psfpix2, NSCALE=nscale2,                             $
                               FILE=fitsdir+"psf.pix."+fstring+".2.fits"
charge_diff2v0, psfpix2, psfcdpix2, DIFFSCALE=double(nscale2) * diff_frac, $
                FILE=fitsdir+"psf.cd.pix."+fstring+".2.fits"
;
charge_diff2v0, psf, psfcd3, DIFFSCALE=double(nscale3) * diff_frac,        $
                             FILE=fitsdir+"psf.cd."+fstring+".3.fits"
pix_response2v0, psf, psfpix3, NSCALE=nscale3,                             $
                               FILE=fitsdir+"psf.pix."+fstring+".3.fits"
charge_diff2v0, psfpix3, psfcdpix3, DIFFSCALE=double(nscale3) * diff_frac, $
                FILE=fitsdir+"psf.cd.pix."+fstring+".3.fits"
;
charge_diff2v0, psf, psfcd4, DIFFSCALE=double(nscale4) * diff_frac,        $
                             FILE=fitsdir+"psf.cd."+fstring+".4.fits"
pix_response2v0, psf, psfpix4, NSCALE=nscale4,                             $
                               FILE=fitsdir+"psf.pix."+fstring+".4.fits"
charge_diff2v0, psfpix4, psfcdpix4, DIFFSCALE=double(nscale4) * diff_frac, $
                FILE=fitsdir+"psf.cd.pix."+fstring+".4.fits"
;
charge_diff2v0, psf, psfcd5, DIFFSCALE=double(nscale5) * diff_frac,        $
                             FILE=fitsdir+"psf.cd."+fstring+".5.fits"
pix_response2v0, psf, psfpix5, NSCALE=nscale5,                             $
                               FILE=fitsdir+"psf.pix."+fstring+".5.fits"
charge_diff2v0, psfpix5, psfcdpix5, DIFFSCALE=double(nscale5) * diff_frac, $
                FILE=fitsdir+"psf.cd.pix."+fstring+".5.fits"
;

psize = size(psf, /DIM)
mtf = shift(fft(psf), psize[0] / 2L, psize[1] / 2L)
mtfcd1 = shift(fft(psfcd1), psize[0] / 2L, psize[1] / 2L)
mtfpix1 = shift(fft(psfpix1), psize[0] / 2L, psize[1] / 2L)
mtfcdpix1 = shift(fft(psfcdpix1), psize[0] / 2L, psize[1] / 2L)

mtfcd2 = shift(fft(psfcd2), psize[0] / 2L, psize[1] / 2L)
mtfpix2 = shift(fft(psfpix2), psize[0] / 2L, psize[1] / 2L)
mtfcdpix2 = shift(fft(psfcdpix2), psize[0] / 2L, psize[1] / 2L)

plot_mtf, mtf, FILE="./plots/jdem_"+fstring+"_mtf",             $
               PSFSCALE=psfscale, SIZE=psize
plot_mtf, mtfcd1, FILE="./plots/jdem_"+fstring+"_mtfcd1",       $
                  PSFSCALE=psfscale, SIZE=psize
plot_mtf, mtfpix1, FILE="./plots/jdem_"+fstring+"_mtfpix1",     $
                   PSFSCALE=psfscale, SIZE=psize
plot_mtf, mtfcdpix1, FILE="./plots/jdem_"+fstring+"_mtfcdpix1", $
                     PSFSCALE=psfscale, SIZE=psize
plot_mtf, mtfcd2, FILE="./plots/jdem_"+fstring+"_mtfcd2",       $
                  PSFSCALE=psfscale, SIZE=psize
plot_mtf, mtfpix2, FILE="./plots/jdem_"+fstring+"_mtfpix2",     $
                   PSFSCALE=psfscale, SIZE=psize
plot_mtf, mtfcdpix2, FILE="./plots/jdem_"+fstring+"_mtfcdpix2", $
                     PSFSCALE=psfscale, SIZE=psize
;
; Then make the images...
;
pwin, 0

im1 = testimage2v0([npix1, npix1],      $
                   CENT=xc,             $
                   XSCALE=gimxscale[0], $
                   YSCALE=gimyscale[0], $
                   PSF=psfcdpix1,          $
                   NSCALE=nscale1,      $
                   NOVER=nover,         $
                   NOISE=noise,         $
                   PLOT=plot,           $
                   TAU=tau,             $
                   DITHER=reform(dithers[0, *]), $
                   FILENAME=fitsdir+"test1."+strtrim(fstring, 2)+".fits")
plt_image, im1, /COL, /FRAME
;stop

im2 = testimage2v0([npix2, npix2],      $
                   CENT=xc,             $
                   XSCALE=gimxscale[1], $
                   YSCALE=gimyscale[1], $
                   PSF=psfcdpix2,          $
                   NSCALE=nscale2,      $
                   NOVER=nover,         $
                   NOISE=noise,         $
                   PLOT=plot,           $
                   TAU=tau,             $
                   DITHER=reform(dithers[1, *]), $
                   FILENAME=fitsdir+"test2."+strtrim(fstring, 2)+".fits")
plt_image, im2, /COL, /FRAME
;stop

im3 = testimage2v0([npix3, npix3],      $
                   CENT=xc,             $
                   XSCALE=gimxscale[2], $
                   YSCALE=gimyscale[2], $
                   PSF=psfcdpix3,          $
                   NSCALE=nscale3,      $
                   NOVER=nover,         $
                   NOISE=noise,         $
                   PLOT=plot,           $
                   TAU=tau,             $
                   DITHER=reform(dithers[2, *]), $
                   FILENAME=fitsdir+"test3."+strtrim(fstring, 2)+".fits")
plt_image, im3, /COL, /FRAME
;stop

im4 = testimage2v0([npix4, npix4],      $
                   CENT=xc,             $
                   XSCALE=gimxscale[3], $
                   YSCALE=gimyscale[3], $
                   PSF=psfcdpix4,          $
                   NSCALE=nscale4,      $
                   NOVER=nover,         $
                   NOISE=noise,         $
                   PLOT=plot,           $
                   TAU=tau,             $
                   DITHER=reform(dithers[3, *]), $
                   FILENAME=fitsdir+"test4."+strtrim(fstring, 2)+".fits")
plt_image, im4, /COL, /FRAME
;stop

im5 = testimage2v0([npix5, npix5],      $
                   CENT=xc,             $
                   XSCALE=gimxscale[4], $
                   YSCALE=gimyscale[4], $
                   PSF=psfcdpix5,          $
                   NSCALE=nscale5,      $
                   NOVER=nover,         $
                   NOISE=noise,         $
                   PLOT=plot,           $
                   TAU=tau,             $
                   DITHER=reform(dithers[4, *]), $
                   FILENAME=fitsdir+"test5."+strtrim(fstring, 2)+".fits")
plt_image, im5, /COL, /FRAME
;stop

imout = testimage2v0([npixout, npixout],  $
                     CENT=xc,             $
                     XSCALE=outxscale,    $
                     YSCALE=outyscale,    $
                     PSF=psfcdpix1,          $
                     NSCALE=nscaleout,    $
                     NOVER=nover,         $
                     NOISE=noise,         $
                     PLOT=plot,           $
                     TAU=tau,             $
                     DITHER=[0.d0, 0.d0], $
                     FILENAME=fitsdir+"testout."+strtrim(fstring, 2)+".fits")

plt_image, imout, /COL, /FRAME
;stop

; then write config files for our test cases

imcom_write_config3v0, CONFIGFILE=configfile+".psfvar",                   $
                       NEXP=5,                                            $
                       GIMFILE=fitsdir+'test'+['1', '2', '3', '4', '5']   $
                              +'.'+strtrim(fstring,2)+'.fits',            $
                       GIMXSCALE=gimxscale,                               $
                       GIMYSCALE=gimyscale,                               $
                       PSFFILE=fitsdir+'psf.cd.pix.'+fstring+'.'          $
                                      +['1', '2', '3', '4', '5']+'.fits', $
                       PSFXSCALE=psfscale,                                $
                       PSFYSCALE=psfscale,                                $
                       ROTANGDEG=replicate(0.d0, 5),                      $
                       NOISE=[1.d0, 1.d0, 1.d0, 1.d0, 1.d0], $ ;noise,    $
                       DITHERS=dithers,                                   $
                       GAMFILE=fitsdir+"psf.cd.pix."+fstring+".1.fits",   $
                       OUTFILE=outfile,                                   $
                       OUTXSCALE=outxscale,                               $
                       OUTYSCALE=outyscale,                               $
                       OUTPOS=[0.d0, 0.d0],                               $
                       NOUT=[npixout, npixout],                           $
                       AFILE=afile,                                       $
                       BFILE=bfile,                                       $
                       QFILE=qfile,                                       $
                       LFILE=lfile,                                       $
                       PFILE=pfile,                                       $
                       TFILE=fitsdir+"/T.test.psfvar.fits",               $
                       SFILE=fitsdir+"/S.test.psfvar.fits",               $
                       UFILE=fitsdir+"/U.test.psfvar.fits",               $
                       KFILE=fitsdir+"/K.test.psfvar.fits"

imcom_write_config3v0, CONFIGFILE=configfile+".psfconst",                 $
                       NEXP=5,                                            $
                       GIMFILE=fitsdir+'test'+['1', '2', '3', '4', '5']   $
                              +'.'+strtrim(fstring,2)+'.fits',            $
                       GIMXSCALE=gimxscale,                               $
                       GIMYSCALE=gimyscale,                               $
                       PSFFILE=fitsdir+'psf.cd.pix.'+fstring+'.'          $
                                      +['1', '1', '1', '1', '1']+'.fits', $
                       PSFXSCALE=psfscale,                                $
                       PSFYSCALE=psfscale,                                $
                       ROTANGDEG=replicate(0.d0, 5),                      $
                       NOISE=[1.d0, 1.d0, 1.d0, 1.d0, 1.d0], $ ;noise,    $
                       DITHERS=dithers,                                   $
                       GAMFILE=fitsdir+"psf.cd.pix."+fstring+".1.fits",   $
                       OUTFILE=outfile,                                   $
                       OUTXSCALE=outxscale,                               $
                       OUTYSCALE=outyscale,                               $
                       OUTPOS=[0.d0, 0.d0],                               $
                       NOUT=[npixout, npixout],                           $
                       AFILE=afile,                                       $
                       BFILE=bfile,                                       $
                       QFILE=qfile,                                       $
                       LFILE=lfile,                                       $
                       PFILE=pfile,                                       $
                       TFILE=fitsdir+"/T.test.psfconst.fits",             $
                       SFILE=fitsdir+"/S.test.psfconst.fits",             $
                       UFILE=fitsdir+"/U.test.psfconst.fits",             $
                       KFILE=fitsdir+"/K.test.psfconst.fits"

imcom_write_config3v0, CONFIGFILE=configfile+".sysNone",                  $
                       NEXP=5,                                            $
                       GIMFILE=fitsdir+'test'+['1', '2', '3', '4', '5']   $
                              +'.'+strtrim(fstring,2)+'.fits',            $
                       GIMXSCALE=gimxscale,                               $
                       GIMYSCALE=gimyscale,                               $
                       PSFFILE=fitsdir+'psf.cd.pix.'+fstring+'.'          $
                                      +['1', '1', '1', '1', '1']+'.fits', $
                       PSFXSCALE=psfscale,                                $
                       PSFYSCALE=psfscale,                                $
                       ROTANGDEG=replicate(0.d0, 5),                      $
                       NOISE=[1.d0, 1.d0, 1.d0, 1.d0, 1.d0], $ ;noise,    $
                       DITHERS=dithers,                                   $
                       GAMFILE=fitsdir+"psf.cd.pix."+fstring+".1.fits",   $
                       OUTFILE=outfile,                                   $
                       OUTXSCALE=outxscale,                               $
                       OUTYSCALE=outyscale,                               $
                       OUTPOS=[0.d0, 0.d0],                               $
                       NOUT=[npixout, npixout],                           $
                       AFILE="None",                                      $
                       BFILE="None",                                      $
                       QFILE="None",                                      $
                       LFILE="None",                                      $
                       PFILE="None",                                      $
                       TFILE=fitsdir+"/T.test.sysNone.fits",              $
                       SFILE=fitsdir+"/S.test.sysNone.fits",              $
                       UFILE=fitsdir+"/U.test.sysNone.fits",              $
                       KFILE=fitsdir+"/K.test.sysNone.fits"

END
