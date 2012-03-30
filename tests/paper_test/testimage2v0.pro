FUNCTION imfunc, x, y, x0, y0, tau
im = 50.d0 * exp(-sqrt((x - x0) * (x - x0) / 1.8d0 $
                     + (y - y0) * (y - y0)) / tau)
im = im + 0.1d0
return, im
END

;---

FUNCTION testimage2v0, npix, CENT=cent,         $
                             DITHER=dither,     $
                             ROTANG=rotang,     $
                             XSCALE=xscale,     $
                             YSCALE=yscale,     $
                             NSCALE=nscale,     $
                             FILENAME=filename, $
                             PLOT=plot,         $
                             NOISE=noise,       $
                             PSF=psf,           $
                             TAU=tau,           $
                             NOVER=nover,       $
                             NPOLY=npoly,       $
                             ERRIM=errim, PSHIRES=pshires

ON_ERROR, 2
if not keyword_set(xscale) then xscale = 1.d0
if not keyword_set(yscale) then yscale = 1.d0
if not keyword_set(dither) then dither = [0.d0, 0.d0]
if not keyword_set(rotang) then rotang = 0B
if not keyword_set(filename) then filename = "scratch.fits"
if not keyword_set(psf) then message, "Must set PSF keyword!"
if not keyword_set(nscale) then nscale = 10L ; factor by which PSF is oversampled relative to input images
if keyword_set(plot) then begin
  window, 0, xs=500, ys=550
endif
if not keyword_set(tau) then tau = 0.25d0
if not keyword_set(npoly) then npoly = 11L
if not keyword_set(nover) then nover = 5L

; note nscale must be xscale / psfxscale = yscale / psfyscale!
dnscale = nscale  ; emphasis the potential floaty-ness of this number!

npsf = size(psf, /DIM)
npad = long(double(npix) * dnscale + 1.d0) + 2L * npsf
nopad = npad * nover
;nim_over = long(double(npix) * dnscale + 1.d0) * nover
;  Note location of +1.d0 and rounding... this is to ensure that
;  nopad = nim_over + 2L * npsf * nover (exactly, always)
;

oxscale = xscale / dnscale / double(nover)
oyscale = yscale / dnscale / double(nover)
make_testxy, nopad, x_opad, y_opad, $
             XSCALE=oxscale,        $
             YSCALE=oyscale,        $
             ROTANG=rotang,         $
             DITHER=dither - (npsf * nover * [oxscale, oyscale])
; also make corresponding rotated *image* coordinates (for later interpolation)
if abs(rotang) gt (2.d0 * !dpi * 1.d-15) then begin
  rotate_xy, x_opad, y_opad, -rotang, xim_opad, yim_opad, /CENTRE
endif else begin
  xim_opad = x_opad
  yim_opad = y_opad
endelse
; Then make image
im_opad = imfunc(x_opad, y_opad, cent[0], cent[1], tau)

psf_pad = dblarr((npad))
psf_pad[0L:npsf[0] - 1L, 0L:npsf[1] - 1L] = psf
psf_pad = temporary(shift(psf_pad, -npsf / 2L))
fpsf_pad = temporary(fft(psf_pad, -1, /DOUBLE) * double(product(npad)))
fpsf_opad = dcomplexarr(nopad)
fpsf_opad[0L:npad[0] - 1L, 0L:npad[1] - 1L] = temporary(shift(fpsf_pad, npad / 2L))
fpsf_opad = temporary(shift(fpsf_opad, -npad / 2L))

fim_opad = temporary(fft(im_opad, -1, /DOUBLE))
imconv_opad = temporary(fft(fim_opad * fpsf_opad, +1, /DOUBLE))
;imconv_over = imconv_opad[0:nim_over[0] - 1L, $
;                          0:nim_over[1] - 1L]

make_testxy, npix, x, y, XSCALE=xscale, YSCALE=yscale, DITHER=dither, $
             ROTANG=rotang
; make corresponding sample locations in "image" coords - applying the
;    same transformations as used to create xim_over etc. above
if abs(rotang) gt (2.d0 * !dpi * 1.d-15) then begin
  rotate_xy, x, y, -rotang, xim, yim, /CENTRE
endif else begin
  xim = x
  yim = y
endelse
im = dblarr(npix)
;errim = im
for i=0L, npix[0] - 1L do begin

  for j=0L, npix[1] - 1L do begin

    im[i, j] = interp_image(imconv_opad, xim_opad, yim_opad, xim[i, j], $
                            yim[i, j], npoly, err)
;    errim[i, j] = err

  endfor

endfor
; add Gaussian Noise
if keyword_set(noise) then im = temporary(im + randomn(iseed, [npix[0], $
                                          npix[1]], /DOUBLE) * noise)
; Write
message, "Writing FITS test image to "+filename, /INFO
fits_write, filename, im
return, im
END
