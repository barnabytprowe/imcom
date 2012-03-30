FUNCTION sinc, x
q = where(abs(x) gt 1.d-14, COMPLEMENT=notq)
out = dblarr(size(x, /DIM))
out[q] = sin(!dpi * x[q]) / x[q] / !dpi
out[notq] = 1.d0
return, out
END


PRO pix_response2v0, psf, psfpix, NSCALE=nscale, FILEOUT=fileout, NPAD=npad
if not keyword_set(nscale) then message, "Set NSCALE keyword!"
s = size(psf, /DIM)
if not keyword_set(npad) then npad = 2L
impad = dblarr(s * npad)
impad[0L:s[0] - 1L, 0L:s[1] - 1L] = psf
fftim = temporary(fft(impad, -1, /DOUBLE))

shapelets_make_dblxarr, s * npad, ux, uy, X0=[0.5d0, 0.5d0]
ux = temporary(ux - double((s[0] * npad) / 2L))
uy = temporary(uy - double((s[1] * npad) / 2L))
ux = temporary(shift(ux, (s[0] * npad + 1L) / 2L) / double(s[0] * npad))
uy = temporary(shift(uy, 0L, (s[1] * npad + 1L) / 2L) / double(s[1] * npad))
sincim = temporary(dcomplex(sinc(ux * double(nscale)) * $
                            sinc(uy * double(nscale))))
inv_ft = temporary(fft(fftim * sincim, +1, /DOUBLE))
im_conv = temporary(abs(inv_ft))
psfpix = temporary(im_conv[0:s[0]-1L, 0:s[1]-1L] / $
                   total(im_conv[0:s[0]-1L, 0:s[1]-1L]))
if keyword_set(fileout) then fits_write, fileout, psfpix / total(psfpix)
END
