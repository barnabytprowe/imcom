PRO charge_diff2v0, psf, psfcd, DIFFSCALE=diffscale, FILEOUT=fileout

if not keyword_set(diffscale) then message, "Set DIFFSCALE keyword (size of diffusion rms in units of PSF image pixels)"
s = size(psf, /DIM)
spad = 2L * s
impad = dblarr(spad)
impad[0L:s[0] -1L, 0L:s[1] -1L] = psf
shapelets_make_dblxarr, spad, x, y, X0=[0.5d0, 0.5d0]
x = shift(x, spad / 2L) - double(s[0])
y = shift(y, spad / 2L) - double(s[1])
gauss = exp(-(x * x + y * y) / 2.d0 / diffscale / diffscale) $
      / 2.d0 / !dpi / diffscale / diffscale
ftgauss = fft(gauss, -1, /DOUBLE) * product(double(spad))
outpad = abs(fft(ftgauss * fft(impad, -1, /DOUBLE), +1, /DOUBLE))
psfcd = outpad[0L:s[0]-1L, 0L:s[1]-1L] / total(outpad[0L:s[0]-1L, 0L:s[1]-1L])
if keyword_set(fileout) then fits_write, fileout, psfcd / total(psfcd)
END
