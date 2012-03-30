PRO make_testxy, npix, x, y, DITHER=dither, $
                             ROTANG=rotang, $
                             XSCALE=xscale, $
                             YSCALE=yscale

if not keyword_set(xscale) then xscale = 1.d0
if not keyword_set(yscale) then yscale = 1.d0
if not keyword_set(dither) then dither = [0.d0, 0.d0]
if keyword_set(rotang) then begin
; Always roll around centre of the image
  shapelets_make_dblxarr, npix, xu, yu
  rotate_xy, xu, yu, rotang, x, y
  x = x + 0.5d0 * double(npix[0]) - 0.5d0
  y = y + 0.5d0 * double(npix[1]) - 0.5d0
endif else begin
  shapelets_make_dblxarr, npix, x, y, X0=[0.5d0, 0.5d0]
endelse
x = x * double(xscale) + double(dither[0])
y = y * double(yscale) + double(dither[1])
END
