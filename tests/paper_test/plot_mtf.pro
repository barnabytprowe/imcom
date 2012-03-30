PRO plot_mtf, mtf, FILESTRING=filestring,       $
                   PSFSCALE=psfscale,           $
                   SIZE_ORIGINAL=size_original, $
                   PLOT=plot

if not keyword_set(psfscale) then psfscale = 1.
dims = size(mtf, /DIM)
if not keyword_set(size_original) then size_original = dims
maxu = float(dims) / float(size_original) / 2.0 / psfscale ; cycles arcsec^-1

if keyword_set(plot) then begin
  window, 0, XS=400, YS=440
  plt_image, float(mtf), /COL, FRAME=[-maxu[0], maxu[0] , -maxu[1], maxu[1]]
  oplot, [0,0], [-1000,1000]
  oplot, [-1000,1000], [0,0]
  window, 1, XS=400, YS=440
  plt_image, imaginary(mtf), /COL, FRAME=[-maxu[0], maxu[0] , -maxu[1], maxu[1]]
  oplot, [0,0], [-1000,1000]
  oplot, [-1000,1000], [0,0] 
  window, 2, XS=400, YS=440
  plt_image, abs(mtf), /COL, FRAME=[-maxu[0], maxu[0] , -maxu[1], maxu[1]]
  oplot, [0,0], [-1000,1000]
  oplot, [-1000,1000], [0,0]
endif

entry_device = !d.name
set_plot, 'PS'
device, FILENAME=filestring+".real.ps", XSIZE=4, YSIZE=4.7, /INCHES, /ENCAP
myplt_image, float(mtf), /COL, $
                         FRAME=[-maxu[0], maxu[0] , -maxu[1], maxu[1]], $
                         /SCALABLE, TITLE="Real(MTF)",          $
                         POS=[0.15,0.12,0.95,0.8],              $
                         CSIZE=0.03, /INVERSE,                  $
                         XTITLE="u!Dx!N [cycles arcsec!U-1!N]", $
                         YTITLE="u!Dy!N [cycles arcsec!U-1!N]" 
oplot, [0,0], [-1000,1000]              
oplot, [-1000,1000], [0,0]
device, /CLOSE_FILE

device, FILENAME=filestring+".imag.ps", XSIZE=4, YSIZE=4.7, /INCHES, /ENCAP
myplt_image, imaginary(mtf), /COL, $
                         FRAME=[-maxu[0], maxu[0] , -maxu[1], maxu[1]], $
                         /SCALABLE, TITLE="Imaginary(MTF)",     $
                         POS=[0.15,0.12,0.95,0.8],              $ 
                         CSIZE=0.03, /INVERSE,                  $
                         XTITLE="u!Dx!N [cycles arcsec!U-1!N]", $
                         YTITLE="u!Dy!N [cycles arcsec!U-1!N]" 
oplot, [0,0], [-1000,1000]              
oplot, [-1000,1000], [0,0]
device, /CLOSE_FILE

device, FILENAME=filestring+".magn.ps", XSIZE=4., YSIZE=4.7, /INCHES, /ENCAP
myplt_image, abs(mtf), /COL, $
                         FRAME=[-maxu[0], maxu[0] , -maxu[1], maxu[1]], $
                         /SCALABLE, TITLE="Magnitude(MTF)",     $
                         POS=[0.15,0.12,0.95,0.8],              $
                         CSIZE=0.03, /INVERSE,                  $
                         XTITLE="u!Dx!N [cycles arcsec!U-1!N]", $
                         YTITLE="u!Dy!N [cycles arcsec!U-1!N]" 
oplot, [0,0], [-1000,1000]              
oplot, [-1000,1000], [0,0]
device, /CLOSE_FILE

set_plot, entry_device

END


