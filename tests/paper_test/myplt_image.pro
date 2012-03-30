pro myplt_image, a, FRAME=frame,COLBAR=colbar,CRAN=cran, $
  TITLE=title,XTITLE=xtitle,YTITLE=ytitle,CTITLE=ctitle, $
  INVERSE=inverse,SCALABLE=scalable,NOERASE=noerase,     $
  CSIZE=csize, CHARSIZE=charsize, POSITION=position,     $
  FSUP=fsup, XGRID=xgrid, YGRID=ygrid, TICKLEN=ticklen, $
  TICKINTERVAL=tickinterval


; Hacked by Barney Rowe to allow position keyword!
; August 2000 - modified by A.R. to allow plotting over the full
;   plotting region
; September 3, 1997 - modified by A.R. to allow active change for
;   the color bar range
; October 14, 1997 - modified by A.R. to allow plotting of a color
;   bar and flexible frame labels
; July 7, 1997 - written by A. Refregier
;
; PURPOSE: plot an image scaled to the current coordinate frame.
; Further plotting (such as contours) can be performed over the
; resulting plot. Optionally, an annotated color bar can be drawn at
; the top. The 'scalable' keyword must be invoked when outputing to a
; postscript file.
; NOTE: !p.region=0 can be used to reset the window after using this
; subroutine
; INPUT: a: image array
; OPTIONAL INPUT: scalable: use scalable pixels which is to be
;                   ps devices (default: nonscalable to be used
;                   with x-term device)
;                 frame: draw a frame with pixel index limits;
;                   to change limits set frame=[x0,x1,y0,y1].
;                   set frame=-1 to supress the frame while keeping
;                   the image within it.
;                 noerase: supress erasing
;                 colbar: draw color bar
;                 cran: range of a for colors (default: [min(a),max(a)])
;                 inverse: invert the color coding
;                 title: title for the plot
;                 x,ytitle: titles for the frame
;                 ctitle: title for the color bar
;		  csize: vertical size of the color bar (0-1, default:.12)
; OUTPUT: plot of the scaled image with, optionally, a coordinate
; frame and a color bar.

; set range of a for colors
if not keyword_set(cran) then cran=[min(a),max(a)]

; draw color bar
if keyword_set(colbar) then begin
  myplt_colbar,cran,scalable=scalable,inverse=inverse,title=ctitle,$
    csize=csize,POSITION=position
  noerase=1
endif

; set up the coordinate frame
nx=n_elements(a(*,0))
ny=n_elements(a(0,*))
plot,[0],[0],/nodata,xstyle=12,ystyle=12,noerase=noerase

; map image to color indices
if not keyword_set(cran) then begin
  if keyword_set(inverse) then b=-a else b=a
  b_ran=[min(b),max(b)]
  b=(b-b_ran(0))/(b_ran(1)-b_ran(0))*(!d.table_size-1)    ; this is equivalent to tvscl
endif else begin
  if keyword_set(inverse) then begin    
    b=(1.-((a>cran(0)<cran(1))-cran(0))/(cran(1)-cran(0)))*(!d.table_size-1)
  endif else begin
    b=((a>cran(0)<cran(1))-cran(0))/(cran(1)-cran(0))*(!d.table_size-1)
  endelse
endelse

; plot image scaled to the frame size
if not keyword_set(scalable) then begin
; non-scalable pixel case (to be used for xterm device)
  if keyword_set(frame) then begin
    px=!x.window*!d.x_vsize     ; position of window in device pixels
    py=!y.window*!d.y_vsize
  endif else begin
    px=!x.region*!d.x_vsize     ; position of window in device pixels
    py=!y.region*!d.y_vsize    
  endelse
  sx=px(1) - px(0) +1         ; desired size of image in pixels
  sy=py(1) - py(0) +1
  tv,congrid(b,sx,sy),px(0),py(0)
endif else begin
; scalable pixel case (to be used for ps device)
  if keyword_set(frame) then begin  
    if not keyword_set(position) then begin
      tv,b,$
        !x.window(0),!y.window(0),$
        xsize=!x.window(1)-!x.window(0),$
        ysize=!y.window(1)-!y.window(0),/norm
     endif else begin 
       tv, b, $
         position[0],position[1],        $
         xsize=position[2]-position[0],  $
         ysize=position[3]-position[1],/NORM
     endelse
  endif else begin
    tv,b,$
      !x.region(0),!y.region(0),$
      xsize=!x.region(1)-!x.region(0),$
      ysize=!y.region(1)-!y.region(0),/norm 
  endelse   
endelse	

; draw frame if requested
if keyword_set(frame) then begin
  fsupress=0
  if n_elements(frame) eq 4 then begin
    xr=[frame(0),frame(1)] & yr=[frame(2),frame(3)]
  endif else begin
    if frame eq -1 then begin
      fsupress=1 
    endif else begin
      xr=[0.,nx] & yr=[0.,ny]
    endelse
  endelse
  if not keyword_set(fsup) and fsupress eq 0 then $
    plot,[0],[0],/nodata,/noerase,/xstyle,/ystyle,$
      xrange=xr,yrange=yr,$
      title=title,xtitle=xtitle,ytitle=ytitle, $
      CHARSIZE=charsize, POSITION=position, TICKLEN=ticklen, $
      XGRID=xgrid, YGRID=ygrid, XTICKINTERVAL=tickinterval, $
      YTICKINTERVAL=tickinterval
endif

end


