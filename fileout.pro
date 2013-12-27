prmin = 0.56
prmax = 0.62

for step=0, (t-1) do begin
  ;; Density
  loadct, 1
  tvlct, r, g, b, /get
  rr = reverse(r)
  gg = reverse(g)
  bb = reverse(b)
  tvlct, rr, gg, bb

  contour, ro(*,*,step), x, y, nlevels=256, /isotropic, /fill, charsize = chsize,$
           title = "Density(t=" + string(step*dt, FORMAT='(F7.3)') + ')',$
           xtitle = "x-dimension", ytitle = "y-dimension",$
           xrange = xrange, yrange = yrange, /downhill

  xx = fltarr(ix/vect_space)
  yy = fltarr(jx/vect_space)
  vxx = fltarr(ix/vect_space, jx/vect_space)
  vyy = fltarr(ix/vect_space, jx/vect_space)
  vix = ix/vect_space - 1
  vjx = jx/vect_space - 1
  for i=0, vix do xx(i) = x(i*vect_space,1)
  for j=0, vjx do yy(j) = y(1,j*vect_space)
  for j=0, vjx do begin
     for i=0, vix do begin
        vxx(i,j) = vx(i*vect_space, j*vect_space, step)
        vyy(i,j) = vy(i*vect_space, j*vect_space, step)
     endfor
  endfor
  loadct, 33
  velovect, vxx, vyy, xx, yy, length=1, /isotropic, /overplot,$
            charsize = chsize, title = "a", xtitle="b", ytitle="c",$
            xrange = xrange, yrange = yrange

  tvlct, rr, gg, bb
  color_bar, ro(*,*,step), [0.82*xsize, 0.10*ysize, 0.87*xsize, 0.94*ysize],$
             /ver, charsize = chsize

  write_png, 'rho' + strtrim(step,2) + '.png', tvrd(/true)

  ;; Pressure

  loadct, 3
  tvlct, r, g, b, /get
  rr = reverse(r)
  gg = reverse(g)
  bb = reverse(b)
  tvlct, rr, gg, bb

  contour, pr(*,*,step), x, y, nlevels=256, /isotropic, /fill, charsize = chsize,$
           title = "Pressure(t=" + string(step*dt, FORMAT='(F7.3)') + ')',$
           xtitle = "x-dimension", ytitle = "y-dimension", zrange = [prmin, prmax]

  color_bar, pr(*,*,step), [0.82*xsize, 0.10*ysize, 0.87*xsize, 0.94*ysize],$
             /ver, charsize = chsize, dmin = prmin, dmax = prmax

  write_png, 'pr' + strtrim(step,2) + '.png', tvrd(/true)
endfor

end
