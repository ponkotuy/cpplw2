
prmin = 0.56
prmax = 0.62

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
end
