  file = 'bOutput'
  ;; Variable
  int = lonarr(3)
  ;; File Read
  openr, 1, file

  readu, 1, int
  print, int
  t = int(0)
  ix = int(1)
  jx = int(2)
  
  arr = fltarr(ix, jx, t*4, /nozero)
  x = fltarr(ix, jx, /nozero)
  y = fltarr(ix, jx, /nozero)
  ro = fltarr(ix, jx, t, /nozero)
  pr = fltarr(ix, jx, t, /nozero)
  vx = fltarr(ix, jx, t, /nozero)
  vy = fltarr(ix, jx, t, /nozero)
  readu, 1, x
  readu, 1, y
  readu, 1, arr
  for i=0, (t-1) do begin
     ro(*,*,i) = arr(*, *, i*4)
     pr(*,*,i) = arr(*, *, i*4 + 1)
     vx(*,*,i) = arr(*, *, i*4 + 2)
     vy(*,*,i) = arr(*, *, i*4 + 3)
  endfor
  close, 1
  ;; Initialize
  dt = 0.512
  xsize = 800
  ysize = 600
  xrange = [min(x), max(x)]
  yrange = [min(y), max(y)]
  chsize = 1.5
  vect_space = 32
  window, xsize=xsize, ysize=ysize
  step = 0
end
