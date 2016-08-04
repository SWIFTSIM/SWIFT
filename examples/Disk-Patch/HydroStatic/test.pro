;
;  test energy / angular momentum conservation of test problem
;

iplot = 1 ; if iplot = 1, make plot of E/Lz conservation, else, simply compare final and initial energy

; set physical constants
@physunits

indir    = './'
basefile = 'Disk-Patch_'

; set properties of potential
uL   = phys.pc                  ; unit of length
uM   = phys.msun                ; unit of mass
uV   = 1d5                      ; unit of velocity

; properties of patch
surface_density = 10.
scale_height    = 100.

; derived units
constG   = 10.^(alog10(phys.g)+alog10(uM)-2d0*alog10(uV)-alog10(uL)) ;
pcentre  = [0.,0.,300.] * pc / uL

;
infile = indir + basefile + '*'
spawn,'ls -1 '+infile,res
nfiles = n_elements(res)


; choose: calculate change of energy and Lz, comparing first and last
; snapshots for all particles, or do so for a subset

; compare all
ifile   = 0
inf     = indir + basefile + strtrim(string(ifile,'(i3.3)'),1) + '.hdf5'
id      = h5rd(inf,'PartType0/ParticleIDs')
nfollow = n_elements(id)

; follow a subset
nfollow  = min(4000, nfollow)   ; number of particles to follow

;
if (iplot eq 1) then begin
   nskip = 1
   nsave = nfiles
endif else begin
   nskip = nfiles - 2
   nsave = 2
endelse

;
lout     = fltarr(nfollow, nsave) ; Lz
xout     = fltarr(nfollow, nsave) ; x
yout     = fltarr(nfollow, nsave) ; y
zout     = fltarr(nfollow, nsave) ; z
vzout    = fltarr(nfollow, nsave) ; z
rout     = fltarr(nfollow, nsave) ; rho
hout     = fltarr(nfollow, nsave) ; h
uout     = fltarr(nfollow, nsave) ; thermal energy
eout     = fltarr(nfollow, nsave) ; energies
ekin     = fltarr(nfollow, nsave)
epot     = fltarr(nfollow, nsave) ; 2 pi G Sigma b ln(cosh(z/b)) + const
tout     = fltarr(nsave)

ifile  = 0
isave = 0
for ifile=0,nfiles-1,nskip do begin
   inf    = indir + basefile + strtrim(string(ifile,'(i3.3)'),1) + '.hdf5'
   time   = h5ra(inf, 'Header','Time')
   p      = h5rd(inf,'PartType0/Coordinates')
   v      = h5rd(inf,'PartType0/Velocities')
   id     = h5rd(inf,'PartType0/ParticleIDs')
   rho    = h5rd(inf,'PartType0/Density')
   h      = h5rd(inf,'PartType0/SmoothingLength')
   utherm = h5rd(inf,'PartType0/InternalEnergy')
   indx   = sort(id)

;  if you want to sort particles by ID
   id     = id[indx]
   rho    = rho[indx]
   utherm = utherm[indx]
   h      = h[indx]
   for ic=0,2 do begin
      tmp = reform(p[ic,*]) & p[ic,*] = tmp[indx]
      tmp = reform(v[ic,*]) & v[ic,*] = tmp[indx]
   endfor

; calculate energy
   dd  = size(p,/dimen) & npart = dd[1]
   ener = fltarr(npart)
   dr   = fltarr(npart) & dv = dr
   for ic=0,2 do dr[*] = dr[*] + (p[ic,*]-pcentre[ic])^2
   for ic=0,2 do dv[*] = dv[*] + v[ic,*]^2
   xout[*,isave] = p[0,0:nfollow-1]-pcentre[0]
   yout[*,isave] = p[1,0:nfollow-1]-pcentre[1]
   zout[*,isave] = p[2,0:nfollow-1]-pcentre[2]
   vzout[*,isave]= v[2,0:nfollow-1]
   rout[*,isave] = rho[0:nfollow-1]
   hout[*,isave] = h[0:nfollow-1]
   uout[*,isave] = utherm[0:nfollow-1]
   Lz  = (p[0,*]-pcentre[0]) * v[1,*] - (p[1,*]-pcentre[1]) * v[0,*]
   dz  = reform(p[2,0:nfollow-1]-pcentre[2])
;   print,'time = ',time,p[0,0],v[0,0],id[0]
   ek   = 0.5 * dv
   ep   = fltarr(nfollow)
   ep   = 2 * !pi * constG * surface_density * scale_height * alog(cosh(abs(dz)/scale_height))
   ener = ek + ep
   tout(isave) = time
   lout[*,isave] = lz[0:nfollow-1]
   eout(*,isave) = ener[0:nfollow-1]
   ekin(*,isave) = ek[0:nfollow-1]
   epot(*,isave) = ep[0:nfollow-1]
   print,format='('' time= '',f7.1,'' E= '',f9.2,'' Lz= '',e9.2)', time,eout[0],lz[0]
   isave = isave + 1
   
endfor

x0 = reform(xout[0,*])
y0 = reform(xout[1,*])
z0 = reform(xout[2,*])

; calculate relative energy change
de    = 0.0 * eout
dl    = 0.0 * lout
nsave = isave
for ifile=1, nsave-1 do de[*,ifile] = (eout[*,ifile]-eout[*,0])/eout[*,0]
for ifile=1, nsave-1 do dl[*,ifile] = (lout[*,ifile] - lout[*,0])/lout[*,0]


; calculate statistics of energy changes
print,' relatve energy change: (per cent) ',minmax(de) * 100.
print,' relative Lz    change: (per cent) ',minmax(dl) * 100.

; plot enery and Lz conservation for some particles
if(iplot eq 1) then begin
   nplot = nfollow

   ; plot density profile
   wset,0
   xr   = [0, 3*scale_height]
   nbins = 100
   zpos  = findgen(nbins)/float(nbins-1) * max(xr)
   dens  = (surface_density/(2.d0*scale_height)) * 1./cosh(zpos/scale_height)^2
   plot,[0],[0],xr=xr,/xs,yr=[0,max(dens)*1.4],/ys,/nodata,xtitle='|z|',ytitle='density'
   oplot,zpos,dens,color=black,thick=3
   oplot,abs(zout[*,1]),rout[*,1],psym=3
   oplot,abs(zout[*,nsave-1]),rout[*,nsave-1],psym=3,color=red

;; ; plot results on energy conservation for some particles
;;    nplot = min(100, nfollow)
;;    win,0
;;    xr = [min(tout), max(tout)]
;;    yr = [-2,2]*1d-2             ; in percent
;;    plot,[0],[0],xr=xr,yr=yr,/xs,/ys,/nodata,xtitle='time',ytitle='dE/E, dL/L (%)'
;;    for i=0,nplot-1 do oplot,tout,de[i,*]
;;    for i=0,nplot-1 do oplot,tout,dl[i,*],color=red
;;    legend,['dE/E','dL/L'],linestyle=[0,0],color=[black,red],box=0,/bottom,/left
;;    screen_to_png,'e-time.png'

;  plot vertical oscillation
   wset,2
   xr = [min(tout), max(tout)]
   yr = [-3,3]*scale_height
   plot,[0],[0],xr=xr,yr=yr,/xs,/ys,/nodata,xtitle='t',ytitle='z(t)'
   color = floor(findgen(nplot)*255/float(nplot))
   for i=0,nplot-1,50 do oplot,tout,zout[i,*],color=color(i)
   screen_to_png,'orbit.png'


;   plot evolution of density for some particles close to disk
   wset, 6
   rr = (surface_density/(2.d0*scale_height)) * 1./cosh(abs(zout)/100.)^2 ; desired density
   gd = where(abs(zout[*,0]) lt 50, ng)
   plot,[0],[0],xr=[0, max(tout)], yr=10.^[-1,1], /xs, /ys, /nodata, /yl
   nplot = min(40, ng)
   color = floor(findgen(nplot)/float(nplot)*255)
   for i=0, nplot-1 do oplot,tout[1:*],rout[gd[i],1:*]/rr[gd[i], 1:*], color=color[i],psym=-4


;; ; make histogram of energy changes at end
;;    win,6
;;    ohist,de,x,y,-0.05,0.05,0.001
;;    plot,x,y,psym=10,xtitle='de (%)'
;;    screen_to_png,'de-hist.png'


endif

end


