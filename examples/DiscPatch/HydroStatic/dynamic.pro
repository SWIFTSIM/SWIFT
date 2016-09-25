;
;  test energy / angular momentum conservation of test problem
;

iplot = 1 ; if iplot = 1, make plot of E/Lz conservation, else, simply compare final and initial energy

; set physical constants
@physunits

indir    = './'
;basefile = 'Disc-Patch-dynamic_'
basefile = 'Disc-Patch_'

; set properties of potential
uL   = phys.pc                  ; unit of length
uM   = phys.msun                ; unit of mass
uV   = 1d5                      ; unit of velocity

; properties of patch
surface_density = 100.          ; surface density of all mass, which generates the gravitational potential
scale_height    = 100.
z_disk          = 200.          ;
fgas            = 0.1           ; gas fraction
gamma           = 5./3.

; derived units
constG   = 10.^(alog10(phys.g)+alog10(uM)-2d0*alog10(uV)-alog10(uL)) ;
pcentre  = [0.,0.,z_disk] * pc / uL
utherm     = !pi * constG * surface_density * scale_height / (gamma-1.)
temp       = (utherm*uV^2)*phys.m_h/phys.kb
soundspeed = sqrt(gamma * (gamma-1.) * utherm)
t_dyn      = sqrt(scale_height / (constG * surface_density))
rho0       = fgas*(surface_density)/(2.*scale_height)
print,' dynamical time = ',t_dyn,' = ',t_dyn*UL/uV/(1d6*phys.yr),' Myr'
print,' thermal energy per unit mass = ',utherm
print,' central density = ',rho0,' = ',rho0*uM/uL^3/m_h,' particles/cm^3'
print,' central temperature = ',temp
lambda = 2 * !pi * phys.G^1.5 * (scale_height*uL)^1.5 * (surface_density * uM/uL^2)^0.5 * phys.m_h^2 / (gamma-1) / fgas
print,' lambda = ',lambda
stop
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


; compute anlytic profile
nbins = 100
zbins = findgen(nbins)/float(nbins-1) * 2 * scale_height
rbins = (surface_density/(2.*scale_height)) / cosh(abs(zbins)/scale_height)^2


; plot analytic profile
wset,0
plot,[0],[0],xr=[0,2*scale_height],yr=[0,max(rbins)],/nodata,xtitle='|z|',ytitle=textoidl('\rho')
oplot,zbins,rbins,color=blue

ifile  = 0
nskip   = nfiles - 1
isave  = 0
nplot  = 8192 ; randomly plot particles
color = floor(findgen(nfiles)/float(nfiles-1)*255)
;for ifile=0,nfiles-1,nskip do begin
tsave  = [0.]
toplot = [1,nfiles-1]
for iplot=0,1 do begin
   ifile  = toplot[iplot]
   inf    = indir + basefile + strtrim(string(ifile,'(i3.3)'),1) + '.hdf5'
   time   = h5ra(inf, 'Header','Time')
   tsave  = [tsave, time]
   print,' time= ',time
   p      = h5rd(inf,'PartType0/Coordinates')
   v      = h5rd(inf,'PartType0/Velocities')
   id     = h5rd(inf,'PartType0/ParticleIDs')
   rho    = h5rd(inf,'PartType0/Density')
   h      = h5rd(inf,'PartType0/SmoothingLength')
   utherm = h5rd(inf,'PartType0/InternalEnergy')
   indx   = sort(id)

; substract disk centre
   for ic=0,2 do p[ic,*]=p[ic,*] - pcentre[ic]


;; ;  if you want to sort particles by ID
;;    id     = id[indx]
;;    rho    = rho[indx]
;;    utherm = utherm[indx]
;;    h      = h[indx]
;;    for ic=0,2 do begin
;;       tmp = reform(p[ic,*]) & p[ic,*] = tmp[indx]
;;       tmp = reform(v[ic,*]) & v[ic,*] = tmp[indx]
;;    endfor
   
   ip = floor(randomu(ifile+1,nplot)*n_elements(rho))
   color = red
   if(ifile eq 1) then begin
      color=black
   endif else begin
      color=red
   endelse
   oplot,abs(p[2,ip]), rho[ip], psym=3, color=color

   isave = isave + 1
   
endfor

; time in units of dynamical time
tsave = tsave[1:*] / t_dyn

label = ['']
for i=0,n_elements(tsave)-1 do label=[label,'time/t_dynamic='+string(tsave[i],format='(f8.0)')]
label = label[1:*]
legend,['analytic',label[0],label[1]],linestyle=[0,0,0],color=[blue,black,red],box=0,/top,/right

; make histograms of particle velocities
xr    = 1d-3 * [-1,1]
bsize = 1.d-5
ohist,v[0,*]/soundspeed,x,vx,xr[0],xr[1],bsize
ohist,v[1,*]/soundspeed,y,vy,xr[0],xr[1],bsize
ohist,v[2,*]/soundspeed,z,vz,xr[0],xr[1],bsize
wset,2
plot,x,vx,psym=10,xtitle='velocity/soundspeed',ytitle='pdf',/nodata,xr=xr,/xs
oplot,x,vx,psym=10,color=black
oplot,y,vy,psym=10,color=blue
oplot,z,vz,psym=10,color=red
legend,['vx/c','vy/c','vz/c'],linestyle=[0,0,0],color=[black,blue,red],box=0,/top,/right
end


