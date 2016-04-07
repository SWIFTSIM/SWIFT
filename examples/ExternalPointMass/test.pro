;
;  test energy / angular momentum conservation of test problem
;
@physunits

indir    = '/gpfs/data/tt/Codes/Swift-git/swiftsim/examples/'
basefile = 'output_'
nfiles   = 657
nfollow  = 100 ; number of particles to follow
eout     = fltarr(nfollow, nfiles)
ekin     = fltarr(nfollow, nfiles)
epot     = fltarr(nfollow, nfiles)
tout     = fltarr(nfiles)
; set properties of potential
uL  = 1e3 * phys.pc             ; unit of length
uM  = phys.msun                 ; unit of mass
uV  = 1d5                       ; unit of velocity

; derived units
constG   = 10.^(alog10(phys.g)+alog10(uM)-2d0*alog10(uV)-alog10(uL)) ;
pcentre  = [50.,50.,50.] * 1d3 * pc / uL
mextern  = 1d10 * msun / uM
;
;
;
ifile  = 0
for ifile=0,nfiles-1 do begin
;for ifile=0,3 do begin
   inf    = indir + basefile + strtrim(string(ifile,'(i3.3)'),1) + '.hdf5'
   time   = h5ra(inf, 'Header','Time')
   p      = h5rd(inf,'PartType1/Coordinates')
   v      = h5rd(inf,'PartType1/Velocities')
   id     = h5rd(inf,'PartType1/ParticleIDs')
   indx   = sort(id)
;
   id     = id[indx]
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
   dr = sqrt(dr)
;   print,'time = ',time,p[0,0],v[0,0],id[0]
   ek   = 0.5 * dv
   ep   = - constG * mextern / dr
   ener = ek + ep
   tout(ifile) = time
   eout(*,ifile) = ener[0:nfollow-1]
   ekin(*,ifile) = ek[0:nfollow-1]
   epot(*,ifile) = ep[0:nfollow-1]
endfor

; calculate relative energy change
de = 0.0 * eout
for ifile=1, nfiles -1 do de[*,ifile] = (eout[*,ifile]-eout[*,0])/eout[*,0]


end


