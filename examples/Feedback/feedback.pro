base = 'Feedback'
inf  = 'Feedback_005.hdf5'

blast  = [5.650488e-01, 5.004371e-01, 5.010494e-01] ; location of blast
pos    = h5rd(inf,'PartType0/Coordinates')
vel    = h5rd(inf,'PartType0/Velocities')
rho    = h5rd(inf,'PartType0/Density')
utherm = h5rd(inf,'PartType0/InternalEnergy')

; shift to centre
for ic=0,2 do pos[ic,*] = pos[ic,*] - blast[ic]

;; distance from centre
dist = fltarr(n_elements(rho))
for ic=0,2 do dist = dist + pos[ic,*]^2
dist = sqrt(dist)

; radial velocity
vr = fltarr(n_elements(rho))
for ic=0,2 do vr = vr + pos[ic,*]*vel[ic,*]
vr = vr / dist

; 
end
