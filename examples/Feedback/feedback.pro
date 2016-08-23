base = 'Feedback'
inf  = 'Feedback_003.hdf5'

blast  = [0.5,0.5,0.5] ; location of blast
pos    = h5rd(inf,'PartType0/Coordinates')
vel    = h5rd(inf,'PartType0/Velocities')
rho    = h5rd(inf,'PartType0/Density')
utherm = h5rd(inf,'PartType0/InternalEnergy')

; shift to centre
for ic=0,2 do pos[ic,*] = pos[ic,*] - blast[ic]

;; distnace from centre
dist = fltarr(n_elements(rho))
for ic=0,2 do dist = dist + pos[ic,*]^2
dist = sqrt(dist)

; radial velocity
vr = fltarr(n_elements(rho))
for ic=0,2 do vr = vr + pos[ic,*]*vel[ic,*]
vr = vr / dist
end
