###############################################################################
# This file is part of SWIFT.
# Copyright (c) 2016 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
#               2018 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
#               2022 Federico Stasyszyn (fstasyszyn@unc.edu.ar)
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
################################################################################

println("Plotting ... ",ARGS[1])

using HDF5
using GRUtils
using Statistics

################################################################################
############################FUNCTIONS###########################################
################################################################################
function read_snap(filename :: String)
    
    pos  = h5read(filename,"PartType0/Coordinates")
    Bfl  = h5read(filename,"PartType0/MagneticFluxDensities") 
    Vel  = h5read(filename,"PartType0/Velocities")
    divB = h5read(filename,"PartType0/MagneticDivergences")
    #Ids  = h5read(filename,"PartType0/ParticleIDs")
    h    = h5read(filename,"PartType0/SmoothingLengths")
    rho  = h5read(filename,"PartType0/Densities")
    u    = h5read(filename,"PartType0/InternalEnergies")
    P    = h5read(filename,"PartType0/Pressures")
    head = h5readattr(filename,"Header")

    print("Reading ",filename," at time: ",head["Time"],"\n")
#    x = pos[1,:]
#    y = pos[2,:]
#    z = pos[3,:]
    bx = Bfl[1,:]
    by = Bfl[2,:]
    bz = Bfl[3,:]
    Vx = Vel[1,:]
    Vy = Vel[2,:]
    Vz = Vel[3,:]
    
    v2 = Vx.*Vx.+Vy.*Vy.+Vz.*Vz
    b2 = by.*by.+by.*by.+bz.*bz
    
    Npart=size(b2,1)
        
#    print("Min x:",minimum(x)," / Max x:",maximum(x),"\n")
#    print("Min y:",minimum(y)," / Max y:",maximum(y),"\n")
#    print("Min z:",minimum(z)," / Max z:",maximum(z),"\n")
#    print("Min h:",minimum(h)," / Max h:",maximum(h),"\n")
#    print("Min B2:",minimum(b2)," / Max B2:",maximum(b2),"\n")
    (Dict(:H => head, :x=>pos, :bfl=> Bfl, :b2=>b2, :v=>Vel, :v2=>v2, 
            :divB=>divB, :rho=>rho, :hsml=>h, #:bet=> bet, :alp=>alp, 
            :Pres=> P, u=>u, :Npart=> Npart))
end

function do_table(data,Nmax)
    #### Data
    x=data[:x][1,:]
    xmin=minimum(x)
    xmax=maximum(x)
    zone=(round(xmax)-round(xmin))/2.
    println("Zones ",zone)
    x = x .- zone # move to 0
    x = 6.0 .* x ./ zone #normalize -1 to 1 just one side of the tube
    #x = x ./ zone #normalize -1 to 1 just one side of the tube
    
    
    rho=data[:rho]
    divB=data[:divB]
    Bx=data[:bfl][1,:]
    By=data[:bfl][2,:]
    Bz=data[:bfl][3,:]
    
    Vx=data[:v][1,:]
    Vy=data[:v][2,:]
    Vz=data[:v][3,:]
    
    divB=data[:hsml] .* data[:divB]./sqrt.(data[:b2])
    
    
    delta=(2.0)/Nmax
    #### Out data
    xx=range(-1.0,1.0,Nmax) # xx[1]>> min  xx[Nmax]>> max
    rrho =ones(Nmax-1,2)
    ddivB=ones(Nmax-1,2)
    BBx=ones(Nmax-1,2)
    BBy=ones(Nmax-1,2)
    BBz=ones(Nmax-1,2)
    VVx=ones(Nmax-1,2)
    VVy=ones(Nmax-1,2)
    VVz=ones(Nmax-1,2)
    Pre=ones(Nmax-1,2)
    
    for ind=1:Nmax-1
        idx=(x .>= xx[ind] .&& x .< xx[ind+1])
        rrho[ind,1] , rrho[ind,2]  = mean(rho[idx])  , std(rho[idx])
        ddivB[ind,1], ddivB[ind,2] = mean(divB[idx]) , std(divB[idx])
        BBx[ind,1]  , BBx[ind,2]   = mean(Bx[idx])   , std(Bx[idx])
        BBy[ind,1]  , BBy[ind,2]   = mean(By[idx])   , std(By[idx])
        BBz[ind,1]  , BBz[ind,2]   = mean(Bz[idx])   , std(Bz[idx])
        
        VVx[ind,1]  , VVx[ind,2]   = mean(Vx[idx])   , std(Vx[idx])
        VVy[ind,1]  , VVy[ind,2]   = mean(Vy[idx])   , std(Vy[idx])
        VVz[ind,1]  , VVz[ind,2]   = mean(Vz[idx])   , std(Vz[idx])
        
        Pre[ind,1]  , Pre[ind,2]   = mean(data[:Pres][idx])   , std(data[:Pres][idx])
    end
    
    xx=range(-1.0+delta,1.0-delta,Nmax-1) # xx[1]>> mean   xx[Nmax]>> mean
    
    (Dict(:x => xx, :rho=>rrho, :divB=>ddivB,
            :Bx=>BBx, :By=>BBy , :Bz=> BBz,
            :Vx=>VVx, :Vy=>VVy , :Vz=> VVz,
            :Pres=>Pre))
end

function exact_1A(time)
     tfac = time/0.1
     xmin = -1.0
     xmax = 1.0
     #### def vars
     #vz     = 0.
     #Bz     = 0.
     #Bxzero = 0.75
     npts   = 12
     MU0_1  = 1.0/sqrt(4.0*pi)
     #### 
     xpts = zeros(npts)
     rho  = zeros(npts)
     pr   = zeros(npts)
     vx   = zeros(npts)
     vy   = zeros(npts)
     vz   = zeros(npts)
     bx   = ones(npts).* 5.0 
     bz   = zeros(npts)
     by   = zeros(npts)
     ####
    
     xpts[1]      = xmin
     xpts[2:3]   .= -0.386*tfac
     xpts[4:5]   .= -0.01*tfac
     xpts[6:7]   .= 0.0505*tfac
     xpts[8:9]   .= 0.12*tfac
     xpts[10:11] .= 0.37*tfac
     xpts[12]     = xmax

     rho[1:2]   .= 1.0
     rho[3:4]   .= 2.6793
     rho[5:6]   .= 2.6713
     rho[7:8]   .= 3.8508
     rho[9:10]  .= 3.7481
     rho[11:12] .= 1.0

     pr[1:2]    .= 20.0
     pr[3:4]    .= 150.98
     pr[5:8]    .= 150.19
     pr[9:10]   .= 143.57
     pr[11:12]  .= 1.0
 
     vx[1:2]    .= 10.0
     vx[3:4]    .= 0.72113
     vx[5:8]    .= 0.72376
     vx[9:10]   .= 0.70505
     vx[11:12]  .= -10.0
 
     vy[1:2]    .= 0.0
     vy[3:4]    .= 0.23139
     vy[5:8]    .= 0.35668
     vy[9:10]   .= -0.3884
     vy[11:12]  .= 0.0

     vz[1:12]	.= 0.0
 
     cte = sqrt(4.0 * pi)

     by[1:2]   .= 1.4105 .* cte
     by[3:4]   .= 3.8389 .* cte
     by[5:8]   .= 4.038  .* cte
     by[9:10]  .= 5.4271 .* cte
     by[11:12] .= 1.4105 .* cte

     bz[1:12]  .= 0.0
    
     (Dict(:xpts=>xpts, :rho=>rho, :pr=>pr,
            :vx=> vx, :vy=>vy, :vz=>vz,
            :bx=> bx, :by=>by, :bz=>bz))
end

function do_6plot(gsnap)
    a=do_table(gsnap,256)
    ex=exact_1A(0.2)
    GRUtils.hold(false)
    subplot(1,1,1)
########## plot1
    subplot(2,3,1)
    plot(a[:x],a[:rho][:,1],ylim=(-0.1,4.1),xlabel="X",ylabel="Rho")
    GRUtils.hold(true)
    errorbar(a[:x],a[:rho][:,1],a[:rho][:,2])
    oplot(ex[:xpts],ex[:rho],"-r")
######## plot2
    subplot(2,3,2)
    plot(a[:x],a[:Bx][:,1],ylim=(-1.1,20.1),xlabel="X",ylabel="Bx/By/Bz")
  #  plot(a[:x],a[:Bx][:,1],xlabel="X",ylabel="Bx/By")
    GRUtils.hold(true)
    errorbar(a[:x],a[:Bx][:,1],a[:Bx][:,2])
    oplot(ex[:xpts],ex[:bx],"-r")
    plot(a[:x],a[:By][:,1])#,xlabel="X",ylabel="By")
    errorbar(a[:x],a[:By][:,1],a[:By][:,2])
    oplot(ex[:xpts],ex[:by],"-g")
    plot(a[:x],a[:Bz][:,1])#,xlabel="X",ylabel="By")
    errorbar(a[:x],a[:Bz][:,1],a[:Bz][:,2])
    oplot(ex[:xpts],ex[:bz],"-b")
####### plot3
    subplot(2,3,3)
    plot(a[:x],a[:Vx][:,1],ylim=(-11.1,11.1),xlabel="X",ylabel="Vx/Vz")
    GRUtils.hold(true)
    errorbar(a[:x],a[:Vx][:,1],a[:Vx][:,2])
    oplot(ex[:xpts],ex[:vx],"-r")
    plot(a[:x],a[:Vz][:,1])
    errorbar(a[:x],a[:Vz][:,1],a[:Vz][:,2])
    oplot(ex[:xpts],ex[:vz],"-g")
###### plot 4
    subplot(2,3,4)
    plot(a[:x],a[:Vy][:,1],xlabel="X",ylabel="Vy",ylim=(-0.5,0.5))
    GRUtils.hold(true)
    errorbar(a[:x],a[:Vy][:,1],a[:Vy][:,2])
    oplot(ex[:xpts],ex[:vy],"-r")
###### plot 6    
    subplot(2,3,5)
    plot(a[:x],a[:Pres][:,1],xlabel="X",ylabel="Pressure",ylim=(0,160.1))
    GRUtils.hold(true)
    errorbar(a[:x],a[:Pres][:,1],a[:Pres][:,2])
    oplot(ex[:xpts],ex[:pr],"-r")
###### plot 6
    subplot(2,3,6)
    plot(a[:x],log10.(abs.(a[:divB][:,1])),xlabel="X",ylabel="Err|divB|", ylim=(-6,0))
    #GRUtils.hold(true)
    #errorbar(a[:x],abs.(a[:divB][:,1]),a[:divB][:,2])
end

################################################################################
############################DO THINGS###########################################
################################################################################
if isfile(ARGS[1])
   gsnap=read_snap(ARGS[1])
   do_6plot(gsnap)
   savefig("Solution.jpeg")
else
   println("File ",ARGS[1]," ***  NOT FOUND   *** ")
end
