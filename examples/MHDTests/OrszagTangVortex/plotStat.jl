##############################################################################
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

using Plots
#using GR
#using Statistics
using DelimitedFiles
################################################################################
############################FUNCTIONS###########################################
################################################################################

################################################################################
############################DO THINGS###########################################
################################################################################
file=ARGS[1]
if isfile(file)
    re=readdlm("stat_BASE.txt", skipstart=107)
    ri=readdlm(file, skipstart=107)
    eTim=re[:,2]
    eb2=re[:,35]
    edivB=re[:,36]
    eHcro=re[:,37]
    eHmag=re[:,38]
      
    iTim=ri[:,2]
    ib2=ri[:,35]
    idivB=ri[:,36]
    iHcro=ri[:,37]
    iHmag=ri[:,38]
        
    lab=file[begin:end-15]
        
    h1=plot([eTim,iTim],[eb2,ib2],ylabel="b2",label=["Norm" lab])
    h2=plot([eTim,iTim],[edivB,idivB],ylabel="DivB",label=["Norm" lab])
    h3=plot([eTim,iTim],[eHcro,iHcro],ylabel="Hcro",label=["Norm" lab])
    h4=plot([eTim,iTim],[eHmag,iHmag],ylab="Hmag",label=["Norm" lab])
    plot(h1,h2,h3,h4,layout=4, xlabel="Time" , size=(800,600))
    savefig("StatOT.png")
else
    println("File ",file," *** NOT FOUND ***")
end

