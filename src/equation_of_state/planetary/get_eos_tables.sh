#!/bin/bash

# Download the tables of the publicly available planetary equations of state
wget http://virgodb.cosma.dur.ac.uk/swift-webstorage/EoS/planetary_HM80_HHe.txt
wget http://virgodb.cosma.dur.ac.uk/swift-webstorage/EoS/planetary_HM80_ice.txt
wget http://virgodb.cosma.dur.ac.uk/swift-webstorage/EoS/planetary_HM80_rock.txt

mv planetary_HM80_HHe.txt ../../../examples/
mv planetary_HM80_ice.txt ../../../examples/
mv planetary_HM80_rock.txt ../../../examples/
