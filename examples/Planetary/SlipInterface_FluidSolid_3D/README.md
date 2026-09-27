./configure --with-hydro=remix --with-equation-of-state=planetary --with-kernel=wendland-C2 --enable-material-strength --with-strength-artificial-stress=basis-indp

../../../swift --hydro --threads=8 --limiter fluid_solid.yml