./configure --with-hydro=remix --with-equation-of-state=planetary --with-kernel=wendland-C2 --enable-material-strength --with-strength-artificial-stress=basis-indp --with-forcing=boundary-particles

../../../swift --hydro --threads=8 --limiter impact.yml
