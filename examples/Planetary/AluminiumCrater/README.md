./configure --with-hydro=remix --with-equation-of-state=planetary --with-kernel=wendland-C2 --enable-material-strength --with-strength-artificial-stress=basis-independent --with-strength-yield-stress=collins --with-strength-yield-stress-weakening-thermal=tanh-melt --with-strength-yield-stress-weakening-density=power-law --with-forcing=boundary-particles

../../../swift --hydro --threads=8 --limiter impact.yml
