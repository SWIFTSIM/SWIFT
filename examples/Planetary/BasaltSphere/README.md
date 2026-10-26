./configure --with-hydro=remix --with-equation-of-state=planetary --with-kernel=wendland-C2 --enable-material-strength --with-strength-artificial-stress=basis-independent --with-strength-yield-stress=benz-asphaug --with-strength-yield-stress-weakening-density=power-law --with-strength-damage-tensile=benz-asphaug

../../../swift --hydro --threads=8 --limiter impact.yml