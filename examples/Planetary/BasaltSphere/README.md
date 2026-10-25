./configure --with-hydro=remix --with-equation-of-state=planetary --with-kernel=wendland-C2 --enable-material-strength --with-strength-artificial-stress=basis-indp --with-strength-yield-stress=benz-asphaug --enable-strength-yield-stress-weakening-density --with-strength-damage=tensile-benz-asphaug --with-strength-max-flaws=34

../../../swift --hydro --threads=8 --limiter impact.yml