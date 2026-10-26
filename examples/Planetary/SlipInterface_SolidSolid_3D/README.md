./configure --with-hydro=remix --with-equation-of-state=planetary --with-kernel=wendland-C2 --enable-material-strength --with-strength-artificial-stress=basis-independent --enable-strength-object-ids

../../../swift --hydro --threads=8 --limiter solid_solid.yml