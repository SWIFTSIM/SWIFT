./configure --with-hydro=remix --with-equation-of-state=planetary --with-kernel=wendland-C2 --enable-material-strength --with-strength-artificial-stress=basis-indp --with-forcing=boundary-particles

../../../swift --hydro --threads=8 --limiter plate.yml

ffmpeg -framerate 60 -i ./images/plate_%04d.png -c:v libx264 -pix_fmt yuv420p plate.mp4