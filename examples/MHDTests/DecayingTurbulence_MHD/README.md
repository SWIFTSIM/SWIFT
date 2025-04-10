Configuration:
./configure --with-spmhd=direct-induction --with-kernel=quintic-spline --with-equation-of-state=isothermal-gas --disable-hand-vec --disable-doxygen-doc

Running:
../../../swift --hydro --threads=4 DecayingTurbulence.yml
