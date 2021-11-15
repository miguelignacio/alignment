cd cvt_plots
make clean; make
cd ../table
make clean; make
cd ../hipo2root/lz4
make clean; make
cd ../hipo4
make clean; make
cp lib/*.a ../lib
cd ..
make clean; make

