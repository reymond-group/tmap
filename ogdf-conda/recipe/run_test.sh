# !/bin/sh
g++ -I${PREFIX}/include -o test test.cc -L${PREFIX}/lib -lOGDF
./test