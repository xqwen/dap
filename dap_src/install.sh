#!/bin/sh

rm -rf bin/ 
mkdir bin
cd bin/
echo -e "Generating MakeFiles ..."
cmake ..
make
echo "Binary created at /bin/dap-g"

