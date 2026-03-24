#!/bin/bash                                                                                                                                                                                              

BASE_DIR="$( cd "$( dirname "${BASH_SOURCE[0]:-${(%):-%x}}" )" && pwd )"
source ${BASE_DIR}/global_vars.sh         

git clone https://github.com/neut-devel/neut.git 
cd neut
git checkout tags/6.1.4
mkdir build; cd build;
cmake ../ -DNEUT_DOWNLOAD_DATA=OFF -DNEUT_WERROR_ENABLED=OFF
make -j 4
make install
source $(pwd)/Linux/bin/setup.NEUT.sh 