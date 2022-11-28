#!/bin/bash




TPL_PATH=$(pwd)
N=32



mkdir $TPL_PATH/src 
mkdir $TPL_PATH/tpl  



sudo apt-get update 
sudo apt-get upgrade
sudo apt-get upgrade
sudo apt-get install build-essential 
sudo apt-get install gfortran -y
sudo apt-get install libblas-dev liblapack-dev  -y
sudo apt-get install libopenmpi-dev openmpi-bin libhdf5-openmpi-dev  -y
sudo apt-get install cmake -y 
sudo apt-get install metis -y
sudo apt-get install libboost-all-dev -y




#----------gmsh 3.0.6-----------------------------------
echo '========================== installing gmsh  ==========================='
cd $TPL_PATH/tars 
tar -C $TPL_PATH/tpl/ -xf gmsh-3.0.6.tar.gz 




#-----------Trilinos---------------------------------------
cd $TPL_PATH/tars 
unzip Trilinos.zip -d $TPL_PATH/src/ 
cd $TPL_PATH/src/Trilinos 
mkdir build && cd build 
cmake \
      -D CMAKE_INSTALL_PREFIX:PATH="$TPL_PATH/tpl/trilinos-12.17" \
      -D CMAKE_BUILD_TYPE:STRING=RELEASE \
      -D Trilinos_ENABLE_EpetraExt=ON \
      -D Trilinos_ENABLE_Epetra=ON \
      -D Trilinos_ENABLE_Belos=ON \
      -D Trilinos_ENABLE_Teuchos=ON \
      -D Trilinos_ENABLE_Ifpack=ON \
      -D CMAKE_CXX_STANDARD=14 \
      -D TPL_ENABLE_MPI=ON \
      -D TPL_ENABLE_Boost=ON \
      -D BUILD_SHARED_LIBS=ON \
      -D Trilinos_ENABLE_ALL_OPTIONAL_PACKAGES=ON \
      -D Trilinos_ENABLE_TESTS=ON \
      -D Trilinos_ENABLE_EXPLICIT_INSTANTIATION=ON \
      .. 
make -j $N  
make install 








echo '#*********************** GALES ********************************' >> ~/.bashrc
echo "TPL_PATH=$TPL_PATH" >> ~/.bashrc  
echo  >> ~/.bashrc 
echo 'export PATH=$TPL_PATH/gales_mesh_preprocessing:$PATH' >> ~/.bashrc 
echo 'export PYTHONPATH=$TPL_PATH/gales_results_postprocessing:$PYTHONPATH' >> ~/.bashrc 
echo 'export PATH=$TPL_PATH/solwcad:$PATH' >> ~/.bashrc 
echo 'export PATH=$TPL_PATH/tpl/gmsh-3.0.6/bin:$PATH'  >> ~/.bashrc
echo 'export PATH=$TPL_PATH/tpl/trilinos-12.17:$PATH'  >> ~/.bashrc
echo '#***********************************************************************' >> ~/.bashrc
echo  >> ~/.bashrc       



source ~/.bashrc
. ~/.bashrc

#--------------solwcad rebuild---------------------------------------
cd $TPL_PATH/solwcad
./build




source ~/.bashrc



