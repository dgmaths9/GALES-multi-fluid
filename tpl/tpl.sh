#!/bin/bash



TPL_PATH=$(pwd)
N=128



mkdir $TPL_PATH/tpl
mkdir $TPL_PATH/src




echo '#*********************** GALES ********************************' >> ~/.bashrc
echo "TPL_PATH=$TPL_PATH" >> ~/.bashrc  
echo  >> ~/.bashrc 








#--------cmake-------------------------------------
echo '========================== installing cmake  ==========================='
cd $TPL_PATH/tars 
tar -C $TPL_PATH/src/ -xf cmake-3.12.3.tar.gz 
cd $TPL_PATH/src/cmake-3.12.3 
mkdir $TPL_PATH/tpl/cmake-3.12.3 
./bootstrap --prefix=$TPL_PATH/tpl/cmake-3.12.3 | tee -a $TPL_PATH/log.txt 
make -j $N | tee -a $TPL_PATH/log.txt
make install | tee -a $TPL_PATH/log.txt
echo 'export PATH=$TPL_PATH/tpl/cmake-3.12.3/bin:$PATH' >> ~/.bashrc 
source ~/.bashrc
. ~/.bashrc







#-------lapack-------------------------------------
echo '========================== installing lapack  ==========================='
cd $TPL_PATH/tars
tar -C $TPL_PATH/src/ -xf lapack-3.9.0.tar.gz
cd $TPL_PATH/src/lapack-3.9.0
mkdir $TPL_PATH/tpl/lapack-3.9.0 
mkdir build 
cd build
cmake -DCMAKE_INSTALL_PREFIX=$TPL_PATH/tpl/lapack-3.9.0 -DCMAKE_BUILD_TYPE=RELEASE -DBUILD_SHARED_LIBS=ON ..  | tee -a $TPL_PATH/log.txt
make -j $N | tee -a $TPL_PATH/log.txt
make install | tee -a $TPL_PATH/log.txt
echo 'export PATH=$TPL_PATH/tpl/lapack-3.9.0/bin:$PATH' >> ~/.bashrc 
source ~/.bashrc
. ~/.bashrc









#------------boost--------------------------------------
echo '========================== installing boost  ==========================='
cd $TPL_PATH/tars 
tar -C $TPL_PATH/src/ -xf boost_1_70_0.tar.gz 
cd $TPL_PATH/src/boost_1_70_0 
mkdir $TPL_PATH/tpl/boost-1.70.0 
./bootstrap.sh --prefix=$TPL_PATH/tpl/boost-1.70.0    | tee -a $TPL_PATH/log.txt
./bjam -j $N install    | tee -a $TPL_PATH/log.txt
echo 'export PATH=$TPL_PATH/tpl/boost-1.70.0/include:$PATH'  >> ~/.bashrc
echo 'export LD_LIBRARY_PATH=$TPL_PATH/tpl/boost-1.70.0/lib:$LD_LIBRARY_PATH'  >> ~/.bashrc
source ~/.bashrc
. ~/.bashrc







#-----------Trilinos---------------------------------------
echo '========================== installing trilinos  ==========================='
cd $TPL_PATH/tars 
unzip Trilinos.zip -d $TPL_PATH/src/ 
cd $TPL_PATH/src/Trilinos 
mkdir $TPL_PATH/tpl/trilinos-12.17 
mkdir build  
cd build 
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
      ..    | tee -a $TPL_PATH/log.txt
make -j $N     | tee -a $TPL_PATH/log.txt
make install    | tee -a $TPL_PATH/log.txt
echo 'export PATH=$TPL_PATH/tpl/trilinos-12.17:$PATH'  >> ~/.bashrc
echo '#***********************************************************************' >> ~/.bashrc
echo  >> ~/.bashrc       
. ~/.bashrc
source ~/.bashrc






. ~/.bashrc
source ~/.bashrc
. ~/.bashrc






