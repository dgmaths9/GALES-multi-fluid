#ifndef GALES_MPI_FUNCTIONS_HPP
#define GALES_MPI_FUNCTIONS_HPP


#include<mpi.h>
#include<vector>



namespace GALES{




   /// This returns the rank of process.
   int get_rank()
   {
      int rank;
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);
      return rank;      
   }



   /// This returns the size of pool.
   int get_size()
   {
      int size;
      MPI_Comm_size(MPI_COMM_WORLD, &size);
      return size;      
   }




  /// This function returns rank, size, num_el_buffer, u and sum.
  /// This is called in fsi coupling classes
  void MPI_NumElbuffer_U_Sum(int send_buffer_size, std::vector<int>& num_el_buffer_of_each_pid, std::vector<int>& u, int& sum)
  {
      int size = get_size();
      num_el_buffer_of_each_pid.resize(size);   
      MPI_Allgather(&send_buffer_size, 1, MPI_INT, &num_el_buffer_of_each_pid[0], 1, MPI_INT, MPI_COMM_WORLD);

      u.resize(size);
      sum = 0;
      for(int i=0; i<size; i++)
      {
        u[i] = sum;
        sum += num_el_buffer_of_each_pid[i];
      }


      //   example for 4 pids 
      //               pid =   0   1    2   3
      //        s_ux_dof s =  ---|----|---|-----
      //    num_el_buffer  =   3   4    3   5
      //    position     u = 0---3----7---10-----   
      //               sum = 15    
  }






}


#endif
