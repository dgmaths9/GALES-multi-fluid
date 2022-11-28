#ifndef HEAT_EQN_T_HPP
#define HEAT_EQN_T_HPP



#include <vector>
#include <map>




namespace GALES{



  /**
       This class computes temperature at the interface of fluid and heat equation domains, to be passed to fluid solver for coupling.  
  */

  


  template<int dim>
  class heat_eqn_T
  {
    
    public:
        
    void T (model<dim>& T_model, std::map<int, double>& gid_T)
    {
      std::vector<double> my_T;
      std::vector<int> my_nd;

      for(const auto& nd : T_model.mesh().bd_nodes())
      {
        if (nd->flag() == 1)
        {
          my_T.push_back(T_model.state().get_dof(nd->first_dof_lid()));
          my_nd.push_back(nd->gid());
        }
      }

      MPI_Barrier(MPI_COMM_WORLD);
      int rank = get_rank();

      const int send_buffer_size = my_nd.size();      
      int sum;
      std::vector<int> num_el_buffer, u;   
      MPI_NumElbuffer_U_Sum(send_buffer_size, num_el_buffer, u, sum);   
      
      /// "num_el_buffer" is the vector of number of T dofs on each pid.
      std::vector<int> all_nd(sum);
      std::vector<double> all_T(sum);

      MPI_Allgatherv(&my_nd[0], num_el_buffer[rank], MPI_INT, &all_nd[0], &num_el_buffer[0], &u[0], MPI_INT, MPI_COMM_WORLD);
      MPI_Allgatherv(&my_T[0], num_el_buffer[rank], MPI_DOUBLE, &all_T[0], &num_el_buffer[0], &u[0], MPI_DOUBLE, MPI_COMM_WORLD);
      /// The above MPI calls collect nd_gid and T on each pid in form of all_nd and all_T, respectively.
      

      /// Some nodes may be shared among different pids; however map stores only unique keys. 
      for(std::size_t i=0; i<all_nd.size(); i++)
        gid_T[all_nd[i]] = all_T[i];
    }

  };


}

#endif
