#ifndef FSI_RESIDUAL_CHECK_HPP
#define FSI_RESIDUAL_CHECK_HPP


#include <vector>
#include "Epetra_FEVector.h"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_StandardCatchMacros.hpp"
#include "Teuchos_RCPDecl.hpp"
#include "BelosConfigDefs.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosEpetraAdapter.hpp"




namespace GALES{


  /**
       This class checks the solid displacement residual at the fluid-solid interface and returns a bool indicating if
       convergence criteria is reached or not.  
  */



  class fsi_residual_check
  {
     public:
     
      template<typename model_type>
      bool execute(model_type& u_model)
      {
        std::vector<double> my_s_u_nrm;
        std::vector<int> my_nd;
        for(const auto& nd : u_model.mesh().bd_nodes())       /// solid boundary mesh 
        {
          if (nd->flag() == 1)
          {
            std::vector<double> s_u;                                                       
            u_model.extract_node_dofs(*nd, s_u);      /// s_u has dofs of current time
            my_nd.push_back(nd->gid());
            
            double u_nrm = 0.0;
            for(double x : s_u) {u_nrm += x*x;}
            u_nrm = sqrt(u_nrm);                                /// this u_nrm is for each fsi node on each process
                    
            my_s_u_nrm.push_back(u_nrm);
          }
        }
        
        MPI_Barrier(MPI_COMM_WORLD);
        int rank = get_rank();

        const int send_buffer_size = my_nd.size();      
        int sum;
        std::vector<int> num_el_buffer, u;   
        MPI_NumElbuffer_U_Sum(send_buffer_size, num_el_buffer, u, sum);   
        
        std::vector<int> all_nd(sum);
        std::vector<double> all_u_nrm(sum);

        /// Here we combine gid and u_nrm for each fsi node on all processes.
        MPI_Allgatherv(&my_nd[0], num_el_buffer[rank], MPI_INT, &all_nd[0], &num_el_buffer[0], &u[0], MPI_INT, MPI_COMM_WORLD);
        MPI_Allgatherv(&my_s_u_nrm[0], num_el_buffer[rank], MPI_DOUBLE, &all_u_nrm[0], &num_el_buffer[0], &u[0], MPI_DOUBLE, MPI_COMM_WORLD);
        /// The above MPI calls collect nd_gid and s_u_nrm on each pid in form of all_nd and all_u_nrm, respectively.
           
        
        /// some nodes may be shared among different pids so we define a map which stores unique key(node) value(u_nrm) pairs.        
        std::map<int, double> nd_s_nrm_map;
        for(std::size_t i=0; i<all_nd.size(); i++)
           nd_s_nrm_map[all_nd[i]] = all_u_nrm[i];
        
        /// Here we compute the "u_nrm" for all interface nodes, including those which lie on other processes.   
        double u_nrm = 0.0;
        for(auto x : nd_s_nrm_map) {u_nrm += x.second*x.second;}
        u_nrm = sqrt(u_nrm);
           
        
        /**
             criterion:
                          ||u||_i - ||u||_{i-1}
                         ------------------------        <  1.e-5
                         sqrt(num of enteries in u) 
        */
                   
        if((u_nrm - s_u_nrm_prev_it_)/sqrt(nd_s_nrm_map.size()) < 1.e-5)       
        {
           return true;
        }
        else
        {
           s_u_nrm_prev_it_ = u_nrm;
           return false;
        }
        return false;
      }

    private:    
      double s_u_nrm_prev_it_ = 0.0;        
  };      


}

#endif
