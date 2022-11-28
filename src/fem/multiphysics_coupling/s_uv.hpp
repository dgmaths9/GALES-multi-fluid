#ifndef S_UV_HPP
#define S_UV_HPP



#include <vector>
#include <map>




namespace GALES{



  /**
       This class computes the solid displacement and velocity at fluid-solid interface to be passed to fluid solver for FSI coupling.  
  */

  


  template<int dim>  
  class s_uv{};



  template<>
  class s_uv<2>
  {
    
    public:
    
    
    void s_u
    (
      model<2>& u_model,
      std::map<int, std::vector<double> >& s_gid_u
    )
    {
      std::vector<double> my_ux, my_uy;
      std::vector<int> my_nd;

      for(const auto& nd : u_model.mesh().bd_nodes())                 ///solid boundary mesh 
      {
        if (nd->flag() == 1)
        {
          std::vector<std::vector<double>> s_u;                                   /// s_u has dofs of current and previos time 
          u_model.extract_node_dofs(*nd, s_u);        
          my_ux.push_back(s_u[0][0]-s_u[1][0]);                          /// s_u[0][0] = s_ux(current time);   s_u[1][0] = s_ux(previous time)
          my_uy.push_back(s_u[0][1]-s_u[1][1]);                          /// s_u[0][1] = s_uy(current time);   s_u[1][1] = s_uy(previous time)
          my_nd.push_back(nd->gid());
        }
      }

      MPI_Barrier(MPI_COMM_WORLD);
      int rank = get_rank();

      const int send_buffer_size = my_nd.size();      
      int sum;
      std::vector<int> num_el_buffer, u;   
      MPI_NumElbuffer_U_Sum(send_buffer_size, num_el_buffer, u, sum);   
      
      /// "num_el_buffer" is the vector of number of my_ux dofs on each pid.
      std::vector<int> all_nd(sum);
      std::vector<double> all_ux(sum), all_uy(sum);

      MPI_Allgatherv(&my_nd[0], num_el_buffer[rank], MPI_INT, &all_nd[0], &num_el_buffer[0], &u[0], MPI_INT, MPI_COMM_WORLD);
      MPI_Allgatherv(&my_ux[0], num_el_buffer[rank], MPI_DOUBLE, &all_ux[0], &num_el_buffer[0], &u[0], MPI_DOUBLE, MPI_COMM_WORLD);
      MPI_Allgatherv(&my_uy[0], num_el_buffer[rank], MPI_DOUBLE, &all_uy[0], &num_el_buffer[0], &u[0], MPI_DOUBLE, MPI_COMM_WORLD);
      /// The above MPI calls collect nd_gid, s_ux and s_uy on each pid in form of all_nd, all_ux and all_uy, respectively.
      

      /// Some nodes may be shared among different pids; however map stores only unique keys. 
      for(std::size_t i=0; i<all_nd.size(); i++)
        s_gid_u[all_nd[i]] = std::vector<double>{all_ux[i], all_uy[i]};
    }





    void s_v
    (
      model<2>& v_model,
      std::map<int, std::vector<double> >& s_gid_v
    )
    {
      std::vector<double> my_vx, my_vy;
      std::vector<int> my_nd;

      for(const auto& nd : v_model.mesh().bd_nodes())
      {
        if (nd->flag() == 1)
        {
          std::vector<double> s_v;                                               ///we get s_v dofs of only current time
          v_model.extract_node_dofs(*nd, s_v);        
          my_vx.push_back(s_v[0]);                                           /// s_v[0] = s_vx(current time);
          my_vy.push_back(s_v[1]);                                           /// s_v[1] = s_vy(current time);
          my_nd.push_back(nd->gid());
        }
      }

      MPI_Barrier(MPI_COMM_WORLD);
      int rank = get_rank();

      const int send_buffer_size = my_nd.size();      
      int sum;
      std::vector<int> num_el_buffer, u;   
      MPI_NumElbuffer_U_Sum(send_buffer_size, num_el_buffer, u, sum);   

      /// "num_el_buffer" is the vector of number of my_vx dofs on each pid.
      std::vector<int> all_nd(sum);
      std::vector<double> all_vx(sum), all_vy(sum);

      MPI_Allgatherv(&my_nd[0], num_el_buffer[rank], MPI_INT, &all_nd[0], &num_el_buffer[0], &u[0], MPI_INT, MPI_COMM_WORLD);
      MPI_Allgatherv(&my_vx[0], num_el_buffer[rank], MPI_DOUBLE, &all_vx[0], &num_el_buffer[0], &u[0], MPI_DOUBLE, MPI_COMM_WORLD);
      MPI_Allgatherv(&my_vy[0], num_el_buffer[rank], MPI_DOUBLE, &all_vy[0], &num_el_buffer[0], &u[0], MPI_DOUBLE, MPI_COMM_WORLD);
      /// The above MPI calls collect nd_gid, s_vx and s_vy on each pid in form of all_nd, all_vx and all_vy, respectively.

      /// Some nodes may be shared among different pids; however map stores only unique keys. 
      for(std::size_t i=0; i<all_nd.size(); i++)
        s_gid_v[all_nd[i]] = std::vector<double>{all_vx[i], all_vy[i]};
    }


  };













  template<>
  class s_uv<3>
  {

    public:
    
    
    void s_u
    (
      model<3>& u_model,
      std::map<int, std::vector<double> >& s_gid_u
    )
    {
      std::vector<double> my_ux, my_uy, my_uz;
      std::vector<int> my_nd;

      for(const auto& nd : u_model.mesh().bd_nodes())
      {
        if (nd->flag() == 1)
        {
          std::vector<std::vector<double>> s_u;
          u_model.extract_node_dofs(*nd, s_u);        
          my_ux.push_back(s_u[0][0]-s_u[1][0]);
          my_uy.push_back(s_u[0][1]-s_u[1][1]);
          my_uz.push_back(s_u[0][2]-s_u[1][2]);
          my_nd.push_back(nd->gid());
        }
      }

      MPI_Barrier(MPI_COMM_WORLD);
      int rank = get_rank();

      const int send_buffer_size = my_nd.size();      
      int sum;
      std::vector<int> num_el_buffer, u;   
      MPI_NumElbuffer_U_Sum(send_buffer_size, num_el_buffer, u, sum);   

      std::vector<int> all_nd(sum);
      std::vector<double> all_ux(sum), all_uy(sum), all_uz(sum);

      MPI_Allgatherv(&my_nd[0], num_el_buffer[rank], MPI_INT, &all_nd[0], &num_el_buffer[0], &u[0], MPI_INT, MPI_COMM_WORLD);
      MPI_Allgatherv(&my_ux[0], num_el_buffer[rank], MPI_DOUBLE, &all_ux[0], &num_el_buffer[0], &u[0], MPI_DOUBLE, MPI_COMM_WORLD);
      MPI_Allgatherv(&my_uy[0], num_el_buffer[rank], MPI_DOUBLE, &all_uy[0], &num_el_buffer[0], &u[0], MPI_DOUBLE, MPI_COMM_WORLD);
      MPI_Allgatherv(&my_uz[0], num_el_buffer[rank], MPI_DOUBLE, &all_uz[0], &num_el_buffer[0], &u[0], MPI_DOUBLE, MPI_COMM_WORLD);

      for(std::size_t i=0; i<all_nd.size(); i++)
        s_gid_u[all_nd[i]] = std::vector<double>{all_ux[i], all_uy[i], all_uz[i]};
    }





    void s_v
    (
      model<3>& v_model,
      std::map<int, std::vector<double> >& s_gid_v
    )
    {
      std::vector<double> my_vx, my_vy, my_vz;
      std::vector<int> my_nd;

      for(const auto& nd : v_model.mesh().bd_nodes())
      {
        if (nd->flag() == 1)
        {
          std::vector<double> s_v;
          v_model.extract_node_dofs(*nd, s_v);        
          my_vx.push_back(s_v[0]);
          my_vy.push_back(s_v[1]);
          my_vz.push_back(s_v[2]);
          my_nd.push_back(nd->gid());
        }
      }

      MPI_Barrier(MPI_COMM_WORLD);
      int rank = get_rank();

      const int send_buffer_size = my_nd.size();      
      int sum;
      std::vector<int> num_el_buffer, u;   
      MPI_NumElbuffer_U_Sum(send_buffer_size, num_el_buffer, u, sum);   

      std::vector<int> all_nd(sum);
      std::vector<double> all_vx(sum), all_vy(sum), all_vz(sum);

      MPI_Allgatherv(&my_nd[0], num_el_buffer[rank], MPI_INT, &all_nd[0], &num_el_buffer[0], &u[0], MPI_INT, MPI_COMM_WORLD);
      MPI_Allgatherv(&my_vx[0], num_el_buffer[rank], MPI_DOUBLE, &all_vx[0], &num_el_buffer[0], &u[0], MPI_DOUBLE, MPI_COMM_WORLD);
      MPI_Allgatherv(&my_vy[0], num_el_buffer[rank], MPI_DOUBLE, &all_vy[0], &num_el_buffer[0], &u[0], MPI_DOUBLE, MPI_COMM_WORLD);
      MPI_Allgatherv(&my_vz[0], num_el_buffer[rank], MPI_DOUBLE, &all_vz[0], &num_el_buffer[0], &u[0], MPI_DOUBLE, MPI_COMM_WORLD);

      for(std::size_t i=0; i<all_nd.size(); i++)
        s_gid_v[all_nd[i]] = std::vector<double>{all_vx[i], all_vy[i], all_vz[i]};
    }


  };



}

#endif
