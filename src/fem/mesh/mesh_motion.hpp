#ifndef _GALES_MESH_MOTION_HPP_
#define _GALES_MESH_MOTION_HPP_



#include "mesh.hpp"



namespace GALES {


  
  template<int dim>  class Mesh_Motion {};




  template<>
  class Mesh_Motion<2>
  {
    using mesh_type = Mesh<2>;
    using vec = std::vector<double>;
    
    public:
           
    
    // This function reads the updated node coordinates of the mesh if mesh is moving and we need to restart at some time instance.  
    // Note that mesh is written in binary. 
    // We use read_bin function to read the file into a vector which is consisted of first x-coord of all nodes; then y-coord of all nodes.

    void read_updated_node_coords(const std::string mesh_file, mesh_type& mesh)
    {
      std::vector<double> v(2*mesh.tot_nodes());
      read_bin(mesh_file,v);

      std::vector<double> nd_x(mesh.tot_nodes()), nd_y(mesh.tot_nodes());
      for(int i=0; i<mesh.tot_nodes(); i++)
      {
        nd_x[i] = v[i];
        nd_y[i] = v[mesh.tot_nodes()+i];
      }

      for(const auto& nd : mesh.nodes())
      {
        nd->set_x(nd_x[nd->gid()]);
        nd->set_y(nd_y[nd->gid()]);
      }    
    }
    
    
    
    
    
    
    // This function updates the mesh coordinates with deformation computed with elastostatic model
    // Each process has a unique list of nodes (host and ghost) and is responsible to
    // update all its nodes including the ghost ones also.
    // If we do not update the ghost node then it will be a problem as some
    // other process on which the same node is the host one will update the node
    // and it will cause inconsistancy.            
        
    template<typename model_type>
    void update_mesh(model_type& e_model, mesh_type& mesh)
    { 
      for (auto& nd : mesh.nodes())
      {
          std::vector<double> dofs_u;
          e_model.extract_node_dofs(*(e_model.mesh().nodes()[nd->lid()]), dofs_u);      
          nd->update_coord(dofs_u);  
      } 
    }
    
    
    
    
    
  
    // This function writes mesh in parallel
    // Mesh is written in binary. We write only the mesh nodes considering that in our implementation 
    // element topology as well as connections do not change upon mesh updation.
    // Only the nodes deform their position. We write first x-coord of all nodes; then y-coord of all nodes; then z-coord of all nodes.
   
    void write_mesh(const IO& parallel_mesh_writer, mesh_type& mesh)
    {
      std::vector<vec> nodes_coord(2, vec());

      for(const auto& nd : mesh.nodes())                /// this collects nodes_coord
      {
        nodes_coord[0].push_back(nd->get_x());
        nodes_coord[1].push_back(nd->get_y());
      }

      parallel_mesh_writer.write2("nd_x", nodes_coord[0]);      /// this writes nodes_coord[0] in file "nd_x" in parallel using map_ 
      parallel_mesh_writer.write2("nd_y", nodes_coord[1]);

      MPI_Barrier(MPI_COMM_WORLD);
      const int rank = get_rank();
      if (rank == 0)
      {
        std::stringstream o_file;
        o_file <<"results/fluid_mesh/"<<time::get().t()<<std::flush;
        std::string s = o_file.str();
       
        std::ofstream of_a(s.c_str(), std::ios_base::binary);
        std::ifstream if_a("nd_x", std::ios_base::binary);
        std::ifstream if_b("nd_y", std::ios_base::binary);

        of_a << if_a.rdbuf() << if_b.rdbuf();            /// this concate the random binary buffer in of_a

        if_a.close();
        if_b.close();
        of_a.close();
        
        remove("nd_x");
        remove("nd_y");
      }
    }

  }; 












  template<>
  class Mesh_Motion<3>
  {
    using mesh_type = Mesh<3>;
    using vec = std::vector<double>;

    public:
           
    
    // This class reads the updated node coordinates of the mesh if mesh is moving and we need to restart at some time instance.  
    // Note that mesh is written in binary. 
    // We use read_bin function to read the file into a vector which is consisted of first x-coord of all nodes; then y-coord of all nodes; then z-coord of all nodes.

    void read_updated_node_coords(const std::string mesh_file, const mesh_type& mesh)
    {
      std::vector<double> v(3*mesh.tot_nodes());
      read_bin(mesh_file,v);

      std::vector<double> nd_x(mesh.tot_nodes()), nd_y(mesh.tot_nodes()), nd_z(mesh.tot_nodes());
      for(int i=0; i<mesh.tot_nodes(); i++)
      {
        nd_x[i] = v[i];
        nd_y[i] = v[mesh.tot_nodes()+i];
        nd_z[i] = v[2*mesh.tot_nodes()+i];
      }

      for(const auto& nd : mesh.nodes())
      {
        nd->set_x(nd_x[nd->gid()]);
        nd->set_y(nd_y[nd->gid()]);
        nd->set_z(nd_z[nd->gid()]);
      }
    }






    // This function updates the mesh coordinates with deformation computed with elastostatic model
    // Each process has a unique list of nodes (host and ghost) and is responsible to
    // update all its nodes including the ghost ones also.
    // If we do not update the ghost node then it will be a problem as some
    // other process on which the same node is the host one will update the node
    // and it will cause inconsistancy.            
        
    template<typename model_type>
    void update_mesh(model_type& e_model, mesh_type& mesh)
    { 
      for (auto& nd : mesh.nodes())
      {
          std::vector<double> dofs_u;
          e_model.extract_node_dofs(*(e_model.mesh().nodes()[nd->lid()]), dofs_u);      
          nd->update_coord(dofs_u);  
      } 
    }






    // This function writes mesh in parallel
    // Mesh is written in binary. We write only the mesh nodes considering that in our implementation 
    // element topology as well as connections do not change upon mesh updation.
    // Only the nodes deform their position. We write first x-coord of all nodes; then y-coord of all nodes; then z-coord of all nodes.

    void write_mesh(const IO& parallel_mesh_writer, mesh_type& mesh)
    {
      std::vector<vec> nodes_coord(3, vec());

      for (const auto& nd : mesh.nodes())
      {
        nodes_coord[0].push_back(nd->get_x());
        nodes_coord[1].push_back(nd->get_y());
        nodes_coord[2].push_back(nd->get_z());
      }

      parallel_mesh_writer.write2("nd_x", nodes_coord[0]);
      parallel_mesh_writer.write2("nd_y", nodes_coord[1]);
      parallel_mesh_writer.write2("nd_z", nodes_coord[2]);


      MPI_Barrier(MPI_COMM_WORLD);
      const int rank = get_rank();
      if (rank == 0)
      {
        std::stringstream o_file;
        o_file <<"results/fluid_mesh/"<<time::get().t()<<std::flush;
        std::string s = o_file.str();
        
        std::ofstream of_a(s.c_str(), std::ios_base::binary);
        std::ifstream if_a("nd_x", std::ios_base::binary);
        std::ifstream if_b("nd_y", std::ios_base::binary);
        std::ifstream if_c("nd_z", std::ios_base::binary);

        of_a << if_a.rdbuf() << if_b.rdbuf()<<if_c.rdbuf();

        if_a.close();
        if_b.close();
        if_c.close();
        of_a.close();

        remove("nd_x");
        remove("nd_y");
        remove("nd_z");
      }
    }

  }; 





}//namespace GALES
#endif
