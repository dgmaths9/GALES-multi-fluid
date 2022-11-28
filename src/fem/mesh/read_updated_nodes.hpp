#ifndef READ_UPDATED_NODES_HPP
#define READ_UPDATED_NODES_HPP



#include <vector>
#include "../base/binary_io_serial.hpp"



namespace GALES{


  /**
     This class reads the updated node coordinates of the mesh if mesh is moving and we need to restart at some time instance.  
     Note that mesh is written in binary. We use read_bin function to read the file into a vector which is consisted of first x-coord of all nodes; then y-coord of all nodes; then z-coord of all nodes.
  */


  template<int dim> class read_updated_node_coords {};

  
  
  template<>
  class read_updated_node_coords<2>
  {
    using mesh_type = Mesh<2>;
   
    public:
    
    void execute(std::string mesh_file, const mesh_type& mesh)
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
  };







  template<>
  class read_updated_node_coords<3>
  {
    using mesh_type = Mesh<3>;

    public:

    void execute(std::string mesh_file, const mesh_type& mesh)
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
  };



}
#endif
