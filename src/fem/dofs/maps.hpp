#ifndef GALES_MAPS_HPP
#define GALES_MAPS_HPP



#include <vector>
#include <algorithm>
#include <memory>
#include "Epetra_Map.h"
#include "Epetra_MpiComm.h"





namespace GALES{



  /**
     This class defines Epetra_Maps used in GALES. Specifically, we define the following three types of maps:
     
     state_map_ ---->  this is based on shared_nd_dof_list, see the function for more details
     dof_map_  ---->  this is based on unshared_nd_dof_list, see the function for more details
     shared_node_map_  ----> this is based on shared_nd_list, see the function for more details
  */
  
  
  class epetra_maps
  {
      public:
       
       template<typename mesh_type>
       explicit epetra_maps(const mesh_type& mesh)
       {
         auto gids = shared_nd_dof_list(mesh);             
         state_map_ = std::make_shared<Epetra_Map>(-1, gids.size(), &(gids[0]), 0, Epetra_MpiComm(MPI_COMM_WORLD));
         set_nd_first_dof_lid(mesh);


         gids = unshared_nd_dof_list(mesh);             
         dof_map_ = std::make_shared<Epetra_Map>(-1, gids.size(), &(gids[0]), 0, Epetra_MpiComm(MPI_COMM_WORLD));     
         
         
         gids = shared_nd_list(mesh);             
         shared_node_map_ = std::make_shared<Epetra_Map>(-1, gids.size(), &(gids[0]), 0, Epetra_MpiComm(MPI_COMM_WORLD));            
       }  
  

       //-------------------------------------------------------------------------------------------------------------------------------------        
       /// Deleting the copy and move constructors - no duplication/transfer in anyway
       epetra_maps(const epetra_maps&) = delete;               //copy constructor
       epetra_maps& operator=(const epetra_maps&) = delete;    //copy assignment operator
       epetra_maps(epetra_maps&&) = delete;                    //move constructor  
       epetra_maps& operator=(epetra_maps&&) = delete;         //move assignment operator 
       //-------------------------------------------------------------------------------------------------------------------------------------


  
       auto state_map()const {return state_map_;}
       auto dof_map()const {return dof_map_;}
       auto shared_node_map()const {return shared_node_map_;}
       
       
       
       
       private:


       /**
           In our implementation, each process reads and computes its own elements. 
           Then we collect nodes belonging to the elements on each process.
           Some nodes (dofs) are shared among inter process boundaries.
           We call these shared nodes (dofs).
           For each node we set its owning process by calling 'set_ownership()' function in mesh_reader.                      
           If a shared node (dof) lies on a process then it is a host node(dof), otherwise it is a ghost node(dof).
           This function collects gids of dofs of nodes (including host and ghost) on each process.
           Thats why we call this  shared_nd_dof_list.
           This is used to create state_map which carries our dofs system throughout the code.
           Note that the dof_state built on the state_map contains both host and ghost dofs.
           By carrying shared dofs on each process we avoid communication when we need to extract dofs of the ghost nodes.
       */
       template<typename mesh_type>
       auto shared_nd_dof_list(const mesh_type& mesh) const
       {
          std::vector<int> gid_list;
          gid_list.reserve( (mesh.nodes().size()) * (mesh.nodes()[0]->nb_dofs()) );
          const int rank = get_rank();
      
          for(const auto& nd : mesh.nodes())
            for(int j=0; j<nd->nb_dofs(); j++)
             gid_list.push_back( nd->first_dof_gid() + j );
             
         sort_unique(gid_list);
         return gid_list;
       }




       /**
          This function sets the lid of the first dof for each node (including host and ghost).                    
          Once we know first_dof_lid of each node, we can easily access the other dofs of node.
       */   
       template<typename mesh_type>
       void set_nd_first_dof_lid(const mesh_type& mesh) const
       {
         for(const auto& nd : mesh.nodes())
         {
           const int fdgid(nd->first_dof_gid());  
           const int fdlid(state_map_->LID(fdgid));
           nd->first_dof_lid(fdlid);
         }
       }




         
       /** 
           This function collects gids of dofs of owned nodes (only host nodes) on each process.
           Although nodes are shared on inter process boundaries, 
           while reading the mesh we identify the host and ghost nodes on each process by setting pid of each mesh node.             
           This is used to create the dof_map which is used to create linear solver and dofMap_to_stateMap_importer in IO class. 
       */ 
       template<typename mesh_type>
       auto unshared_nd_dof_list(const mesh_type& mesh) const
       {
          std::vector<int> gid_list;
          gid_list.reserve( (mesh.nodes().size()) * (mesh.nodes()[0]->nb_dofs()) );
          const int rank = get_rank();
      
          for(const auto& nd : mesh.nodes())
           if(nd->pid()==rank) 
            for(int j=0; j<nd->nb_dofs(); j++)
             gid_list.push_back( nd->first_dof_gid() + j );

         sort_unique(gid_list);
         return gid_list;
       }





       /**
           This function collects gids of nodes on each process.
           "gid_list" contains gids on nodes, some of them are shared on other processes.
           This is used to create shared_node_map which is used for parallel printing of mesh nodes and secondary dofs.
       */
       template<typename mesh_type>
       auto shared_nd_list(const mesh_type& mesh) const 
       {
         std::vector<int> gid_list;
         
         for(const auto& nd : mesh.nodes())
           gid_list.push_back(nd->gid());

         sort_unique(gid_list);
         return gid_list;
       }


       
       
       std::shared_ptr<Epetra_Map> state_map_;
       std::shared_ptr<Epetra_Map> dof_map_;
       std::shared_ptr<Epetra_Map> shared_node_map_;    
  };




}



#endif
