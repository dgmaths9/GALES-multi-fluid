#ifndef _GALES_MESH_HPP_
#define _GALES_MESH_HPP_

#include <functional>
#include <vector>
#include <limits>


#include "element.hpp"
#include "node.hpp"
#include "quad.hpp"



namespace GALES {


  /**
      This class defines the mesh.

    The mesh file we read already contains the mesh partition information.
   
    
   
    el          0    1    2    3                   nd          0    1    2    3    4
    PID         0    0    1    1                   PID         0    0    1    1    1
             o----o---)o(---o----o                            o----o---)o(---o----o
    el GID      0    1    2    3                   nd GID      0    1    2    3    4
    el LID(0)   0    1                             nd LID(0)   0    1    2
    el LID(1)             0    1                   nd LID(1)             0    1    2   
    
    
    Each process reads and stores only those elements which are owned by it.      
    An element has a unique GID and a unique PID(disjoint partition). Element LID is defined as the unique local counter number on each process. 
    Therefore two different elements lying on two different processes can have the same LID.
    
    Similarly each process stores and reads the nodes owned by it. 
    A node has a unique GID. some nodes lie on inter process boundaries (e.g. see node 2 above).
    In this case, we define the owner PID of the node based on the largest integer number of the shared PIDs.
    Node LID is defined as the unique local counter number on each process.  Two different nodes lying on two different processes can have the same LID.
    
    Since we do most of the operations through LIDs, each process is responsible for executing the operation on its nodes and elements. 
    This applies even on the shared nodes. For example to do mesh update nodes with GIDs 0,1,2 will be updated by PID 0 and nodes with GIDs 2,3,4 will be
    updated by PID 1. Note that Node 2 is stored as two different objects on respective PIDs:  LID 2 on PID 0 and LID 0 on PID 1.  
    So, it will not be updated twice.
  */



  template<int dim>
  class Mesh 
  {
        
    public:
                
        using point_t = point<dim>;
        using node_t = node<dim>;
        using element_t = element<dim>;

                
        Mesh(std::string mesh_name, int nd_nb_dofs, const read_setup& setup) 
        { 
           mesh_reader(mesh_name, nd_nb_dofs); 
        }
        
        
        //-------------------------------------------------------------------------------------------------------------------------------------        
        /// Deleting the copy and move constructors - no duplication/transfer in anyway
        Mesh(const Mesh&) = delete;               //copy constructor
        Mesh& operator=(const Mesh&) = delete;    //copy assignment operator
        Mesh(Mesh&&) = delete;                    //move constructor  
        Mesh& operator=(Mesh&&) = delete;         //move assignment operator 
        //-------------------------------------------------------------------------------------------------------------------------------------
        
        
        
        
        
        
        //------------------------------------------------------------------------------------------------------------------
        void mesh_reader(std::string& mesh_name, int nd_nb_dofs)
        {
          std::ifstream infile(mesh_name);                              // opening the mesh file to read
          if(!infile.is_open())
          {
             Error("unable to open and read mesh file");
          }
      
    
          read_header(infile);
    
    
          double el_read_start = MPI_Wtime();
          std::vector<int> my_nodes;
          read_elements(infile, my_nodes);
          print_only_pid<0>(std::cerr)<<"El reading took: "<< MPI_Wtime()-el_read_start<<" s\n";      
    
    
    
    
          double side_read_start = MPI_Wtime();
          read_sides(infile);
          print_only_pid<0>(std::cerr)<<"Side reading took: "<< MPI_Wtime()-side_read_start<<" s\n";      
    
    
          
          infile.clear();
          infile.seekg(std::ios_base::beg);           ///rewind the file to read nodes
    
    
    
    
          double nd_read_start = MPI_Wtime();
          read_nodes(infile, my_nodes, nd_nb_dofs);
          print_only_pid<0>(std::cerr)<<"Nd reading took: "<< MPI_Wtime()-nd_read_start<<" s\n";
    
    
          infile.close();                                           // closing the mesh file after reading
    
    
          /// With this function each element has el_node_vec_ of its connected nodes; nodes know their lid and gid on respective pid
          double set_element_node_links_start = MPI_Wtime();
          set_element_node_links(my_nodes);
          print_only_pid<0>(std::cerr)<<"seting Element node links took: "<< MPI_Wtime()-set_element_node_links_start<<" s\n";
    
    
          /// This functions reads the info about side_nodes, side_flag etc for each element 
          double el_side_info_start = MPI_Wtime();
          el_side_info();
          print_only_pid<0>(std::cerr)<<"El side info setting took: "<< MPI_Wtime()-el_side_info_start<<" s\n";
            


          /// This functions computes the data at quadrature points for each element 
          double el_quad_data_start = MPI_Wtime();
          compute_el_quad_data();
          print_only_pid<0>(std::cerr)<<"El quad data setting took: "<< MPI_Wtime()-el_quad_data_start<<" s\n";
        }
        //------------------------------------------------------------------------------------------------------------------
            
            
            
            
            
    
        //-----------------------------------------------------------------------------------------------------------------
        void read_header(std::ifstream& infile)
        {             
          std::vector<std::string> result;
          
          read_one_line( infile, result);   ///MESH! 2D            or         MESH! 3D
    
          read_one_line( infile, result);   ///nodes  <n>   
          tot_nodes(std::stoi(result[1]));    /// this sets total number of nodes of mesh
          print_only_pid<0>(std::cerr)<< "nodes:  "<< tot_nodes() << "\n" ;
    
          read_one_line( infile, result);  ///elements  <n>
          tot_elements(std::stoi(result[1]));  /// this sets total number of elements of mesh
          print_only_pid<0>(std::cerr)<< "elements:  "<<tot_elements() << "\n";
          
          read_one_line( infile, result);   ///sides  <n>
          tot_sides(std::stoi(result[1]));  /// this sets total number of bd sides of mesh
          print_only_pid<0>(std::cerr)<< "sides:  "<<tot_sides() << "\n";
    
          for(int i=0; i<dim; i++)  ///X 0.0 1.0;      Y 0.0 1.0       Z 
          read_one_line( infile, result);   /// This reads  mesh dimension
        }
        //-----------------------------------------------------------------------------------------------------------------
        
        
    
    
    
    
    
    
    
    
    
        //-----------------------------------------------------------------------------------------------------------------
        /// Reads the elements owned by the processor.
        /// Elements are stored with complete id, gid node list and boundary information.
        /// Each element has gids of the connected nodes in the form of "node_gid_vec_".    
        void read_elements(std::ifstream& infile, std::vector<int>& my_nodes)
        {
          const int rank = get_rank();
          const int size = get_size();
          const int avg_num_els_per_pid = (int) (tot_elements()/size);      
          const int avg_num_nds_per_pid = (int) (tot_nodes()/size);      
          elements().reserve(avg_num_els_per_pid);
          my_nodes.reserve(avg_num_nds_per_pid);
          
          int count(0);
          std::vector<std::string> result;
    
          for(int i=0; i<tot_nodes(); i++)     /// skip all lines starting with Node
             skip_one_line(infile);
          
          for(int i=0; i<tot_elements(); i++)  /// start reading elements
          {
            read_one_line(infile, result);          ///for quad element:     Element 0 0 4 3 13 16 12 1
            const int pid = std::stoi(result[2]);
                      
            ///  Only elements on the current rank are read.
            if(pid == rank)
            {
              const int gid = std::stoi(result[1]); 
              const int num_nodes = std::stoi(result[3]); 
      
              auto el =  std::make_shared<element<dim>>(num_nodes); 
              el->pid(pid);
              el->gid(gid);
              el->lid(count);
                
              std::vector<int> node_gid_vec(num_nodes);          
              for(int j=0; j<num_nodes; j++)
              {
                const int nd_gid = std::stoi(result[4+j]);
                node_gid_vec[j] = nd_gid;            
                my_nodes.push_back(nd_gid);
              }
              el->node_gid_vec(std::move(node_gid_vec));            
      
              const int flag = std::stoi(result[4+num_nodes] );         
              el->on_boundary(flag>0);
              elements().push_back(std::move(el));
      
              count++;
            }
          }
          
       
          /// Fill "bd_elements_" container  
          for(const auto& el : elements())
            if(el->on_boundary())
              bd_elements().push_back(el);    
          
    
          sort_unique(my_nodes);
          /// my_nodes contains a sorted and non repeating list of gids of nodes defined on each pid.
          /// some node gids are shared among multiple pids due to sharing of same node on multiple pids.
        }
        //-----------------------------------------------------------------------------------------------------------------
    
    
    
        
        
        
        
       //-----------------------------------------------------------------------------------------------------------------
       /// Reads the sides owned by the process.
       void read_sides(std::ifstream& infile)
       {
          std::vector<std::string> result;
          const int rank = get_rank();
    
          for(int i=0; i<tot_sides(); i++)        /// start reading sides
          {
            read_one_line(infile, result);         ///for Side of quad and tri:      Side 0 2 2 0 4 5
            const int pid = std::stoi(result[2]);
                      
            ///  Only sides of elements on current rank are read.
            if(pid == rank)
            { 
               const int num_side_nds_ = std::stoi(result[3]);

               std::vector<int> side(num_side_nds_);
               for(int j=0; j<num_side_nds_; j++) 
                 side[j] = std::stoi(result[4+j]);

               sides().push_back(std::move(side));   
                  
               side_flag().push_back(std::stoi(result[4+num_side_nds_]));
            }
          }                 
        }    
        //-----------------------------------------------------------------------------------------------------------------
                      
    
    
        
        
        
        
        
        
    
        //-----------------------------------------------------------------------------------------------------------------
        /// Reads the nodes owned by each processor.
        /// reading nodes: complete id, coordinates;  Each processor reads its own nodes defined in my_nodes
        /// Each node of the mesh on each processor knows its lid and gid
        
        void read_nodes(std::ifstream& infile, const std::vector<int>& my_nodes, int nd_nb_dofs)
        {
          nodes().resize(my_nodes.size());
          int count(0);
          std::vector<std::string> result;
          auto it(my_nodes.begin());
          
          while(it != my_nodes.end()) 
          {
            read_one_line(infile, result);
            if(result[0]!="Node") continue;
            const int gid(stoi(result[1]));
           
            if(*it==gid)
            {
              auto nd = std::make_shared<node<dim>>();        /// nd_pid is set below
              nd->gid(gid);
              nd->lid(count);                          ///count is node_lid on each pid
              nd->nb_dofs(nd_nb_dofs);             
           
              nd->set_x(std::stod(result[2]));
              nd->set_y(std::stod(result[3]));
              if(dim == 3)
                nd->set_z(std::stod(result[4]));

              const int flag(std::stoi(result[dim+2]));
              nd->on_boundary(flag>0);          
              nd->flag(flag);
                                   
              nodes()[count] = std::move(nd);
    
              it++;
              count++;
            }        
          }
          
    
          /// Fill "bd_nodes_" container  
          /// we should reserve the length of   bd_nodes()  vector by some good estimate!!!!
          for(const auto& nd : nodes())
            if(nd->on_boundary())
              bd_nodes().push_back(nd);         
              
          
          
          /// compute number of fsi nodes    
          count = 0;
          for(const auto& nd : bd_nodes())
           if (nd->flag() == 1)
            count++;
          num_fsi_nodes(count);
                 
    
          /// Next we set the ownership of each node i.e. to which process it belongs.
          /// Owner is set according to the highest pid number.
          /// If node 11 is shared on pid 3 and 4 then we set owner = 4
          /// We should set it according to the nd pid file generated by Metis
          std::vector<int> nd_gids(nodes().size());
          count = 0;
          for(const auto& nd : nodes())
          {
              nd_gids[count] = nd->gid();
              count++;
          }
          
          const int n_nds = tot_nodes();
          const int rank = get_rank();
          std::vector<int> ownership(n_nds,-1);
          for(auto i : nd_gids)
            ownership[i] = rank;
      
          std::vector<int> ownership_tmp(n_nds);
          MPI_Allreduce(&ownership[0], &ownership_tmp[0], n_nds, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
      
          for(auto& nd : nodes())
            nd->pid(ownership_tmp[nd->gid()]);            
        }
        //-----------------------------------------------------------------------------------------------------------------
     
    
      
    
       
       
    
    
    
      //-----------------------------------------------------------------------------------------------------------------
      /// Define link between each element and its nodes
      void set_element_node_links(const std::vector<int>& my_nodes)
      {
        for(auto& el : elements())
        {    
          for(int i=0; i<el->nb_nodes(); i++)
          {
            const int nd_gid = el->node_gid_vec()[i];
            auto it = find(begin(my_nodes), end(my_nodes), nd_gid);
            const int nd_lid = it - my_nodes.begin();              /// this is corresponding to mesh_nodes on pid  
            auto& nd (nodes()[nd_lid]);
            el->el_node_vec(nd);                            /// fill the el_node_vec_ for each element with nodes having complete lid and gid
            nd->node_el_vec(el);                            /// this fill node_el_vec for each node with elements on local pid        
          }
          //now we know the nodes belonging to the element so we can compute the bound which requires node coordinates info
          el->compute_bounds();          
        }  
      }
        //-----------------------------------------------------------------------------------------------------------------
       
        
        
    
    
    
        //-----------------------------------------------------------------------------------------------------------------
        ///  This functions reads the info about side_nodes, side_flag etc for each element
        void el_side_info()
        {
          for(auto& el : elements())
          {
            if(el->on_boundary())
            {
              const auto& nds = el->el_node_vec();        
              for(int i=0; i<el->nb_sides(); i++)
              {
                 std::vector<int> side_nds;
                 for(int j=0; j<el->nb_side_nds(); j++)
                 {
                    side_nds.push_back(nds[el->side(i,j)]->gid());
                 }
                 el->set_side_nodes(std::move(side_nds));                                      
    
                 int sd_flag = 0;
                 int k = 0;
                 for(const auto& j : el->side(i))
                 {
                   if(nds[j]->on_boundary())   
                     k++;
                 }  
                 if(k == el->nb_side_nds())
                 {
                    for(std::size_t j=0; j<num_my_sides(); j++) 
                      if(are_equal(side_nds, side(j)))
                        sd_flag = side_flag(j);             
                 }                                          
                 el->set_side_flag(sd_flag);
    
              
                 bool side_on_bd = false;
                 if(sd_flag>0) side_on_bd = true;
                 el->set_side_on_boundary(side_on_bd);              
              }
            }                
          } 
          
          /// compute number of fsi sides    
          int count = 0;
          for(const auto& el : bd_elements())
           for(int i=0; i<el->nb_sides(); i++)
            if(el->is_side_on_boundary(i) &&  el->side_flag(i) == 1)
              count++;
          num_fsi_sides(count);
             
        }
        //-----------------------------------------------------------------------------------------------------------------
        
        






    /// This function computes the data at quadrature points for each element 
    //----------------------------------------------------------------------------------------------------------------
    void compute_el_quad_data()
    {
      for(auto& el : elements())
      {
        for(int i=0; i<el->nb_gp(); i++)
        {
           auto q = std::make_shared<quad>(*el, el->gp(i));
           el->quad_i()[i] = std::move(q);
        }        
      }

      for(auto& el : bd_elements())
      {
          for(int i=0; i<el->nb_sides(); i++)
          {
             std::vector<std::shared_ptr<quad>> v(el->nb_side_gp(i), nullptr);
             if(el->is_side_on_boundary(i))
             {
                for(int j=0; j<el->nb_side_gp(i); j++)
                {
                  auto q = std::make_shared<quad>(*el, el->side_gp(i,j), el->side_gp(i,j));
                  v[j] = std::move(q);
                }
             }
             el->quad_b()[i] = v;
          }
      }
    }
    //-----------------------------------------------------------------------------------------------------------------
         
        
        
        
        



     
        
                
        static int dimension(){ return dim; }                               /// it returns the mesh dimension                
                
        void tot_elements(int n){tot_elements_=n;}                          /// set total number of mesh elements
        int tot_elements()const{return tot_elements_;}                      /// get total number of mesh elements

        auto& elements() { return elements_; }                              /// it returns the vector of shared_ptr to elements on current pid
        const auto& elements()const { return elements_; }                   /// it returns the vector of shared_ptr to elements on current pid (const version)

        int num_my_elements()const{return elements_.size();}                /// get number of mesh elements on current pid
        
        auto& bd_elements() { return bd_elements_; }                        /// it returns the vector of shared_ptr to boundary elements on current pid
        const auto& bd_elements()const { return bd_elements_; }             /// it returns the vector of shared_ptr to boundary elements on current pid (const version)        

        int num_my_bd_elements()const{return bd_elements_.size();}          /// get number of boundary elements on current pid

        int nb_el_nodes()const {return elements_[0]->nb_nodes();}           /// get number of element nodes (for a single element)
        int nb_el_sides()const {return elements_[0]->nb_sides();}           /// get number of element sides (for a single element)
        int nb_el_side_nodes()const {return elements_[0]->nb_side_nds();}           /// get number of nodes of a side of an element (for a single element)





        void tot_nodes(int n){tot_nodes_=n;}                                /// set total number of mesh nodes
        int tot_nodes()const{return tot_nodes_;}                            /// get total number of mesh nodes

        auto& nodes() { return nodes_; }                                    /// it returns the vector of shared_ptr to nodes on current pid
        const auto& nodes()const { return nodes_; }                         /// it returns the vector of shared_ptr to nodes on current pid (const version) 

        int num_my_nodes()const{return nodes_.size();}                      /// get number of mesh nodes on current pid

        auto& bd_nodes() { return bd_nodes_; }                              /// it returns the vector of shared_ptr to boundary nodes on current pid
        const auto& bd_nodes()const { return bd_nodes_; }                   /// it returns the vector of shared_ptr to boundary nodes on current pid (const version)

        int num_my_bd_nodes()const{return bd_nodes_.size();}                /// get number of boundary nodes on current pid

        void num_fsi_nodes(int i){num_fsi_nodes_ = i;}                      /// set number of fsi nodes on current pid
        int num_fsi_nodes()const{return num_fsi_nodes_;}                    /// get number of fsi nodes on current pid

        int nb_nd_dofs()const {return nodes_[0]->nb_dofs();}                /// get number of node dofs (for a single node)





        void tot_sides(int n){tot_sides_=n;}                                /// set total number of boundary sides
        int tot_sides()const{return tot_sides_;}                            /// get total number of boundary sides
        
        auto& sides() { return sides_; }                                    /// it returns sides_ on current pid
        const auto& sides()const { return sides_; }                         /// it returns sides_ on current pid (const version)
        
        auto& side(int i) { return sides_[i]; }                             /// it returns i th side from sides_ vector
        const auto& side(int i)const { return sides_[i]; }                  /// it returns i th side from sides_ vector(const version)

        int num_my_sides()const{return sides_.size();}                      /// get number of mesh sides on current pid
        
        auto& side_flag() { return side_flag_; }                            /// it returns side_flag_ on current pid
        const auto& side_flag()const { return side_flag_; }                 /// it returns side_flag_ on current pid(const version)

        auto& side_flag(int i) { return side_flag_[i]; }                    /// it returns i th flag from side_flag_ vector
        const auto& side_flag(int i)const { return side_flag_[i]; }         /// it returns i th flag from side_flag_ vector(const version)

        void num_fsi_sides(int i){num_fsi_sides_ = i;}                      /// set number of fsi sides on current pid
        int num_fsi_sides()const{return num_fsi_sides_;}                    /// get number of fsi sides on current pid


        
        
  private:

        int tot_elements_;                                                  /// total number of mesh elements 
        std::vector<std::shared_ptr<element_t>> elements_;                  /// vector of shared_ptr to elements on current pid                           
        int num_my_elements_;                                               /// number of mesh elements on current pid
        std::vector<std::shared_ptr<element_t>> bd_elements_;               /// vector of shared_ptr to boundary elements on current pid
        int num_my_bd_elements_;                                            /// number of boundary elements on current pid

        int tot_nodes_;                                                     /// total number of mesh nodes
        std::vector<std::shared_ptr<node_t>> nodes_;                        /// vector of shared_ptr to nodes on current pid  
        int num_my_nodes_;                                                  /// number of mesh nodes on current pid
        std::vector<std::shared_ptr<node_t>> bd_nodes_;                     /// vector of shared_ptr to boundary nodes on current pid  
        int num_my_bd_nodes_;                                               /// number of boundary nodes on current pid
        
        int num_fsi_nodes_;
        int num_fsi_sides_;
     
        int tot_sides_;                                                     /// total number of boundary sides
        std::vector<std::vector<int>> sides_;                               /// vector of sides where each side is vector of composing nodes gids
        std::vector<int> side_flag_;                                        /// vector of side flags         
  };
    
    



}
//namespace GALES
#endif


