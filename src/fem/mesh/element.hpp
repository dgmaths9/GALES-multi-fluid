#ifndef _GALES_ELEMENT_HPP_
#define _GALES_ELEMENT_HPP_

#include <vector>
#include "node.hpp"



namespace GALES {


/**
  A generic mesh element.
  This class implements a generic mesh elements: it can be composed by any number of nodes; also, it can have any number of adjacent elements.
*/


  //-----------------------------------------------------------------------------------------------------------------------   
  //  this is forward decleration of element class
  class quad;    
  //-----------------------------------------------------------------------------------------------------------------------







  //-----------------------------------------------------------------------------------------------------------------------   
  // This is generic template for node<dim>
  template<int dim>
  class element{};
  //-----------------------------------------------------------------------------------------------------------------------





  //-----------------------------------------------------------------------------------------------------------------------
  // This class is tempate specialized of node<dim> with dim = 2 
  template<>
  class element<2>
  {
    public:

    element(int nb_nodes): nb_nodes_(nb_nodes) {gauss_data();}
    
    
    //-------------------------------------------------------------------------------------------------------------------------------------        
    /// Deleting the copy and move constructors - no duplication/transfer in anyway
    element(const element&) = delete;               //copy constructor
    element& operator=(const element&) = delete;    //copy assignment operator
    element(element&&) = delete;                    //move constructor  
    element& operator=(element&&) = delete;         //move assignment operator 
    //-------------------------------------------------------------------------------------------------------------------------------------






  //-----------------------------------------------------------------------------------------------------------------------   
    void gauss_data()
    {

      /**
                          |eta
                          |
                          | 
                          |               
        	    3 ___2|____2
        	     |    |    |
        	     |    |____|____________xi
        	    3|         |1
        	     |_________|
                     0    0    1
      */

      if(nb_nodes_ == 4)
      {
        const double g_minus = -1.0/sqrt(3.0);      // -0.57735026918962573105886804115
        const double g_plus = 1.0/sqrt(3.0);       //  0.57735026918962573105886804115;

        ///------------------------------interior-------------------------------------------
        /// gauss points for element interior 
        gp_.resize(4);
        gp_[0] = point<2>(g_minus,g_minus);
        gp_[1] = point<2>(g_plus,g_minus);
        gp_[2] = point<2>(g_plus,g_plus);
        gp_[3] = point<2>(g_minus,g_plus);
        ///---------------------------------------------------------------------------------
        


        ///------------------------------boundary-------------------------------------------
        /// The order of nodes is anticlockwise looking from outside of the cell.
        nb_sides_ = 4;
        nb_side_nds_ = 2;
        
        side_.resize(4);
        side_[0] = {0,1};  ///bottom --->  eta = -1
        side_[1] = {1,2};  ///right  --->  xi = 1
        side_[2] = {2,3};  ///top    --->  eta = 1
        side_[3] = {3,0};  ///left   --->  xi = -1
        
        /// To facilitate fsi fluid traction computation, gauss points are ordered in accordance to side nodes. 
        side_gp_.resize(4);
        side_gp_[0] = {point<2>(g_minus,-1.0), point<2>(g_plus,-1.0)};  ///bottom
        side_gp_[1] = {point<2>(1.0,g_minus),  point<2>(1.0,g_plus)};   ///right
        side_gp_[2] = {point<2>(g_plus,1.0),  point<2>(g_minus,1.0)};   ///top
        side_gp_[3] = {point<2>(-1.0,g_plus), point<2>(-1.0,g_minus)};  ///left
        ///---------------------------------------------------------------------------------
      }                        

      /**
                     |eta
                     | 
                     |
                     |
                    2|
        	     |\
        	     |  \ 0
        	    1|    \
        	     |______\_________________
                     0   2    1               xi
      */  



      else if(nb_nodes_ == 3)
      {
        ///------------------------------interior-------------------------------------------
        gp_.resize(3);
        gp_[0] = point<2>(1.0/6.0, 1.0/6.0);
        gp_[1] = point<2>(4.0/6.0, 1.0/6.0);
        gp_[2] = point<2>(1.0/6.0, 4.0/6.0);
        ///---------------------------------------------------------------------------------



        ///------------------------------boundary-------------------------------------------
        nb_sides_ = 3;
        nb_side_nds_ = 2;
        
        side_.resize(3);
        /// the order of nodes is anticlockwise looking from outside of the cell
        side_[0] = {1,2};  /// xi + eta = 1
        side_[1] = {2,0};  /// xi = 0
        side_[2] = {0,1};  /// eta = 0


        /// xi, eta for master triangle edges vary in range [0-1];
        /// To compute gauss quadrature rule on line integral from [-1,1], we define a mapping 
        /// f(t) = (1-t)/2  which maps -1--->1 and 1---->0:   [-1,1]---->[0,1] with determinant of transformation is 1/2
        /// therefore the gauss points are 0.5(1-1/sqrt(3)) and 0.5*(1+1/sqrt(3)) instead of 1/sqrt(3) and -1/sqrt(3)
                
        const double g_minus = 0.5*(1.0 - 1.0/sqrt(3.0));
        const double g_plus = 0.5*(1.0 + 1.0/sqrt(3.0));         
               
        /// to facilitate fsi fluid traction computation, gauss points are ordered in accordance to side nodes 
        side_gp_.resize(3);
        side_gp_[0] = { point<2>(g_plus,1.0-g_plus), point<2>(g_minus,1.0-g_minus) };  /// x+y=1
        side_gp_[1] = { point<2>(0.0,g_plus), point<2>(0.0,g_minus) };  /// x=0;  
        side_gp_[2] = { point<2>(g_minus,0.0), point<2>(g_plus,0.0) };  /// y=0;
        ///---------------------------------------------------------------------------------
      }
      
      for(int i=0; i<nb_gp(); i++)
        quad_i_.push_back(nullptr);

      for(int i=0; i<nb_sides(); i++)
      {
        std::vector<std::shared_ptr<quad>> v(nb_side_gp(i), nullptr);
        quad_b_.push_back(v);
      }  
    }    
  //-----------------------------------------------------------------------------------------------------------------------   








  //-----------------------------------------------------------------------------------------------------------------------   
    void compute_bounds()  // this function computes the element bounds [min_x, max_x, min_y, max_y]
    {
       std::vector<double> x, y;
       x.reserve(nb_nodes_);
       y.reserve(nb_nodes_);
       for(const auto& nd : el_node_vec_)
       {
         x.push_back(nd->get_x());
         y.push_back(nd->get_y());
       }  
       auto minmax_x = std::minmax_element(x.begin(), x.end());
       auto minmax_y = std::minmax_element(y.begin(), y.end());
       
       min_x_ = *minmax_x.first;       max_x_ = *minmax_x.second;       
       min_y_ = *minmax_y.first;       max_y_ = *minmax_y.second;
    }
  //-----------------------------------------------------------------------------------------------------------------------   





    double min_x()const {return min_x_;}
    double max_x()const {return max_x_;}
    double min_y()const {return min_y_;}
    double max_y()const {return max_y_;}
       
                
    void pid(int i){pid_=i;}                                                                  /// set pid of the node
    int pid()const{return pid_;}                                                              /// get pid of the node  
    
    void gid(int i){gid_=i;}                                                                  /// set gid of the node
    int gid()const{return gid_;}                                                              /// get gid of the node

    void lid(int i){lid_=i;}                                                                  /// set lid of the node
    int lid()const{return lid_;}                                                              /// get lid of the node
  
    void node_gid_vec(const std::vector<int>& node_gid_vec) { node_gid_vec_ = std::move(node_gid_vec);}                                /// set vector of node gids 
    auto node_gid_vec()const { return node_gid_vec_; }                                                                                 /// get vector of node gids

    int nb_nodes()const { return nb_nodes_; }                                                 /// get the number of nodes of the element
    
    void on_boundary(bool a) { on_boundary_ = a; }                                            /// set the flag if el is on boundary or not
    bool on_boundary()const { return on_boundary_; }                                          /// get the flag if el is on boundary or not                       



    void set_side_nodes(const std::vector<int>& v) {side_nodes_.push_back(v);}                /// set nodes for a side of the element
    auto side_nodes(int sd_idx)const {return side_nodes_[sd_idx];}                            /// get nodes for a(sd_idx) side of the element
    auto side_nodes(int sd_idx, int nd_idx)const {return side_nodes_[sd_idx][nd_idx];}        /// get node for a(sd_idx and nd_indx) side of the element

    void set_side_flag(int side_flag) {side_flag_.push_back(side_flag);}                      /// set side flag for a side of the element
    int side_flag(int sd_idx)const {return side_flag_[sd_idx];}                               /// get side flag for a(sd_idx) side of the element
    
    void set_side_on_boundary(bool flag) {is_side_on_boundary_.push_back(flag);}              /// set vector of side flags that tells which sides are on boundary
    bool is_side_on_boundary(int sd_idx)const {return is_side_on_boundary_[sd_idx];}          /// get a flag that tells if a(sd_idx) side is on boundary or not
             


    void el_node_vec(const std::shared_ptr<node<2>>& nd) { el_node_vec_.push_back(nd);}       /// set vector of nodes
    auto el_node_vec() const {return el_node_vec_;}                                           /// get vector of nodes 

    int nd_gid(int i) const{return el_node_vec_[i]->gid();}                                   /// get gid of ith node 
    int nd_lid(int i) const{return el_node_vec_[i]->lid();}                                   /// get lid of ith node                

    void set_x(int nd_index, double x) {el_node_vec_[nd_index]->set_x(x);}                    /// set x-coord of the node
    void set_y(int nd_index, double y) {el_node_vec_[nd_index]->set_y(y);}                    /// set y-coord of the node

    double get_x(int nd_index) const {return el_node_vec_[nd_index]->get_x();}                /// get x-coord of the node
    double get_y(int nd_index) const {return el_node_vec_[nd_index]->get_y();}                /// get y-coord of the node
    
             
    
    auto gp()const {return gp_;}
    auto gp(int i)const {return gp_[i];}    
    int nb_gp()const{return gp_.size();}
    
    auto side()const {return side_;}
    auto side(int side_index)const {return side_[side_index];}
    auto side(int side_index, int node_index)const {return side_[side_index][node_index];}    
    
    auto side_gp()const {return side_gp_;}
    auto side_gp(int side_index)const {return side_gp_[side_index];}
    auto side_gp(int side_index, int gp_index)const {return side_gp_[side_index][gp_index];}    
    int nb_side_gp(int side_index) const {return side_gp_[side_index].size();}
    
    int nb_sides()const {return nb_sides_;}
    int nb_side_nds()const {return nb_side_nds_;}
    
    
    
    auto& quad_i() { return quad_i_; }                        
    const auto& quad_i()const { return quad_i_; }                    

    auto& quad_b(){return quad_b_;}
    const auto& quad_b()const {return quad_b_;} 
           
    auto& quad_b(int side_index){return quad_b_[side_index];}    
    const auto& quad_b(int side_index)const{return quad_b_[side_index];} 
       
    auto& quad_b(int side_index, int gp_index){return quad_b_[side_index][gp_index];}    
    const auto& quad_b(int side_index, int gp_index)const{return quad_b_[side_index][gp_index];}    
    

  private:

    int pid_;                                              /// process id of the item
    int gid_;                                              /// global id of the item among all processes
    int lid_;                                              /// local id of the item on the owner process
    bool on_boundary_;                                     /// flag: if el is on bundary
    int nb_nodes_;                                         /// number of nodes with which el is composed 
    std::vector<int> node_gid_vec_;                        /// vector of node gids belong to the element
    std::vector<std::vector<int>> side_nodes_;             /// vector of node_gids for each side of the el
    std::vector<int> side_flag_;                           /// vector of side_flag for each side of the el
    std::vector<bool> is_side_on_boundary_;                /// vector of flag(if side is on bd) for each side of the el
    std::vector<std::shared_ptr<node<2>>> el_node_vec_;    /// vector of nodes belong to the element
    double min_x_, max_x_, min_y_, max_y_;

    std::vector<point<2>> gp_;                          /// Vector of gauss points for integration of element interior
    std::vector<std::vector<int>> side_;                  /// Vector of sides of an element
    std::vector<std::vector<point<2>>> side_gp_;        /// Vector of gauss points for integration of each side
    int nb_sides_;                                         /// Number of sides of an element
    int nb_side_nds_;                                     /// Number of nodes of side of an element
    double W_in_;                                         /// Weight of gauss points for integration of element interior 
    double W_bd_;                                         /// Weight of gauss points for integration of element face

    std::vector<std::shared_ptr<quad>> quad_i_;                            ///quad data computed at each gp for interior integration  
    std::vector<std::vector<std::shared_ptr<quad>>> quad_b_;               ///quad data computed at each gp for boundary integration  

    
  };
  //-----------------------------------------------------------------------------------------------------------------------



   












  //-----------------------------------------------------------------------------------------------------------------------
  // This class is tempate specialized of node<dim> with dim = 3 
  template<>
  class element<3>
  {
    public:
                
    element(int nb_nodes): nb_nodes_(nb_nodes) {gauss_data();}
    
    
    //-------------------------------------------------------------------------------------------------------------------------------------        
    /// Deleting the copy and move constructors - no duplication/transfer in anyway
    element(const element&) = delete;               //copy constructor
    element& operator=(const element&) = delete;    //copy assignment operator
    element(element&&) = delete;                    //move constructor  
    element& operator=(element&&) = delete;         //move assignment operator 
    //-------------------------------------------------------------------------------------------------------------------------------------







    //-------------------------------------------------------------------------------------------------------------------------------------        
    void gauss_data()
    {

      /*  
                           |eta
                           |     
                           |     
   \                3------|------2                       
    \               |\     |      |\           
     \              | \    |      | \          
      \             |  \   |      |  \         
       \            |   7=============6        
        \           |   || +------|--||--------> xi   
         \          0---||--\-----1  ||        
          \          \  ||   \     \ ||        
           \          \ ||    \     \||        
            \          \||     \     ||        
             \          4=============5        
              \                  \      
               \                  \ zeta 
                \                         
      */


      if(nb_nodes_ == 8)
      {
        const double g_minus = -1.0/sqrt(3.0);
        const double g_plus = 1.0/sqrt(3.0);

        gp_.resize(8);
        gp_[0] = point<3>(g_minus,g_minus,g_minus);        
        gp_[1] = point<3>(g_plus,g_minus,g_minus);
        gp_[2] = point<3>(g_plus,g_plus,g_minus);
        gp_[3] = point<3>(g_minus,g_plus,g_minus);
        gp_[4] = point<3>(g_minus,g_minus,g_plus);
        gp_[5] = point<3>(g_plus,g_minus,g_plus);
        gp_[6] = point<3>(g_plus,g_plus,g_plus);
        gp_[7] = point<3>(g_minus,g_plus,g_plus);
        
        nb_sides_ = 6;
        nb_side_nds_ = 4;        
        side_.resize(6);
        /// the order of nodes is anticlockwise looking from outside of the cell
        side_[0] = {4,5,6,7};  ///front   --->   zeta = 1
        side_[1] = {0,3,2,1};  ///back    --->   zeta = -1
        side_[2] = {0,4,7,3};  ///left    --->   xi = -1
        side_[3] = {1,2,6,5};  ///right   --->   xi = 1
        side_[4] = {2,3,7,6};  ///top     --->   eta = 1 
        side_[5] = {1,5,4,0};  ///bottom  --->   eta = -1
        
        /// to facilitate fsi fluid traction computation, gauss points are ordered in accordance to side nodes 
        side_gp_.resize(6);
        side_gp_[0] = {point<3>(g_minus,g_minus,1.0),  point<3>(g_plus,g_minus,1.0),  point<3>(g_plus,g_plus,1.0),  point<3>(g_minus,g_plus,1.0)};   ///front                
        side_gp_[1] = {point<3>(g_minus,g_minus,-1.0), point<3>(g_minus,g_plus,-1.0), point<3>(g_plus,g_plus,-1.0), point<3>(g_plus,g_minus,-1.0)};  ///back                
        side_gp_[2] = {point<3>(-1.0,g_minus,g_minus), point<3>(-1.0,g_minus,g_plus), point<3>(-1.0,g_plus,g_plus), point<3>(-1.0,g_plus,g_minus)};  ///left                        
        side_gp_[3] = {point<3>(1.0,g_minus,g_minus),  point<3>(1.0,g_plus,g_minus),  point<3>(1.0,g_plus,g_plus),  point<3>(1.0,g_minus,g_plus)};   ///right                        
        side_gp_[4] = {point<3>(g_plus,1.0,g_minus),   point<3>(g_minus,1.0,g_minus), point<3>(g_minus,1.0,g_plus), point<3>(g_plus,1.0,g_plus)};    ///top                
        side_gp_[5] = {point<3>(g_plus,-1.0,g_minus),  point<3>(g_plus,-1.0,g_plus),  point<3>(g_minus,-1.0,g_plus), point<3>(g_minus,-1.0,g_minus)}; ///bottom
      }                        


      /*
                           |eta
                           |
                           | 
                           |2
                          /|\,
                         / |  \,
                        /  |    \,
                       /   |      \,
                      /    |        \,
                     /     |__________\______________xi
                    /    ,/0          /1
                   /   ,/         / ^
                  /  ,/       / ^
                 / ,/    / ^
                /,/, / ^
              3/// ^ 
              /     
             /
            /zeta

      */ 

      else if(nb_nodes_ == 4 )
      {
        gp_.resize(4);
        gp_[0] = point<3>(0.1381966011250105, 0.1381966011250105, 0.1381966011250105);
        gp_[1] = point<3>(0.5854101966249685, 0.1381966011250105, 0.1381966011250105);
        gp_[2] = point<3>(0.1381966011250105, 0.5854101966249685, 0.1381966011250105);
        gp_[3] = point<3>(0.1381966011250105, 0.1381966011250105, 0.5854101966249685);

        nb_sides_ = 4;
        nb_side_nds_ = 3;
        
        side_.resize(4);
        /// the order of nodes is anticlockwise looking from outside of the cell
        side_[0] = {1,2,3};  /// xi + eta + zeta = 1
        side_[1] = {0,3,2};  /// xi = 0
        side_[2] = {0,1,3};  /// eta = 0
        side_[3] = {0,2,1};  /// zeta = 0

        side_gp_.resize(4);
        side_gp_[0] = {point<3>(4.0/6.0, 1.0/6.0, 1.0/6.0), point<3>(1.0/6.0, 4.0/6.0, 1.0/6.0), point<3>(1.0/6.0, 1.0/6.0, 4.0/6.0)};  /// xi + eta + zeta = 1
        side_gp_[1] = {point<3>(0.0, 1.0/6.0, 1.0/6.0), point<3>(0.0, 1.0/6.0, 4.0/6.0), point<3>(0.0, 4.0/6.0, 1.0/6.0)};  /// xi = 0
        side_gp_[2] = {point<3>(1.0/6.0, 0.0, 1.0/6.0), point<3>(4.0/6.0, 0.0, 1.0/6.0), point<3>(1.0/6.0, 0.0, 4.0/6.0)};  /// eta = 0
        side_gp_[3] = {point<3>(1.0/6.0, 1.0/6.0, 0.0), point<3>(1.0/6.0, 4.0/6.0, 0.0), point_3d(4.0/6.0, 1.0/6.0, 0.0)};  /// zeta = 0
      }
      
      for(int i=0; i<nb_gp(); i++)
        quad_i_.push_back(nullptr);

      for(int i=0; i<nb_sides(); i++)
      {
        std::vector<std::shared_ptr<quad>> v(nb_side_gp(i), nullptr);
        quad_b_.push_back(v);
      }  
    }    
    //-------------------------------------------------------------------------------------------------------------------------------------        








    //-------------------------------------------------------------------------------------------------------------------------------------        
    void compute_bounds()  // this function computes the element bounds [min_x, max_x, min_y, max_y, min_z, max_z]
    {
       std::vector<double> x, y, z;
       x.reserve(nb_nodes_);
       y.reserve(nb_nodes_);
       z.reserve(nb_nodes_);
       for(const auto& nd : el_node_vec_)
       {
         x.push_back(nd->get_x());
         y.push_back(nd->get_y());
         z.push_back(nd->get_z());
       }  
       auto minmax_x = std::minmax_element(x.begin(), x.end());
       auto minmax_y = std::minmax_element(y.begin(), y.end());
       auto minmax_z = std::minmax_element(z.begin(), z.end());
       
       min_x_ = *minmax_x.first;       max_x_ = *minmax_x.second;       
       min_y_ = *minmax_y.first;       max_y_ = *minmax_y.second;
       min_z_ = *minmax_z.first;       max_z_ = *minmax_z.second;  
    }
    //-------------------------------------------------------------------------------------------------------------------------------------        





    double min_x()const {return min_x_;}
    double max_x()const {return max_x_;}
    double min_y()const {return min_y_;}
    double max_y()const {return max_y_;}
    double min_z()const {return min_z_;}
    double max_z()const {return max_z_;}

    
    void pid(int i){pid_=i;}                                                                  /// set pid of the node
    int pid()const{return pid_;}                                                              /// get pid of the node  
    
    void gid(int i){gid_=i;}                                                                  /// set gid of the node
    int gid()const{return gid_;}                                                              /// get gid of the node

    void lid(int i){lid_=i;}                                                                  /// set lid of the node
    int lid()const{return lid_;}                                                              /// get lid of the node
  
    void node_gid_vec(const std::vector<int>& node_gid_vec) { node_gid_vec_ = std::move(node_gid_vec);}                                /// set vector of node gids 
    auto node_gid_vec()const { return node_gid_vec_; }                                                                                 /// get vector of node gids

    int nb_nodes()const { return nb_nodes_; }                                                 /// get the number of nodes of the element
    
    void on_boundary(bool a) { on_boundary_ = a; }                                            /// set the flag if el is on boundary or not
    bool on_boundary()const { return on_boundary_; }                                          /// get the flag if el is on boundary or not                       



    void set_side_nodes(const std::vector<int>& v) {side_nodes_.push_back(v);}                /// set nodes for a side of the element
    auto side_nodes(int sd_idx)const {return side_nodes_[sd_idx];}                            /// get nodes for a(sd_idx) side of the element
    auto side_nodes(int sd_idx, int nd_idx)const {return side_nodes_[sd_idx][nd_idx];}        /// get node for a(sd_idx and nd_indx) side of the element

    void set_side_flag(int side_flag) {side_flag_.push_back(side_flag);}                      /// set side flag for a side of the element
    int side_flag(int sd_idx)const {return side_flag_[sd_idx];}                               /// get side flag for a(sd_idx) side of the element
    
    void set_side_on_boundary(bool flag) {is_side_on_boundary_.push_back(flag);}              /// set vector of side flags that tells which sides are on boundary
    bool is_side_on_boundary(int sd_idx)const {return is_side_on_boundary_[sd_idx];}          /// get a flag that tells if a(sd_idx) side is on boundary or not
             


    void el_node_vec(const std::shared_ptr<node<3>>& nd) { el_node_vec_.push_back(nd);}       /// set vector of nodes
    auto el_node_vec() const {return el_node_vec_;}                                           /// get vector of nodes 

    int nd_gid(int i) const{return el_node_vec_[i]->gid();}                                   /// get gid of ith node 
    int nd_lid(int i) const{return el_node_vec_[i]->lid();}                                   /// get lid of ith node                
    
    void set_x(int nd_index, double x) {el_node_vec_[nd_index]->set_x(x);}                    /// set x-coord of the node
    void set_y(int nd_index, double y) {el_node_vec_[nd_index]->set_y(y);}                    /// set y-coord of the node
    void set_z(int nd_index, double z) {el_node_vec_[nd_index]->set_z(z);}                    /// set z-coord of the node

    double get_x(int nd_index) const {return el_node_vec_[nd_index]->get_x();}                /// get x-coord of the node
    double get_y(int nd_index) const {return el_node_vec_[nd_index]->get_y();}                /// get y-coord of the node
    double get_z(int nd_index) const {return el_node_vec_[nd_index]->get_z();}                /// get y-coord of the node
        
    

    

    auto gp()const {return gp_;}
    auto gp(int i)const {return gp_[i];}
    int nb_gp()const{return gp_.size();}
    
    auto side()const {return side_;}
    auto side(int side_index)const {return side_[side_index];}
    auto side(int side_index, int node_index)const {return side_[side_index][node_index];}  
      
    auto side_gp()const {return side_gp_;}
    auto side_gp(int side_index)const {return side_gp_[side_index];}
    auto side_gp(int side_index, int gp_index)const {return side_gp_[side_index][gp_index];}
    int nb_side_gp(int side_index) const {return side_gp_[side_index].size();}
    
    int nb_sides()const {return nb_sides_;}
    int nb_side_nds()const {return nb_side_nds_;}
 
       
    
    auto& quad_i() { return quad_i_; }                        
    const auto& quad_i()const { return quad_i_; }                     

    auto& quad_b(){return quad_b_;}
    const auto& quad_b()const {return quad_b_;}        
    
    auto& quad_b(int side_index){return quad_b_[side_index];}   
    const auto& quad_b(int side_index)const{return quad_b_[side_index];}    
    
    auto& quad_b(int side_index, int gp_index){return quad_b_[side_index][gp_index];}    
    const auto& quad_b(int side_index, int gp_index)const{return quad_b_[side_index][gp_index];}    
    
    

  private:
  
    int pid_;                                              /// process id of the item
    int gid_;                                              /// global id of the item among all processes
    int lid_;                                              /// local id of the item on the owner process
    bool on_boundary_;                                     /// flag: if el is on bundary
    int nb_nodes_;                                         /// number of nodes with which el is composed 
    std::vector<int> node_gid_vec_;                        /// vector of node gids belong to the element
    std::vector<std::vector<int>> side_nodes_;             /// vector of node_gids for each side of the el
    std::vector<int> side_flag_;                           /// vector of side_flag for each side of the el
    std::vector<bool> is_side_on_boundary_;                /// vector of flag(if side is on bd) for each side of the el
    std::vector<std::shared_ptr<node<3>>> el_node_vec_;    /// vector of nodes belong to the element
    double min_x_, max_x_, min_y_, max_y_, min_z_, max_z_;

    std::vector<point<3>> gp_;                          /// Vector of gauss points for integration of element interior
    std::vector<std::vector<int>> side_;                  /// Vector of sides of an element
    std::vector<std::vector<point<3>>> side_gp_;        /// Vector of gauss points for integration of each side
    int nb_sides_;                                         /// Number of sides of an element
    int nb_side_nds_;                                     /// Number of nodes of side of an element

    std::vector<std::shared_ptr<quad>> quad_i_;                            ///quad data computed at each gp for interior integration  
    std::vector<std::vector<std::shared_ptr<quad>>> quad_b_;               ///quad data computed at each gp for boundary integration  
    
  };
  //-----------------------------------------------------------------------------------------------------------------------
  






}//namespace GALES
#endif
