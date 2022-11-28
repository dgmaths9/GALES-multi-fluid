#ifndef _GALES_NODE_HPP_
#define _GALES_NODE_HPP_

#include <vector>
#include <tuple>
#include "point.hpp"
#include "element.hpp"



namespace GALES {



   /**
       This class contains all info for an element node such as coordinates, boundary_flag, number of dofs it carries, etc.   
   */
      
   
   //-----------------------------------------------------------------------------------------------------------------------   
   //  this is forward decleration of element class
   template<int dim>
   class element;
   //-----------------------------------------------------------------------------------------------------------------------
   
   
   //-----------------------------------------------------------------------------------------------------------------------   
   // This is generic template for node<dim>
   template<int dim>
   class node{};
   //-----------------------------------------------------------------------------------------------------------------------
   



   

   //-----------------------------------------------------------------------------------------------------------------------
   // This class is tempate specialized of node<dim> with dim = 2 
   template<>
   class node<2>
   {
     public:

         //-------------------------------------------------------------------------------------------------------------------------------------        
         node() = default;   //default constructor         
         node(const point<2>& pt) : point_(pt){} // custom constructor (this is used in pressure profiles)        
         /// Deleting the copy and move constructors - no duplication/transfer in anyway
         node(const node&) = delete;               //copy constructor
         node& operator=(const node&) = delete;    //copy assignment operator
         node(node&&) = delete;                    //move constructor  
         node& operator=(node&&) = delete;         //move assignment operator 
         //-------------------------------------------------------------------------------------------------------------------------------------
     
         void pid(int i){pid_=i;}                                                                /// set pid of the node
         int pid()const{return pid_;}                                                            /// get pid of the node  
       
         void gid(int i){gid_=i;}                                                                /// set gid of the node
         int gid()const{return gid_;}                                                            /// get gid of the node
   
         void lid(int i){lid_=i;}                                                                /// set lid of the node
         int lid()const{return lid_;}                                                            /// get lid of the node
   
         void on_boundary(bool a) { on_boundary_=a;}                                             /// set whether the node is on boundary
         bool on_boundary() const { return on_boundary_;}                                        /// returns whether the node is on boundary
   
         void flag(int flag) {flag_ = flag;}                                                     /// set node flag for location;   interior: 0;    boundary: non zero;   fluid-solid interface: 1
         int flag() const {return flag_;}                                                        /// returns the flag
         
         void nb_dofs(int nb_dofs) {nb_dofs_ = nb_dofs;}                                         /// set the number of dofs for node            
         int nb_dofs()const { return nb_dofs_; }                                                 /// returns the number of dofs for node
   
         int first_dof_gid()const {return gid_*nb_dofs_;}                                       /// This returns first_dof_gid of node

         void first_dof_lid(int lid) {first_dof_lid_ = lid;}                                     /// This is called from maps.hpp to set the first_dof_lid
         int first_dof_lid()const {return first_dof_lid_;}                                       /// This is called from dof_extractor

         void set_x(double x) { point_.set_x(x);}                                  /// set x-coord of the node
         void set_y(double y) { point_.set_y(y);}                                  /// set y-coord of the node
         void set_z(double z) { /* do nothing */; }                                                  
   
         double get_x() const {return point_.get_x(); }                            /// get x-coord of the node
         double get_y() const {return point_.get_y(); }                            /// get y-coord of the node

         void coord(const point<2>& pt) {point_ = pt;}                                /// set nd coordinates in form of point
         const point<2>& coord()const { return point_;}                               /// get nd coordinates in form of point
        
         static int dimension() { return 2; }                                      /// returns the node dimension
   
         void node_el_vec(const std::shared_ptr<element<2>>& el) { node_el_vec_.push_back(el);}     /// fill the node_el_vec_
         auto node_el_vec() const {return node_el_vec_;}                                         /// returns the node_el_vec_

         ///---------------- These are for mesh motion -------------------------------------------------------------------------
         void clear_updated_coord() { updated_node_coords_=false;}                               /// it indicates that nd is free to be updated
         bool is_updated_coord() { return updated_node_coords_;}                                 /// it returns the value of the flag
         void update_coord(const std::vector<double>& v) 
         {
           set_x(get_x() + v[0]);           // update the x-coord
           set_y(get_y() + v[1]);           // update the y-coord            
           updated_node_coords_ = true;   // set the flag which indicates that the node coords are updated
         }
         ///---------------- ---------------------------------------------------------------------------------------------------

   
     
     private:
         int pid_;                                              /// process id of the item
         int gid_;                                              /// global id of the item among all processes
         int lid_;                                              /// local id of the item on the owner process
         bool on_boundary_;                                     /// flag to check if node is on boundary or not  
         int flag_;                                             /// flag for location: interior(0) and boundary(non zero)
         int nb_dofs_;                                          /// number of dofs carried by the node
         int first_dof_lid_;                                    /// the lid of the first dof of the node
         bool updated_node_coords_ = false;                             /// flag to check if node coord are updated or not in case of mesh motion
         point<2> point_;                                            /// coordinates of the node
         std::vector<std::shared_ptr<element<2>>> node_el_vec_;    /// vector of elements containing the node    
   };
   //-----------------------------------------------------------------------------------------------------------------------



   




   //-----------------------------------------------------------------------------------------------------------------------
   // This class is tempate specialized of node<dim> with dim = 3 
   template<>
   class node<3>
   {
     public:
     
         //-------------------------------------------------------------------------------------------------------------------------------------        
         node() = default;   //default constructor
         node(const point<3>& pt) : point_(pt){} // custom constructor  (this is used in pressure profiles)          
         /// Deleting the copy and move constructors - no duplication/transfer in anyway
         node(const node&) = delete;               //copy constructor
         node& operator=(const node&) = delete;    //copy assignment operator
         node(node&&) = delete;                    //move constructor  
         node& operator=(node&&) = delete;         //move assignment operator 
         //-------------------------------------------------------------------------------------------------------------------------------------

         void pid(int i){pid_=i;}                                                                /// set pid of the node
         int pid()const{return pid_;}                                                            /// get pid of the node  
       
         void gid(int i){gid_=i;}                                                                /// set gid of the node
         int gid()const{return gid_;}                                                            /// get gid of the node
   
         void lid(int i){lid_=i;}                                                                /// set lid of the node
         int lid()const{return lid_;}                                                            /// get lid of the node
   
         void on_boundary(bool a) { on_boundary_=a;}                                             /// set whether the node is on boundary
         bool on_boundary() const { return on_boundary_;}                                        /// returns whether the node is on boundary
   
         void flag(int flag) {flag_ = flag;}                                                     /// set node flag for location;   interior: 0;    boundary: non zero;   fluid-solid interface: 1
         int flag() const {return flag_;}                                                        /// returns the flag
         
         void nb_dofs(int nb_dofs) {nb_dofs_ = nb_dofs;}                                         /// set the number of dofs for node            
         int nb_dofs()const { return nb_dofs_; }                                                 /// returns the number of dofs for node
   
         int first_dof_gid()const {return gid()*nb_dofs_;}                                       /// This returns first_dof_gid of node
         void first_dof_lid(int lid) {first_dof_lid_ = lid;}                                     /// This is called from maps.hpp to set the first_dof_lid
         int first_dof_lid()const {return first_dof_lid_;}                                       /// This is called from dof_extractor
   
         void set_x(double x) { point_.set_x(x);}                                  /// set x-coord of the node
         void set_y(double y) { point_.set_y(y);}                                  /// set y-coord of the node
         void set_z(double z) { point_.set_z(z);}                                  /// set z-coord of the node
   
         double get_x() const {return point_.get_x(); }                                              /// get x-coord of the node
         double get_y() const {return point_.get_y(); }                                              /// get y-coord of the node
         double get_z() const {return point_.get_z(); }                                              /// get z-coord of the node
        
         void coord(const point<3>& pt) {point_ = pt;}                                /// set nd coordinates in form of point
         const point<3>& coord()const { return point_;}                               /// get nd coordinates in form of point

         static int dimension() { return 3; }                                                    /// returns the node dimension
   
         void node_el_vec(const std::shared_ptr<element<3>>& el) { node_el_vec_.push_back(el);}     /// fill the node_el_vec_
         auto node_el_vec() const {return node_el_vec_;}                                         /// returns the node_el_vec_

         ///---------------- These are for mesh motion -------------------------------------------------------------------------
         void update_coord(const std::vector<double>& v) 
         {
           set_x(get_x() + v[0]);           // update the x-coord
           set_y(get_y() + v[1]);           // update the y-coord            
           set_z(get_z() + v[2]);           // update the z-coord            
         }
         ///---------------- ---------------------------------------------------------------------------------------------------

     
     private:
         int pid_;                                              /// process id of the item
         int gid_;                                              /// global id of the item among all processes
         int lid_;                                              /// local id of the item on the owner process
         bool on_boundary_;                                     /// flag to check if node is on boundary or not  
         int flag_;                                             /// flag for location: interior(0) and boundary(non zero)
         int nb_dofs_;                                          /// number of dofs carried by the node
         int first_dof_lid_;                                    /// the lid of the first dof of the node
         point<3> point_;                                            /// coordinates of the node
         std::vector<std::shared_ptr<element<3>>> node_el_vec_;    /// vector of elements containing the node    
   };
   //-----------------------------------------------------------------------------------------------------------------------
   
   
   


} //namespace GALES
#endif

