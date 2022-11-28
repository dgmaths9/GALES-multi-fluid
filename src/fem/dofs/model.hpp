#ifndef GALES_MODEL_HPP
#define GALES_MODEL_HPP


#include "dof_state.hpp"
#include "../mesh/mesh.hpp"


namespace GALES {


 /**
    This is simply a container class which carries mesh, dofs state, maps and setup together.
    so instead of passing all of them individually we pass the model and then from there get the stuff.
    
    In future we will also put here the initial and bounary conditions as well and hence carry all info with the object created from this class. 
 */


  template<int dim>
  class model 
  {
    
    public:

    using mesh_t = Mesh<dim>;


    model(mesh_t& m, dof_state& s, epetra_maps& maps, read_setup& setup):  mesh_(m), state_(s), maps_(maps), setup_(setup){}
    
    //-------------------------------------------------------------------------------------------------------------------------------------        
    /// Deleting the copy and move constructors - no duplication/transfer in anyway
    model(const model&) = delete;               //copy constructor
    model& operator=(const model&) = delete;    //copy assignment operator
    model(model&&) = delete;                    //move constructor  
    model& operator=(model&&) = delete;         //move assignment operator 
    //-------------------------------------------------------------------------------------------------------------------------------------


    const auto& mesh()const {return mesh_; }
    auto& mesh() {return mesh_;}

    const auto& state()const { return state_; }
    auto& state() {return state_; }
    
    const auto& maps()const { return maps_; }
    auto& maps() {return maps_; }
    
    const auto& setup()const { return setup_; }
    auto& setup() {return setup_; }
    
    
    

    ///---------------------------------------------------------------------------------------------
    /// This function extracts dofs of a node for state slot 0.
    void extract_node_dofs(const node<dim>& nd, std::vector<double>& dofs) const
    {
      dofs.resize(nd.nb_dofs());
      auto start = state_.dofs(0).begin() + nd.first_dof_lid();          
      std::copy(start, start+nd.nb_dofs(), dofs.begin());
    }
    ///---------------------------------------------------------------------------------------------





    ///---------------------------------------------------------------------------------------------
    /// This function extracts dofs of a node for all state slots.
    void extract_node_dofs(const node<dim>& nd, std::vector<std::vector<double>>& dofs) const
    {
      dofs.resize(state_.num_slot());
      
      for(int i=0; i<state_.num_slot(); i++)    
        dofs[i].resize(nd.nb_dofs());

      for(int i=0; i<state_.num_slot(); i++)
      {
        auto start = state_.dofs(i).begin() + nd.first_dof_lid();          
        std::copy(start, start+nd.nb_dofs(), dofs[i].begin());
      }          
    }
    ///---------------------------------------------------------------------------------------------


    




    ///---------------------------------------------------------------------------------------------
    /// This function extracts dofs of nodes of an element for state slot 0.     
    void extract_element_dofs(const element<dim>& el, boost::numeric::ublas::vector<double>& dofs) const
    {
      dofs.resize(mesh_.nb_el_nodes() * mesh_.nb_nd_dofs());          //resizing the container according to the number of dofs

      int j=0; 
      for(const auto& nd : el.el_node_vec())
      {
         std::vector<double> nd_dofs;
         extract_node_dofs(*nd, nd_dofs);
         std::copy(nd_dofs.begin(), nd_dofs.end(), dofs.begin()+j);
         j += nd->nb_dofs();
      }
    }
    ///---------------------------------------------------------------------------------------------






    ///---------------------------------------------------------------------------------------------
    /// This function extracts dofs of nodes of an element for all state slots. 
    void extract_element_dofs(const element<dim>& el, std::vector<boost::numeric::ublas::vector<double>>& dofs) const 
    {
      dofs.resize(state_.num_slot());            
      
      for(int i=0; i<state_.num_slot(); i++)  
        dofs[i].resize(mesh_.nb_el_nodes() * mesh_.nb_nd_dofs());          //resizing the container according to the number of dofs

      int j=0; 
      for(const auto& nd : el.el_node_vec())
      {
         std::vector<std::vector<double>> nd_dofs;
         extract_node_dofs(*nd, nd_dofs);
         for(int i=0; i<state_.num_slot(); i++)
           std::copy(nd_dofs[i].begin(), nd_dofs[i].end(), dofs[i].begin()+j);    
         j += nd->nb_dofs();
      }
    }
    ///---------------------------------------------------------------------------------------------










    ///---------------------------------------------------------------------------------------------
    /// This function extracts history dofs of an element and is used for viscoelasticity. 
    void extract_element_PP_history_dofs(const element<dim>& el, int nb_Maxwell_el,  boost::numeric::ublas::vector<double>& PP_history) const
    {
        const int el_lid(el.lid());        

        int j;        
        if(dim==2) j = el.nb_nodes()*nb_Maxwell_el*3;
        else if(dim==3) j = el.nb_nodes()*nb_Maxwell_el*6;

        PP_history.resize(j);
        
        //extraction 
        for(int i=el_lid*j, k=0; i<el_lid*j+j; i++)
        {
          PP_history[k] = state_.dofs(2)[i];
          k++;
        }
    }
    ///---------------------------------------------------------------------------------------------








    ///---------------------------------------------------------------------------------------------
    /// This function trimms the dofs between lower and upper bounds according to the dof index. 
    /// for example to trim Y_dofs(p, vx, vy, vz, T, Y) between 0 and 1 we call 
    ///   dof_treamer(5, 0.0, 1.0)
    void dof_trimmer(int dof_index, double l_bound, double u_bound)
    {      
      for(const auto& nd : mesh_.nodes())
      {
        auto dof = state_.get_dof(nd->first_dof_lid() + dof_index);
        auto new_value = std::max(l_bound, std::min(dof, u_bound) );
        state_.set_dof(nd->first_dof_lid() + dof_index, new_value);
      }  
    }
    ///---------------------------------------------------------------------------------------------





   
    
    ///---------------------------------------------------------------------------------------------
    // This function returns the dofs according to the dof index of each node in form of std::vector<double>
    // For example this can be used to collect only pressure dofs or Temperature dofs etc.     
    auto get_dofs_vec(int dof_index)const 
    {
      std::vector<double> dofs;
      dofs.reserve(mesh_.nodes().size());
      
      for(const auto& nd : mesh_.nodes())
        dofs.push_back(state_.get_dof(nd->first_dof_lid() + dof_index));
       
      return dofs;
    }    
    ///---------------------------------------------------------------------------------------------
   
    


  private:
    mesh_t& mesh_;
    dof_state& state_;
    epetra_maps& maps_;
    read_setup& setup_;
  };

}//namespace GALES
#endif

