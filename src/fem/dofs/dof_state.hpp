#ifndef GALES_DOF_STATE_HPP
#define GALES_DOF_STATE_HPP




namespace GALES 
{


  /** 
    This class defines the structure of dofs for parallel computation.
    "state_" is a vector of dofs for current and previous time/times. The time slots are defined as current=0, previous=1, ...    
    num_my_elements  = state_map_.NumMyElements() 
  */  

  class dof_state 
  {

    using vec = boost::numeric::ublas::vector<double>;

    public:

    dof_state(int num_my_elements, int n_slot):  n_slot_(n_slot), state_(n_slot)
    {
      for(int i=0; i<n_slot_; i++)
	state_[i] = std::make_unique<vec>(num_my_elements, 0.0);
    }


    //-------------------------------------------------------------------------------------------------------------------------------------        
    /// Deleting the copy and move constructors - no duplication/transfer in anyway
    dof_state(const dof_state&) = delete;               //copy constructor
    dof_state& operator=(const dof_state&) = delete;    //copy assignment operator
    dof_state(dof_state&&) = delete;                    //move constructor  
    dof_state& operator=(dof_state&&) = delete;         //move assignment operator 
    //-------------------------------------------------------------------------------------------------------------------------------------

    
    auto num_slot()const { return n_slot_; }       
    
    // this function returns the dofs values of state_ vec for given state slot index
    auto& dofs(int slot_index)const {return *state_[slot_index];}    
            
    // get dof value from state_ corresponding to given lid and state slot index 
    auto get_dof(int lid, int slot_index = 0)const {return (*state_[slot_index])[lid];}

    // set dof value of state_ equal to v  corresponding to given lid and state slot index 
    void set_dof(int lid, double d, int slot_index = 0) { (*state_[slot_index])[lid] = d;} 


  private:
  
    int n_slot_;
    std::vector<std::unique_ptr<vec>> state_;
  

  };

}//namespace GALES

#endif
