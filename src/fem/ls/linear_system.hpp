#ifndef GALES_LINEAR_SYSTEM_
#define GALES_LINEAR_SYSTEM_



#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"
#include "Epetra_Map.h"
#include "Epetra_MpiComm.h"
#include "Epetra_FECrsGraph.h"
#include "Epetra_FECrsMatrix.h"
#include "Epetra_FEVector.h"



namespace GALES {


  /**
      This class defines the linear system Ax=b in GALES. We define only A and b in this class as x is defined and computed in the solver class
      A ---> pointer to Epetra_FECrsMatrix (global matrix)
      b ---> pointer to Epetra_FEvector  (global rhs)

      In the constructor, we create matrix and rhs of the linear system.
      The matrix is created by Epetra_FECrsGraph for the performance gain and the rhs is created using dof_map.
      Dof map is a non shared distributed Epetra_Map constructed from gids of dofs of owned nodes on each process and is defined in the maps class.  
      
      Note that the most efficient but least flexible fill method is to create the Epetra_FECrsMatrix with a constant graph. 
      That is, to create a Epetra_FECrsGraph separately, fill it, call its FillComplete() method, then pass the graph to the Epetra_FECrsMatrix constructor. 
      This completely constrains the structure of the CrsMatrix. We can only set or modify values in the matrix, not the structure. 
      This means we can not call InsertGlobalValues() or InsertMyValues(), only the Replace*Values(), SumInto*Values, Scale(), and PutScalar() methods.  
  */


  class linear_system 
  {
    public:
    
    template<typename model_type>
    explicit linear_system(model_type& model)
    {
      const auto& mesh(model.mesh());
      const auto& maps(model.maps());         

      /// graph construction
      auto graph = std::make_shared<Epetra_FECrsGraph>(Copy, *(maps.dof_map()), 1);                   
      for(const auto& el: mesh.elements())
      {
        const int num_el_nds = el->nb_nodes();
        const int num_nd_dofs = el->el_node_vec()[0]->nb_dofs();
        const int n = num_el_nds*num_nd_dofs;
        std::vector<int> dof_gid_index(n);        
        for(int i=0; i<num_el_nds; i++) 
        {
	  const int nd_lid = el->nd_lid(i);
  	  const auto& nd(mesh.nodes()[nd_lid]);
	  const int fd_gid = nd->first_dof_gid();
	  for(int k=0; k<num_nd_dofs; k++)
	    dof_gid_index[num_nd_dofs*i + k] = fd_gid + k;
        }      
        graph->InsertGlobalIndices(n, &dof_gid_index[0], n, &dof_gid_index[0]);   
      }        
      graph->GlobalAssemble();      
      matrix(std::make_shared<Epetra_FECrsMatrix>(Copy, *graph));        /// creating A_ with graph
      rhs(std::make_shared<Epetra_FEVector>(*(maps.dof_map())));         /// creating b_ with dof_map  
    }



    //-------------------------------------------------------------------------------------------------------------------------------------        
    /// Deleting the copy and move constructors - no duplication/transfer in anyway
    linear_system(const linear_system&) = delete;               //copy constructor
    linear_system& operator=(const linear_system&) = delete;    //copy assignment operator
    linear_system(linear_system&&) = delete;                    //move constructor  
    linear_system& operator=(linear_system&&) = delete;         //move assignment operator 
    //-------------------------------------------------------------------------------------------------------------------------------------




    void matrix(std::shared_ptr<Epetra_FECrsMatrix> m){A_ = m;}
    auto matrix() {return A_;}
    const auto matrix()const {return A_;}

    void rhs(std::shared_ptr<Epetra_FEVector> v){b_ = v;}
    auto rhs() {return b_;}
    const auto rhs()const {return b_;}
    
        
    void clear()                              /// This function fills A_ and b_ with zeros
    {
       matrix()->PutScalar(0.0);
       rhs()->PutScalar(0.0);
    }
    
    void assembled()
    {
       matrix()->GlobalAssemble();    
       rhs()->GlobalAssemble(); 
    }
    
    private:
    std::shared_ptr<Epetra_FECrsMatrix> A_ = nullptr;
    std::shared_ptr<Epetra_FEVector> b_ = nullptr;
  };

  
}
#endif
