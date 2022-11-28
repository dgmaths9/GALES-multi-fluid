#ifndef _GALES_SCATTER_HPP_
#define _GALES_SCATTER_HPP_


#include <vector>




namespace GALES {


  /**
      This class scatters the element level computations into global matrix and vector.  
  */


  template<typename dofs_ic_bc_type>
  class scatter
  {

    public:
    
    explicit scatter(dofs_ic_bc_type& dofs_ic_bc):  dofs_ic_bc_(dofs_ic_bc){}    


    //-------------------------------------------------------------------------------------------------------------------------------------        
    /// Deleting the copy and move constructors - no duplication/transfer in anyway
    scatter(const scatter&) = delete;               //copy constructor
    scatter& operator=(const scatter&) = delete;    //copy assignment operator
    scatter(scatter&&) = delete;                    //move constructor  
    scatter& operator=(scatter&&) = delete;         //move assignment operator 
    //-------------------------------------------------------------------------------------------------------------------------------------


    template <typename element_type, typename mesh_type, typename lp_type, typename m_type, typename v_type>
    void execute(const element_type& el, const mesh_type& mesh, lp_type& lp, m_type& m, v_type& r) const 
    { 
      const int n_item(m.size1());                          /// number of dofs in local matrix
      std::vector<int> dof_gid_index(n_item);
      std::vector<std::pair<bool,double>> to_be_blocked(n_item, std::make_pair(false,0.0));


      /// creating the dof_gid_index and dof constrained info for the dofs belonging to the el 
      int count(0);
      for (int i=0; i<el.nb_nodes(); i++) 
      {
	const int nd_lid = el.nd_lid(i);
	const auto& nd(mesh.nodes()[nd_lid]);
	const int fd_gid = nd->first_dof_gid();
	for(int k_i=0; k_i<nd->nb_dofs(); k_i++)
	{
	  dof_gid_index[count] = fd_gid + k_i;
          to_be_blocked[count] =  dofs_ic_bc_.dof_constraint(dof_gid_index[count], *nd) ;          
	  count++;
	}
      }




     /// Applying dirichlet bc at the local matrix.
     /// The diagonal value in matrix is not 1.0 instead is set to the average value of the matrix diagonal for better condition number.

     double avg_val(1.0);
     for(int i=0; i<n_item; ++i)
       avg_val += m(i,i);
     avg_val /= n_item;

     for(int i=0; i<n_item; ++i)   // i is for row
       if(to_be_blocked[i].first)
       { 
         for(int j=0; j<n_item; ++j)   // j is for column
           m(i,j) = 0.0;                  
           
         m(i,i) = avg_val;     
         r[i] = avg_val*to_be_blocked[i].second;
       }




     /// To pass matrix m into lp.matrix() with just one call to SumIntoGlobalValues function, we convert m(2d structure) to a m_vec(1d structure);
     std::vector<double> m_vec(n_item*n_item);
     for(int i=0; i<n_item; ++i)  
     {
      for(int j=0; j<n_item; ++j)  
      {
        m_vec[i*n_item + j] = m(i,j);
      }
     }   



      /// Here we fill the global matrix and RHS vector.
      lp.rhs()->SumIntoGlobalValues(n_item, &dof_gid_index[0], &r[0]);    
      lp.matrix()->SumIntoGlobalValues(n_item, &dof_gid_index[0], n_item, &dof_gid_index[0], &m_vec[0], Epetra_FECrsMatrix::ROW_MAJOR);           
    }


      private:
      dofs_ic_bc_type& dofs_ic_bc_;
  };

  }//namespace GALES
#endif


