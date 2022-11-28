#ifndef FLUID_ASSEMBLY_HPP_
#define FLUID_ASSEMBLY_HPP_


#include "../../../src/fem/fem.hpp"



namespace GALES {





  template<typename integral_type, typename scatter_type, int dim>
  class fluid_assembly
  {
      using model_type = model<dim>;
      using vec = boost::numeric::ublas::vector<double>;
      using mat = boost::numeric::ublas::matrix<double>;

    public:
    
    fluid_assembly
    (
      integral_type& integral,
      scatter_type& sct,
      model_type& f_model,
      linear_system& lp      
    ):
      integral_(integral),
      scatter_(sct),
      f_model_(f_model),
      lp_(lp)
    {}

        
        
    double execute()
    {
      double start(MPI_Wtime());      
      
      const auto& mesh( f_model_.mesh() );
      std::vector<vec>dofs_fluid(f_model_.state().num_slot());

      lp_.clear();
      for(const auto& el: mesh.elements())
      {
	f_model_.extract_element_dofs(*el, dofs_fluid);

        mat m(dofs_fluid[0].size(), dofs_fluid[0].size(), 0.0);
 	vec r(dofs_fluid[0].size(), 0.0);   

	integral_.execute_i(*el, dofs_fluid, m, r);
	scatter_.execute(*el, mesh, lp_, m, r);

	if(el->on_boundary())
	{
	  r.clear();
	  m.clear();
 	  integral_.execute_b(*el, dofs_fluid, m, r);
	  scatter_.execute(*el, mesh, lp_, m, r);	  
	}
      }
      
      lp_.assembled();       
      return MPI_Wtime()-start;                  
    }
    
  

   private:
    integral_type& integral_;
    scatter_type& scatter_;
    model_type& f_model_;
    linear_system& lp_;    
  };




}//namespace GALES
#endif
