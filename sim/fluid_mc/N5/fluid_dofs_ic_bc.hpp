#ifndef _GALES_FLUID_DOFS_IC_BC_HPP_
#define _GALES_FLUID_DOFS_IC_BC_HPP_


#include "../../../src/fem/fem.hpp"



namespace GALES {


  template<typename ic_bc_type, int dim>
  class fluid_dofs_ic_bc{};






  template<typename ic_bc_type>
  class fluid_dofs_ic_bc <ic_bc_type, 2>
  {
      using model_type = model<2>;
  
      public:
      
      fluid_dofs_ic_bc(model_type& f_model, ic_bc_type& ic_bc, int nb_comp)
      :  
      state_(f_model.state()), mesh_(f_model.mesh()), ic_bc_(ic_bc), nb_comp_(nb_comp)
      {
        ic();
      }



    void ic()
    {
      vec Y(nb_comp_-1,0.0);
      for(const auto& nd : mesh_.nodes())
      {
        const int first_dof_lid (nd->first_dof_lid());
        state_.set_dof(first_dof_lid, ic_bc_.initial_p(*nd));
        state_.set_dof(first_dof_lid+1, ic_bc_.initial_vx(*nd));
        state_.set_dof(first_dof_lid+2, ic_bc_.initial_vy(*nd));
	state_.set_dof(first_dof_lid+3, ic_bc_.initial_T(*nd));

        ic_bc_.initial_Y(*nd,Y);
        for(int i=0; i<nb_comp_-1; i++)
          state_.set_dof(first_dof_lid+4+i, Y[i]);
      }
    }


      void dirichlet_bc()
      {
        std::pair<bool,double> result;

        for(const auto& nd : mesh_.nodes())
        {
          const int first_dof_lid (nd->first_dof_lid());

          result= ic_bc_.dirichlet_p(*nd);
          if (result.first)          state_.set_dof(first_dof_lid, result.second);

          result= ic_bc_.dirichlet_vx(*nd);
          if (result.first)          state_.set_dof(first_dof_lid+1, result.second);

          result= ic_bc_.dirichlet_vy(*nd);
          if (result.first)          state_.set_dof(first_dof_lid+2, result.second);

          result= ic_bc_.dirichlet_T(*nd);
          if (result.first)          state_.set_dof(first_dof_lid+3, result.second);

          for(int i=0; i<nb_comp_-1; i++)
          {
            result = ic_bc_.dirichlet_Y(*nd, i);
            if(result.first)         state_.set_dof(first_dof_lid+4+i, result.second);
          }    
        }
      }




    template<typename nd_type>
    auto dof_constraint(int dof, const nd_type &nd)const  
    {
      const int nb_dofs = 4 + nb_comp_-1;
      const int n = dof%nb_dofs;
      std::pair<bool,double> result(false, 0.0);
      std::pair<bool,double> value(false, 0.0);

      switch(n)
      {
          case 0:
              result = ic_bc_.dirichlet_p(nd);
              if(result.first == true) value = std::make_pair(true, 0.0);
              break;
          case 1:
              result = ic_bc_.dirichlet_vx(nd);
              if(result.first == true) value = std::make_pair(true, 0.0);
              break;
          case 2:
              result = ic_bc_.dirichlet_vy(nd);
              if(result.first == true) value = std::make_pair(true, 0.0);
              break;
          case 3:
              result = ic_bc_.dirichlet_T(nd);
              if(result.first == true) value = std::make_pair(true, 0.0);
              break;
          default:
              result = ic_bc_.dirichlet_Y(nd,n%4);
              if(result.first == true) value = std::make_pair(true, 0.0);
              break;
      }
      return value;
    }
 

    private:
    dof_state& state_;
    Mesh<2>& mesh_;
    ic_bc_type& ic_bc_;
    int nb_comp_;

  };









  template<typename ic_bc_type>
  class fluid_dofs_ic_bc <ic_bc_type, 3>
  {
      using model_type = model<3>;
  
      public:
      
      fluid_dofs_ic_bc(model_type& f_model, ic_bc_type& ic_bc, int nb_comp)
      :  
      state_(f_model.state()), mesh_(f_model.mesh()), ic_bc_(ic_bc), nb_comp_(nb_comp)
      {
        ic();
      }


    void ic()
    {
      vec Y(nb_comp_-1,0.0);
      for(const auto& nd : mesh_.nodes())
      {
        const int first_dof_lid (nd->first_dof_lid());
        state_.set_dof(first_dof_lid, ic_bc_.initial_p(*nd));
        state_.set_dof(first_dof_lid+1, ic_bc_.initial_vx(*nd));
        state_.set_dof(first_dof_lid+2, ic_bc_.initial_vy(*nd));
        state_.set_dof(first_dof_lid+3, ic_bc_.initial_vz(*nd));
	state_.set_dof(first_dof_lid+4, ic_bc_.initial_T(*nd));

        ic_bc_.initial_Y(*nd,Y);
        for(int i=0; i<nb_comp_-1; i++)
          state_.set_dof(first_dof_lid+5+i, Y[i]);
      }
    }

   
      void dirichlet_bc()
      {
        std::pair<bool,double> result;

        for(const auto& nd : mesh_.nodes())
        {
          const int first_dof_lid (nd->first_dof_lid());

          result= ic_bc_.dirichlet_p(*nd);
          if (result.first)          state_.set_dof(first_dof_lid, result.second);

          result= ic_bc_.dirichlet_vx(*nd);
          if (result.first)          state_.set_dof(first_dof_lid+1, result.second);

          result= ic_bc_.dirichlet_vy(*nd);
          if (result.first)          state_.set_dof(first_dof_lid+2, result.second);

          result= ic_bc_.dirichlet_vz(*nd);
          if (result.first)          state_.set_dof(first_dof_lid+3, result.second);

          result= ic_bc_.dirichlet_T(*nd);
          if (result.first)          state_.set_dof(first_dof_lid+4, result.second);

          for(int i=0; i<nb_comp_-1; i++)
          {
            result = ic_bc_.dirichlet_Y(*nd, i);
            if(result.first)         state_.set_dof(first_dof_lid+5+i, result.second);
          }    
        }
      }



    template<typename nd_type>
    auto dof_constraint(int dof, const nd_type &nd)const  
    {
      const int nb_dofs = 5 + nb_comp_-1;
      const int n = dof%nb_dofs;
      std::pair<bool,double> result(false, 0.0);
      std::pair<bool,double> value(false, 0.0);

      switch(n)
      {
          case 0:
              result = ic_bc_.dirichlet_p(nd);
              if(result.first == true) value = std::make_pair(true, 0.0);
              break;
          case 1:
              result = ic_bc_.dirichlet_vx(nd);
              if(result.first == true) value = std::make_pair(true, 0.0);
              break;
          case 2:
              result = ic_bc_.dirichlet_vy(nd);
              if(result.first == true) value = std::make_pair(true, 0.0);
              break;
          case 3:
              result = ic_bc_.dirichlet_vz(nd);
              if(result.first == true) value = std::make_pair(true, 0.0);
              break;
          case 4:
              result = ic_bc_.dirichlet_T(nd);
              if(result.first == true) value = std::make_pair(true, 0.0);
              break;
          default:
              result = ic_bc_.dirichlet_Y(nd,n%5);
              if(result.first == true) value = std::make_pair(true, 0.0);
      }
      return value;
    }



    private:
    dof_state& state_;
    Mesh<3>& mesh_;
    ic_bc_type& ic_bc_;
    int nb_comp_;

  };


} /* namespace GALES */

#endif

