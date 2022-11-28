#ifndef __FLUID_IC_BC_HPP
#define __FLUID_IC_BC_HPP



#include "../../../src/fem/fem.hpp"



namespace GALES {




  template<int dim>
  class fluid_ic_bc 
  {

    using nd_type = node<dim>;
    using point_type = point<dim>;
    using vec = boost::numeric::ublas::vector<double>;


  public :
 
    fluid_ic_bc():
      automatic_correction_(false),
      central_pressure_profile_(true),
      x_dependent_pressure_profile_(false),
      p0_(70.e6),
      Y_min_(0.0), 
      Y_max_(1.0-Y_min_),      
      p_ref_(p0_), v1_ref_(0.0), v2_ref_(0.0), v3_ref_(0.0), T_ref_(1300.0)   // these are for dc_2006      
      {
        max_ = -3000.0;
        min_ = -9000.0;
        h_ = 1.0;
        total_steps_ = (int)((max_-min_)/h_);     
      }




    void body_force(vec& gravity, const point<2>& position)
    {
      gravity[0] = 0.0;
      gravity[1] = -9.81;
    }


    void body_force(vec& gravity,  const point<3>& position)
    {
      gravity[0] = 0.0;
      gravity[1] = 0.0;
      gravity[2] = 0.0;
    }




   
    //-------------------- IC ----------------------------------------

    double initial_p(const nd_type &nd) const { return p0_; }
    double initial_vx(const nd_type &nd) const { return 0.0; }
    double initial_vy(const nd_type &nd) const { return 0.0; }
    double initial_vz(const nd_type &nd) const { return 0.0; }
    double initial_T(const nd_type &nd)const  { return 1300.0; } 

    void initial_Y(const nd_type &nd, vec &Y)const 
    {
       double y = nd.get_y();
       if (y <= -3420.0)  Y[0] = Y_max_;
       else   Y[0] = Y_min_;
    }
 	





  // --------------------Dirichlet BC------------------------------------

    auto dirichlet_p(const nd_type &nd) const 
    {
      return std::make_pair(false,0.0);
    }


    auto dirichlet_vx(const nd_type &nd) const 
    { 
      if( nd.flag() ==  5)     return std::make_pair(true,0.0); 
	
      return std::make_pair(false,0.0);
    }


    auto dirichlet_vy(const nd_type &nd) 
    {
	if( nd.flag() ==  5)     return std::make_pair(true,0.0); 
	  
      return std::make_pair(false,0.0);
    }


    auto dirichlet_vz(const nd_type &nd) 
    {
      return std::make_pair(false,0.0);
    }


    auto dirichlet_T(const nd_type &nd) const 
    {
      return std::make_pair(true, 1300.0);
    }
 

    // comp_index is the index of component. e.g. for 3 component mixture we have first 2 components as dofs; Y[0] and Y[1]
    // 0---first component;   1----second component
    auto dirichlet_Y(const nd_type &nd, int comp_index)const
    {
//       if( nd.flag() == 2 && comp_index == 0 )          return std::make_pair(true,0.01);

      return std::make_pair(false,0.0);
    }





  //-------------------------Neumann BC--------------------------------------------------

    auto neumann_tau11(const std::vector<int>& bd_nodes, int side_flag) 
    { 
      return std::make_pair(false,0.0);  
    }

    auto neumann_tau12(const std::vector<int>& bd_nodes, int side_flag) 
    { 
      return std::make_pair(false,0.0);
    }

    auto neumann_tau13(const std::vector<int>& bd_nodes, int side_flag) 
    {  
      return std::make_pair(false,0.0);
    }
  
    auto neumann_tau22(const std::vector<int>& bd_nodes, int side_flag) 
    { 
      return std::make_pair(false,0.0); 
    }
                
    auto neumann_tau23(const std::vector<int>& bd_nodes, int side_flag) 
    { 
      return std::make_pair(false,0.0);
    }

    auto neumann_tau33(const std::vector<int>& bd_nodes, int side_flag) 
    {  
      return std::make_pair(false,0.0); 
    }
                
    auto neumann_q1(const std::vector<int>& bd_nodes, int side_flag)  
    {   
      if(side_flag == 5)     return std::make_pair(true,0.0);

      return std::make_pair(false,0.0); 
    }

    auto neumann_q2(const std::vector<int>& bd_nodes, int side_flag) 
    { 
      if(side_flag == 5)     return std::make_pair(true,0.0);

      return std::make_pair(false,0.0); 
    }

    auto neumann_q3(const std::vector<int>& bd_nodes, int side_flag) 
    { 
      return std::make_pair(false,0.0); 
    }

    auto neumann_J1(int comp_index, const std::vector<int>& bd_nodes, int side_flag) 
    { 
      if(side_flag == 5 && comp_index==0)     return std::make_pair(true,0.0);

      return std::make_pair(false,0.0); 
    }

    auto neumann_J2(int comp_index, const std::vector<int>& bd_nodes, int side_flag) 
    { 
      if(side_flag == 5 && comp_index==0)     return std::make_pair(true,0.0);

      return std::make_pair(false,0.0); 
    }

    auto neumann_J3(int comp_index, const std::vector<int>& bd_nodes, int side_flag) 
    { 
      return std::make_pair(false,0.0); 
    }

    // rho*v1 towards outwards normal 
    auto mass_flux1(const std::vector<int>& bd_nodes, int side_flag)
    {
      return std::make_pair(false,0.0); 
    }

    // rho*v2 towards outwards normal
    auto mass_flux2(const std::vector<int>& bd_nodes, int side_flag)
    {
      return std::make_pair(false,0.0); 
    }

    // rho*v3 towards outwards normal
    auto mass_flux3(const std::vector<int>& bd_nodes, int side_flag)
    {
      return std::make_pair(false,0.0); 
    }


  public:
    bool automatic_correction_;
    bool central_pressure_profile_;
    bool x_dependent_pressure_profile_;
    double p0_;
    double Y_min_, Y_max_;    
    double p_ref_, v1_ref_, v2_ref_, v3_ref_, T_ref_;
    double max_, min_;
    int total_steps_;
    double h_;


  };
  
} //namespace
#endif



