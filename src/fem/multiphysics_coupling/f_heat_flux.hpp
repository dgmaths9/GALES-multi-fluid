#ifndef F_HEAT_FLUX_HPP
#define F_HEAT_FLUX_HPP



#include <numeric>
#include <algorithm>
#include <vector>
#include <map>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/tuple/tuple.hpp>




namespace GALES{



  /**
       This class computes the fluid heat flux at fluid-solid inetrface to be passed to heat equation solver for coupling.
       The tractions are first computed at gauss nodes of each side and then averaged on the nodes.
  */





  template<int dim> class f_heat_flux{};



  template<>
  class f_heat_flux<2>
  {
    public:

    void heat_flux
    (
      model<2>& f_model,
      fluid_properties& props, 
      std::map<int, std::vector<double> >& f_nd_heat_flux
    )
    {
        std::vector<int> my_nd;
        std::vector<double> my_q1, my_q2;

        const int nb_comp = props.nb_comp();
        vec dofs_fluid;


        /// We first collect the nd gids and the heat flux vectors at gauss nodes of each interface side at each process
        for(const auto& el : f_model.mesh().bd_elements())
        {
          f_model.extract_element_dofs(*el, dofs_fluid);
          const int nb_dofs_f = dofs_fluid.size()/el->nb_nodes();
          vec I(nb_dofs_f,0.0);
          vec dI_dx(nb_dofs_f,0.0), dI_dy(nb_dofs_f,0.0);
          double dT_dx(0.0), dT_dy(0.0);
                
          for(int i=0; i<el->nb_sides(); i++)
          {
            if(el->is_side_on_boundary(i) && el->side_flag(i) == 1)
            {
                for(int j=0; j<el->side_gp(i).size(); j++)               
                {
                  auto quad_ptr = std::make_unique<quad>(*el, el->side_gp(i,j), el->side_gp(i,j));
                  
                  quad_ptr->interpolate(dofs_fluid, I);
                  quad_ptr->x_derivative_interpolate(dofs_fluid, dI_dx);
                  quad_ptr->y_derivative_interpolate(dofs_fluid, dI_dy);
                  const double p(I[0]);                  

                  const double shear_rate =  compute_shear_rate(dI_dx, dI_dy);     
                  if(nb_comp==1)
                  {
                     if(I.size()==3)
                     {
                         props.properties(p, shear_rate);               // fluid_sc_isothermal
                         // dT_dx = 0.0;   dT_dy = 0.0
                     }                         
                     else if (I.size()==4)
                     {
                         props.properties(p, I[I.size()-1], shear_rate);   // fluid_sc
                         dT_dx = dI_dx[I.size()-1];
                         dT_dy = dI_dy[I.size()-1];
                     }   
                  }
                  else if(nb_comp>1)
                  {
                    boost::numeric::ublas::vector<double> Y(nb_comp, 0.0);
                    for(int i=0; i<nb_comp-1; i++)
                      Y[i] = I[I.size()-(nb_comp-1) +i];
                    Y[nb_comp-1] = 1.0 - std::accumulate(&Y[0], &Y[nb_comp-1], 0.0);
                    
                    if(I.size()-(nb_comp-1)==3)
                    {
                       props.properties(p, Y, shear_rate);           //fluid_mc_isothermal
                       // dT_dx = 0.0;   dT_dy = 0.0
                    } 
                    else if (I.size()-(nb_comp-1)==4) 
                    {
                       props.properties(p, I[I.size()-(nb_comp-1)-1], Y, shear_rate);  // fluid_mc
                       dT_dx = dI_dx[I.size()-(nb_comp-1)-1];
                       dT_dy = dI_dy[I.size()-(nb_comp-1)-1];
                    }
                  }                 
                  
                  const auto kappa = props.kappa_;
                  const auto JN = quad_ptr->JN();

                  my_nd.push_back(el->side_nodes(i,j));
                  my_q1.push_back(-kappa*dT_dx*JN[0]);
                  my_q2.push_back(-kappa*dT_dy*JN[1]);
                }
            }
          }
        }

        ///This function collects the nds and tractions from each pid and average over the repeated nodes, and fill f_nd_tr
        all_pid_avg_over_repeated_nodes(my_nd, my_q1, my_q2, f_nd_heat_flux);
    }

 };










  template<>
  class f_heat_flux<3>
  {
    public:

    void heat_flux
    (
      model<3>& f_model,
      fluid_properties& props, 
      std::map<int, std::vector<double> >& f_nd_heat_flux
    )
    {
        std::vector<int> my_nd;
        std::vector<double> my_q1, my_q2, my_q3;

        const int nb_comp = props.nb_comp();
        vec dofs_fluid;


        /// We first collect the nd gids and the traction vectors at gauss nodes of each fsi side at each process
        for(const auto& el : f_model.mesh().bd_elements())
        {
          f_model.extract_element_dofs(*el, dofs_fluid);
          const int nb_dofs_f = dofs_fluid.size()/el->nb_nodes();
          vec I(nb_dofs_f,0.0);
          vec dI_dx(nb_dofs_f,0.0), dI_dy(nb_dofs_f,0.0), dI_dz(nb_dofs_f,0.0);
          double dT_dx(0.0), dT_dy(0.0), dT_dz(0.0);
                
          for(int i=0; i<el->nb_sides(); i++)
          {
            if(el->is_side_on_boundary(i) && el->side_flag(i) == 1)
            {
                for(int j=0; j<el->side_gp(i).size(); j++)               
                {
                  auto quad_ptr = std::make_unique<quad>(*el, el->side_gp(i,j), el->side_gp(i,j));

                  quad_ptr->interpolate(dofs_fluid, I);
                  quad_ptr->x_derivative_interpolate(dofs_fluid, dI_dx);
                  quad_ptr->y_derivative_interpolate(dofs_fluid, dI_dy);
                  quad_ptr->z_derivative_interpolate(dofs_fluid, dI_dz);
                  const double p(I[0]);
                  
                  const double shear_rate =  compute_shear_rate(dI_dx, dI_dy, dI_dz);     
                  if(nb_comp==1)
                  {
                     if(I.size()==4)
                     {
                         props.properties(p, shear_rate);               // fluid_sc_isothermal
                         // dT_dx = 0.0;   dT_dy = 0.0;   dT_dz = 0.0
                     }                         
                     else if (I.size()==5)
                     {
                         props.properties(p, I[I.size()-1], shear_rate);   // fluid_sc
                         dT_dx = dI_dx[I.size()-1];
                         dT_dy = dI_dy[I.size()-1];
                         dT_dz = dI_dz[I.size()-1];
                     }   
                  }
                  else if(nb_comp>1)
                  {
                    boost::numeric::ublas::vector<double> Y(nb_comp, 0.0);
                    for(int i=0; i<nb_comp-1; i++)
                      Y[i] = I[I.size()-(nb_comp-1) +i];
                    Y[nb_comp-1] = 1.0 - std::accumulate(&Y[0], &Y[nb_comp-1], 0.0);
                    
                    if(I.size()-(nb_comp-1)==4)
                    {
                       props.properties(p, Y, shear_rate);           //fluid_mc_isothermal
                       // dT_dx = 0.0;   dT_dy = 0.0;   dT_dz = 0.0
                    } 
                    else if (I.size()-(nb_comp-1)==5) 
                    {
                       props.properties(p, I[I.size()-(nb_comp-1)-1], Y, shear_rate);  // fluid_mc
                       dT_dx = dI_dx[I.size()-(nb_comp-1)-1];
                       dT_dy = dI_dy[I.size()-(nb_comp-1)-1];
                       dT_dz = dI_dz[I.size()-(nb_comp-1)-1];
                    }
                  }                 
                  
                  const auto kappa = props.kappa_;
                  const auto JN = quad_ptr->JN();

                  my_nd.push_back(el->side_nodes(i,j));
                  my_q1.push_back(-kappa*dT_dx*JN[0]);
                  my_q2.push_back(-kappa*dT_dy*JN[1]);
                  my_q3.push_back(-kappa*dT_dz*JN[2]);
                }
            }
          }
        }

        ///This function collects the nds and tractions from each pid and average over the repeated nodes, and fill f_nd_tr
        all_pid_avg_over_repeated_nodes(my_nd, my_q1, my_q2, my_q3, f_nd_heat_flux);
    }

 };



}

#endif
