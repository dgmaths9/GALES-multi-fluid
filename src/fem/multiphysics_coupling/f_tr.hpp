#ifndef F_TR_HPP
#define F_TR_HPP



#include <numeric>
#include <algorithm>
#include <vector>
#include <map>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/tuple/tuple.hpp>




namespace GALES{



  /**
       This class computes the fluid traction forces at fluid-solid inetrface to be passed to solid solver for FSI coupling.
       The tractions are first computed at gauss nodes of each side and then averaged on the nodes.
  */





  template<int dim> class f_tr{};


  template<>
  class f_tr<2>
  {
    public:

    void tr
    (
      model<2>& f_model,
      fluid_properties& props, 
      std::map<int, std::vector<double> >& f_nd_tr
    )
    {
        std::vector<int> my_nd;
        std::vector<double> my_tx, my_ty;

        const int nb_comp = props.nb_comp();
        vec dofs_fluid;


        /// We first collect the nd gids and the traction vectors at gauss nodes of each fsi side at each process
        for(const auto& el : f_model.mesh().bd_elements())
        {
          f_model.extract_element_dofs(*el, dofs_fluid);
          const int nb_dofs_f = dofs_fluid.size()/el->nb_nodes();
          vec I(nb_dofs_f,0.0);
          vec dI_dx(nb_dofs_f,0.0), dI_dy(nb_dofs_f,0.0);
                
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
                     if(I.size()==3)   props.properties(p, shear_rate);               // fluid_sc_isothermal                      
                     else if (I.size()==4)   props.properties(p, I[I.size()-1], shear_rate);   // fluid_sc
                  }
                  else if(nb_comp>1)
                  {
                    boost::numeric::ublas::vector<double> Y(nb_comp, 0.0);
                    for(int i=0; i<nb_comp-1; i++)
                      Y[i] = I[I.size()-(nb_comp-1) +i];
                    Y[nb_comp-1] = 1.0 - std::accumulate(&Y[0], &Y[nb_comp-1], 0.0);
                    
                    if(I.size()-(nb_comp-1)==3) props.properties(p, Y, shear_rate);           //fluid_mc_isothermal
                    else if (I.size()-(nb_comp-1)==4) props.properties(p, I[I.size()-(nb_comp-1)-1], Y, shear_rate);  // fluid_mc
                  }                 
                  
                  const double mu = props.mu_;
                  const double lambda = -2.0/3*mu;
             
                  const double sigma11 = lambda*(dI_dx[1] + dI_dy[2]) + 2*mu*dI_dx[1] - p;
                  const double sigma22 = lambda*(dI_dx[1] + dI_dy[2]) + 2*mu*dI_dy[2] - p;
                  const double sigma12 = mu*(dI_dy[1] + dI_dx[2]);     

                  const auto JN = quad_ptr->JN();

                  my_nd.push_back(el->side_nodes(i,j));
                  my_tx.push_back(sigma11*JN[0] + sigma12*JN[1]);
                  my_ty.push_back(sigma12*JN[0] + sigma22*JN[1]);
                }
            }
          }
        }

        ///This function collects the nds and tractions from each pid and average over the repeated nodes, and fill f_nd_tr
        all_pid_avg_over_repeated_nodes(my_nd, my_tx, my_ty, f_nd_tr);
    }

 };









  template<>
  class f_tr<3>
  {
    public:

    void tr
    (
      model<3>& f_model,
      fluid_properties& props, 
      std::map<int, std::vector<double> >& f_nd_tr
    )
    {
        std::vector<int> my_nd;
        std::vector<double> my_tx, my_ty, my_tz;

        const int nb_comp = props.nb_comp();
        vec dofs_fluid;


        /// We first collect the nd gids and the traction vectors at gauss nodes of each fsi side at each process
        for(const auto& el : f_model.mesh().bd_elements())
        {
          f_model.extract_element_dofs(*el, dofs_fluid);
          const int nb_dofs_f = dofs_fluid.size()/el->nb_nodes();
          vec I(nb_dofs_f,0.0);
          vec dI_dx(nb_dofs_f,0.0), dI_dy(nb_dofs_f,0.0), dI_dz(nb_dofs_f,0.0);
                
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
                     if(I.size()==4)   props.properties(p, shear_rate);               // fluid_sc_isothermal                      
                     else if (I.size()==5)   props.properties(p, I[I.size()-1], shear_rate);   // fluid_sc
                  }
                  else if(nb_comp>1)
                  {
                    boost::numeric::ublas::vector<double> Y(nb_comp, 0.0);
                    for(int i=0; i<nb_comp-1; i++)
                      Y[i] = I[I.size()-(nb_comp-1) +i];
                    Y[nb_comp-1] = 1.0 - std::accumulate(&Y[0], &Y[nb_comp-1], 0.0);
                    
                    if(I.size()-(nb_comp-1)==4) props.properties(p, Y, shear_rate);           //fluid_mc_isothermal
                    else if (I.size()-(nb_comp-1)==5) props.properties(p, I[I.size()-(nb_comp-1)-1], Y, shear_rate);  // fluid_mc
                  }                 
                  
                  const double mu = props.mu_;
                  const double lambda = -2.0/3*mu;
             
                  const double sigma11 = lambda*(dI_dx[1] + dI_dy[2] + dI_dz[3]) + 2*mu*dI_dx[1] - p;
                  const double sigma22 = lambda*(dI_dx[1] + dI_dy[2] + dI_dz[3]) + 2*mu*dI_dy[2] - p;
                  const double sigma33 = lambda*(dI_dx[1] + dI_dy[2] + dI_dz[3]) + 2*mu*dI_dz[3] - p;     
                  const double sigma12 = mu*(dI_dy[1] + dI_dx[2]);
                  const double sigma23 = mu*(dI_dz[2] + dI_dy[3]);
                  const double sigma13 = mu*(dI_dx[3] + dI_dz[1]);   

                  const auto JN = quad_ptr->JN();

                  my_nd.push_back(el->side_nodes(i,j));
                  my_tx.push_back(sigma11*JN[0] + sigma12*JN[1] + sigma13*JN[2]);
                  my_ty.push_back(sigma12*JN[0] + sigma22*JN[1] + sigma23*JN[2]);
                  my_tz.push_back(sigma13*JN[0] + sigma23*JN[1] + sigma33*JN[2]);
                }
            }
          }
        }

        ///This function collects the nds and tractions from each pid and average over the repeated nodes, and fill f_nd_tr
        all_pid_avg_over_repeated_nodes(my_nd, my_tx, my_ty, my_tz, f_nd_tr);
    }

 };



}

#endif
