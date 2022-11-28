#ifndef _GALES_ICO_F_DRAG_LIFT_HPP_
#define _GALES_ICO_F_DRAG_LIFT_HPP_



#include <numeric>
#include <algorithm>
#include <vector>
#include <map>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/tuple/tuple.hpp>




namespace GALES {


  /**
      This file computed drag and list coefficients for ico solver  
  */



  template<int dim> class ico_f_drag_lift{};


  template<>
  class ico_f_drag_lift<2>
  {
    using vec = boost::numeric::ublas::vector<double>;
    point<2> dummy_;   
    
    public:
    
    void drag_lift
    (
      model<2>& f_model,
      fluid_properties& props, 
      read_setup& setup,
      double& C_d,
      double& C_l
    )
    {        
        const auto& mesh(f_model.mesh());
        const double nb_el_nodes = f_model.mesh().nb_el_nodes();
        vec dofs_fluid;
        
        vec I(4,0.0);
        vec dI_dx(4,0.0), dI_dy(4,0.0);
        vec ux_e(nb_el_nodes,0.0), uy_e(nb_el_nodes,0.0);

        double fx(0.0), fy(0.0);
        
        for(const auto& el : mesh.bd_elements())
        {
          f_model.extract_element_dofs(*el, dofs_fluid);
          
          for(int i=0; i<el->nb_sides(); i++)
          {
            if(el->is_side_on_boundary(i))
            {
              auto side_flag = el->side_flag(i);
              
              if(side_flag==setup.pp_bd_flag())
              {
                for(const auto& quad_ptr: el->quad_b(i))
                {                  
                  quad_ptr->interpolate(dofs_fluid, I);
                  quad_ptr->x_derivative_interpolate(dofs_fluid, dI_dx);
                  quad_ptr->y_derivative_interpolate(dofs_fluid, dI_dy);

                  const double dv1_1(dI_dx[1]), dv2_1(dI_dx[2]), dv1_2(dI_dy[1]), dv2_2(dI_dy[2]);
                  const double shear_rate =  sqrt(std::pow(dv1_1,2) + std::pow(dv2_2,2) + 2.0*std::pow(dv1_2 + dv2_1,2));
                  const double p = I[0];
                  const double T = I[3];
                                    
                  props.properties(p, T, shear_rate);
                  const double mu = props.mu_;
                  const double lambda = -2.0/3*mu; 
                   
                  const double sigma11 = lambda*(dv1_1 + dv2_2) + 2*mu*dv1_1 - p;
                  const double sigma22 = lambda*(dv1_1 + dv2_2) + 2*mu*dv2_2 - p;
                  const double sigma12 = mu*(dv1_2 + dv2_1);
                                    
                  const auto JN = quad_ptr->JN();
                  fx += sigma11*JN[0] + sigma12*JN[1];
                  fy += sigma12*JN[0] + sigma22*JN[1];                  
                }
              }  
            }                
          }
        }
        
        double Fx, Fy;
        MPI_Reduce(&fx, &Fx, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(&fy, &Fy, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        
        const double rho_inf = setup.pp_rho_inf();
        const double v_inf = setup.pp_v_inf();
        const double D = setup.pp_cylinder_diameter();
                
        C_d = Fx/(0.5*rho_inf*v_inf*v_inf*D);
        C_l = Fy/(0.5*rho_inf*v_inf*v_inf*D);
        
    }    


  };
 
 
 
 
  
}

#endif
