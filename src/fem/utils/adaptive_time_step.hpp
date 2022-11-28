#ifndef ADAPTIVE_TIME_HPP
#define ADAPTIVE_TIME_HPP



#include <vector>
#include "boost/numeric/ublas/vector.hpp"


namespace GALES{


        
    /**
        This file defines the adaptive time step criterion for fluid flow solvers
    */
    
    
        
    template<int dim>
    struct adaptive_time_criteria
    {

      using vec = boost::numeric::ublas::vector<double>;
      using mat = boost::numeric::ublas::matrix<double>;

     //------------------------------- criterion for sc ----------------------------------------------------------------      
      
      void f_sc(model<dim>& f_model, fluid_properties& props, read_setup& setup, int nb_dofs)
      {
         const double CFL = setup.CFL();
         const int nb_el_nodes = f_model.mesh().nb_el_nodes();
         
         const auto& mesh(f_model.mesh());
         vec dofs_fluid;
         point<dim> pt;
         vec I(nb_dofs, 0.0), dI_dx(nb_dofs, 0.0), dI_dy(nb_dofs, 0.0), dI_dz(nb_dofs, 0.0);
         std::vector<double> dt_vec;
         const double tol = std::numeric_limits<double>::epsilon();
         double shear_rate(0.0);           
         
         for(const auto& el : mesh.elements())         
         {
           f_model.extract_element_dofs(*el, dofs_fluid);

           std::vector<double> dt_gp_vec;
           dt_gp_vec.reserve(el->nb_gp());
           for(auto& quad_ptr : el->quad_i())       // here we loop over each gauss point of the element and compute dt        
           {             
             quad_ptr->interpolate(dofs_fluid, I);
             quad_ptr->x_derivative_interpolate(dofs_fluid, dI_dx);
             quad_ptr->y_derivative_interpolate(dofs_fluid, dI_dy);
             if(dim==3) quad_ptr->z_derivative_interpolate(dofs_fluid, dI_dz);
             auto del_v = compute_del_v(dI_dx, dI_dy, dI_dz);

             const double p = I[0];           
             vec v(dim);
             for(int i=0; i<dim; i++)
              v[i] = I[1+i];
             const double T = I[dim+1]; 
             
             props.properties(p,T,shear_rate);
             const double rho = props.rho_;
             const double mu = props.mu_;
             const double kappa = props.kappa_;
             const double c = props.sound_speed_;
             const double cv = props.cv_;
             const double v_nrm = boost::numeric::ublas::norm_2(v);
                        
             if(v_nrm < tol)
             {
               dt_gp_vec.push_back(1.2*time::get().delta_t());	
             }
             else
             {
               double h = quad_ptr->compute_h(v, del_v, pt);
               double v_tau = v_nrm;
               const double Mach = v_nrm/c;
               if(Mach > 0.3)
               {
                 if(dim==2) v_tau = sqrt(v_nrm*v_nrm + 1.5*c*c + c*sqrt(16.0*v_nrm*v_nrm + c*c));
                 else v_tau = sqrt(v_nrm*v_nrm + 2.0*c*c + c*sqrt(4.0*v_nrm*v_nrm + c*c));
               }
  
               const double use = 2.0*std::max(2.0*mu/rho, kappa/(rho*cv))/(h*h) + v_tau/h;                                                                 
               dt_gp_vec.push_back(CFL/use);
             }
           }
           auto it = min_element(std::begin(dt_gp_vec), std::end(dt_gp_vec));  
           dt_vec.push_back(*it);                  // here we compute the minimum dt on all gauss points of the element and put that in  dt_vec
         }                 
         set_dt(setup, dt_vec); 
      }
      
      //-----------------------------------------------------------------------------------------------      









     //------------------------------- criterion for sc_isothermal ----------------------------------------------------------------      
      
      void f_sc_isothermal(model<dim>& f_model, fluid_properties& props, read_setup& setup, int nb_dofs)
      {
         const double CFL = setup.CFL();
         const int nb_el_nodes = f_model.mesh().nb_el_nodes();
         
         const auto& mesh(f_model.mesh());
         vec dofs_fluid;
         point<dim> pt;
         vec I(nb_dofs, 0.0), dI_dx(nb_dofs, 0.0), dI_dy(nb_dofs, 0.0), dI_dz(nb_dofs, 0.0);
         std::vector<double> dt_vec;
         const double tol = std::numeric_limits<double>::epsilon();
         double shear_rate(0.0);           
         
         for(const auto& el : mesh.elements())         
         {
           f_model.extract_element_dofs(*el, dofs_fluid);

           std::vector<double> dt_gp_vec;
           dt_gp_vec.reserve(el->nb_gp());
           for(auto& quad_ptr : el->quad_i())       // here we loop over each gauss point of the element and compute dt        
           {             
             quad_ptr->interpolate(dofs_fluid, I);
             quad_ptr->x_derivative_interpolate(dofs_fluid, dI_dx);
             quad_ptr->y_derivative_interpolate(dofs_fluid, dI_dy);
             if(dim==3) quad_ptr->z_derivative_interpolate(dofs_fluid, dI_dz);
             auto del_v = compute_del_v(dI_dx, dI_dy, dI_dz);

             const double p = I[0];           
             vec v(dim);
             for(int i=0; i<dim; i++)
              v[i] = I[1+i];           
             
             props.properties(p,shear_rate);
             const double rho = props.rho_;
             const double mu = props.mu_;
             const double c = props.sound_speed_;
             const double v_nrm = boost::numeric::ublas::norm_2(v);
                        
             if(v_nrm < tol)
             {
               dt_gp_vec.push_back(1.2*time::get().delta_t());	
             }
             else
             {
               double h = quad_ptr->compute_h(v, del_v, pt);
               double v_tau = v_nrm;
               const double Mach = v_nrm/c;
               if(Mach > 0.3)
               {
                 if(dim==2) v_tau = sqrt(v_nrm*v_nrm + 1.5*c*c + c*sqrt(16.0*v_nrm*v_nrm + c*c));
                 else v_tau = sqrt(v_nrm*v_nrm + 2.0*c*c + c*sqrt(4.0*v_nrm*v_nrm + c*c));
               }
  
               const double use = 4.0*mu/(rho*h*h) + v_tau/h;                                                                 
               dt_gp_vec.push_back(CFL/use);
             }
           }
           auto it = min_element(std::begin(dt_gp_vec), std::end(dt_gp_vec));  
           dt_vec.push_back(*it);                  // here we compute the minimum dt on all gauss points of the element and put that in  dt_vec
         }                 
        
         set_dt(setup, dt_vec); 
      }
      
      //-----------------------------------------------------------------------------------------------      









     //------------------------------- criterion for MC ----------------------------------------------------------------      
      
      void f_mc
      (
         model<dim>& f_model, fluid_properties& props, read_setup& setup, int nb_dofs
      )
      {
         const double CFL = setup.CFL();
         const int nb_el_nodes = f_model.mesh().nb_el_nodes();
         
         const auto& mesh(f_model.mesh());
         vec dofs_fluid;
         point<dim> pt;
         vec I(nb_dofs), dI_dx(nb_dofs), dI_dy(nb_dofs), dI_dz(nb_dofs);
         std::vector<double> dt_vec;
         const double tol = std::numeric_limits<double>::epsilon();
         double shear_rate(0.0);           
         
         for(const auto& el : mesh.elements())         
         {
           f_model.extract_element_dofs(*el, dofs_fluid);

           std::vector<double> dt_gp_vec;
           dt_gp_vec.reserve(el->nb_gp());
           for(auto& quad_ptr : el->quad_i())       // here we loop over each gauss point of the element and compute dt        
           {             
             quad_ptr->interpolate(dofs_fluid, I);
             quad_ptr->x_derivative_interpolate(dofs_fluid, dI_dx);
             quad_ptr->y_derivative_interpolate(dofs_fluid, dI_dy);
             if(dim==3) quad_ptr->z_derivative_interpolate(dofs_fluid, dI_dz);
             auto del_v = compute_del_v(dI_dx, dI_dy, dI_dz);


             const double p = I[0];           
             vec v(dim, 0.0);
             for(int i=0; i<dim; i++)
               v[i] = I[1+i];
             const double T = I[dim+1]; 
             vec Y(props.nb_comp(), 0.0);
             for (int i=0; i<props.nb_comp()-1; i++)
               Y[i] = I[dim+2+i];
             Y[props.nb_comp()-1] = 1.0- std::accumulate(&Y[0],&Y[props.nb_comp()-1],0.0);           
             
             
              
             props.properties(p,T,Y,shear_rate);
             const double rho = props.rho_;
             const double mu = props.mu_;
             const double kappa = props.kappa_;
             const double c = props.sound_speed_;
             const double cv = props.cv_;
             const double v_nrm = boost::numeric::ublas::norm_2(v);           
  
  
  
             if(v_nrm < tol)
             {
               dt_gp_vec.push_back(1.2*time::get().delta_t());	
             }
             else
             {
               double h = quad_ptr->compute_h(v, del_v, pt);
               double v_tau = v_nrm;
               const double Mach = v_nrm/c;
               if(Mach > 0.3)
               {
                 if(dim==2) v_tau = sqrt(v_nrm*v_nrm + 1.5*c*c + c*sqrt(16.0*v_nrm*v_nrm + c*c));
                 else v_tau = sqrt(v_nrm*v_nrm + 2.0*c*c + c*sqrt(4.0*v_nrm*v_nrm + c*c));
               }
  
               const double use = 2.0*std::max(2.0*mu/rho, kappa/(rho*cv))/(h*h) + v_tau/h;                                                                 
               dt_gp_vec.push_back(CFL/use);                          
             }
           }
           auto it = min_element(std::begin(dt_gp_vec), std::end(dt_gp_vec));  
           dt_vec.push_back(*it);                  // here we compute the minimum dt on all gauss points of the element and put that in  dt_vec
         }                 

         set_dt(setup, dt_vec); 
      }
      
      //-----------------------------------------------------------------------------------------------      








     //------------------------------- criterion for MC_isothermal ----------------------------------------------------------------      
      
      void f_mc_isothermal
      (
         model<dim>& f_model, fluid_properties& props, read_setup& setup, int nb_dofs
      )
      {
         const double CFL = setup.CFL();
         const int nb_el_nodes = f_model.mesh().nb_el_nodes();
         
         const auto& mesh(f_model.mesh());
         vec dofs_fluid;
         point<dim> pt;
         vec I(nb_dofs), dI_dx(nb_dofs), dI_dy(nb_dofs), dI_dz(nb_dofs);
         std::vector<double> dt_vec;
         const double tol = std::numeric_limits<double>::epsilon();
         double shear_rate(0.0);           
         
         for(const auto& el : mesh.elements())         
         {
           f_model.extract_element_dofs(*el, dofs_fluid);

           std::vector<double> dt_gp_vec;
           dt_gp_vec.reserve(el->nb_gp());
           for(auto& quad_ptr : el->quad_i())       // here we loop over each gauss point of the element and compute dt        
           {             
             quad_ptr->interpolate(dofs_fluid, I);
             quad_ptr->x_derivative_interpolate(dofs_fluid, dI_dx);
             quad_ptr->y_derivative_interpolate(dofs_fluid, dI_dy);
             if(dim==3) quad_ptr->z_derivative_interpolate(dofs_fluid, dI_dz);
             auto del_v = compute_del_v(dI_dx, dI_dy, dI_dz);


             const double p = I[0];           
             vec v(dim, 0.0);
             for(int i=0; i<dim; i++)
               v[i] = I[1+i];
             vec Y(props.nb_comp(), 0.0);
             for (int i=0; i<props.nb_comp()-1; i++)
               Y[i] = I[dim+1+i];
             Y[props.nb_comp()-1] = 1.0- std::accumulate(&Y[0],&Y[props.nb_comp()-1],0.0);           
             
             
              
             props.properties(p,Y,shear_rate);
             const double rho = props.rho_;
             const double mu = props.mu_;
             const double c = props.sound_speed_;
             const double v_nrm = boost::numeric::ublas::norm_2(v);           
  
  
  
             if(v_nrm < tol)
             {
               dt_gp_vec.push_back(1.2*time::get().delta_t());	
             }
             else
             {
               double h = quad_ptr->compute_h(v, del_v, pt);
               double v_tau = v_nrm;
               const double Mach = v_nrm/c;
               if(Mach > 0.3)
               {
                 if(dim==2) v_tau = sqrt(v_nrm*v_nrm + 1.5*c*c + c*sqrt(16.0*v_nrm*v_nrm + c*c));
                 else v_tau = sqrt(v_nrm*v_nrm + 2.0*c*c + c*sqrt(4.0*v_nrm*v_nrm + c*c));
               }
  
               const double use = 4.0*mu/(rho*h*h) + v_tau/h;                                                                 
               dt_gp_vec.push_back(CFL/use);                          
             }
           }
           auto it = min_element(std::begin(dt_gp_vec), std::end(dt_gp_vec));  
           dt_vec.push_back(*it);                  // here we compute the minimum dt on all gauss points of the element and put that in  dt_vec
         }                 

         set_dt(setup, dt_vec); 
      }
      
      //-----------------------------------------------------------------------------------------------      
      





      mat compute_del_v(const vec& dI_dx, const vec& dI_dy, const vec& dI_dz)
      {
        mat del_v(dim,dim,0.0);             
        if(dim==2)
        {
          for(int i=0; i<2; i++)
          {
             del_v(i,0) = dI_dx[1+i];       del_v(i,1) = dI_dy[1+i];     
          }
        }  
        else
        {
          for(int i=0; i<3; i++)
          {
             del_v(i,0) = dI_dx[1+i];       del_v(i,1) = dI_dy[1+i];     del_v(i,2) = dI_dz[1+i];
          }
        }
        return del_v;
      }



      
      
      void set_dt(read_setup& setup, const std::vector<double>& dt_vec)
      {
         auto it = min_element(std::begin(dt_vec), std::end(dt_vec));  
         double dt_pid = *it;                  
         double dt_min(0.0); 
         MPI_Allreduce(&dt_pid, &dt_min, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

         auto dt = time::get().delta_t();         
         if(dt_min > (1.01*dt)) dt = 1.01*dt;            //Bound maximum increase in dt by 1% of previous dt value
         else dt = dt_min;
                                  
         dt = std::max(std::min(dt, setup.dt_max()), setup.dt_min());     
         time::get().delta_t(dt);      
      }

 };
    
    

}    
    

#endif


