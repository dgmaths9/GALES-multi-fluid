#ifndef _FLUID_2D_HPP_
#define _FLUID_2D_HPP_



#include "../../../src/fem/fem.hpp"




namespace GALES
{


    

  template<typename ic_bc_type, int dim>
  class fluid{};



  template<typename ic_bc_type>
  class fluid<ic_bc_type, 2>
  {
    using element_type = element<2>;
    using vec = boost::numeric::ublas::vector<double>;
    using mat = boost::numeric::ublas::matrix<double>;

    public :

    fluid(ic_bc_type& ic_bc, fluid_properties& props, model<2>& model)
    :
    nb_el_nodes_(model.mesh().nb_el_nodes()),
    nb_dof_fluid_(4 + props.nb_comp()-1),
    nb_comp_(props.nb_comp()),

    ic_bc_(ic_bc),
    props_(props),
    setup_(model.setup()),

    JNxW_(2,0.0),
    JxW_(0.0),
    dt_(0.0),

    P_orig_fluid_(nb_dof_fluid_*nb_el_nodes_),
    I_orig_fluid_(nb_dof_fluid_*nb_el_nodes_),
    P_(nb_dof_fluid_,0.0),
    I_(nb_dof_fluid_,0.0),
    dI_dx_(nb_dof_fluid_,0.0),
    dI_dy_(nb_dof_fluid_,0.0),

    alpha_(0.0),
    beta_(0.0),
    rho_(0.0),
    mu_(0.0),
    cp_(0.0),
    cv_(0.0),
    kappa_(0.0),
    sound_speed_(0.0),

    rho_comp_(nb_comp_,0.0),
    internal_energy_comp_(nb_comp_,0.0),
    chemical_diffusivity_comp_(nb_comp_,0.0),

    U_P_(nb_dof_fluid_,0.0),
    U_(nb_dof_fluid_,0.0),
    F1_adv_(nb_dof_fluid_,0.0),
    F2_adv_(nb_dof_fluid_,0.0),
    F1_dif_(nb_dof_fluid_,0.0),
    F2_dif_(nb_dof_fluid_,0.0),
    S_(nb_dof_fluid_,0.0),

    DUDY_(nb_dof_fluid_,nb_dof_fluid_,0.0),
    A1_(nb_dof_fluid_,nb_dof_fluid_,0.0),
    A2_(nb_dof_fluid_,nb_dof_fluid_,0.0),
    K11_(nb_dof_fluid_,nb_dof_fluid_,0.0),
    K12_(nb_dof_fluid_,nb_dof_fluid_,0.0),
    K21_(nb_dof_fluid_,nb_dof_fluid_,0.0),
    K22_(nb_dof_fluid_,nb_dof_fluid_,0.0),
    S0_(nb_dof_fluid_,nb_dof_fluid_,0.0),

    tau_(nb_dof_fluid_,nb_dof_fluid_,0.0),
    shock_matrix_(nb_dof_fluid_,nb_dof_fluid_,0.0),
    u_ref_inv_(nb_dof_fluid_, nb_dof_fluid_, 0.0),

    p_indx_(0), vx_indx_(1), vy_indx_(2), T_indx_(3), Y_indx_(nb_comp_),
    rho_indx_(0), rho_vx_indx_(1), rho_vy_indx_(2), rho_E_indx_(3), rho_Y_indx_(nb_comp_),

    r_(nb_dof_fluid_*nb_el_nodes_,0.0),
    m_(nb_dof_fluid_*nb_el_nodes_, nb_dof_fluid_*nb_el_nodes_,0.0),
    gravity_(2,0.0),

    tol_(std::numeric_limits<double>::epsilon())
    {
      for(int i=0; i<nb_comp_-1; i++)
      {
        Y_indx_[i] = 4 + i;
        rho_Y_indx_[i] = 4 + i;
      }
      
      ic_bc_.body_force(gravity_, dummy_);
      if(setup_.dc_2006())  u_ref_inv();   
    }





     void assign(const vec &v, double &p, double &vx, double &vy, double &T, vec &Y)
     {
        p =  v[p_indx_];
        vx = v[vx_indx_];
        vy = v[vy_indx_];
        T =  v[T_indx_];

        for(int i=0; i<nb_comp_-1;i++)
          Y[i] = v[Y_indx_[i]];
        Y[nb_comp_-1] = 1.0 - std::accumulate(&Y[0],&Y[nb_comp_-1],0.0);
     }



     void u_ref_inv()
     {
         vec u_ref(nb_dof_fluid_,0.0);
         vec Y_ref(nb_comp_, 0.5);
         props_.properties(ic_bc_.p_ref_, ic_bc_.T_ref_, Y_ref, 0.0);         
         U(ic_bc_.p_ref_, ic_bc_.v1_ref_, ic_bc_.v2_ref_, ic_bc_.T_ref_, Y_ref, u_ref);
 
         for(int i=0; i<nb_dof_fluid_; i++)
         {
             if(u_ref[i] < tol_)  u_ref_inv_(i,i) = 0.0;
             else u_ref_inv_(i,i) = 1.0/u_ref[i];
         }
     }
 
 
 
     void get_dofs(const std::vector<vec>& dofs_fluid)
     {
         dt_= time::get().delta_t();
 
         P_orig_fluid_ = dofs_fluid[1];
         I_orig_fluid_ = dofs_fluid[0];
     }
 
 
 
 
     void get_properties()
     {
         rho_ = props_.rho_;
         mu_ = props_.mu_;
         kappa_ = props_.kappa_;
         sound_speed_ = props_.sound_speed_;
         cv_ = props_.cv_;
         cp_ = props_.cp_;
         alpha_ = props_.alpha_;
         beta_ = props_.beta_;

        rho_comp_ = props_.rho_comp_;
        chemical_diffusivity_comp_ = props_.chemical_diffusivity_comp_;
        internal_energy_comp_ = props_.internal_energy_comp_;
        
        if(setup_.print_props())   props_.print_props_mc();
     }
 

 
 
     void execute_i(const element_type& el, const std::vector<vec>& dofs_fluid, mat& m, vec& r)
     {
         get_dofs(dofs_fluid);
         r_.clear();
         m_.clear();

          for(const auto& quad_ptr : el.quad_i())
          {
              quad_ptr_ = quad_ptr;
              JxW_ = quad_ptr_->JxW();

           double p,vx,vy,T;
           vec Y(nb_comp_);
           quad_ptr_->interpolate(P_orig_fluid_,P_);
           assign(P_,p,vx,vy,T,Y);
           props_.properties(p, T, Y, 0.0);           
           U(p,vx,vy,T,Y, U_P_);

           quad_ptr_->interpolate(I_orig_fluid_,I_);
           quad_ptr_->x_derivative_interpolate(I_orig_fluid_,dI_dx_);
           quad_ptr_->y_derivative_interpolate(I_orig_fluid_,dI_dy_);
           assign(I_,p,vx,vy,T,Y);
           const double shear_rate = compute_shear_rate(dI_dx_, dI_dy_);
           props_.properties(p, T, Y, shear_rate);           
           get_properties();
           flux_matrices_and_vectors(p,vx,vy,T,Y);

           RM_i(el, p,vx,vy,T,Y);
         }
         m = m_;
         r = -r_;
     }






     void execute_b(const element_type& el, const std::vector<vec>& dofs_fluid, mat& m, vec& r)
     {
         get_dofs(dofs_fluid);
         r_.clear();
         m_.clear();

         for(int i=0; i<el.nb_sides(); i++)
         {
           if(el.is_side_on_boundary(i))
           {
             auto bd_nodes = el.side_nodes(i);
             auto side_flag = el.side_flag(i);
             
             for(const auto& quad_ptr: el.quad_b(i))
             {
               quad_ptr_ = quad_ptr;               
               JNxW_ = quad_ptr_->JNxW();

                quad_ptr_->interpolate(I_orig_fluid_,I_);
                quad_ptr_->x_derivative_interpolate(I_orig_fluid_,dI_dx_);
                quad_ptr_->y_derivative_interpolate(I_orig_fluid_,dI_dy_);

                double p,vx,vy,T;
                vec Y(nb_comp_);
                assign(I_,p,vx,vy,T,Y);
                const double shear_rate = compute_shear_rate(dI_dx_, dI_dy_);
                props_.properties(p, T, Y, shear_rate);           
                get_properties();
                flux_matrices_and_vectors(p,vx,vy,T,Y,bd_nodes, side_flag);

                RM_b();
             }
           }
         }
         m = m_;
         r = -r_;
     }





     void RM_i(const element_type& el, double p, double vx, double vy, double T, const vec &Y)
     {
         for (int b=0; b<nb_el_nodes_; b++)
         for (int jb=0; jb<nb_dof_fluid_; jb++)
         {
             r_[nb_dof_fluid_*b+jb] += quad_ptr_->sh(b)*(U_[jb]*JxW_ - U_P_[jb]*JxW_)/dt_;
             r_[nb_dof_fluid_*b+jb] -= quad_ptr_->dsh_dx(b)*F1_adv_[jb]*JxW_;
             r_[nb_dof_fluid_*b+jb] -= quad_ptr_->dsh_dy(b)*F2_adv_[jb]*JxW_;
             r_[nb_dof_fluid_*b+jb] += quad_ptr_->dsh_dx(b)*F1_dif_[jb]*JxW_;
             r_[nb_dof_fluid_*b+jb] += quad_ptr_->dsh_dy(b)*F2_dif_[jb]*JxW_;
             r_[nb_dof_fluid_*b+jb] -= quad_ptr_->sh(b)*S_[jb]*JxW_;
         }

         for(int b=0; b<nb_el_nodes_; b++)
         for(int jb=0; jb<nb_dof_fluid_; jb++)
         for(int a=0; a<nb_el_nodes_; a++)
         for(int ja=0; ja<nb_dof_fluid_; ja++)
         {
             m_(nb_dof_fluid_*b+jb,nb_dof_fluid_*a+ja) +=  quad_ptr_->sh(b) * DUDY_(jb,ja) * quad_ptr_->sh(a) * JxW_/dt_;
             m_(nb_dof_fluid_*b+jb,nb_dof_fluid_*a+ja) -=  quad_ptr_->dsh_dx(b) * A1_(jb,ja) *quad_ptr_->sh(a) * JxW_;
             m_(nb_dof_fluid_*b+jb,nb_dof_fluid_*a+ja) -=  quad_ptr_->dsh_dy(b) * A2_(jb,ja) *quad_ptr_->sh(a) * JxW_;
             m_(nb_dof_fluid_*b+jb,nb_dof_fluid_*a+ja) +=  quad_ptr_->dsh_dx(b) * K11_(jb,ja) * quad_ptr_->dsh_dx(a) * JxW_;
             m_(nb_dof_fluid_*b+jb,nb_dof_fluid_*a+ja) +=  quad_ptr_->dsh_dx(b) * K12_(jb,ja) * quad_ptr_->dsh_dy(a) * JxW_;
             m_(nb_dof_fluid_*b+jb,nb_dof_fluid_*a+ja) +=  quad_ptr_->dsh_dy(b) * K21_(jb,ja) * quad_ptr_->dsh_dx(a) * JxW_;
             m_(nb_dof_fluid_*b+jb,nb_dof_fluid_*a+ja) +=  quad_ptr_->dsh_dy(b) * K22_(jb,ja) * quad_ptr_->dsh_dy(a) * JxW_;
             m_(nb_dof_fluid_*b+jb,nb_dof_fluid_*a+ja) -=  quad_ptr_->sh(b) * S0_(jb,ja) * quad_ptr_->sh(a) * JxW_;
         }

         //--------------------------------------------R_LeastSquares--------------------------------------------------------------------------
         tau_matrix(p,vx,vy,T,Y);

         const vec Res_vec = prod(DUDY_,(I_-P_)/dt_) + prod(A1_,dI_dx_) + prod(A2_,dI_dy_) - S_;
         for(int b=0; b<nb_el_nodes_; b++)
         {
             const mat first_part = A1_*quad_ptr_->dsh_dx(b) + A2_*quad_ptr_->dsh_dy(b) - S0_*quad_ptr_->sh(b);
             const mat first_part_tau = prod(first_part, tau_);
             const vec use = prod(first_part_tau, Res_vec);
             for(int jb=0; jb<nb_dof_fluid_; jb++)
             r_[nb_dof_fluid_*b+jb] += use[jb]*JxW_;
         }

         //--------------------------------------------M_LeastSquares--------------------------------------------------------------------------
         for(int a=0; a<nb_el_nodes_; a++)
         {
             const mat Res_mat = DUDY_*quad_ptr_->sh(a)/dt_ + A1_*quad_ptr_->dsh_dx(a) + A2_*quad_ptr_->dsh_dy(a) - S0_*quad_ptr_->sh(a);
             for(int b=0; b<nb_el_nodes_; b++)
             {
                 const mat first_part = A1_*quad_ptr_->dsh_dx(b) + A2_*quad_ptr_->dsh_dy(b) - S0_*quad_ptr_->sh(b);
                 const mat first_part_tau = prod(first_part, tau_);
                 const mat use = prod(first_part_tau, Res_mat);
                 for(int jb=0; jb<nb_dof_fluid_; jb++)
                 for(int ja=0; ja<nb_dof_fluid_; ja++)
                 m_(nb_dof_fluid_*b+jb, nb_dof_fluid_*a+ja) += use(jb,ja)*JxW_;
             }
         }


         //-----------------dc------------------------
         dc(p,vx,vy,T,Y);

         if(setup_.dc_2006())
         {
             const vec app_1 = prod(shock_matrix_,dI_dx_);
             const vec app_2 = prod(shock_matrix_,dI_dy_);

             for(int b=0;b<nb_el_nodes_;b++)
             for(int jb=0;jb<nb_dof_fluid_;jb++)
             {
                 r_[nb_dof_fluid_*b+jb] += quad_ptr_->dsh_dx(b)*app_1[jb]*JxW_;
                 r_[nb_dof_fluid_*b+jb] += quad_ptr_->dsh_dy(b)*app_2[jb]*JxW_;
             }

             for(int b=0;b<nb_el_nodes_;b++)
             for(int jb=0;jb<nb_dof_fluid_;jb++)
             for(int a=0;a<nb_el_nodes_;a++)
             for(int ja=0; ja<nb_dof_fluid_;ja++)
             {
                 m_(nb_dof_fluid_*b+jb,nb_dof_fluid_*a+ja) += quad_ptr_->dsh_dx(b)*shock_matrix_(jb,ja)*quad_ptr_->dsh_dx(a) * JxW_;
                 m_(nb_dof_fluid_*b+jb,nb_dof_fluid_*a+ja) += quad_ptr_->dsh_dy(b)*shock_matrix_(jb,ja)*quad_ptr_->dsh_dy(a) * JxW_;
             }
         }

     }







     void RM_b()
     {
         for(int b=0;b<nb_el_nodes_;b++)
         for(int j=0;j<nb_dof_fluid_;j++)
         {
             r_[nb_dof_fluid_*b+j] += quad_ptr_->sh(b)*F1_adv_[j]*JNxW_[0];
             r_[nb_dof_fluid_*b+j] += quad_ptr_->sh(b)*F2_adv_[j]*JNxW_[1];
             r_[nb_dof_fluid_*b+j] -= quad_ptr_->sh(b)*F1_dif_[j]*JNxW_[0];
             r_[nb_dof_fluid_*b+j] -= quad_ptr_->sh(b)*F2_dif_[j]*JNxW_[1];
         }

         for(int b=0; b<nb_el_nodes_; b++)
         for(int jb=0;jb<nb_dof_fluid_;jb++)
         for(int a=0; a<nb_el_nodes_; a++)
         for(int ja=0;ja<nb_dof_fluid_;ja++)
         {
             m_(nb_dof_fluid_*b+jb,nb_dof_fluid_*a+ja) += quad_ptr_->sh(b) * A1_(jb,ja) * JNxW_[0] * quad_ptr_->sh(a);
             m_(nb_dof_fluid_*b+jb,nb_dof_fluid_*a+ja) += quad_ptr_->sh(b) * A2_(jb,ja) * JNxW_[1] * quad_ptr_->sh(a);
             m_(nb_dof_fluid_*b+jb,nb_dof_fluid_*a+ja) -= quad_ptr_->sh(b) * K11_(jb,ja) * JNxW_[0] * quad_ptr_->dsh_dx(a);
             m_(nb_dof_fluid_*b+jb,nb_dof_fluid_*a+ja) -= quad_ptr_->sh(b) * K12_(jb,ja) * JNxW_[0] * quad_ptr_->dsh_dy(a);
             m_(nb_dof_fluid_*b+jb,nb_dof_fluid_*a+ja) -= quad_ptr_->sh(b) * K21_(jb,ja) * JNxW_[1] * quad_ptr_->dsh_dx(a);
             m_(nb_dof_fluid_*b+jb,nb_dof_fluid_*a+ja) -= quad_ptr_->sh(b) * K22_(jb,ja) * JNxW_[1] * quad_ptr_->dsh_dy(a);
         }
     }






     //this is called to generate U_P_, U_ and U_ref_
     void U(double p,  double vx,  double vy, double T, const vec &Y, vec& u)
     {
         double rho =  props_.rho_;
         double cv = props_.cv_;
         double E = cv*T + 0.5*(vx*vx+vy*vy);

         u[rho_indx_] = rho;
         u[rho_vx_indx_] = rho*vx;
         u[rho_vy_indx_] = rho*vy;
         u[rho_E_indx_] = rho*E;

         for(int i=0; i<nb_comp_-1; i++)
           u[rho_Y_indx_[i]] = rho*Y[i];
     }





     void flux_matrices_and_vectors(double p, double vx, double vy, double T, const vec &Y)
     {
         const double E = cv_*T + 0.5*(vx*vx + vy*vy);
         const double drho_dp = rho_*beta_;
         const double drho_dT = -rho_*alpha_;
         const double de_dp =  (beta_*p - alpha_*T)/rho_;
         const double de_dT =  cp_ - p*alpha_/rho_;
         
         vec drho_dY(nb_comp_-1,0.0), dE_dY(nb_comp_-1,0.0);
         for(int i=0; i<nb_comp_-1; i++)
         {
           drho_dY[i] = -rho_*rho_* (1./rho_comp_[i]-1/rho_comp_[nb_comp_-1]);
           dE_dY[i] = internal_energy_comp_[i]-internal_energy_comp_[nb_comp_-1];
         }  

         //------------------------------------------------------------------------
         DUDY_.clear();
         DUDY_(rho_indx_,p_indx_) = drho_dp;
         DUDY_(rho_indx_,T_indx_) = drho_dT;
         for(int c=0; c<nb_comp_-1; c++)
           DUDY_(rho_indx_, Y_indx_[c]) = drho_dY[c];            

         DUDY_(rho_vx_indx_,p_indx_) = drho_dp*vx;
         DUDY_(rho_vx_indx_,vx_indx_) = rho_;
         DUDY_(rho_vx_indx_,T_indx_) = drho_dT*vx;
         for(int c=0; c<nb_comp_-1; c++)
           DUDY_(rho_vx_indx_, Y_indx_[c]) = vx*drho_dY[c];

         DUDY_(rho_vy_indx_,p_indx_) = drho_dp*vy;
         DUDY_(rho_vy_indx_,vy_indx_) = rho_;
         DUDY_(rho_vy_indx_,T_indx_) = drho_dT*vy;
         for(int c=0; c<nb_comp_-1; c++)
           DUDY_(rho_vy_indx_, Y_indx_[c]) = vy*drho_dY[c];

         DUDY_(rho_E_indx_,p_indx_) = drho_dp*E + rho_*de_dp;
         DUDY_(rho_E_indx_,vx_indx_) = rho_*vx;
         DUDY_(rho_E_indx_,vy_indx_) = rho_*vy;
         DUDY_(rho_E_indx_,T_indx_) = drho_dT*E + rho_*de_dT;
         for(int c=0; c<nb_comp_-1; c++)
           DUDY_(rho_E_indx_, Y_indx_[c]) = E*drho_dY[c] + rho_*dE_dY[c];

         for(int i=0; i<nb_comp_-1; i++)
         {
           DUDY_(rho_Y_indx_[i], p_indx_) = drho_dp*Y[i];
           DUDY_(rho_Y_indx_[i], T_indx_) = drho_dT*Y[i];
           for(int c=0; c<nb_comp_-1; c++)
             DUDY_(rho_Y_indx_[i], Y_indx_[c]) = Y[i]*drho_dY[c];
         }

         for(int i=0; i<nb_comp_-1; i++)
           DUDY_(rho_Y_indx_[i], Y_indx_[i]) += rho_;


         //------------------------------------------------------------------------
         A1_.clear();
         A1_(rho_indx_,p_indx_)  = drho_dp*vx;
         A1_(rho_indx_,vx_indx_) = rho_;
         A1_(rho_indx_,T_indx_)  = drho_dT*vx;
         for(int c=0; c<nb_comp_-1; c++)
           A1_(rho_indx_, Y_indx_[c]) = drho_dY[c]*vx;

         A1_(rho_vx_indx_,p_indx_)  = drho_dp*vx*vx+1.0;
         A1_(rho_vx_indx_,vx_indx_) = 2.0*rho_*vx;
         A1_(rho_vx_indx_,T_indx_)  = drho_dT*vx*vx;
         for(int c=0; c<nb_comp_-1; c++)
           A1_(rho_vx_indx_, Y_indx_[c]) = drho_dY[c]*vx*vx;

         A1_(rho_vy_indx_,p_indx_)  = drho_dp*vy*vx;
         A1_(rho_vy_indx_,vx_indx_) = rho_*vy;
         A1_(rho_vy_indx_,vy_indx_) = rho_*vx;
         A1_(rho_vy_indx_,T_indx_)  = drho_dT*vy*vx;
         for(int c=0; c<nb_comp_-1; c++)
           A1_(rho_vy_indx_, Y_indx_[c]) = drho_dY[c]*vy*vx;

         A1_(rho_E_indx_,p_indx_)  = vx + vx*(drho_dp*E + rho_*de_dp);
         A1_(rho_E_indx_,vx_indx_) = p + rho_*E + rho_*vx*vx;
         A1_(rho_E_indx_,vy_indx_) = rho_*vx*vy;
         A1_(rho_E_indx_,T_indx_)  = vx*(drho_dT*E + rho_*de_dT);
         for(int c=0; c<nb_comp_-1; c++)
           A1_(rho_E_indx_, Y_indx_[c]) = (E*drho_dY[c] + rho_*dE_dY[c])*vx;

         for(int i=0; i<nb_comp_-1; i++)
         {
           A1_(rho_Y_indx_[i], p_indx_)  = drho_dp*Y[i]*vx;
           A1_(rho_Y_indx_[i], vx_indx_) = rho_*Y[i];
           A1_(rho_Y_indx_[i], T_indx_)  = drho_dT*Y[i]*vx;
           for(int c=0; c<nb_comp_-1; c++)
             A1_(rho_Y_indx_[i], Y_indx_[c]) = Y[i]*drho_dY[c]*vx;
         }

         for(int i=0; i<nb_comp_-1; i++)
           A1_(rho_Y_indx_[i], Y_indx_[i]) += rho_*vx;



         //------------------------------------------------------------------------
         A2_.clear();
         A2_(rho_indx_,p_indx_)  = drho_dp*vy;
         A2_(rho_indx_,vy_indx_) = rho_;
         A2_(rho_indx_,T_indx_)  = drho_dT*vy;
         for(int c=0; c<nb_comp_-1; c++)
           A2_(rho_indx_, Y_indx_[c]) = drho_dY[c]*vy;

         A2_(rho_vx_indx_,p_indx_)  = drho_dp*vx*vy;
         A2_(rho_vx_indx_,vx_indx_) = rho_*vy;
         A2_(rho_vx_indx_,vy_indx_) = rho_*vx;
         A2_(rho_vx_indx_,T_indx_)  = drho_dT*vx*vy;
         for(int c=0; c<nb_comp_-1; c++)
           A2_(rho_vx_indx_, Y_indx_[c]) = drho_dY[c]*vx*vy;

         A2_(rho_vy_indx_,p_indx_)  = drho_dp*vy*vy+1.0;
         A2_(rho_vy_indx_,vy_indx_) = 2.0*rho_*vy;
         A2_(rho_vy_indx_,T_indx_)  = drho_dT*vy*vy;
         for(int c=0; c<nb_comp_-1; c++)
           A2_(rho_vy_indx_, Y_indx_[c]) = drho_dY[c]*vy*vy;

         A2_(rho_E_indx_,p_indx_)  = vy + vy*(drho_dp*E + rho_*de_dp);
         A2_(rho_E_indx_,vx_indx_) = rho_*vx*vy;
         A2_(rho_E_indx_,vy_indx_) = p + rho_*E + rho_*vy*vy;
         A2_(rho_E_indx_,T_indx_)  = vy*(drho_dT*E + rho_*de_dT);
         for(int c=0; c<nb_comp_-1; c++)
           A2_(rho_E_indx_, Y_indx_[c]) = (E*drho_dY[c] + rho_*dE_dY[c])*vy;

         for(int i=0; i<nb_comp_-1; i++)
         {
           A2_(rho_Y_indx_[i], p_indx_)  = drho_dp*Y[i]*vy;
           A2_(rho_Y_indx_[i], vy_indx_) = rho_*Y[i];
           A2_(rho_Y_indx_[i], T_indx_)  = drho_dT*Y[i]*vy;
           for(int c=0; c<nb_comp_-1; c++)
             A2_(rho_Y_indx_[i], Y_indx_[c]) = Y[i]*drho_dY[c]*vy;
         }

         for(int i=0; i<nb_comp_-1; i++)
           A2_(rho_Y_indx_[i], Y_indx_[i]) += rho_*vy;
       
       
       


         //------------------------------------------------------------------------
         K11_.clear(); K12_.clear(); 
         K21_.clear(); K22_.clear(); 

         vec H(nb_comp_,0.0);
         for(int i=0; i<nb_comp_; i++)
           H[i] = internal_energy_comp_[i] + p/rho_comp_[i];

         const double lambda = -2.0*mu_/3.0;
         auto D = chemical_diffusivity_comp_;

         K11_(rho_vx_indx_,vx_indx_) = lambda + 2*mu_;
         K11_(rho_vy_indx_,vy_indx_) = mu_;
         K11_(rho_E_indx_,vx_indx_) = (lambda + 2*mu_)*vx;
         K11_(rho_E_indx_,vy_indx_) = mu_*vy;
         K11_(rho_E_indx_,T_indx_) =  kappa_;
         for(int i=0; i<nb_comp_-1; i++)
         {
           K11_(rho_Y_indx_[i], rho_Y_indx_[i]) = rho_*D[i];
           K11_(rho_E_indx_, rho_Y_indx_[i]) = rho_*D[i]*H[i];
         }  

         K12_(rho_vx_indx_,vy_indx_) = lambda;
         K12_(rho_vy_indx_,vx_indx_) =  mu_;
         K12_(rho_E_indx_,vx_indx_) = mu_*vy;
         K12_(rho_E_indx_,vy_indx_) = lambda*vx;

         K21_(rho_vx_indx_,vy_indx_) = mu_;
         K21_(rho_vy_indx_,vx_indx_) = lambda;
         K21_(rho_E_indx_,vx_indx_) = lambda*vy;
         K21_(rho_E_indx_,vy_indx_) = mu_*vx;

         K22_(rho_vx_indx_,vx_indx_) = mu_;
         K22_(rho_vy_indx_,vy_indx_) = lambda + 2*mu_;
         K22_(rho_E_indx_,vx_indx_) = mu_*vx;
         K22_(rho_E_indx_,vy_indx_) = (lambda + 2*mu_)*vy;
         K22_(rho_E_indx_,T_indx_) =  kappa_;
         for(int i=0; i<nb_comp_-1; i++)
         {
           K22_(rho_Y_indx_[i], rho_Y_indx_[i]) = rho_*D[i];
           K22_(rho_E_indx_, rho_Y_indx_[i]) = rho_*D[i]*H[i];
         }  

         S0_.clear();
         S0_(rho_E_indx_,vx_indx_) = rho_*gravity_[0];
         S0_(rho_E_indx_,vy_indx_) = rho_*gravity_[1];

         U(p,vx,vy,T,Y,U_);
         
                  

         F1_adv_[rho_indx_] = rho_*vx;
         F1_adv_[rho_vx_indx_] = rho_*vx*vx + p;
         F1_adv_[rho_vy_indx_] = rho_*vx*vy;
         F1_adv_[rho_E_indx_] = vx*(rho_*E+p);

         F2_adv_[rho_indx_] = rho_*vy;
         F2_adv_[rho_vx_indx_] = rho_*vy*vx;
         F2_adv_[rho_vy_indx_] = rho_*vy*vy + p;
         F2_adv_[rho_E_indx_] = vy*(rho_*E+p);

         for(int i=0; i<nb_comp_-1; i++)
         {
           F1_adv_[rho_Y_indx_[i]] = rho_*vx*Y[i];
           F2_adv_[rho_Y_indx_[i]] = rho_*vy*Y[i];
         }

         S_[rho_vx_indx_] = rho_*gravity_[0];
         S_[rho_vy_indx_] = rho_*gravity_[1];
         S_[rho_E_indx_] = rho_*(gravity_[0]*vx + gravity_[1]*vy);

         F1_dif_ =  prod(K11_,dI_dx_) + prod(K12_,dI_dy_);
         F2_dif_ =  prod(K21_,dI_dx_) + prod(K22_,dI_dy_);
     }






     void flux_matrices_and_vectors(double p, double vx, double vy, double T, const vec &Y, const std::vector<int>& bd_nodes, int side_flag)
     {
         flux_matrices_and_vectors(p,vx,vy,T,Y);

         const double E = cv_*T + 0.5*(vx*vx + vy*vy);
         const double drho_dp = rho_*beta_;
         const double drho_dT = -rho_*alpha_;
         const double de_dp =  (beta_*p - alpha_*T)/rho_;
         const double de_dT =  cp_ - p*alpha_/rho_;

         vec drho_dY(nb_comp_-1,0.0), dE_dY(nb_comp_-1,0.0);
         for(int i=0; i<nb_comp_-1; i++)
         {
           drho_dY[i] = -rho_*rho_* (1./rho_comp_[i]-1/rho_comp_[nb_comp_-1]);
           dE_dY[i] = internal_energy_comp_[i]-internal_energy_comp_[nb_comp_-1];
         }  

         if(ic_bc_.mass_flux1(bd_nodes, side_flag).first)
         {
            const double mass_flux1 = ic_bc_.mass_flux1(bd_nodes, side_flag).second;
         
            U_.clear();
            U_[rho_indx_] = rho_;
            U_[rho_vx_indx_] = mass_flux1;
            U_[rho_vy_indx_] = rho_*vy;
            U_[rho_E_indx_] = rho_*cv_*T + mass_flux1*vx/2 + rho_*vy*vy/2;
            for(int i=0; i<nb_comp_-1; i++)
              U_[rho_Y_indx_[i]] = rho_*Y[i];
        
            DUDY_.clear();
            DUDY_(rho_indx_,p_indx_) = drho_dp;
            DUDY_(rho_indx_,T_indx_) = drho_dT;
            for(int c=0; c<nb_comp_-1; c++)
              DUDY_(rho_indx_, Y_indx_[c]) = drho_dY[c];            

            DUDY_(rho_vy_indx_,p_indx_) = drho_dp*vy;
            DUDY_(rho_vy_indx_,vy_indx_) = rho_;
            DUDY_(rho_vy_indx_,T_indx_) = drho_dT*vy;
            for(int c=0; c<nb_comp_-1; c++)
              DUDY_(rho_vy_indx_, Y_indx_[c]) = vy*drho_dY[c];

            DUDY_(rho_E_indx_,p_indx_) = drho_dp*(cv_*T + 0.5*vy*vy) + rho_*de_dp;
            DUDY_(rho_E_indx_,vx_indx_) = mass_flux1/2;
            DUDY_(rho_E_indx_,vy_indx_) = rho_*vy;
            DUDY_(rho_E_indx_,T_indx_) = drho_dT*(cv_*T + 0.5*vy*vy) + rho_*de_dT;
            for(int c=0; c<nb_comp_-1; c++)
              DUDY_(rho_E_indx_, Y_indx_[c]) = drho_dY[c]*(cv_*T + 0.5*vy*vy)  + rho_*dE_dY[c];

            for(int i=0; i<nb_comp_-1; i++)
            {
              DUDY_(rho_Y_indx_[i], p_indx_) = drho_dp*Y[i];
              DUDY_(rho_Y_indx_[i], T_indx_) = drho_dT*Y[i];
              for(int c=0; c<nb_comp_-1; c++)
                DUDY_(rho_Y_indx_[i], Y_indx_[c]) = Y[i]*drho_dY[c];
            }

            for(int i=0; i<nb_comp_-1; i++)
              DUDY_(rho_Y_indx_[i], Y_indx_[i]) += rho_;


            F1_adv_.clear();
            F1_adv_[rho_indx_] = mass_flux1;
            F1_adv_[rho_vx_indx_] = mass_flux1*vx + p;
            F1_adv_[rho_vy_indx_] = mass_flux1*vy;
            F1_adv_[rho_E_indx_] = mass_flux1*E + p*vx;         
            for(int i=0; i<nb_comp_-1; i++)
              F1_adv_[rho_Y_indx_[i]] = mass_flux1*Y[i];

            A1_.clear();
            A1_(rho_vx_indx_,p_indx_)  = 1.0;
            A1_(rho_vx_indx_,vx_indx_) = mass_flux1;
            A1_(rho_vy_indx_,vy_indx_) = mass_flux1;
            A1_(rho_E_indx_,p_indx_)  = vx + mass_flux1*de_dp;
            A1_(rho_E_indx_,vx_indx_) = p + mass_flux1*vx;
            A1_(rho_E_indx_,vy_indx_) = mass_flux1*vy;
            A1_(rho_E_indx_,T_indx_)  = mass_flux1*de_dT;
            for(int c=0; c<nb_comp_-1; c++)
              A1_(rho_E_indx_, Y_indx_[c]) = mass_flux1*dE_dY[c];

            for(int i=0; i<nb_comp_-1; i++)
              A1_(rho_Y_indx_[i], Y_indx_[i]) = mass_flux1;
         }    


         if(ic_bc_.mass_flux2(bd_nodes, side_flag).first)
         {
            const double mass_flux2 = ic_bc_.mass_flux2(bd_nodes, side_flag).second;
         
            U_.clear();
            U_[rho_indx_] = rho_;
            U_[rho_vx_indx_] = rho_*vx;
            U_[rho_vy_indx_] = mass_flux2;
            U_[rho_E_indx_] = rho_*cv_*T + mass_flux2*vy/2 + rho_*vx*vx/2;
            for(int i=0; i<nb_comp_-1; i++)
              U_[rho_Y_indx_[i]] = rho_*Y[i];
        
            DUDY_.clear();
            DUDY_(rho_indx_,p_indx_) = drho_dp;
            DUDY_(rho_indx_,T_indx_) = drho_dT;
            for(int c=0; c<nb_comp_-1; c++)
              DUDY_(rho_indx_, Y_indx_[c]) = drho_dY[c];            

            DUDY_(rho_vx_indx_,p_indx_) = drho_dp*vx;
            DUDY_(rho_vx_indx_,vx_indx_) = rho_;
            DUDY_(rho_vx_indx_,T_indx_) = drho_dT*vx;
            for(int c=0; c<nb_comp_-1; c++)
              DUDY_(rho_vx_indx_, Y_indx_[c]) = vx*drho_dY[c];

            DUDY_(rho_E_indx_,p_indx_) = drho_dp*(cv_*T + 0.5*vx*vx) + rho_*de_dp;
            DUDY_(rho_E_indx_,vx_indx_) = rho_*vx;
            DUDY_(rho_E_indx_,vy_indx_) = mass_flux2/2;
            DUDY_(rho_E_indx_,T_indx_) = drho_dT*(cv_*T + 0.5*vx*vx) + rho_*de_dT;
            for(int c=0; c<nb_comp_-1; c++)
              DUDY_(rho_E_indx_, Y_indx_[c]) = drho_dY[c]*(cv_*T + 0.5*vx*vx)  + rho_*dE_dY[c];

            for(int i=0; i<nb_comp_-1; i++)
            {
              DUDY_(rho_Y_indx_[i], p_indx_) = drho_dp*Y[i];
              DUDY_(rho_Y_indx_[i], T_indx_) = drho_dT*Y[i];
              for(int c=0; c<nb_comp_-1; c++)
                DUDY_(rho_Y_indx_[i], Y_indx_[c]) = Y[i]*drho_dY[c];
            }

            for(int i=0; i<nb_comp_-1; i++)
              DUDY_(rho_Y_indx_[i], Y_indx_[i]) += rho_;


            F2_adv_.clear();
            F2_adv_[rho_indx_] = mass_flux2;
            F2_adv_[rho_vx_indx_] = mass_flux2*vx;
            F2_adv_[rho_vy_indx_] = mass_flux2*vy + p;
            F2_adv_[rho_E_indx_] = mass_flux2*E + p*vy;         
            for(int i=0; i<nb_comp_-1; i++)
              F2_adv_[rho_Y_indx_[i]] = mass_flux2*Y[i];

            A2_.clear();
            A2_(rho_vx_indx_,vx_indx_) = mass_flux2;
            A2_(rho_vy_indx_,p_indx_)  = 1.0;
            A2_(rho_vy_indx_,vy_indx_) = mass_flux2;
            A2_(rho_E_indx_,p_indx_)  = vy + mass_flux2*de_dp;
            A2_(rho_E_indx_,vx_indx_) = mass_flux2*vx;
            A2_(rho_E_indx_,vy_indx_) = p + mass_flux2*vy;
            A2_(rho_E_indx_,T_indx_)  = mass_flux2*de_dT;
            for(int c=0; c<nb_comp_-1; c++)
              A2_(rho_E_indx_, Y_indx_[c]) = mass_flux2*dE_dY[c];

            for(int i=0; i<nb_comp_-1; i++)
              A2_(rho_Y_indx_[i], Y_indx_[i]) = mass_flux2;
         }    


         const double vx_x(dI_dx_[1]), vx_y(dI_dy_[1]);
         const double vy_x(dI_dx_[2]), vy_y(dI_dy_[2]);
         const double T_x(dI_dx_[3]),  T_y(dI_dy_[3]);
         vec Y_x(nb_comp_-1,0.0), Y_y(nb_comp_-1,0.0);    
         for(int i=0; i<nb_comp_-1; i++)
         {
            Y_x[i] = dI_dx_[Y_indx_[i]];     Y_y[i] = dI_dy_[Y_indx_[i]]; 
         }

         const double lambda = -2.0*mu_/3.0;
         double tau11 = lambda*(vx_x + vy_y) + 2.0*mu_*vx_x;
         double tau22 = lambda*(vx_x + vy_y) + 2.0*mu_*vy_y;
         double tau12 = mu_*(vx_y + vy_x);
         double q1 = -kappa_*T_x;
         double q2 = -kappa_*T_y;
         vec J1(nb_comp_-1,0.0), J2(nb_comp_-1,0.0);
         auto d = chemical_diffusivity_comp_;
         for(int i=0; i<nb_comp_-1; i++)
         {
           J1[i] = rho_*d[i]*Y_x[i];    J2[i] = rho_*d[i]*Y_y[i];        
         }  
 

         if(ic_bc_.neumann_tau11(bd_nodes, side_flag).first)    tau11 = ic_bc_.neumann_tau11(bd_nodes, side_flag).second;
         if(ic_bc_.neumann_tau22(bd_nodes, side_flag).first)    tau22 = ic_bc_.neumann_tau22(bd_nodes, side_flag).second;
         if(ic_bc_.neumann_tau12(bd_nodes, side_flag).first)    tau12 = ic_bc_.neumann_tau12(bd_nodes, side_flag).second;
         if(ic_bc_.neumann_q1(bd_nodes, side_flag).first)          q1 = ic_bc_.neumann_q1(bd_nodes, side_flag).second;
         if(ic_bc_.neumann_q2(bd_nodes, side_flag).first)          q2 = ic_bc_.neumann_q2(bd_nodes, side_flag).second;
         
         for(int i=0; i<nb_comp_-1; i++)
         {
           if(ic_bc_.neumann_J1(i,bd_nodes, side_flag).first)          J1[i] = ic_bc_.neumann_J1(i,bd_nodes, side_flag).second;
           if(ic_bc_.neumann_J2(i,bd_nodes, side_flag).first)          J2[i] = ic_bc_.neumann_J2(i,bd_nodes, side_flag).second;
         }

         vec H(nb_comp_-1,0.0);
         for(int i=0; i<nb_comp_-1; i++)
           H[i] = internal_energy_comp_[i] + p/rho_comp_[i];

         double H_J1(0.0), H_J2(0.0);
         for(int i=0; i<nb_comp_-1; i++)
         {
           H_J1 += H[i]*J1[i];        H_J2 += H[i]*J2[i];    
         }

         F1_dif_[rho_vx_indx_] = tau11;
         F1_dif_[rho_vy_indx_] = tau12;
         F1_dif_[rho_E_indx_] = tau11*vx + tau12*vy - q1 - H_J1;

         F2_dif_[rho_vx_indx_] = tau12;
         F2_dif_[rho_vy_indx_] = tau22;
         F2_dif_[rho_E_indx_] = tau12*vx + tau22*vy - q2 - H_J2;

         for(int i=0; i<nb_comp_-1; i++)
         {
           F1_dif_[rho_Y_indx_[i]] = J1[i];
           F2_dif_[rho_Y_indx_[i]] = J2[i];
         }

         if(ic_bc_.neumann_tau11(bd_nodes, side_flag).first)
         {
           K11_(rho_vx_indx_, vx_indx_) = 0.0;
           K12_(rho_vx_indx_, vy_indx_) = 0.0;
           K11_(rho_E_indx_, vx_indx_) = 0.0;
           K12_(rho_E_indx_, vy_indx_) = 0.0;
         }

         if(ic_bc_.neumann_tau22(bd_nodes, side_flag).first)
         {
           K21_(rho_vy_indx_, vx_indx_) = 0.0;
           K22_(rho_vy_indx_, vy_indx_) = 0.0;
           K21_(rho_E_indx_, vx_indx_) = 0.0;
           K22_(rho_E_indx_, vy_indx_) = 0.0;
         }

         if(ic_bc_.neumann_tau12(bd_nodes, side_flag).first)
         {
           K11_(rho_vy_indx_, vy_indx_) = 0.0;
           K12_(rho_vy_indx_, vx_indx_) = 0.0;
           K21_(rho_vx_indx_, vy_indx_) = 0.0;
           K22_(rho_vx_indx_, vx_indx_) = 0.0;
           K11_(rho_E_indx_, vy_indx_) = 0.0;
           K12_(rho_E_indx_, vx_indx_) = 0.0;
           K21_(rho_E_indx_, vy_indx_) = 0.0;
           K22_(rho_E_indx_, vx_indx_) = 0.0;
         }

         if(ic_bc_.neumann_q1(bd_nodes, side_flag).first)
           K11_(rho_E_indx_, T_indx_) =  0.0;

         if(ic_bc_.neumann_q2(bd_nodes, side_flag).first)
           K22_(rho_E_indx_, T_indx_) =  0.0;      

         for(int i=0; i<nb_comp_-1; i++)
         {
            if(ic_bc_.neumann_J1(i, bd_nodes, side_flag).first)
            {         
              K11_(rho_Y_indx_[i], rho_Y_indx_[i]) = 0.0;
              K11_(rho_E_indx_, rho_Y_indx_[i]) = 0.0;
            }

            if(ic_bc_.neumann_J2(i, bd_nodes, side_flag).first)
            {         
              K22_(rho_Y_indx_[i], rho_Y_indx_[i]) = 0.0;
              K22_(rho_E_indx_, rho_Y_indx_[i]) = 0.0;
            }
         }
     }








     void tau_matrix(double p, double vx,  double vy, double T, const vec &Y)
     {
         tau_.clear();

         const double v_nrm = sqrt(vx*vx + vy*vy);
         if(v_nrm < tol_) return;

         const double c(sound_speed_);
         vec v(2,0.0), alpha(2,0.0);
         v[0] = vx;    v[1] = vy;  

         mat del_v(2,2,0.0);
         for(int i=0; i<2; i++)
         {
             del_v(i,0) = dI_dx_[1+i];       del_v(i,1) = dI_dy_[1+i];     
         }
         auto h = quad_ptr_->compute_h(v, del_v, dummy_);

         if(setup_.tau_non_diag_comp_2001())
         {
             double lambda_sqr = alpha[0]*(v_nrm*v_nrm + c*c) + (alpha[1]) *0.5*c*c;
             lambda_sqr += 0.5*c*sqrt(c*c*alpha[1]*alpha[1] + 16*alpha[0]*alpha[0]*v_nrm*v_nrm);

             if(!setup_.steady_state()) lambda_sqr += 4/(dt_*dt_);
             const double lambda = 1.0/sqrt(lambda_sqr) + h/(2*v_nrm);

             double use(0.0);
             if(setup_.steady_state())  use = lambda;
             else  use = std::min(dt_*0.5, lambda);

             tau_(rho_indx_, rho_indx_) = use;
             tau_(rho_vx_indx_, rho_vx_indx_) = std::min(use, rho_*h*h/(mu_*12));
             tau_(rho_vy_indx_, rho_vy_indx_) = tau_(rho_vx_indx_, rho_vx_indx_);
             tau_(rho_E_indx_, rho_E_indx_) = std::min(use, rho_*h*h*cv_/(kappa_*12));

             for(int i=0; i<nb_comp_-1; i++)
               tau_(rho_Y_indx_[i], rho_Y_indx_[i]) = std::min(use, h*h/(chemical_diffusivity_comp_[i]*12));

             const mat DYDU = inverse(DUDY_);
             tau_= prod(DYDU,tau_);

             // tau correction for weak compressibility (if M < 0.1)
             if(v_nrm/c < 0.1)
             {
               const auto G = quad_ptr_->G();
               tau_(rho_indx_, rho_indx_) = pow( pow(tau_(0,0),-1) + pow(rho_*tau_(1,1)*(G(0,0)+G(1,1)),-1), -1); 
               tau_(rho_indx_, rho_vx_indx_) = 0.0;
               tau_(rho_indx_, rho_vy_indx_) = 0.0;
               tau_(rho_indx_, rho_E_indx_) = 0.0;
               for(int i=0; i<nb_comp_-1; i++)
                 tau_(rho_indx_, rho_Y_indx_[i]) = 0.0;
             }
         }


         else if(setup_.tau_non_diag_comp_2019())
         {
             double use(0.0);
             if(setup_.steady_state())  use = h/(2.0*(v_nrm + c)) + h/(2.0*v_nrm);
             else  use = std::min(dt_*0.5, h/(2.0*(v_nrm + c)) + h/(2.0*v_nrm));

             tau_(rho_indx_, rho_indx_) = use;
             tau_(rho_vx_indx_, rho_vx_indx_) = std::min(use, h*h*rho_/(12.0*mu_));
             tau_(rho_vy_indx_, rho_vy_indx_) = tau_(rho_vx_indx_, rho_vx_indx_);
             tau_(rho_E_indx_, rho_E_indx_) = std::min(use, h*h*rho_*cv_/(12.0*kappa_));

             for(int i=0; i<nb_comp_-1; i++)
               tau_(rho_Y_indx_[i], rho_Y_indx_[i]) = std::min(use, h*h/(12.0*chemical_diffusivity_comp_[i]));

             const mat DYDU = inverse(DUDY_);
             tau_= prod(DYDU,tau_);            

             // tau correction for weak compressibility (if M < 0.1)
             if(v_nrm/c < 0.1)
             {
               const auto G = quad_ptr_->G();
               tau_(rho_indx_, rho_indx_) = pow( pow(tau_(0,0),-1) + pow(rho_*tau_(1,1)*(G(0,0)+G(1,1)),-1), -1); 
               tau_(rho_indx_, rho_vx_indx_) = 0.0;
               tau_(rho_indx_, rho_vy_indx_) = 0.0;
               tau_(rho_indx_, rho_E_indx_) = 0.0;
               for(int i=0; i<nb_comp_-1; i++)
                 tau_(rho_indx_, rho_Y_indx_[i]) = 0.0;
             }
         }


         else if(setup_.tau_diag_incomp_2007())
         {
          const auto G = quad_ptr_->G();
          double use(0.0);
          for(int i=0; i<2; i++)
           for(int j=0; j<2; j++)
            use += G(i,j)*G(i,j);
           
          using namespace boost::numeric::ublas; 
          double tau_m = inner_prod(v, prod(G, v)) + 9.0*use*std::pow(mu_/rho_, 2);
          if(!setup_.steady_state()) tau_m += 4.0/(dt_*dt_);
          tau_m = 1.0/sqrt(tau_m); 

          const auto g = quad_ptr_->g();
          const double tau_c = 1.0/(rho_*tau_m*inner_prod(g,g));

          double tau_e = inner_prod(v, prod(G, v)) + 9.0*use*std::pow(kappa_/(rho_*cv_), 2);
          if(!setup_.steady_state()) tau_e += 4.0/(dt_*dt_);
          tau_e = 1.0/sqrt(tau_e); 
          
          tau_(rho_indx_, rho_indx_) = tau_c;
          tau_(rho_vx_indx_, rho_vx_indx_) = tau_m/rho_;
          tau_(rho_vy_indx_, rho_vy_indx_) = tau_m/rho_;
          tau_(rho_E_indx_, rho_E_indx_) = tau_e/(rho_*cv_);

          for(int i=0; i<nb_comp_-1; i++)
            tau_(rho_Y_indx_[i], rho_Y_indx_[i]) = tau_m;

         }


         else if(setup_.tau_diag_2014())
         {
            auto A0 = diag_mat_of_mat(DUDY_);
            auto A1 = diag_mat_of_mat(A1_);
            auto A2 = diag_mat_of_mat(A2_);
            auto K11 = diag_mat_of_mat(K11_);
            auto K12 = diag_mat_of_mat(K12_);
            auto K21 = diag_mat_of_mat(K21_);
            auto K22 = diag_mat_of_mat(K22_);
                        
            const auto G = quad_ptr_->G();
            mat tau_t_inv_sqr = 4.0/(dt_*dt_)*prod(A0,A0); 
            mat tau_a_inv_sqr = G(0,0)*prod(A1,A1) + G(0,1)*prod(A1,A2) + G(1,0)*prod(A2,A1) + G(1,1)*prod(A2,A2);
            mat tau_d_inv_sqr = G(0,0)*G(0,0)*prod(K11,K11) + G(0,1)*G(0,1)*prod(K12,K12) + G(1,0)*G(1,0)*prod(K21,K21) + G(1,1)*G(1,1)*prod(K22,K22);
            
            vec tau_t_inv(nb_dof_fluid_, 0.0), tau_a_inv(nb_dof_fluid_, 0.0), tau_d_inv(nb_dof_fluid_, 0.0);
            for(int i=0;i<nb_dof_fluid_;i++)
            {
              tau_t_inv[i] = sqrt(tau_t_inv_sqr(i,i));
              tau_a_inv[i] = sqrt(tau_a_inv_sqr(i,i));
              tau_d_inv[i] = sqrt(tau_d_inv_sqr(i,i));
            }
                              
            for(int i=0;i<nb_dof_fluid_;i++)
               tau_(i,i) = std::pow(tau_t_inv[i] + tau_a_inv[i] + tau_d_inv[i] + 1.e-8, -1);

            tau_(0,0) = std::pow(tau_(0,0) + rho_*tau_(1,1)*(G(0,0)+G(1,1)), -1);            
         }
     }




     void dc(double p, double vx, double vy, double T, const vec &Y)
     {
         if (setup_.dc_2006())
         {
             shock_matrix_.clear();
             const vec dU_dx = prod(DUDY_, dI_dx_);
             const vec dU_dy = prod(DUDY_, dI_dy_);
             const double del_rho_nrm = sqrt(pow(dU_dx[rho_indx_],2) + pow(dU_dy[rho_indx_],2));
             if(del_rho_nrm < tol_) return;
             vec J(2,0.0);
             J[0] = dU_dx[rho_indx_]/del_rho_nrm;
             J[1] = dU_dy[rho_indx_]/del_rho_nrm;
             double h(0.0);
             for(int i=0; i<nb_el_nodes_; i++)
             {
                 h += abs(J[0]*quad_ptr_->dsh_dx(i) + J[1]*quad_ptr_->dsh_dy(i) );
             }
             if(h > tol_) h = 1.0/h;
             else return;

             vec Z = prod(A1_,dI_dx_) + prod(A2_,dI_dy_) - S_;
//             if(!setup_.steady_state()) Z += prod(DUDY_,(I_-P_)/dt_);
             const double u_ref_inv_Z_nrm = norm_2(prod(u_ref_inv_,Z));
             const double u_ref_inv_U_nrm = norm_2(prod(u_ref_inv_,U_));

             const double use1 = pow(norm_2(prod(u_ref_inv_, dU_dx)),2) + pow(norm_2(prod(u_ref_inv_, dU_dy)),2);
             const double a = setup_.dc_sharp();
             double nu_shoc(0.0);
             if(a == 1.0) nu_shoc = u_ref_inv_Z_nrm*h*pow(use1, -0.5);
             else if(a == 2.0) nu_shoc = u_ref_inv_Z_nrm*h*h/u_ref_inv_U_nrm;
             else if(a == 1.5) nu_shoc = u_ref_inv_Z_nrm*h*(pow(use1, -0.5) + h/u_ref_inv_U_nrm)/2.0;

             mat K(nb_dof_fluid_,nb_dof_fluid_,0.0);
             for(int i=0; i<nb_dof_fluid_; i++)
             K(i,i) = nu_shoc;

             shock_matrix_ = setup_.dc_scale_fact()*prod(K, DUDY_);
         }
     }




 private:

     int nb_el_nodes_;
     int nb_dof_fluid_;
     int nb_comp_;

     ic_bc_type& ic_bc_;
     fluid_properties& props_;
     read_setup& setup_;

     double JxW_; 
     vec JNxW_;
     double dt_;


     vec P_orig_fluid_;
     vec I_orig_fluid_;
     vec P_;
     vec I_;
     vec dI_dx_;
     vec dI_dy_;

     double rho_;
     double mu_;
     double cp_, cv_;
     double kappa_;
     double sound_speed_;
     double beta_, alpha_;

    vec rho_comp_;
    vec internal_energy_comp_;
    vec chemical_diffusivity_comp_;

    std::vector<int> Y_indx_, rho_Y_indx_;

     vec U_P_;
     vec U_;
     vec F1_adv_;
     vec F2_adv_;
     vec F1_dif_;
     vec F2_dif_;
     vec S_;

     mat DUDY_;
     mat A1_,A2_;
     mat K11_,K12_,K21_,K22_;
     mat S0_;

     mat tau_;
     mat shock_matrix_;
     mat u_ref_inv_;

     int p_indx_, vx_indx_, vy_indx_, T_indx_;
     int rho_indx_, rho_vx_indx_, rho_vy_indx_, rho_E_indx_;

     vec gravity_;

     vec  r_;
     mat  m_;

     double tol_;
     std::shared_ptr<quad> quad_ptr_ = nullptr;
     point<2> dummy_;


  };



}/* namespace GALES */
#endif
