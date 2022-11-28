#ifndef _GALES_QUAD_HPP_
#define _GALES_QUAD_HPP_

#include <cmath>
#include <limits> 

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

#include "../base/read_setup.hpp"
#include "point.hpp"
#include "node.hpp"


namespace GALES {



  /**
       This class defines all functionality for the local parametric element such as
       - shape of local element 
       - shape functions, their local and global derivatives,  
       - mapping, jacobian of mapping
       - normal vectors for sides
       - data for numerical integration - gauss points, weights 
       - interpolate, x_derivative_interpolate, y_derivative_interpolate, z_derivative_interpolate functions for a vector
  */


  class quad 
  {

  public :

    using vec = boost::numeric::ublas::vector<double>;
    using mat = boost::numeric::ublas::matrix<double>;
    


    ///--------------------------------------------------------2D--------------------------------------------------------------------------
    ///------constructor for interior-------------------------------
    quad(const element<2>& el, const point<2>& pt)
    :
    nb_el_nodes_(el.nb_nodes()),
    sh_(el.nb_nodes()),
    dsh_dxi_(el.nb_nodes()), dsh_deta_(el.nb_nodes()), dsh_dzeta_(el.nb_nodes()),
    dsh_dx_(el.nb_nodes()), dsh_dy_(el.nb_nodes()), dsh_dz_(el.nb_nodes()),
    g_up_(2,2,0.0),
    g_(2, 0.0),
    G_(2,2,0.0)        
    {
       sh_and_local_derivatives(pt);
       calc(el);
       g_up_g_G(el);       
    }    
    //-------------------------------------------------------------------------------------------------------------------------------------
    





    ///--------------------------------------------------------2D--------------------------------------------------------------------------
    ///------ constructor for boundary-------------------------------
    quad(const element<2>& el, const point<2>& pt,  const point<2>& dummy)  
    :
    nb_el_nodes_(el.nb_nodes()),
    sh_(el.nb_nodes()),
    dsh_dxi_(el.nb_nodes()), dsh_deta_(el.nb_nodes()), dsh_dzeta_(el.nb_nodes()),
    dsh_dx_(el.nb_nodes()), dsh_dy_(el.nb_nodes()), dsh_dz_(el.nb_nodes()),
    N_(2,0.0), 
    JN_(2,0.0)
    {
       sh_and_local_derivatives(pt);
       calc(el);
       JN(el);       
    }    
    //-------------------------------------------------------------------------------------------------------------------------------------







    ///--------------------------------------------------------3D--------------------------------------------------------------------------
    ///------constructor for interior-------------------------------
    quad(const element<3>& el, const point<3>& pt)  // for interior
    :
    nb_el_nodes_(el.nb_nodes()),
    sh_(el.nb_nodes()),
    dsh_dxi_(el.nb_nodes()), dsh_deta_(el.nb_nodes()), dsh_dzeta_(el.nb_nodes()),
    dsh_dx_(el.nb_nodes()), dsh_dy_(el.nb_nodes()), dsh_dz_(el.nb_nodes()),
    g_up_(3,3,0.0),
    g_(3, 0.0),
    G_(3,3,0.0)        
    {
       sh_and_local_derivatives(pt);
       calc(el);
       g_up_g_G(el);       
    }    
    //-------------------------------------------------------------------------------------------------------------------------------------
    





    ///--------------------------------------------------------3D--------------------------------------------------------------------------
    ///------constructor for boundary-------------------------------
    quad(const element<3>& el, const point<3>& pt,  const point<3>& dummy)  
    :
    nb_el_nodes_(el.nb_nodes()),
    sh_(el.nb_nodes()),
    dsh_dxi_(el.nb_nodes()), dsh_deta_(el.nb_nodes()), dsh_dzeta_(el.nb_nodes()),
    dsh_dx_(el.nb_nodes()), dsh_dy_(el.nb_nodes()), dsh_dz_(el.nb_nodes()),
    N_(3,0.0), 
    JN_(3,0.0)
    {
       sh_and_local_derivatives(pt);
       calc(el);
       JN(el);       
    }    
    //-------------------------------------------------------------------------------------------------------------------------------------






    //-------------------------------------------------------------------------------------------------------------------------------------        
    /// Deleting the copy and move constructors - no duplication/transfer in anyway
    quad(const quad&) = delete;               //copy constructor
    quad& operator=(const quad&) = delete;    //copy assignment operator
    quad(quad&&) = delete;                    //move constructor  
    quad& operator=(quad&&) = delete;         //move assignment operator 
    //-------------------------------------------------------------------------------------------------------------------------------------


    
  



    ///--------------------------------------------------------2D--------------------------------------------------------------------------
    ///-----shape functions and their derivatives w.r.t local coordinates--------------------------------
    void sh_and_local_derivatives(const point<2>& pt)
    {
      xi_ = pt.get_x();
      eta_ = pt.get_y();
      
      if(nb_el_nodes_ == 4)
      {
        sh_[0] = 0.25*(1.0-xi_)*(1.0-eta_);
        sh_[1] = 0.25*(1.0+xi_)*(1.0-eta_);
        sh_[2] = 0.25*(1.0+xi_)*(1.0+eta_);
        sh_[3] = 0.25*(1.0-xi_)*(1.0+eta_);

        dsh_dxi_[0] = -0.25*(1.-eta_);         dsh_deta_[0] = -0.25*(1.-xi_);
        dsh_dxi_[1] =  0.25*(1.-eta_);         dsh_deta_[1] = -0.25*(1.+xi_); 
        dsh_dxi_[2] =  0.25*(1.+eta_);         dsh_deta_[2] =  0.25*(1.+xi_);
        dsh_dxi_[3] = -0.25*(1.+eta_);         dsh_deta_[3] =  0.25*(1.-xi_);

        W_in_ = 1.0;                               
        W_bd_ = 1.0;                               
      }
      else if(nb_el_nodes_ == 3)
      {
        sh_[0] = 1.0-xi_-eta_;
        sh_[1] = xi_;
        sh_[2] = eta_;

        dsh_dxi_[0] = -1.0;         dsh_deta_[0] = -1.0;
        dsh_dxi_[1] = 1.0;          dsh_deta_[1] = 0.0;
        dsh_dxi_[2] = 0.0;          dsh_deta_[2] = 1.0;                           

        W_in_ = 1.0/6.0;

        /// to compute the integrals we need to multiply det*weight
        /// according to gauss integration on [-1,1] weight is one 
        /// since we use additional transformation: det([x,y]-->[-1,1]) = det([x,y]-->[0,1]) * det([0,1]-->[-1,1])
        /// det([0,1]-->[-1,1]) = 1/2, and weight([-1,1]) = 1; 
        /// for implementation we set weight = 1/2 instead of 1, which accounts for the transformation defined above

        W_bd_ = 1.0/2.0;
      }        
    }
    //----------------------------------------------------------------------------------------------------





    ///--------------------------------------------------------2D--------------------------------------------------------------------------
    void calc(const element<2>& el)
    {
      ///clearing previous computed values of variables for computation at gauss point of the element
      dx_dxi_ = 0.0;     dx_deta_ = 0.0;  
      dy_dxi_ = 0.0;     dy_deta_ = 0.0;    
      x_ = 0.0;          y_ = 0.0; 

       /// mapping
      for(int i=0; i<nb_el_nodes_; i++)  
      {
        x_ += sh_[i] * el.get_x(i);
        y_ += sh_[i] * el.get_y(i);
        dx_dxi_   += dsh_dxi_[i]   * el.get_x(i);
        dx_deta_  += dsh_deta_[i]  * el.get_x(i);
        dy_dxi_   += dsh_dxi_[i]   * el.get_y(i);
        dy_deta_  += dsh_deta_[i]  * el.get_y(i);
      }

      J_= dx_dxi_*dy_deta_ - dy_dxi_*dx_deta_;
      if(J_ <= 0.0) 
      {
        Error("Jacobian is not positive - either element is too much stretched or element nodes are not anticlockwise");
      }

      dxi_dx_ =   (1/J_)*(dy_deta_);
      deta_dx_ =  (1/J_)*(-dy_dxi_);
      dxi_dy_ =   (1/J_)*(-dx_deta_);
      deta_dy_ =  (1/J_)*(dx_dxi_);

      //derivatives of shape functions w.r.t. global coordinates
      for(int i=0; i<nb_el_nodes_; i++)  
      {
        dsh_dx_[i] = dsh_dxi_[i]*dxi_dx_ + dsh_deta_[i]*deta_dx_;
        dsh_dy_[i] = dsh_dxi_[i]*dxi_dy_ + dsh_deta_[i]*deta_dy_;
      }
    } 
    //-------------------------------------------------------------------------------------------------------------------------------------










  
    ///--------------------------------------------------------2D--------------------------------------------------------------------------
    void JN(const element<2>& el)
    {    
      /// area vectors(JN_), jacobian(J_) and normal vectors(N_) 
      if(nb_el_nodes_ == 4)
      {
        if(eta_ == 1.0 )
        { 
           JN_[0] = -dy_dxi_;   
           JN_[1] = dx_dxi_;
        }
        else if(eta_ == -1.0)
        { 
           JN_[0] = dy_dxi_;    
           JN_[1] = -dx_dxi_;
        }
        else if(xi_ == 1.0)  
        { 
           JN_[0] = dy_deta_;   
           JN_[1] = -dx_deta_;
        }
        else if(xi_ == -1.0) 
        { 
           JN_[0] = -dy_deta_;  
           JN_[1] = dx_deta_;
        }
      }
      else if(nb_el_nodes_ == 3)
      {
        if(xi_ == 0.0)
        {
          JN_[0] = -dy_deta_;  
          JN_[1] = dx_deta_;
        }
        else if(eta_ == 0.0)
        {
           JN_[0] = dy_dxi_;    
           JN_[1] = -dx_dxi_;
        }
        else
        {
           JN_[0] = dy_deta_-dy_dxi_;  
           JN_[1] = dx_dxi_-dx_deta_;
        }
      }
      
      J_bd_ = sqrt(JN_[0]*JN_[0] + JN_[1]*JN_[1]);
      N_ = JN_/J_bd_;
    }
    //-------------------------------------------------------------------------------------------------------------------------------------
  







    ///--------------------------------------------------------2D--------------------------------------------------------------------------
    /// Computation of element metric tensors
    /// These are used to cmpute the element size along the streamline direction
    void g_up_g_G(const element<2>& el)
    {
      /*  _         _      _         _ T
         |           |    |           |
         |  d(x,y)   |    |  d(x,y)   |  
         |  ------   | *  |  ------   | 
         | d(xi,eta) |    | d(xi,eta) |
         |_         _|    |_         _|
      */
      g_up_(0,0) = dx_dxi_*dx_dxi_ + dx_deta_*dx_deta_;
      g_up_(0,1) = dx_dxi_*dy_dxi_ + dx_deta_*dy_deta_;
      g_up_(1,0) = g_up_(0,1);
      g_up_(1,1) = dy_dxi_*dy_dxi_ + dy_deta_*dy_deta_;
      

      /*           _          _     
                  |            |    
                  |  d(xi,eta) |      
        trace of  |  --------  |    
                  |   d(x,y)   |    
                  |_          _|    
      */
      
      g_[0] = dxi_dx_ + deta_dx_;
      g_[1] = dxi_dy_ + deta_dy_;
      



      /*  _          _ T    _         _ 
         |            |    |           |
         | d(xi,eta)  |    | d(xi,eta) |  
         |  ------    | *  |  ------   | 
         |  d(x,y)    |    |  d(x,y)   |
         |_          _|    |_         _|
      */
      G_(0,0) = dxi_dx_*dxi_dx_ + deta_dx_*deta_dx_;
      G_(1,1) = dxi_dy_*dxi_dy_ + deta_dy_*deta_dy_;
      G_(0,1) = dxi_dx_*dxi_dy_ + deta_dx_*deta_dy_;
      G_(1,0) = G_(0,1);
    }
    //-------------------------------------------------------------------------------------------------------------------------------------








    ///--------------------------------------------------------2D--------------------------------------------------------------------------
    double compute_h(const vec &v, const mat &del_v, const point<2>& dummy) const
    {
      vec alpha(2,0.0);   //The unit vector (alpha[i]) calculation in direction of the streamline system of reference
         
      using namespace boost::numeric::ublas;
      const double nrm_v = norm_2(v);
      vec s0(2,0.0);
      for(int i=0; i<2; i++)
        s0[i] = v[i]/nrm_v;
      alpha[0] = inner_prod(s0, prod(G_, s0));

      const vec acc = prod(del_v, v);    /// convective acc = v.del(v)
      vec s1 = acc - inner_prod(acc,s0)*s0;
      const double nrm_s1 = norm_2(s1);
      if(nrm_s1 > std::numeric_limits<double>::epsilon() )
      {
        for(int i=0; i<2; i++)
          s1[i] /= nrm_s1;
      }
      alpha[1] = inner_prod(s1, prod(G_, s1));
         
      double h(0.0);
      if(alpha[0] != 0.0)  h += 1.0/alpha[0];
      if(alpha[1] != 0.0)  h += 1.0/alpha[1];
      h = 2*sqrt(abs(h));    

      return h;
    }
    //-------------------------------------------------------------------------------------------------------------------------------------

  

    



    ///--------------------------------------------------------2D--------------------------------------------------------------------------
    /// For mesh moving problems this returns the jacobian of mapping from local el to global el before mesh motion
    double det_low(const element<2>& el, const vec& P_ux_e, const vec& P_uy_e)
    {
       double dx_dxi(0.0), dx_deta(0.0);
       double dy_dxi(0.0), dy_deta(0.0);

       for (int i=0; i<nb_el_nodes_; i++)
       {
         dx_dxi   += dsh_dxi_[i]   * (el.get_x(i) - P_ux_e[i]);
         dx_deta  += dsh_deta_[i]  * (el.get_x(i) - P_ux_e[i]);
         dy_dxi   += dsh_dxi_[i]   * (el.get_y(i) - P_uy_e[i]);
         dy_deta  += dsh_deta_[i]  * (el.get_y(i) - P_uy_e[i]);
       }
      double det_low = dx_dxi*dy_deta - dy_dxi*dx_deta;
      return det_low;
    }  
    //-------------------------------------------------------------------------------------------------------------------------------------







    ///--------------------------------------------------------2D--------------------------------------------------------------------------
    /// interpolated mesh velocity
    vec mesh_v(const vec& P_ux_e, const vec& P_uy_e, double delta_t)
    {
      vec mesh_u(2,0.0);
      interpolate(P_ux_e, mesh_u[0]);
      interpolate(P_uy_e, mesh_u[1]);
      return mesh_u/delta_t;
    }
    //-------------------------------------------------------------------------------------------------------------------------------------





    
    
    

































    ///--------------------------------------------------------3D--------------------------------------------------------------------------
    ///-----shape functions and their derivatives w.r.t local coordinates--------------------------------
    void sh_and_local_derivatives(const point<3>& pt)
    {
      xi_ = pt.get_x();
      eta_ = pt.get_y();
      zeta_ = pt.get_z();
      
      if(nb_el_nodes_ == 8)
      {
        sh_[0] = 0.125*(1.0-xi_)*(1.0-eta_)*(1.0-zeta_);
        sh_[1] = 0.125*(1.0+xi_)*(1.0-eta_)*(1.0-zeta_);
        sh_[2] = 0.125*(1.0+xi_)*(1.0+eta_)*(1.0-zeta_);
        sh_[3] = 0.125*(1.0-xi_)*(1.0+eta_)*(1.0-zeta_);
        sh_[4] = 0.125*(1.0-xi_)*(1.0-eta_)*(1.0+zeta_);
        sh_[5] = 0.125*(1.0+xi_)*(1.0-eta_)*(1.0+zeta_);
        sh_[6] = 0.125*(1.0+xi_)*(1.0+eta_)*(1.0+zeta_);
        sh_[7] = 0.125*(1.0-xi_)*(1.0+eta_)*(1.0+zeta_);

        dsh_dxi_[0] = -0.125*(1.0-eta_)*(1.0-zeta_);             dsh_deta_[0] = -0.125*(1.0-xi_)*(1.0-zeta_);             dsh_dzeta_[0] = -0.125*(1.0-xi_)*(1.0-eta_);
        dsh_dxi_[1] =  0.125*(1.0-eta_)*(1.0-zeta_);             dsh_deta_[1] = -0.125*(1.0+xi_)*(1.0-zeta_);             dsh_dzeta_[1] = -0.125*(1.0+xi_)*(1.0-eta_);
        dsh_dxi_[2] =  0.125*(1.0+eta_)*(1.0-zeta_);             dsh_deta_[2] =  0.125*(1.0+xi_)*(1.0-zeta_);             dsh_dzeta_[2] = -0.125*(1.0+xi_)*(1.0+eta_);
        dsh_dxi_[3] = -0.125*(1.0+eta_)*(1.0-zeta_);             dsh_deta_[3] =  0.125*(1.0-xi_)*(1.0-zeta_);             dsh_dzeta_[3] = -0.125*(1.0-xi_)*(1.0+eta_);
        dsh_dxi_[4] = -0.125*(1.0-eta_)*(1.0+zeta_);             dsh_deta_[4] = -0.125*(1.0-xi_)*(1.0+zeta_);             dsh_dzeta_[4] =  0.125*(1.0-xi_)*(1.0-eta_);
        dsh_dxi_[5] =  0.125*(1.0-eta_)*(1.0+zeta_);             dsh_deta_[5] = -0.125*(1.0+xi_)*(1.0+zeta_);             dsh_dzeta_[5] =  0.125*(1.0+xi_)*(1.0-eta_);
        dsh_dxi_[6] =  0.125*(1.0+eta_)*(1.0+zeta_);             dsh_deta_[6] =  0.125*(1.0+xi_)*(1.0+zeta_);             dsh_dzeta_[6] =  0.125*(1.0+xi_)*(1.0+eta_);
        dsh_dxi_[7] = -0.125*(1.0+eta_)*(1.0+zeta_);             dsh_deta_[7] =  0.125*(1.0-xi_)*(1.0+zeta_);             dsh_dzeta_[7] =  0.125*(1.0-xi_)*(1.0+eta_);


        W_in_ = 1.0;                               
        W_bd_ = 1.0;                               
      }
      else if(nb_el_nodes_ == 4)
      {
        sh_[0] = 1.0-xi_-eta_-zeta_;
        sh_[1] = xi_;
        sh_[2] = eta_;
        sh_[3] = zeta_;

        dsh_dxi_[0] = -1.0;        dsh_deta_[0] = -1.0;          dsh_dzeta_[0] = -1.0;
        dsh_dxi_[1] = 1.0;         dsh_deta_[1] = 0.0;           dsh_dzeta_[1] = 0.0;
        dsh_dxi_[2] = 0.0;         dsh_deta_[2] = 1.0;           dsh_dzeta_[2] = 0.0;
        dsh_dxi_[3] = 0.0;         dsh_deta_[3] = 0.0;           dsh_dzeta_[3] = 1.0;

        W_in_ = 1.0/24.0;
        W_bd_ = 1.0/6.0;
      } 
    }    
    ///-----------------------------------------------------------------------------------------------





    ///--------------------------------------------------------3D--------------------------------------------------------------------------
    void calc(const element<3>& el)
    {
      ///clearing previous computed values of variables for computation at gauss point of the element
      dx_dxi_ = 0.0;     dx_deta_ = 0.0;      dx_dzeta_ = 0.0;
      dy_dxi_ = 0.0;     dy_deta_ = 0.0;      dy_dzeta_ = 0.0;
      dz_dxi_ = 0.0;     dz_deta_ = 0.0;      dz_dzeta_ = 0.0;
      x_ = 0.0;          y_ = 0.0;            z_ = 0.0;  

      /// mapping
      for(int i=0; i<nb_el_nodes_; i++)  
      {
        x_ += sh_[i] * el.get_x(i);
        y_ += sh_[i] * el.get_y(i);
        z_ += sh_[i] * el.get_z(i);
        dx_dxi_   += dsh_dxi_[i]   * el.get_x(i);
        dx_deta_  += dsh_deta_[i]  * el.get_x(i);
        dx_dzeta_ += dsh_dzeta_[i] * el.get_x(i);
        dy_dxi_   += dsh_dxi_[i]   * el.get_y(i);
        dy_deta_  += dsh_deta_[i]  * el.get_y(i);
        dy_dzeta_ += dsh_dzeta_[i] * el.get_y(i);
        dz_dxi_   += dsh_dxi_[i]   * el.get_z(i);
        dz_deta_  += dsh_deta_[i]  * el.get_z(i);
        dz_dzeta_ += dsh_dzeta_[i] * el.get_z(i);
      }

      J_= dx_dxi_*(dy_deta_*dz_dzeta_-dz_deta_*dy_dzeta_) - dy_dxi_*(dx_deta_*dz_dzeta_-dz_deta_*dx_dzeta_) +  dz_dxi_*(dx_deta_*dy_dzeta_-dy_deta_*dx_dzeta_);
      if(J_ <= 0.0) 
      {
        Error("Jacobian is not positive - either element is too much stretched or element nodes are not anticlockwise");
      }

      dxi_dx_ =   (1/J_)*( dy_deta_*dz_dzeta_ - dz_deta_*dy_dzeta_);
      deta_dx_ =  (1/J_)*(-dy_dxi_*dz_dzeta_  + dz_dxi_*dy_dzeta_);
      dzeta_dx_ = (1/J_)*( dy_dxi_*dz_deta_   - dz_dxi_*dy_deta_);
      dxi_dy_ =   (1/J_)*(-dx_deta_*dz_dzeta_ + dz_deta_*dx_dzeta_);
      deta_dy_ =  (1/J_)*( dx_dxi_*dz_dzeta_  - dz_dxi_*dx_dzeta_);
      dzeta_dy_ = (1/J_)*(-dx_dxi_*dz_deta_   + dz_dxi_*dx_deta_);
      dxi_dz_ =   (1/J_)*( dx_deta_*dy_dzeta_ - dy_deta_*dx_dzeta_);
      deta_dz_ =  (1/J_)*(-dx_dxi_*dy_dzeta_  + dy_dxi_*dx_dzeta_);
      dzeta_dz_ = (1/J_)*( dx_dxi_*dy_deta_   - dy_dxi_*dx_deta_);

      //derivatives of shape functions w.r.t. global coordinates
      for(int i=0; i<nb_el_nodes_; i++)  
      {
        dsh_dx_[i] = dsh_dxi_[i]*dxi_dx_ + dsh_deta_[i]*deta_dx_ + dsh_dzeta_[i]*dzeta_dx_;
        dsh_dy_[i] = dsh_dxi_[i]*dxi_dy_ + dsh_deta_[i]*deta_dy_ + dsh_dzeta_[i]*dzeta_dy_;
        dsh_dz_[i] = dsh_dxi_[i]*dxi_dz_ + dsh_deta_[i]*deta_dz_ + dsh_dzeta_[i]*dzeta_dz_;
      }
    }
    //-------------------------------------------------------------------------------------------------------------------------------------







    ///--------------------------------------------------------3D--------------------------------------------------------------------------
    /// area vectors(JN_), jacobian(J_) and normal vectors(N_) 
    void JN(const element<3>& el)
    {
      if(nb_el_nodes_ == 8)
      {
        if(zeta_ == 1.0)
        {
          JN_[0] = dy_dxi_*dz_deta_-dz_dxi_*dy_deta_;
          JN_[1] = dz_dxi_*dx_deta_-dx_dxi_*dz_deta_;
          JN_[2] = dx_dxi_*dy_deta_-dy_dxi_*dx_deta_;
        }
        else if(zeta_ == -1.0)
        {
          JN_[0] = -dy_dxi_*dz_deta_+dz_dxi_*dy_deta_;
          JN_[1] = -dz_dxi_*dx_deta_+dx_dxi_*dz_deta_;
          JN_[2] = -dx_dxi_*dy_deta_+dy_dxi_*dx_deta_;
        }
        else if(eta_ == 1.0 )
        {
          JN_[0] = dy_dzeta_*dz_dxi_-dy_dxi_*dz_dzeta_;
          JN_[1] = dz_dzeta_*dx_dxi_-dz_dxi_*dx_dzeta_;
          JN_[2] = dx_dzeta_*dy_dxi_-dx_dxi_*dy_dzeta_;
        }
        else if(eta_ == -1.0)
        {
          JN_[0] = -dy_dzeta_*dz_dxi_+dy_dxi_*dz_dzeta_;
          JN_[1] = -dz_dzeta_*dx_dxi_+dz_dxi_*dx_dzeta_;
          JN_[2] = -dx_dzeta_*dy_dxi_+dx_dxi_*dy_dzeta_;
        }
        else if(xi_ == 1.0 )
        {
          JN_[0] = dy_deta_*dz_dzeta_-dy_dzeta_*dz_deta_;
          JN_[1] = dz_deta_*dx_dzeta_-dz_dzeta_*dx_deta_;
          JN_[2] = dx_deta_*dy_dzeta_-dx_dzeta_*dy_deta_;
        }
        else if(xi_ == -1.0)
        {
          JN_[0] = -dy_deta_*dz_dzeta_+dy_dzeta_*dz_deta_;
          JN_[1] = -dz_deta_*dx_dzeta_+dz_dzeta_*dx_deta_;
          JN_[2] = -dx_deta_*dy_dzeta_+dx_dzeta_*dy_deta_;
        }
      }

      else if(nb_el_nodes_ == 4)
      {
        if(xi_ == 0.0)
        {
          JN_[0] =  dz_deta_*dy_dzeta_ - dy_deta_*dz_dzeta_;
          JN_[1] =  dx_deta_*dz_dzeta_ - dz_deta_*dx_dzeta_;  
          JN_[2] =  dy_deta_*dx_dzeta_ - dx_deta_*dy_dzeta_;             
        }
        else if(eta_ == 0.0)
        {
          JN_[0] =  dz_dzeta_*dy_dxi_ - dy_dzeta_*dz_dxi_;
          JN_[1] =  dz_dxi_*dx_dzeta_ - dx_dxi_*dz_dzeta_;  
          JN_[2] =  dy_dzeta_*dx_dxi_ - dx_dzeta_*dy_dxi_;   
        }   
        else if(zeta_ == 0.0)
        {
          JN_[0] =  dz_dxi_*dy_deta_ - dy_dxi_*dz_deta_;
          JN_[1] =  dx_dxi_*dz_deta_ - dx_deta_*dz_dxi_;   
          JN_[2] =  dy_dxi_*dx_deta_ - dx_dxi_*dy_deta_;
        }             
        else
        {
          JN_[0] =  (dz_dzeta_-dz_dxi_)*(dy_deta_-dy_dxi_) - (dy_dzeta_-dy_dxi_)*(dz_deta_-dz_dxi_); 
          JN_[1] =  (dx_dzeta_-dx_dxi_)*(dz_deta_-dz_dxi_) - (dx_deta_-dx_dxi_)*(dz_dzeta_-dz_dxi_);  
          JN_[2] =  (dx_deta_-dx_dxi_)*(dy_dzeta_-dy_dxi_) - (dx_dzeta_-dx_dxi_)*(dy_deta_-dy_dxi_);   
        }   
      }  

      J_bd_ = sqrt(JN_[0]*JN_[0] + JN_[1]*JN_[1] + JN_[2]*JN_[2]);
      N_ = JN_/J_bd_;
    }
    //-------------------------------------------------------------------------------------------------------------------------------------



  
  


    ///--------------------------------------------------------3D--------------------------------------------------------------------------
    /// Computation of element metric tensors
    /// These are used to cmpute the element size along the streamline direction
    void g_up_g_G(const element<3>& el)
    {
      g_up_(0,0) = dx_dxi_*dx_dxi_ + dx_deta_*dx_deta_ + dx_dzeta_*dx_dzeta_;
      g_up_(0,1) = dx_dxi_*dy_dxi_ + dx_deta_*dy_deta_ + dx_dzeta_*dy_dzeta_;
      g_up_(0,2) = dx_dxi_*dz_dxi_ + dx_deta_*dz_deta_ + dx_dzeta_*dz_dzeta_;
      g_up_(1,0) = g_up_(0,1);
      g_up_(1,1) = dy_dxi_*dy_dxi_ + dy_deta_*dy_deta_ + dy_dzeta_*dy_dzeta_;
      g_up_(1,2) = dy_dxi_*dz_dxi_ + dy_deta_*dz_deta_ + dy_dzeta_*dz_dzeta_;
      g_up_(2,0) = g_up_(0,2);
      g_up_(2,1) = g_up_(1,2);
      g_up_(2,2) = dz_dxi_*dz_dxi_ + dz_deta_*dz_deta_ + dz_dzeta_*dz_dzeta_;

      g_[0] = dxi_dx_ + deta_dx_ + dzeta_dx_;
      g_[1] = dxi_dy_ + deta_dy_ + dzeta_dy_;
      g_[2] = dxi_dz_ + deta_dz_ + dzeta_dz_;

      G_(0,0) = dxi_dx_*dxi_dx_ + deta_dx_*deta_dx_ + dzeta_dx_*dzeta_dx_;
      G_(1,1) = dxi_dy_*dxi_dy_ + deta_dy_*deta_dy_ + dzeta_dy_*dzeta_dy_;
      G_(2,2) = dxi_dz_*dxi_dz_ + deta_dz_*deta_dz_ + dzeta_dz_*dzeta_dz_;
      G_(0,1) = dxi_dx_*dxi_dy_ + deta_dx_*deta_dy_ + dzeta_dx_*dzeta_dy_;
      G_(1,2) = dxi_dy_*dxi_dz_ + deta_dy_*deta_dz_ + dzeta_dy_*dzeta_dz_;
      G_(0,2) = dxi_dz_*dxi_dx_ + deta_dz_*deta_dx_ + dzeta_dz_*dzeta_dx_;
      G_(1,0) = G_(0,1);
      G_(2,1) = G_(1,2);
      G_(2,0) = G_(0,2);
    }
    //-------------------------------------------------------------------------------------------------------------------------------------







    ///--------------------------------------------------------3D--------------------------------------------------------------------------
    double compute_h(const vec &v, const mat &del_v, const point<3>& dummy) const
    {
      vec alpha(3,0.0);
      using namespace boost::numeric::ublas;
      const double nrm_v = norm_2(v);
      vec s0(3,0.0);
      for(int i=0; i<3; i++)
        s0[i] = v[i]/nrm_v;
      alpha[0] = inner_prod(s0, prod(G_, s0));

      const vec acc = prod(del_v, v);    /// convective acc = v.del(v)
      vec s1 = acc - inner_prod(acc,s0)*s0;
      const double nrm_s1 = norm_2(s1);
      if(nrm_s1 > std::numeric_limits<double>::epsilon() )
      {
        for(int i=0; i<3; i++)
          s1[i] /= nrm_s1;
      }
      alpha[1] = inner_prod(s1, prod(G_, s1));

      vec s2(3,0.0);
      s2[0] = s0[1]*s1[2] - s0[2]*s1[1];
      s2[1] = s1[0]*s0[2] - s1[2]*s0[0];
      s2[2] = s0[0]*s1[1] - s0[1]*s1[0];
      const double nrm_s2 = norm_2(s2);
      if(nrm_s2 > std::numeric_limits<double>::epsilon() )
      {
        for(int i=0; i<3; i++)
          s2[i] /= nrm_s2;
      }
      alpha[2] = inner_prod(s2, prod(G_, s2));

      double h(0.0);
      if(alpha[0] != 0.0)  h += 1.0/alpha[0];
      if(alpha[1] != 0.0)  h += 1.0/alpha[1];
      if(alpha[2] != 0.0)  h += 1.0/alpha[2];
      h = 2*sqrt(abs(h));    
      return h;
    }
    //-------------------------------------------------------------------------------------------------------------------------------------




     

    ///--------------------------------------------------------3D--------------------------------------------------------------------------
    /// for mesh moving problems this returns the jacobian of mapping from local el to global el before mesh motion
    double det_low(const element<3>& el, const vec& P_ux_e, const vec& P_uy_e, const vec& P_uz_e)
    {
       double dx_dxi(0.0), dx_deta(0.0), dx_dzeta(0.0);
       double dy_dxi(0.0), dy_deta(0.0), dy_dzeta(0.0);
       double dz_dxi(0.0), dz_deta(0.0), dz_dzeta(0.0);

       for (int i=0; i<nb_el_nodes_; i++)
       {
         dx_dxi   += dsh_dxi_[i]   * (el.get_x(i) - P_ux_e[i]);
         dx_deta  += dsh_deta_[i]  * (el.get_x(i) - P_ux_e[i]);
         dx_dzeta += dsh_dzeta_[i] * (el.get_x(i) - P_ux_e[i]);

         dy_dxi   += dsh_dxi_[i]   * (el.get_y(i) - P_uy_e[i]);
         dy_deta  += dsh_deta_[i]  * (el.get_y(i) - P_uy_e[i]);
         dy_dzeta += dsh_dzeta_[i] * (el.get_y(i) - P_uy_e[i]);

         dz_dxi   += dsh_dxi_[i]   * (el.get_z(i) - P_uz_e[i]);
         dz_deta  += dsh_deta_[i]  * (el.get_z(i) - P_uz_e[i]);
         dz_dzeta += dsh_dzeta_[i] * (el.get_z(i) - P_uz_e[i]);
       }
      double det_low = dx_dxi*(dy_deta*dz_dzeta-dz_deta*dy_dzeta) -  dy_dxi*(dx_deta*dz_dzeta-dz_deta*dx_dzeta) + dz_dxi*(dx_deta*dy_dzeta-dy_deta*dx_dzeta);
      return det_low;
    }  
    //-------------------------------------------------------------------------------------------------------------------------------------







    ///--------------------------------------------------------3D--------------------------------------------------------------------------
    /// interpolated mesh velocity
    vec mesh_v(const vec& P_ux_e, const vec& P_uy_e, const vec& P_uz_e, double delta_t)
    {
      vec mesh_u(3,0.0);
      interpolate(P_ux_e, mesh_u[0]);
      interpolate(P_uy_e, mesh_u[1]);
      interpolate(P_uz_e, mesh_u[2]);
      return mesh_u/delta_t;
    }
    //-------------------------------------------------------------------------------------------------------------------------------------




   










    //-----------------------------------------------------------2D/3D----------------------------------------------------------
    /// This interpolates a vector "v" into a double "a".
    void interpolate(const vec& v, double &a) const
    {
      a = 0.0;
      for (int i=0; i<nb_el_nodes_; i++)
	a += sh_[i] * v[i];
    }


    /// This interpolates a vector "v" into a vector "a".
    void interpolate(const vec& v, vec &a)  const
    {
      std::fill(a.begin(),a.end(), 0.0);
      int nb_dof (a.size());
      for (int i=0; i<nb_el_nodes_; i++)
        for (int j=0; j<nb_dof; j++) 
  	  a[j] += sh_[i] * v[nb_dof*i+j];
    }


    /// This interpolates the x-derivative of a vector "v" into a double "a".
    void x_derivative_interpolate(const vec& v, double &a) const 
    {
      a = 0.0;
      for (int i=0; i<nb_el_nodes_; i++)
	a += dsh_dx_[i] * v[i];
    }


    /// This interpolates the x-derivative of a vector "v" into a vector "a".
    void x_derivative_interpolate(const vec& v,  vec& a) const 
    {
      std::fill(a.begin(),a.end(), 0.0);
      int nb_dof (a.size());
      for (int i=0; i<nb_el_nodes_; i++)
	for (int j=0; j<nb_dof; j++) 
	  a[j] += dsh_dx_[i]* v[nb_dof*i+j];
    }


    /// This interpolates the y-derivative of a vector "v" into a double "a".
    void y_derivative_interpolate(const vec&v, double &a) const 
    {
      a = 0.0;
      for (int i=0; i<nb_el_nodes_; i++)
	a += dsh_dy_[i] * v[i];
    }



    /// This interpolates the y-derivative of a vector "v" into a vector "a".
    void y_derivative_interpolate(const vec& v, vec& a) const 
    {
      std::fill(a.begin(),a.end(), 0.0);
      int nb_dof (a.size());
      for (int i=0; i<nb_el_nodes_; i++)
	for (int j=0; j<nb_dof; j++) 
	  a[j] += dsh_dy_[i]* v[nb_dof*i+j];
    }


    /// This interpolates the z-derivative of a vector "v" into a double "a".
    void z_derivative_interpolate(const vec& v, double &a) const 
    {
      a = 0.0;
      for (int i=0; i<nb_el_nodes_; i++)
	a += dsh_dz_[i] * v[i];
    }


    /// This interpolates the z-derivative of a vector "v" into a vector "a".
    void z_derivative_interpolate(const vec& v, vec &a) const 
    {
      std::fill(a.begin(),a.end(), 0.0);
      int nb_dof (a.size());
      for (int i=0; i<nb_el_nodes_; i++)
	for (int j=0; j<nb_dof; j++) 
	  a[j] += dsh_dz_[i]* v[nb_dof*i+j];
    }
    
   
   
    /// This function returns the normal component of velocity vector
    auto normal_component_of_v(const vec& v) const
    {            
        return boost::numeric::ublas::inner_prod(v, N_); 
    }    


    /// This function returns the normal component of velocity * normal vector
    /// note that v = v_normal * unit_normal + v_tangential * unit_tangential
    auto normal_v(const vec& v) const
    {            
        return normal_component_of_v(v)*N_; 
    }    


    /// This function returns the tangential component of velocity * unit_tangential vector
    auto tangential_v(const vec& v) const
    {            
        return v - normal_v(v); 
    }    








    ///--------------------------------------------------------inspector functions-------------------------------------------
    auto nb_el_nodes()const {return nb_el_nodes_;}
    auto x()const {return x_;}
    auto y()const {return y_;}
    auto z()const {return z_;}
    auto sh(int i)const {return sh_[i];}
    auto sh()const {return sh_;}
    auto dsh_dxi(int i)const {return dsh_dxi_[i];}
    auto dsh_dxi()const {return dsh_dxi_;}
    auto dsh_deta(int i)const {return dsh_deta_[i];}
    auto dsh_deta()const {return dsh_deta_;}
    auto dsh_dzeta(int i)const {return dsh_dzeta_[i];}
    auto dsh_dzeta()const {return dsh_dzeta_;}    
    auto dsh_dx(int i)const {return dsh_dx_[i];}
    auto dsh_dx()const {return dsh_dx_;}
    auto dsh_dy(int i)const {return dsh_dy_[i];}
    auto dsh_dy()const {return dsh_dy_;}
    auto dsh_dz(int i)const {return dsh_dz_[i];}
    auto dsh_dz()const {return dsh_dz_;}
    auto J()const {return J_;}
    auto J_bd()const {return J_bd_;}
    auto W_in()const {return W_in_;}
    auto W_bd()const {return W_bd_;}
    auto JxW()const {return J_*W_in_;}
    auto N(int i)const {return N_[i];}
    auto N()const {return N_;}
    auto JN(int i)const {return JN_[i];}
    auto JN()const {return JN_;}
    auto JNxW(int i)const {return JN_[i]*W_bd_;}
    auto JNxW()const {return JN_*W_bd_;}
    auto g()const {return g_;}
    auto G()const {return G_;}
    auto g_up()const {return g_up_;}









  private:
  
    int nb_el_nodes_;                                     /// Number of nodes of element (3 for triangle, 4 for tetrahedral, etc.)
    vec sh_;                                              /// Vector of shape functions
    vec dsh_dxi_, dsh_deta_, dsh_dzeta_;                  /// Vectors of derivatives of shape functions w.r.t. isoparametric coordinate system (xi, eta, zeta)
    vec dsh_dx_, dsh_dy_, dsh_dz_;                        /// Vectors of derivatives of shape functions w.r.t. global coordinate system (x, y, z)
    double J_, J_bd_;                                     /// Jacobian of the transformation from (xi, eta, zeta) ----> (x, y, z)
    vec N_;                                               /// Outward normal vector on the face of the element 
    vec JN_;                                              /// Outward area normal vector on the face of the element
    mat g_up_;                                            /// element metric tensor
    vec g_;                                               /// element metric tensor
    mat G_;                                               /// element metric tensor
    double xi_, eta_, zeta_;                              /// Isoparametric coordinate system axes
    double x_, y_, z_;                                    /// Point coordinates corresponding to gauss point (xi, eta, zeta)
    double dx_dxi_, dx_deta_, dx_dzeta_;                  /// Partial derivate of   x   w.r.t.   (xi, eta, zeta)
    double dy_dxi_, dy_deta_, dy_dzeta_;                  /// Partial derivate of   y   w.r.t.   (xi, eta, zeta)
    double dz_dxi_, dz_deta_, dz_dzeta_;                  /// Partial derivate of   z   w.r.t.   (xi, eta, zeta)
    double dxi_dx_, deta_dx_, dzeta_dx_;                  /// Partial derivate of   xi   w.r.t.   (x, y, z)
    double dxi_dy_, deta_dy_, dzeta_dy_;                  /// Partial derivate of   eta   w.r.t.   (x, y, z)
    double dxi_dz_, deta_dz_, dzeta_dz_;                  /// Partial derivate of   zeta   w.r.t.   (x, y, z)    
    double W_in_;                                         /// Weight of gauss points for integration of element interior 
    double W_bd_;                                         /// Weight of gauss points for integration of element face
  };
   

} 

#endif

