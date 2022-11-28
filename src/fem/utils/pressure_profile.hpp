#ifndef PP_HPP
#define PP_HPP



#include <vector>
#include "boost/numeric/ublas/vector.hpp"



namespace GALES{


  /**
       This file computes hydrostatic pressure profile   
  */ 



  typedef boost::numeric::ublas::vector<double> v_type;




  /// This function interpolates the value of pressure from p_vec depending on y position in y_vec
  double p_interpolation(const v_type& y_vec, const v_type& p_vec, double y)
  {   
    std::size_t i;
    if( *(y_vec.begin()) < *(y_vec.end()-1))   i = std::lower_bound(y_vec.begin(), y_vec.end(), y) - y_vec.begin();
    else  i = std::lower_bound(y_vec.begin(), y_vec.end(), y, std::greater<double>()) - y_vec.begin();
         
    if(i>0)  i--;
    if(i >= y_vec.size()-1)
    {
         Error("index i goes out of the vector scope");
    }
    
    const double f = (y-y_vec[i])/(y_vec[i+1]-y_vec[i]);
    const double p = f*p_vec[i+1]+(1.0-f)*p_vec[i];
    return p;
  }





  /// This class computes hydrostatic pressure profile for monocoponent fluid
  template<typename ic_bc_type, int dim>
  class pressure_profile_sc
  {
      using model_type = model<dim>;
      using nd_type = node<dim>;

      public:
   
      pressure_profile_sc(ic_bc_type& ic_bc, model_type& f_model, fluid_properties& props): 
      ic_bc_(ic_bc),
      state_(f_model.state()),  
      props_(props),
      mesh_(f_model.mesh()),
      total_steps_(ic_bc.total_steps_),
      h_(ic_bc.h_),
      max_(ic_bc.max_)
      {                         
         if(ic_bc_.automatic_correction_)  return;                                          /// automatic correction means doing nothing
         v_type y_vec, p_vec;         
         if(ic_bc_.central_pressure_profile_) correct(y_vec, p_vec, dummy_);                ///This is when pressure is same at a horizontal level but very only vertically
                  
         for(const auto& nd : mesh_.nodes())
         {
             if(ic_bc_.x_dependent_pressure_profile_)                                      ///This is when pressure varies at a horizontal and vertical levels
             {
                 const double x(nd->get_x());
                 correct(x, y_vec, p_vec, dummy_);
             }

            const double y(nd->get_y());
            const double p_interpolated = p_interpolation(y_vec, p_vec, y);
            state_.set_dof(nd->first_dof_lid(), p_interpolated);
         }         
      }

      
      private:

      void correct(v_type& y_vec, v_type& p_vec, const point_2d& unused)
      {
        p_vec.resize(total_steps_+1), y_vec.resize(total_steps_+1);        
        v_type T(total_steps_+1,0.0), T1(total_steps_+1,0.0);
  
        for(int i=0; i<=total_steps_; i++)
        {
          y_vec[i]= max_ - h_*i;
          point_2d b(0.0,y_vec[i]);
          nd_type nd1(b);
          T[i] = ic_bc_.initial_T(nd1);

          auto y1 = y_vec[i] - h_*0.5;
          point_2d c(0.0,y1);
          nd_type nd2(c);          
          T1[i] = ic_bc_.initial_T(nd2);
        }
        Runge_Kutta(T,T1,p_vec);
      }
  
  
      void correct(v_type& y_vec, v_type& p_vec, const point_3d& unused)
      {
        p_vec.resize(total_steps_+1), y_vec.resize(total_steps_+1);
        v_type T(total_steps_+1,0.0), T1(total_steps_+1,0.0);
  
        for(int i=0; i<=total_steps_; i++)
        {
          y_vec[i]= max_ - h_*i;
          point_3d b(0.0,y_vec[i],0.0);
          nd_type nd1(b);
          T[i] = ic_bc_.initial_T(nd1);

          auto y1 = y_vec[i] - h_*0.5;
          point_3d c(0.0,y1,0.0);
          nd_type nd2(c);
          T1[i] = ic_bc_.initial_T(nd2);
        }
        Runge_Kutta(T,T1,p_vec);
      }
  
  
      void correct(double x, v_type& y_vec, v_type& p_vec, const point_2d& unused)
      {
        p_vec.resize(total_steps_+1), y_vec.resize(total_steps_+1);
        v_type T(total_steps_+1,0.0), T1(total_steps_+1,0.0);
  
        for(int i=0; i<=total_steps_; i++)
        {
          y_vec[i]= max_ - h_*i;
          point_2d b(x,y_vec[i]);
          nd_type nd1(b);
          T[i] = ic_bc_.initial_T(nd1);

          auto y1 = y_vec[i] - h_*0.5;
          point_2d c(x,y1);
          nd_type nd2(c);
          T1[i] = ic_bc_.initial_T(nd2);
        }
        Runge_Kutta(T,T1,p_vec);
      }
  
  
      void correct(double x, v_type& y_vec, v_type& p_vec, const point_3d& unused)
      {
        p_vec.resize(total_steps_+1), y_vec.resize(total_steps_+1);
        v_type T(total_steps_+1,0.0), T1(total_steps_+1,0.0);
  
        for(int i=0; i<=total_steps_; i++)
        {
          y_vec[i]= max_ - h_*i;
          point_3d b(x,y_vec[i],0.0);
          nd_type nd1(b);
          T[i] = ic_bc_.initial_T(nd1);

          auto y1 = y_vec[i] - h_*0.5;
          point_3d c(x,y1,0.0);
          nd_type nd2(c);
          T1[i] = ic_bc_.initial_T(nd2);
        }
        Runge_Kutta(T,T1,p_vec);
      }
  
  
      void Runge_Kutta(const v_type& T, const v_type& T1, v_type& p_vec)
      {
        double m1,m2,m3,m4;
        p_vec[0] = ic_bc_.p0_;
        for(int i=0; i<total_steps_; i++)
        {
          props_.properties(p_vec[i],T[i],0.0);
          m1 = 9.81*props_.rho_;
          props_.properties(p_vec[i]+0.5*h_*m1,T1[i],0.0);
          m2 = 9.81*props_.rho_;
          props_.properties(p_vec[i]+0.5*h_*m2,T1[i],0.0);
          m3 = 9.81*props_.rho_;
          props_.properties(p_vec[i]+h_*m3,T[i+1],0.0);
          m4 = 9.81*props_.rho_;
          p_vec[i+1]=p_vec[i]+h_/6.*(m1+2.0*(m2+m3)+m4);
        }
      }
            
      ic_bc_type& ic_bc_;
      dof_state& state_;
      fluid_properties& props_;
      Mesh<dim>& mesh_;
      int total_steps_;
      double h_;
      double max_;
      point<dim> dummy_;               
   };






















  /// This class computes hydrostatic pressure profile for monocoponent isothermal fluid
  template<typename ic_bc_type, int dim>
  class pressure_profile_sc_isothermal
  {
      using model_type = model<dim>;

      public:
   
      pressure_profile_sc_isothermal(ic_bc_type& ic_bc, model_type& f_model, fluid_properties& props): 
      ic_bc_(ic_bc),
      state_(f_model.state()),  
      props_(props),
      mesh_(f_model.mesh()),
      total_steps_(ic_bc.total_steps_),
      h_(ic_bc.h_),
      max_(ic_bc.max_)
      {                         
         if(ic_bc_.automatic_correction_)  return;                                          /// automatic correction means doing nothing
         v_type y_vec, p_vec;         
         if(ic_bc_.central_pressure_profile_) correct(y_vec, p_vec, dummy_);                ///This is when pressure is same at a horizontal level but very only vertically
                  
         for(const auto& nd : mesh_.nodes())
         {
             if(ic_bc_.x_dependent_pressure_profile_)                                      ///This is when pressure varies at a horizontal and vertical levels
             {
                 const double x(nd->get_x());
                 correct(x, y_vec, p_vec, dummy_);
             }

            const double y(nd->get_y());
            const double p_interpolated = p_interpolation(y_vec, p_vec, y);
            state_.set_dof(nd->first_dof_lid(), p_interpolated);
         }         
      }

      
      private:

      void correct(v_type& y_vec, v_type& p_vec, const point<dim>& unused)
      {
        p_vec.resize(total_steps_+1); 
        y_vec.resize(total_steps_+1);        
  
        for(int i=0; i<=total_steps_; i++)
          y_vec[i]= max_ - h_*i;

        Runge_Kutta(p_vec);
      }
  
  
      void correct(double x, v_type& y_vec, v_type& p_vec, const point<dim>& unused)
      {
        p_vec.resize(total_steps_+1); 
        y_vec.resize(total_steps_+1);        
  
        for(int i=0; i<=total_steps_; i++)
          y_vec[i]= max_ - h_*i;

        Runge_Kutta(p_vec);
      }


      void Runge_Kutta(v_type& p_vec)
      {
        double m1,m2,m3,m4;
        p_vec[0] = ic_bc_.p0_;
        for(int i=0; i<total_steps_; i++)
        {
          props_.properties(p_vec[i],0.0);
          m1 = 9.81*props_.rho_;
          props_.properties(p_vec[i]+0.5*h_*m1,0.0);
          m2 = 9.81*props_.rho_;
          props_.properties(p_vec[i]+0.5*h_*m2,0.0);
          m3 = 9.81*props_.rho_;
          props_.properties(p_vec[i]+h_*m3,0.0);
          m4 = 9.81*props_.rho_;
          p_vec[i+1]=p_vec[i]+h_/6.*(m1+2.0*(m2+m3)+m4);
        }
      }
            
      ic_bc_type& ic_bc_;
      dof_state& state_;
      fluid_properties& props_;
      Mesh<dim>& mesh_;
      int total_steps_;
      double h_;
      double max_;
      point<dim> dummy_;               
   };


















  /// This class computes hydrostatic pressure profile for multicoponent fluid
  template<typename ic_bc_type, int dim>
  class pressure_profile_mc
  {
      using model_type = model<dim>;
      using nd_type = node<dim>;

      public:
   
      pressure_profile_mc(ic_bc_type& ic_bc, model_type& f_model, fluid_properties& props): 
      ic_bc_(ic_bc),
      state_(f_model.state()),  
      props_(props),
      mesh_(f_model.mesh()),
      total_steps_(ic_bc.total_steps_),
      h_(ic_bc.h_),
      max_(ic_bc.max_)
      {              
         if(ic_bc_.automatic_correction_)  return;                    /// automatic correction means doing nothing

         v_type y_vec, p_vec;         
         if(ic_bc_.central_pressure_profile_)   correct(y_vec, p_vec, dummy_);    ///This is when pressure is same at a horizontal level but very only vertically
                  
         for(const auto& nd : mesh_.nodes())
         {
             if(ic_bc_.x_dependent_pressure_profile_)                    ///This is when pressure varies at a horizontal and vertical levels
             {
                 const double x(nd->get_x());
                 correct(x, y_vec, p_vec, dummy_);
             }

            const double y(nd->get_y());
            const double p_interpolated = p_interpolation(y_vec, p_vec, y);
            state_.set_dof(nd->first_dof_lid(), p_interpolated);
         }
      }

      
      private:

      void correct(v_type& y_vec, v_type& p_vec, const point_2d& unused)
      {
        p_vec.resize(total_steps_+1), y_vec.resize(total_steps_+1);        
        v_type T(total_steps_+1,0.0), T1(total_steps_+1,0.0);
        std::vector<v_type> A(total_steps_+1), A1(total_steps_+1);
  
        for(int i=0; i<=total_steps_; i++)
        {
          y_vec[i]= max_ - h_*i;
          point_2d b(0.0, y_vec[i]);
          nd_type nd1(b);
          T[i] = ic_bc_.initial_T(nd1);
          A[i].resize(props_.nb_comp_);
          ic_bc_.initial_Y(nd1, A[i]);          
          A[i][props_.nb_comp_-1] = 1.0- std::accumulate(&A[i][0], &A[i][props_.nb_comp_-1], 0.0);

          auto y1 = y_vec[i] - h_*0.5;
          point_2d c(0.0, y1);
          nd_type nd2(c);
          T1[i] = ic_bc_.initial_T(nd2);
          A1[i].resize(props_.nb_comp_);
          ic_bc_.initial_Y(nd2, A1[i]);
          A1[i][props_.nb_comp_-1] = 1.0- std::accumulate(&A1[i][0], &A1[i][props_.nb_comp_-1], 0.0);
        }        
        
        
        Runge_Kutta(T, T1, A, A1, p_vec);
      }
  
  
      void correct(v_type& y_vec, v_type& p_vec, const point_3d& unused)
      {
        p_vec.resize(total_steps_+1), y_vec.resize(total_steps_+1);
        v_type T(total_steps_+1,0.0), T1(total_steps_+1,0.0);
        std::vector<v_type> A(total_steps_+1), A1(total_steps_+1);
  
        for(int i=0; i<=total_steps_; i++)
        {
          y_vec[i]= max_ - h_*i;
          point_3d b(0.0,y_vec[i],0.0);
          nd_type nd1(b);
          T[i] = ic_bc_.initial_T(nd1);
          A[i].resize(props_.nb_comp_);
          ic_bc_.initial_Y(nd1,A[i]);
          A[i][props_.nb_comp_-1] = 1.0- std::accumulate(&A[i][0],&A[i][props_.nb_comp_-1],0.0);

          auto y1 = y_vec[i] - h_*0.5;
          point_3d c(0.0,y1,0.0);
          nd_type nd2(c);
          T1[i] = ic_bc_.initial_T(nd2);
          A1[i].resize(props_.nb_comp_);
          ic_bc_.initial_Y(nd2,A1[i]);
          A1[i][props_.nb_comp_-1] = 1.0- std::accumulate(&A1[i][0],&A1[i][props_.nb_comp_-1],0.0);
        }
        Runge_Kutta(T,T1,A,A1,p_vec);
      }
  
  
      void correct(double x, v_type& y_vec, v_type& p_vec, const point_2d& unused)
      {
        p_vec.resize(total_steps_+1), y_vec.resize(total_steps_+1);
        v_type T(total_steps_+1,0.0), T1(total_steps_+1,0.0);
        std::vector<v_type> A(total_steps_+1), A1(total_steps_+1);
  
        for(int i=0; i<=total_steps_; i++)
        {
          y_vec[i]= max_ - h_*i;
          point_2d b(x,y_vec[i]);
          nd_type nd1(b);
          T[i] = ic_bc_.initial_T(nd1);
          A[i].resize(props_.nb_comp_);
          ic_bc_.initial_Y(nd1,A[i]);
          A[i][props_.nb_comp_-1] = 1.0- std::accumulate(&A[i][0],&A[i][props_.nb_comp_-1],0.0);

          auto y1 = y_vec[i] - h_*0.5;
          point_2d c(x,y1);
          nd_type nd2(c);
          T1[i] = ic_bc_.initial_T(nd2);
          A1[i].resize(props_.nb_comp_);
          ic_bc_.initial_Y(nd2,A1[i]);
          A1[i][props_.nb_comp_-1] = 1.0- std::accumulate(&A1[i][0],&A1[i][props_.nb_comp_-1],0.0);
        }
        Runge_Kutta(T,T1,A,A1,p_vec);
      }
  
  
      void correct(double x, v_type& y_vec, v_type& p_vec, const point_3d& unused)
      {
        p_vec.resize(total_steps_+1), y_vec.resize(total_steps_+1);
        v_type T(total_steps_+1,0.0), T1(total_steps_+1,0.0);
        std::vector<v_type> A(total_steps_+1), A1(total_steps_+1);
  
        for(int i=0; i<=total_steps_; i++)
        {
          y_vec[i]= max_ - h_*i;
          point_3d b(x,y_vec[i],0.0);
          nd_type nd1(b);
          T[i] = ic_bc_.initial_T(nd1);
          A[i].resize(props_.nb_comp_);
          ic_bc_.initial_Y(nd1,A[i]);
          A[i][props_.nb_comp_-1] = 1.- std::accumulate(&A[i][0],&A[i][props_.nb_comp_-1],0.0);

          auto y1 = y_vec[i] - h_*0.5;
          point_3d c(x,y1,0.0);
          nd_type nd2(c);
          T1[i] = ic_bc_.initial_T(nd2);
          A1[i].resize(props_.nb_comp_);
          ic_bc_.initial_Y(nd2,A1[i]);
          A1[i][props_.nb_comp_-1] = 1.- std::accumulate(&A1[i][0],&A1[i][props_.nb_comp_-1],0.0);
        }
        Runge_Kutta(T,T1,A,A1,p_vec);
      }
  
  
      void Runge_Kutta(const v_type& T, const v_type& T1, const std::vector<v_type>& A, const std::vector<v_type>& A1, v_type& p_vec)
      {
        double m1,m2,m3,m4;
        p_vec[0] = ic_bc_.p0_;
 
        for(int i=0; i<total_steps_; i++)
        {
          props_.properties(p_vec[i],T[i],A[i],0.0);
          m1 = 9.81*props_.rho_;
          props_.properties(p_vec[i]+0.5*h_*m1,T1[i],A1[i],0.0);
          m2 = 9.81*props_.rho_;
          props_.properties(p_vec[i]+0.5*h_*m2,T1[i],A1[i],0.0);
          m3 = 9.81*props_.rho_;          
          props_.properties(p_vec[i]+h_*m3,T[i+1],A[i+1],0.0);
          m4 = 9.81*props_.rho_;          
          p_vec[i+1] = p_vec[i]+h_/6.*(m1+2.0*(m2+m3)+m4);
        }
      }
  
      ic_bc_type& ic_bc_;
      dof_state& state_;
      fluid_properties& props_;
      Mesh<dim>& mesh_;
      int total_steps_;
      double h_;
      double max_;
      point<dim> dummy_;               
              
   };









  /// This class computes hydrostatic pressure profile for multicoponent isothermal fluid
  template<typename ic_bc_type, int dim>
  class pressure_profile_mc_isothermal
  {
      using model_type = model<dim>;
      using nd_type = node<dim>;

      public:
   
      pressure_profile_mc_isothermal(ic_bc_type& ic_bc, model_type& f_model, fluid_properties& props): 
      ic_bc_(ic_bc),
      state_(f_model.state()),  
      props_(props),
      mesh_(f_model.mesh()),
      total_steps_(ic_bc.total_steps_),
      h_(ic_bc.h_),
      max_(ic_bc.max_)
      {              
         if(ic_bc_.automatic_correction_)  return;                    /// automatic correction means doing nothing

         v_type y_vec, p_vec;         
         if(ic_bc_.central_pressure_profile_)   correct(y_vec, p_vec, dummy_);    ///This is when pressure is same at a horizontal level but very only vertically
                  
         for(const auto& nd : mesh_.nodes())
         {
             if(ic_bc_.x_dependent_pressure_profile_)                    ///This is when pressure varies at a horizontal and vertical levels
             {
                 const double x(nd->get_x());
                 correct(x, y_vec, p_vec, dummy_);
             }

            const double y(nd->get_y());
            const double p_interpolated = p_interpolation(y_vec, p_vec, y);
            state_.set_dof(nd->first_dof_lid(), p_interpolated);
         }
      }

      
      private:

      void correct(v_type& y_vec, v_type& p_vec, const point_2d& unused)
      {
        p_vec.resize(total_steps_+1), y_vec.resize(total_steps_+1);        
        std::vector<v_type> A(total_steps_+1), A1(total_steps_+1);
  
        for(int i=0; i<=total_steps_; i++)
        {
          y_vec[i]= max_ - h_*i;
          point_2d b(0.0, y_vec[i]);
          nd_type nd1(b);
          A[i].resize(props_.nb_comp_);
          ic_bc_.initial_Y(nd1, A[i]);          
          A[i][props_.nb_comp_-1] = 1.0- std::accumulate(&A[i][0], &A[i][props_.nb_comp_-1], 0.0);

          auto y1 = y_vec[i] - h_*0.5;
          point_2d c(0.0, y1);
          nd_type nd2(c);
          A1[i].resize(props_.nb_comp_);
          ic_bc_.initial_Y(nd2, A1[i]);
          A1[i][props_.nb_comp_-1] = 1.0- std::accumulate(&A1[i][0], &A1[i][props_.nb_comp_-1], 0.0);
        }        
        
        
        Runge_Kutta(A, A1, p_vec);
      }
  
  
      void correct(v_type& y_vec, v_type& p_vec, const point_3d& unused)
      {
        p_vec.resize(total_steps_+1), y_vec.resize(total_steps_+1);
        std::vector<v_type> A(total_steps_+1), A1(total_steps_+1);
  
        for(int i=0; i<=total_steps_; i++)
        {
          y_vec[i]= max_ - h_*i;
          point_3d b(0.0,y_vec[i],0.0);
          nd_type nd1(b);
          A[i].resize(props_.nb_comp_);
          ic_bc_.initial_Y(nd1,A[i]);
          A[i][props_.nb_comp_-1] = 1.0- std::accumulate(&A[i][0],&A[i][props_.nb_comp_-1],0.0);

          auto y1 = y_vec[i] - h_*0.5;
          point_3d c(0.0,y1,0.0);
          nd_type nd2(c);
          A1[i].resize(props_.nb_comp_);
          ic_bc_.initial_Y(nd2,A1[i]);
          A1[i][props_.nb_comp_-1] = 1.0- std::accumulate(&A1[i][0],&A1[i][props_.nb_comp_-1],0.0);
        }
        Runge_Kutta(A,A1,p_vec);
      }
  
  
      void correct(double x, v_type& y_vec, v_type& p_vec, const point_2d& unused)
      {
        p_vec.resize(total_steps_+1), y_vec.resize(total_steps_+1);
        std::vector<v_type> A(total_steps_+1), A1(total_steps_+1);
  
        for(int i=0; i<=total_steps_; i++)
        {
          y_vec[i]= max_ - h_*i;
          point_2d b(x,y_vec[i]);
          nd_type nd1(b);
          A[i].resize(props_.nb_comp_);
          ic_bc_.initial_Y(nd1,A[i]);
          A[i][props_.nb_comp_-1] = 1.0- std::accumulate(&A[i][0],&A[i][props_.nb_comp_-1],0.0);

          auto y1 = y_vec[i] - h_*0.5;
          point_2d c(x,y1);
          nd_type nd2(c);
          A1[i].resize(props_.nb_comp_);
          ic_bc_.initial_Y(nd2,A1[i]);
          A1[i][props_.nb_comp_-1] = 1.0- std::accumulate(&A1[i][0],&A1[i][props_.nb_comp_-1],0.0);
        }
        Runge_Kutta(A,A1,p_vec);
      }
  
  
      void correct(double x, v_type& y_vec, v_type& p_vec, const point_3d& unused)
      {
        p_vec.resize(total_steps_+1), y_vec.resize(total_steps_+1);
        std::vector<v_type> A(total_steps_+1), A1(total_steps_+1);
  
        for(int i=0; i<=total_steps_; i++)
        {
          y_vec[i]= max_ - h_*i;
          point_3d b(x,y_vec[i],0.0);
          nd_type nd1(b);
          A[i].resize(props_.nb_comp_);
          ic_bc_.initial_Y(nd1,A[i]);
          A[i][props_.nb_comp_-1] = 1.- std::accumulate(&A[i][0],&A[i][props_.nb_comp_-1],0.0);

          auto y1 = y_vec[i] - h_*0.5;
          point_3d c(x,y1,0.0);
          nd_type nd2(c);
          A1[i].resize(props_.nb_comp_);
          ic_bc_.initial_Y(nd2,A1[i]);
          A1[i][props_.nb_comp_-1] = 1.- std::accumulate(&A1[i][0],&A1[i][props_.nb_comp_-1],0.0);
        }
        Runge_Kutta(A,A1,p_vec);
      }
  
  
      void Runge_Kutta(const std::vector<v_type>& A, const std::vector<v_type>& A1, v_type& p_vec)
      {
        double m1,m2,m3,m4;
        p_vec[0] = ic_bc_.p0_;
 
        for(int i=0; i<total_steps_; i++)
        {
          props_.properties(p_vec[i],A[i],0.0);
          m1 = 9.81*props_.rho_;
          props_.properties(p_vec[i]+0.5*h_*m1,A1[i],0.0);
          m2 = 9.81*props_.rho_;
          props_.properties(p_vec[i]+0.5*h_*m2,A1[i],0.0);
          m3 = 9.81*props_.rho_;          
          props_.properties(p_vec[i]+h_*m3,A[i+1],0.0);
          m4 = 9.81*props_.rho_;          
          p_vec[i+1] = p_vec[i]+h_/6.*(m1+2.0*(m2+m3)+m4);
        }
      }
  
      ic_bc_type& ic_bc_;
      dof_state& state_;
      fluid_properties& props_;
      Mesh<dim>& mesh_;
      int total_steps_;
      double h_;
      double max_;
      point<dim> dummy_;               
              
   };





    
}    
    
#endif    

