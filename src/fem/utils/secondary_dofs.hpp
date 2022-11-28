#ifndef SECONDARY_DOFS_HPP
#define SECONDARY_DOFS_HPP




#include "boost/numeric/ublas/vector.hpp"
#include "boost/tuple/tuple.hpp"
#include <vector>
#include "../base/read_setup.hpp"
#include "../base/io.hpp"



namespace GALES{

  
  
  /**
      This class is to compute density volume fraction etc for magma simulations as postprocessing in pp_2d.hpp. 
      The following must be set in setup.txt file to compute the required properties.
      
      pp_start_time        0
      pp_final_time        10000
      pp_delta_t           500

      pp_rho               T
      pp_mu                T
      pp_cp                T
      pp_cv                T
      pp_alpha             T
      pp_beta              T
      pp_sound_speed       T
      pp_kappa             T
      pp_molar_mass        T
      pp_vf_g              T      
  */
  
  
  template<int dim>
  class secondary_dofs 
  {
     using vec = std::vector<double>;

    public:

     secondary_dofs(Epetra_Map& map, fluid_properties& props, read_setup& setup)
     : 
     setup_(setup), nb_comp_(props.nb_comp()), sec_io_(map)
     {
        create_dirs();
     }



     void execute(model<dim>& f_model, fluid_properties& props) 
     {      
        const auto& mesh(f_model.mesh());
        const int l(props.l_), g(props.g_);
        
        
        /// vectors of secondary variables
        vec rho, mu, cp, cv, alpha, beta, kappa, sound_speed;
        vec vf_g_on_total, vf_h2o_g_on_total, vf_co2_g_on_total; //, vf_h2o_g_on_g_in_total, vf_co2_g_on_g_in_total;
        vec wf_g_on_total, wf_h2o_g_on_total, wf_co2_g_on_total; //, wf_h2o_g_on_g_in_total, wf_co2_g_on_g_in_total;
        std::vector<vec> vf_g_comp(nb_comp_, std::vector<double>());                
        std::vector<vec> vf_comp_vec(nb_comp_, std::vector<double>());                
        
        
        /// loop over mesh nodes to compute secondary variables and fill in vectors        
        for(const auto& nd : mesh.nodes()) 
        {
           const int first_dof_lid (nd->first_dof_lid());      
           boost::numeric::ublas::vector<double> Y(nb_comp_, 0.0);	

           const double p(f_model.state().get_dof(first_dof_lid));           
           const double T(f_model.state().get_dof(first_dof_lid+dim+1));          
           double sum(0.0);
           for(int i=0; i<nb_comp_-1; i++)
           {
                  Y[i] = f_model.state().get_dof(first_dof_lid+dim+2+i);
                  sum += Y[i];
           }            
           Y[nb_comp_-1] = 1.0 - sum;   
           props.properties(p, T, Y, 0.0);

           
           ///start filling in vectors 
           if(setup_.pp_rho())           rho.push_back(props.rho_);
           if(setup_.pp_mu())            mu.push_back(props.mu_);
           if(setup_.pp_cp())            cp.push_back(props.cp_);
           if(setup_.pp_cv())            cv.push_back(props.cv_);
           if(setup_.pp_alpha())         alpha.push_back(props.alpha_);
           if(setup_.pp_beta())          beta.push_back(props.beta_);
           if(setup_.pp_kappa())         kappa.push_back(props.kappa_);
           if(setup_.pp_sound_speed())   sound_speed.push_back(props.sound_speed_);


           ///-------------------vf wf computation----------------------------
           if(setup_.pp_vf_g())         
           {
             double vf_g_on_tot(0.0);              //   vol(g)/vol(total magma)
             double vf_h2o_g_on_tot(0.0);          //   vol(h2o_g)/vol(total magma)
             double vf_co2_g_on_tot(0.0);          //   vol(co2_g)/vol(total magma)
//             double vf_h2o_g_on_g_in_tot(0.0);     //   vol(h2o_g)/vol(g) in total magma
//             double vf_co2_g_on_g_in_tot(0.0);     //   vol(co2_g)/vol(g) in total magma
             double wf_g_on_tot(0.0);              //   mass(g)/mass(total magma)
             double wf_h2o_g_on_tot(0.0);          //   mass(h2o_g)/mass(total magma)
             double wf_co2_g_on_tot(0.0);          //   mass(co2_g)/mass(total magma)
//             double wf_h2o_g_on_g_in_tot(0.0);     //   mass(h2o_g)/mass(g) in total magma
//             double wf_co2_g_on_g_in_tot(0.0);     //   mass(co2_g)/mass(g) in total magma
 
             for (int c=0; c<nb_comp_; c++)
             {
               double vf_comp = Y[c]*props.rho_/props.mixture_ptrs_[c]->rho_;                  //   vol(comp)/vol(mix)
               vf_comp_vec[c].push_back(vf_comp);
               vf_g_comp[c].push_back(props.mixture_ptrs_[c]->vf_ph_[g]*vf_comp);          //   vol(g)/vol(comp)
              
               vf_g_on_tot += (props.mixture_ptrs_[c]->vf_ph_[g])*vf_comp;            
               vf_h2o_g_on_tot += (props.mixture_ptrs_[c]->vf_[3])*vf_comp;
               vf_co2_g_on_tot += (props.mixture_ptrs_[c]->vf_[4])*vf_comp;            
//               vf_h2o_g_on_g_in_tot += vf_h2o_g_on_tot/vf_g_on_tot;
//               vf_co2_g_on_g_in_tot += vf_co2_g_on_tot/vf_g_on_tot;
 
               wf_g_on_tot += (props.mixture_ptrs_[c]->wf_ph_[g])*Y[c];
               wf_h2o_g_on_tot += (props.mixture_ptrs_[c]->wf_[3])*Y[c];
               wf_co2_g_on_tot += (props.mixture_ptrs_[c]->wf_[4])*Y[c];            
//               wf_h2o_g_on_g_in_tot += wf_h2o_g_on_tot/wf_g_on_tot;
//               wf_co2_g_on_g_in_tot += wf_co2_g_on_tot/wf_g_on_tot;
             }
          
             
             vf_g_on_total.push_back(vf_g_on_tot);             
             vf_h2o_g_on_total.push_back(vf_h2o_g_on_tot);             
             vf_co2_g_on_total.push_back(vf_co2_g_on_tot);             
//             vf_h2o_g_on_g_in_total.push_back(vf_h2o_g_on_g_in_tot);             
//             vf_co2_g_on_g_in_total.push_back(vf_co2_g_on_g_in_tot);             
             wf_g_on_total.push_back(wf_g_on_tot);             
             wf_h2o_g_on_total.push_back(wf_h2o_g_on_tot);             
             wf_co2_g_on_total.push_back(wf_co2_g_on_tot);             
//             wf_h2o_g_on_g_in_total.push_back(wf_h2o_g_on_g_in_tot);             
//             wf_co2_g_on_g_in_total.push_back(wf_co2_g_on_g_in_tot);
           }                                            
        }  


           ///-----------------print all vecs-------------------------------------  
           if(setup_.pp_rho())              print_sec_dofs("rho", rho);        
           if(setup_.pp_mu())               print_sec_dofs("mu", mu);        
           if(setup_.pp_cp())               print_sec_dofs("cp", cp);        
           if(setup_.pp_cv())               print_sec_dofs("cv", cv);        
           if(setup_.pp_alpha())            print_sec_dofs("alpha", alpha);
           if(setup_.pp_beta())             print_sec_dofs("beta", beta);        
           if(setup_.pp_kappa())            print_sec_dofs("kappa", kappa);
           if(setup_.pp_sound_speed())      print_sec_dofs("sound_speed", sound_speed);
           if(setup_.pp_vf_g())
           {
             for (int c=0; c<nb_comp_; c++)
             {
               const std::string s = "vf_g_comp_" + std::to_string(c);
               print_sec_dofs(s, vf_g_comp[c]);  
               
               const std::string s1 = "vf_comp_" + std::to_string(c);               
               print_sec_dofs(s1, vf_comp_vec[c]);                 
             }
             
             print_sec_dofs("vf_g_on_total", vf_g_on_total);
             print_sec_dofs("vf_h2o_g_on_total", vf_h2o_g_on_total);
      	     print_sec_dofs("vf_co2_g_on_total", vf_co2_g_on_total);
//    	     print_sec_dofs("vf_h2o_g_on_g_in_total", vf_h2o_g_on_g_in_total);
//    	     print_sec_dofs("vf_co2_g_on_g_in_total", vf_co2_g_on_g_in_total);
    	     print_sec_dofs("wf_g_on_total", wf_g_on_total);
    	     print_sec_dofs("wf_h2o_g_on_total", wf_h2o_g_on_total);
    	     print_sec_dofs("wf_co2_g_on_total", wf_co2_g_on_total);
//    	     print_sec_dofs("wf_h2o_g_on_g_in_total", wf_h2o_g_on_g_in_total);
//    	     print_sec_dofs("wf_co2_g_on_g_in_total", wf_co2_g_on_g_in_total);
           }
              
     }
     




     /// This function writes the vector "v" in parallel in directory "results/sec_dofs/{s}_{time}"
     void print_sec_dofs(std::string s, const vec& v)
     {
       std::stringstream ss;
       ss<<"results/sec_dofs/"<<s<<"/"<<time::get().t();
       sec_io_.write2(ss.str(), v);  	
     }
  




     /// This function geerate directories with names of properties in "results/sec_dofs/"
     void create_dirs()
     {
//       auto rank = get_rank();
//       if(rank==0)
//       {
        std::vector<std::string> v;
        if(setup_.pp_rho()) v.push_back("rho");
        if(setup_.pp_mu()) v.push_back("mu");
        if(setup_.pp_cp()) v.push_back("cp");
        if(setup_.pp_cv()) v.push_back("cv");
        if(setup_.pp_alpha()) v.push_back("alpha");
        if(setup_.pp_beta()) v.push_back("beta");
        if(setup_.pp_sound_speed()) v.push_back("sound_speed");
        if(setup_.pp_kappa()) v.push_back("kappa");

	if(setup_.pp_vf_g()) 
        {
          for(int i=0; i<nb_comp_; i++)
          {
            std::string s = "vf_g_comp_" + std::to_string(i);
            v.push_back(s);

            std::string s1 = "vf_comp_" + std::to_string(i);
            v.push_back(s1);
          }  
          
          v.push_back("vf_g_on_total");
          v.push_back("vf_h2o_g_on_total");
          v.push_back("vf_co2_g_on_total");
//          v.push_back("vf_h2o_g_on_g_in_total");
//          v.push_back("vf_co2_g_on_g_in_total");
          v.push_back("wf_g_on_total");
          v.push_back("wf_h2o_g_on_total");
          v.push_back("wf_co2_g_on_total");
//          v.push_back("wf_h2o_g_on_g_in_total");
//          v.push_back("wf_co2_g_on_g_in_total");
        }
 
        for(auto a : v)
          make_dirs("results/sec_dofs/" + a);
//       } 
     }




     private:     
     read_setup& setup_;
     int nb_comp_;
     IO sec_io_;  

  }; 































  template<int dim>
  class secondary_dofs_isothermal 
  {
     using vec = std::vector<double>;

    public:

     secondary_dofs_isothermal(Epetra_Map& map, fluid_properties& props, read_setup& setup)
     : 
     setup_(setup), nb_comp_(props.nb_comp()), sec_io_(map)
     {
        create_dirs();
     }



     void execute(model<dim>& f_model, fluid_properties& props) 
     {      
        const auto& mesh(f_model.mesh());
        const int l(props.l_), g(props.g_);
        
        
        /// vectors of secondary variables
        vec rho, mu, beta, sound_speed;
        vec vf_g_on_total, vf_h2o_g_on_total, vf_co2_g_on_total; //, vf_h2o_g_on_g_in_total, vf_co2_g_on_g_in_total;
        vec wf_g_on_total, wf_h2o_g_on_total, wf_co2_g_on_total; //, wf_h2o_g_on_g_in_total, wf_co2_g_on_g_in_total;
        std::vector<vec> vf_g_comp(nb_comp_, std::vector<double>());                
        std::vector<vec> vf_comp_vec(nb_comp_, std::vector<double>());                
        
        
        /// loop over mesh nodes to compute secondary variables and fill in vectors        
        for(const auto& nd : mesh.nodes()) 
        {
           const int first_dof_lid (nd->first_dof_lid());      
           boost::numeric::ublas::vector<double> Y(nb_comp_, 0.0);	

           const double p(f_model.state().get_dof(first_dof_lid));           
           double sum(0.0);
           for(int i=0; i<nb_comp_-1; i++)
           {
                  Y[i] = f_model.state().get_dof(first_dof_lid+dim+1+i);
                  sum += Y[i];
           }            
           Y[nb_comp_-1] = 1.0 - sum;   
           props.properties(p, Y, 0.0);

           
           ///start filling in vectors 
           if(setup_.pp_rho())           rho.push_back(props.rho_);
           if(setup_.pp_mu())            mu.push_back(props.mu_);
           if(setup_.pp_beta())          beta.push_back(props.beta_);
           if(setup_.pp_sound_speed())   sound_speed.push_back(props.sound_speed_);


           ///-------------------vf wf computation----------------------------
           if(setup_.pp_vf_g())         
           {
             double vf_g_on_tot(0.0);              //   vol(g)/vol(total magma)
             double vf_h2o_g_on_tot(0.0);          //   vol(h2o_g)/vol(total magma)
             double vf_co2_g_on_tot(0.0);          //   vol(co2_g)/vol(total magma)
//             double vf_h2o_g_on_g_in_tot(0.0);     //   vol(h2o_g)/vol(g) in total magma
//             double vf_co2_g_on_g_in_tot(0.0);     //   vol(co2_g)/vol(g) in total magma
             double wf_g_on_tot(0.0);              //   mass(g)/mass(total magma)
             double wf_h2o_g_on_tot(0.0);          //   mass(h2o_g)/mass(total magma)
             double wf_co2_g_on_tot(0.0);          //   mass(co2_g)/mass(total magma)
//             double wf_h2o_g_on_g_in_tot(0.0);     //   mass(h2o_g)/mass(g) in total magma
//             double wf_co2_g_on_g_in_tot(0.0);     //   mass(co2_g)/mass(g) in total magma
 
             for (int c=0; c<nb_comp_; c++)
             {
               double vf_comp = Y[c]*props.rho_/props.mixture_ptrs_[c]->rho_;                  //   vol(comp)/vol(mix)
               vf_comp_vec[c].push_back(vf_comp);
               vf_g_comp[c].push_back(props.mixture_ptrs_[c]->vf_ph_[g]*vf_comp);          //   vol(g)/vol(comp)
              
               vf_g_on_tot += (props.mixture_ptrs_[c]->vf_ph_[g])*vf_comp;            
               vf_h2o_g_on_tot += (props.mixture_ptrs_[c]->vf_[3])*vf_comp;
               vf_co2_g_on_tot += (props.mixture_ptrs_[c]->vf_[4])*vf_comp;            
//               vf_h2o_g_on_g_in_tot += vf_h2o_g_on_tot/vf_g_on_tot;
//               vf_co2_g_on_g_in_tot += vf_co2_g_on_tot/vf_g_on_tot;
 
               wf_g_on_tot += (props.mixture_ptrs_[c]->wf_ph_[g])*Y[c];
               wf_h2o_g_on_tot += (props.mixture_ptrs_[c]->wf_[3])*Y[c];
               wf_co2_g_on_tot += (props.mixture_ptrs_[c]->wf_[4])*Y[c];            
//               wf_h2o_g_on_g_in_tot += wf_h2o_g_on_tot/wf_g_on_tot;
//               wf_co2_g_on_g_in_tot += wf_co2_g_on_tot/wf_g_on_tot;
             }
          
             
             vf_g_on_total.push_back(vf_g_on_tot);             
             vf_h2o_g_on_total.push_back(vf_h2o_g_on_tot);             
             vf_co2_g_on_total.push_back(vf_co2_g_on_tot);             
//             vf_h2o_g_on_g_in_total.push_back(vf_h2o_g_on_g_in_tot);             
//             vf_co2_g_on_g_in_total.push_back(vf_co2_g_on_g_in_tot);             
             wf_g_on_total.push_back(wf_g_on_tot);             
             wf_h2o_g_on_total.push_back(wf_h2o_g_on_tot);             
             wf_co2_g_on_total.push_back(wf_co2_g_on_tot);             
//             wf_h2o_g_on_g_in_total.push_back(wf_h2o_g_on_g_in_tot);             
//             wf_co2_g_on_g_in_total.push_back(wf_co2_g_on_g_in_tot);
           }                                            
        }  


           ///-----------------print all vecs-------------------------------------  
           if(setup_.pp_rho())              print_sec_dofs("rho", rho);        
           if(setup_.pp_mu())               print_sec_dofs("mu", mu);        
           if(setup_.pp_beta())             print_sec_dofs("beta", beta);        
           if(setup_.pp_sound_speed())      print_sec_dofs("sound_speed", sound_speed);
           if(setup_.pp_vf_g())
           {
             for (int c=0; c<nb_comp_; c++)
             {
               const std::string s = "vf_g_comp_" + std::to_string(c);
               print_sec_dofs(s, vf_g_comp[c]);  
               
               const std::string s1 = "vf_comp_" + std::to_string(c);               
               print_sec_dofs(s1, vf_comp_vec[c]);                 
             }
             
             print_sec_dofs("vf_g_on_total", vf_g_on_total);
             print_sec_dofs("vf_h2o_g_on_total", vf_h2o_g_on_total);
      	     print_sec_dofs("vf_co2_g_on_total", vf_co2_g_on_total);
//    	     print_sec_dofs("vf_h2o_g_on_g_in_total", vf_h2o_g_on_g_in_total);
//    	     print_sec_dofs("vf_co2_g_on_g_in_total", vf_co2_g_on_g_in_total);
    	     print_sec_dofs("wf_g_on_total", wf_g_on_total);
    	     print_sec_dofs("wf_h2o_g_on_total", wf_h2o_g_on_total);
    	     print_sec_dofs("wf_co2_g_on_total", wf_co2_g_on_total);
//    	     print_sec_dofs("wf_h2o_g_on_g_in_total", wf_h2o_g_on_g_in_total);
//    	     print_sec_dofs("wf_co2_g_on_g_in_total", wf_co2_g_on_g_in_total);
           }
              
     }
     




     /// This function writes the vector "v" in parallel in directory "results/sec_dofs/{s}_{time}"
     void print_sec_dofs(std::string s, const vec& v)
     {
       std::stringstream ss;
       ss<<"results/sec_dofs/"<<s<<"/"<<time::get().t();
       sec_io_.write2(ss.str(), v);  	
     }
  




     /// This function geerate directories with names of properties in "results/sec_dofs/"
     void create_dirs()
     {
//       auto rank = get_rank();
//       if(rank==0)
//       {
        std::vector<std::string> v;
        if(setup_.pp_rho()) v.push_back("rho");
        if(setup_.pp_mu()) v.push_back("mu");
        if(setup_.pp_beta()) v.push_back("beta");
        if(setup_.pp_sound_speed()) v.push_back("sound_speed");

	if(setup_.pp_vf_g()) 
        {
          for(int i=0; i<nb_comp_; i++)
          {
            std::string s = "vf_g_comp_" + std::to_string(i);
            v.push_back(s);

            std::string s1 = "vf_comp_" + std::to_string(i);
            v.push_back(s1);
          }  
          
          v.push_back("vf_g_on_total");
          v.push_back("vf_h2o_g_on_total");
          v.push_back("vf_co2_g_on_total");
//          v.push_back("vf_h2o_g_on_g_in_total");
//          v.push_back("vf_co2_g_on_g_in_total");
          v.push_back("wf_g_on_total");
          v.push_back("wf_h2o_g_on_total");
          v.push_back("wf_co2_g_on_total");
//          v.push_back("wf_h2o_g_on_g_in_total");
//          v.push_back("wf_co2_g_on_g_in_total");
        }
 
        for(auto a : v)
          make_dirs("results/sec_dofs/" + a);
//       } 
     }




     private:     
     read_setup& setup_;
     int nb_comp_;
     IO sec_io_;  

  }; 









}
#endif
