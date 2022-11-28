#ifndef READ_SETUP_HPP
#define READ_SETUP_HPP


#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>

#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/trim.hpp>



namespace GALES{



  class read_setup
  {      
   
   public:
   
   
    
   read_setup()                
   {
      read_data();
      set_vars();
   }    
      

    
    //-------------------------------------------------------------------------------------------------------------------------------------        
    /// Deleting the copy and move constructors - no duplication/transfer in anyway
    read_setup(const read_setup&) = delete;               //copy constructor
    read_setup& operator=(const read_setup&) = delete;    //copy assignment operator
    read_setup(read_setup&&) = delete;                    //move constructor  
    read_setup& operator=(read_setup&&) = delete;         //move assignment operator 
    //-------------------------------------------------------------------------------------------------------------------------------------
    
    
    
    /// This function reads data from "setup.txt" file and put it in a map<string, string>
    void read_data()
    {  
      std::ifstream file("setup.txt");
      if(!file.is_open())
      {
         Error("unable to open and read  'setup.txt'");
      }
      if(file.is_open())
      {
        std::string s;
        std::vector<std::string> split_result;
        while(getline(file, s))
        { 
          boost::trim(s);       
          if(!s.empty())
          {
               boost::split(split_result, s, boost::is_any_of(" "), boost::token_compress_on);
               data_[split_result[0]] = split_result[1];
          }          
       }
       file.close();       
     }
   }      








  
   /// This function reads data from the map and assigns variables of this class. Data type conversion is done by stoi and stod functions.
   void set_vars()
   {
   
     for(auto s : data_)
     {
       /// general parameters  
       if (s.first == "mesh_file")            mesh_file_ = "input/" + s.second;          
       if (s.first == "fluid_mesh_file")      fluid_mesh_file_ = "input/" + s.second;          
       if (s.first == "solid_mesh_file")      solid_mesh_file_ = "input/" + s.second;
       if (s.first == "heat_eq_mesh_file")    heat_eq_mesh_file_ = "input/" + s.second;
       if (s.first == "adv_diff_mesh_file")    adv_diff_mesh_file_ = "input/" + s.second;

       if (s.first == "delta_t")              delta_t_ = stod(s.second);
       if (s.first == "final_time")           end_time_ = stod(s.second);

       if (s.first == "restart")
          if(s.second == "T" || s.second == "True" || s.second == "TRUE" || s.second == "t"  || s.second == "true")      restart_ = true;

       if (s.first == "restart_time")         restart_time_ = stod(s.second);
       if (s.first == "n_max_it")             n_max_it_ = stoi(s.second);
       if (s.first == "print_freq")           print_freq_ = stoi(s.second);
 
       if (s.first == "precision")            precision_ = stoi(s.second);

       if (s.first == "t_dependent_dirichlet_bc")
          if(s.second == "T" || s.second == "True" || s.second == "TRUE" || s.second == "t"  || s.second == "true")            t_dependent_dirichlet_bc_ = true;

       if (s.first == "print_props")
          if(s.second == "T" || s.second == "True" || s.second == "TRUE" || s.second == "t"  || s.second == "true")            print_props_ = true;

       
       
       /// linear solver parameters                      
       if (s.first == "ls_solver")
       {
          ///---------------------Belos----------------------------------------
          if(s.second == "GMRES")                  ls_solver_ = "GMRES";         
          else if(s.second == "Flexible_GMRES")    ls_solver_ = "Flexible GMRES";
          else if(s.second == "Block")             ls_solver_ = "Block CG";
          else if(s.second == "PseudoBlockCG")     ls_solver_ = "PseudoBlockCG";
          else if(s.second == "Recycling_GMRES")   ls_solver_ = "Recycling GMRES";
          else if(s.second == "Recycling_CG")      ls_solver_ = "Recycling CG";
          else if(s.second == "MINRES")            ls_solver_ = "MINRES";
          else if(s.second == "LSQR")              ls_solver_ = "LSQR";
          else if(s.second == "TFQMR")             ls_solver_ = "TFQMR";
          else if(s.second == "Pseudoblock")       ls_solver_ = "Pseudoblock";
          else if(s.second == "GmresPoly")         ls_solver_ = "GmresPoly";
          else if(s.second == "CGPoly")            ls_solver_ = "CGPoly";

          ///---------------------Amesos----------------------------------------
          else if(s.second == "Klu")          ls_solver_ = "Klu";
          else if(s.second == "Lapack")       ls_solver_ = "Lapack";
          else if(s.second == "Umfpack")      ls_solver_ = "Umfpack";
          else if(s.second == "Pardiso")      ls_solver_ = "Pardiso";
          else if(s.second == "Taucs")        ls_solver_ = "Taucs";
          else if(s.second == "Superlu")      ls_solver_ = "Superlu";
          else if(s.second == "Superludist")  ls_solver_ = "Superludist";
          else if(s.second == "Mumps")        ls_solver_ = "Mumps";
          else if(s.second == "Dscpack")      ls_solver_ = "Dscpack";

          else ls_solver_ = "Flexible GMRES";
       } 
       if (s.first == "ls_precond")            ls_precond_ = s.second;
       if (s.first == "ls_rel_res_tol")        ls_rel_res_tol_ = stod(s.second);
       if (s.first == "ls_maxsubspace")        ls_maxsubspace_ = stoi(s.second);
       if (s.first == "ls_maxrestarts")        ls_maxrestarts_ = stoi(s.second);
       if (s.first == "ls_maxiters")           ls_maxiters_ = stoi(s.second);
       if (s.first == "ls_fill")               ls_fill_ = stoi(s.second);
       if (s.first == "ls_details")
          if(s.second == "T" || s.second == "True" || s.second == "TRUE" || s.second == "t"  || s.second == "true")            ls_details_ = true;




       /// adaptive time step parameters
       if (s.first == "adaptive_time_step")
          if(s.second == "T" || s.second == "True" || s.second == "TRUE" || s.second == "t"  || s.second == "true")            adaptive_time_step_ = true;

       if (s.first == "dt_min")         dt_min_ = stod(s.second);          
       if (s.first == "dt_max")         dt_max_ = stod(s.second);          
       if (s.first == "N")              N_ = stoi(s.second);          
       if (s.first == "CFL")            CFL_ = stod(s.second);




       /// fluid parameters
       if (s.first == "mesh_stiffness")     mesh_stiffness_ = stod(s.second); 

       if (s.first == "steady_state") 
          if(s.second == "T" || s.second == "True" || s.second == "TRUE" || s.second == "t"  || s.second == "true")          steady_state_ = true;

       if (s.first == "tau_diag_2014") 
          if(s.second == "T" || s.second == "True" || s.second == "TRUE" || s.second == "t"  || s.second == "true")            tau_diag_2014_ = true;

       if (s.first == "tau_diag_incomp_2007") 
          if(s.second == "T" || s.second == "True" || s.second == "TRUE" || s.second == "t"  || s.second == "true")            tau_diag_incomp_2007_ = true;

       if (s.first == "tau_non_diag_comp_2001") 
          if(s.second == "T" || s.second == "True" || s.second == "TRUE" || s.second == "t"  || s.second == "true")             tau_non_diag_comp_2001_ = true;

       if (s.first == "tau_non_diag_comp_2019") 
          if(s.second == "T" || s.second == "True" || s.second == "TRUE" || s.second == "t"  || s.second == "true")             tau_non_diag_comp_2019_ = true;

       if (s.first == "dc_1998")
          if(s.second == "T" || s.second == "True" || s.second == "TRUE" || s.second == "t"  || s.second == "true")             dc_1998_ = true;

       if (s.first == "dc_2006")
          if(s.second == "T" || s.second == "True" || s.second == "TRUE" || s.second == "t"  || s.second == "true")             dc_2006_ = true;

       if (s.first == "dc_sharp")                dc_sharp_ = stod(s.second);
       if (s.first == "dc_scale_fact")           dc_scale_fact_ = stod(s.second);





       /// solid parameters      
       if (s.first == "s_rho_inf")           s_rho_inf_ = stod(s.second);

       if (s.first == "constant_v")
          if(s.second == "T" || s.second == "True" || s.second == "TRUE" || s.second == "t"  || s.second == "true")             constant_v_ = true;

       if (s.first == "zero_a") 
          if(s.second == "T" || s.second == "True" || s.second == "TRUE" || s.second == "t"  || s.second == "true")             zero_a_ = true;





       /// post processing time parameters
       if (s.first == "pp_start_time")         pp_start_time_ = stod(s.second);          
       if (s.first == "pp_final_time")         pp_final_time_ = stod(s.second);          
       if (s.first == "pp_delta_t")            pp_delta_t_ = stod(s.second);          

       
       
       /// post processing secondary dofs
       if (s.first == "pp_rho") 
          if(s.second == "T" || s.second == "True" || s.second == "TRUE" || s.second == "t"  || s.second == "true")            pp_rho_ = true;

       if (s.first == "pp_mu") 
          if(s.second == "T" || s.second == "True" || s.second == "TRUE" || s.second == "t"  || s.second == "true")            pp_mu_ = true;

       if (s.first == "pp_cp") 
          if(s.second == "T" || s.second == "True" || s.second == "TRUE" || s.second == "t"  || s.second == "true")            pp_cp_ = true;

       if (s.first == "pp_cv") 
          if(s.second == "T" || s.second == "True" || s.second == "TRUE" || s.second == "t"  || s.second == "true")            pp_cv_ = true;

       if (s.first == "pp_alpha") 
          if(s.second == "T" || s.second == "True" || s.second == "TRUE" || s.second == "t"  || s.second == "true")            pp_alpha_ = true;

       if (s.first == "pp_beta") 
          if(s.second == "T" || s.second == "True" || s.second == "TRUE" || s.second == "t"  || s.second == "true")            pp_beta_ = true;

       if (s.first == "pp_sound_speed") 
          if(s.second == "T" || s.second == "True" || s.second == "TRUE" || s.second == "t"  || s.second == "true")            pp_sound_speed_ = true;

       if (s.first == "pp_kappa")
          if(s.second == "T" || s.second == "True" || s.second == "TRUE" || s.second == "t"  || s.second == "true")            pp_kappa_ = true;

       if (s.first == "pp_vf_g") 
          if(s.second == "T" || s.second == "True" || s.second == "TRUE" || s.second == "t"  || s.second == "true")            pp_vf_g_ = true;





       /// post processing drag lift
       if (s.first == "pp_rho_inf")         pp_rho_inf_ = stod(s.second);          
       if (s.first == "pp_v_inf")           pp_v_inf_ = stod(s.second);          
       if (s.first == "pp_bd_flag")         pp_bd_flag_ = stod(s.second);          
       if (s.first == "pp_cylinder_diameter")         pp_cylinder_diameter_ = stod(s.second);          
     }
   }  


 
 
 
 







     
    //--------------------------------------------------------inspector functions---------------------------------------------------------------------
    /// general parameters      
    auto mesh_file()const {return mesh_file_;}
    auto fluid_mesh_file()const {return fluid_mesh_file_;}
    auto solid_mesh_file()const {return solid_mesh_file_;}
    auto heat_eq_mesh_file()const {return heat_eq_mesh_file_;}
    auto adv_diff_mesh_file()const {return adv_diff_mesh_file_;}

    auto delta_t()const {return delta_t_;}
    auto end_time()const {return end_time_;}
    auto restart()const {return restart_;}
    auto restart_time()const {return restart_time_;}
    auto n_max_it()const {return n_max_it_;}
    auto print_freq()const {return print_freq_;}
    auto precision()const {return precision_;}
    auto t_dependent_dirichlet_bc()const {return t_dependent_dirichlet_bc_;}
    auto print_props()const {return print_props_;}


    /// linear solver parameters
    auto ls_solver()const {return ls_solver_;}
    auto ls_precond()const {return ls_precond_;}
    auto ls_rel_res_tol()const {return ls_rel_res_tol_;}
    auto ls_maxsubspace()const {return ls_maxsubspace_;}
    auto ls_maxrestarts()const {return ls_maxrestarts_;}
    auto ls_maxiters()const {return ls_maxiters_;}
    auto ls_fill()const {return ls_fill_;}
    auto ls_details()const{return ls_details_;}


    /// adaptive time step parameters
    auto adaptive_time_step()const {return adaptive_time_step_;}
    auto dt_min()const {return dt_min_;}
    auto dt_max()const {return dt_max_;}
    auto N()const {return N_;}
    auto CFL()const {return CFL_;}


    /// fluid parameters
    auto mesh_stiffness()const {return mesh_stiffness_;}
    auto steady_state()const {return steady_state_;}
    auto tau_diag_2014()const {return tau_diag_2014_;}
    auto tau_diag_incomp_2007()const {return tau_diag_incomp_2007_;}
    auto tau_non_diag_comp_2001()const {return tau_non_diag_comp_2001_;}
    auto tau_non_diag_comp_2019()const {return tau_non_diag_comp_2019_;}
    auto dc_1998()const {return dc_1998_;}
    auto dc_2006()const {return dc_2006_;}
    auto dc_sharp()const {return dc_sharp_;}
    auto dc_scale_fact()const {return dc_scale_fact_;}


    /// solid parameters    
    auto s_rho_inf()const {return s_rho_inf_;}
    auto constant_v()const {return constant_v_;}
    auto zero_a()const {return zero_a_;}


    ///for time of post processing 
    auto pp_start_time()const {return pp_start_time_;}
    auto pp_final_time()const {return pp_final_time_;}
    auto pp_delta_t()const {return pp_delta_t_;}


    ///for post processing of secondary dofs
    auto pp_rho()const {return pp_rho_;}
    auto pp_mu()const {return pp_mu_;}
    auto pp_cp()const {return pp_cp_;}
    auto pp_cv()const {return pp_cv_;}
    auto pp_alpha()const {return pp_alpha_;}
    auto pp_beta()const {return pp_beta_;}
    auto pp_sound_speed()const {return pp_sound_speed_;}
    auto pp_kappa()const {return pp_kappa_;}
    auto pp_vf_g()const {return pp_vf_g_;}
    
    
    ///for post processing of drag and lift of cylinder
    auto pp_rho_inf()const {return pp_rho_inf_;}
    auto pp_v_inf()const {return pp_v_inf_;}
    auto pp_bd_flag()const {return pp_bd_flag_;}
    auto pp_cylinder_diameter()const {return pp_cylinder_diameter_;}
    //---------------------------------------------------------------------------------------------------------------------------------------------


    






    
    



    private:
 
    //------------------------------------------------------------------- class member variables ----------------------------------------------------  
    /// this is map of data
    std::map<std::string, std::string> data_;


    /// general parameters  
    std::string mesh_file_ = "";
    std::string fluid_mesh_file_ = "";
    std::string solid_mesh_file_ = "";
    std::string heat_eq_mesh_file_ = "";
    std::string adv_diff_mesh_file_ = "";
    
    double delta_t_      = 0.002;
    double end_time_     = 100000.0;
    bool restart_        = false;
    double restart_time_ = 0.0;
    int n_max_it_        = 3;
    int print_freq_      = 10;
    int precision_       = 6;
    bool t_dependent_dirichlet_bc_ = true;
    bool print_props_ = false;    


    /// linear solver parameters
    std::string ls_solver_    = "Flexible GMRES";
    std::string ls_precond_   = "ILU";
    double ls_rel_res_tol_ = 1.e-9;
    int ls_maxsubspace_  = 200;
    int ls_maxrestarts_  = 5;
    int ls_maxiters_     = 3000;
    int ls_fill_         = 1;
    bool ls_details_     = false;


    /// adaptive time step parameters
    bool adaptive_time_step_ = false;
    double dt_min_       = 0.002;
    double dt_max_       = 1.0;    
    int N_               = 10;    
    double CFL_          = 1.0;        

        
    /// fluid parameters
    double mesh_stiffness_ = 0.0;
    bool steady_state_             = false;
    bool tau_diag_incomp_2007_     = false;
    bool tau_diag_2014_            = false;
    bool tau_non_diag_comp_2001_   = false;
    bool tau_non_diag_comp_2019_   = false;
    bool dc_1998_                  = false;
    bool dc_2006_                  = false;
    double dc_sharp_      = 1.0;
    double dc_scale_fact_ = 1.0;


    /// solid parameters    
    double s_rho_inf_     = 0.5;    
    bool constant_v_              = false;
    bool zero_a_                  = false;


    ///for time of post processing 
    double pp_start_time_         = 0.0;
    double pp_final_time_         = 100000.0;
    double pp_delta_t_            = 10.0;
    
    
    ///for post processing of secondary dofs
    bool pp_rho_                  = false;
    bool pp_mu_                   = false;
    bool pp_cp_                   = false;
    bool pp_cv_                   = false;
    bool pp_alpha_                = false;
    bool pp_beta_                 = false;
    bool pp_sound_speed_          = false;
    bool pp_kappa_                = false;
    bool pp_vf_g_                 = false;
    
    
    ///for post processing of drag and lift of cylinder
    double pp_rho_inf_            = 1.0;
    double pp_v_inf_              = 1.0;
    int pp_bd_flag_               = 5;
    double pp_cylinder_diameter_  = 1.0;
    //-------------------------------------------------------------------------------------------------------------------------
  };


}


#endif

