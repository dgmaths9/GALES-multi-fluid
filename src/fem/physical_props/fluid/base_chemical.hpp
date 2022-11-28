#ifndef GALES_BASE_CHEMICAL_HPP
#define GALES_BASE_CHEMICAL_HPP



#include<vector>
#include<memory>


namespace GALES {


  /**
       base chemical class
       it is supposed to be used as base class when we define other chemicals
  */



class base_chemical
{

 public:


  //-------------------------------------------------------------------------------------------------------------------------------------        
  base_chemical() = default;               //constructor
  /// Deleting the copy and move constructors - no duplication/transfer in anyway
  base_chemical(const base_chemical&) = delete;               //copy constructor
  base_chemical& operator=(const base_chemical&) = delete;    //copy assignment operator
  base_chemical(base_chemical&&) = delete;                    //move constructor  
  base_chemical& operator=(base_chemical&&) = delete;         //move assignment operator 
  //-------------------------------------------------------------------------------------------------------------------------------------


  virtual void properties(double p, double T, double strain_rate) {}
  virtual void properties(double p, double strain_rate) {}

  virtual std::vector<double> oxides_wf()
  { 
     std::vector<double> v(0);
     return v;
  }
  
  virtual std::vector<double> oxides_molar_mass()
  {
     std::vector<double> v(0);
     return v;  	
  }



  /// this is used in fluid_sc_isothermal to print properties for debugging purpose
  void print_props_sc_isothermal()
  {
       print_in_file("properties", "rho", rho_);
       print_in_file("properties", "mu", mu_);
       print_in_file("properties", "sound_speed", sound_speed_);
       print_in_file("properties", "beta", beta_);
       print_in_file("properties", " ");
  }



  /// this is used in fluid_sc to print properties for debugging purpose
  void print_props_sc()
  {
       print_in_file("properties", "rho", rho_);
       print_in_file("properties", "mu", mu_);
       print_in_file("properties", "kappa", kappa_);
       print_in_file("properties", "sound_speed", sound_speed_);
       print_in_file("properties", "alpha", alpha_);
       print_in_file("properties", "beta", beta_);
       print_in_file("properties", "cv", cv_);
       print_in_file("properties", "cp", cp_);
       print_in_file("properties", " ");
  }


  int l_ = 0;
  int g_ = 1;
  double rho_ = 0.0;
  double cv_ = 0.0;
  double cp_ = 0.0;
  double kappa_ = 0.0;
  double chemical_diffusivity_ = 1.e-10;
  double mu_ = 0.0;
  double sound_speed_ = 0.0;
  double alpha_ = 0.0;
  double beta_ = 0.0;
};



}


#endif
