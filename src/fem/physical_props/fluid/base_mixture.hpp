#ifndef GALES_BASE_MIXTURE_HPP
#define GALES_BASE_MIXTURE_HPP


#include <vector>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

#include "base_chemical.hpp"



namespace GALES {



  /**
       base mixture class
       it is supposed to be used as base class when we define other mixtures
  */




class base_mixture : public base_chemical
{

  using vec = boost::numeric::ublas::vector<double>;
  

  public:

  //-------------------------------------------------------------------------------------------------------------------------------------        
  base_mixture() = default;               //constructor
  /// Deleting the copy and move constructors - no duplication/transfer in anyway
  base_mixture(const base_mixture&) = delete;               //copy constructor
  base_mixture& operator=(const base_mixture&) = delete;    //copy assignment operator
  base_mixture(base_mixture&&) = delete;                    //move constructor  
  base_mixture& operator=(base_mixture&&) = delete;         //move assignment operator 
  //-------------------------------------------------------------------------------------------------------------------------------------

     
  virtual void properties(double p, double T, double strain_rate) {}
  virtual void properties(double p, double strain_rate) {}
  virtual void properties(double p, double T, const vec &Y, double strain_rate){}
  virtual void properties(double p, const vec &Y, double strain_rate){}



  void clear()
  {
    rho_ = 0.0;
    cv_ = 0.0;
    cp_ = 0.0;
    kappa_ = 0.0;
    chemical_diffusivity_ = 0.0;
    mu_ = 0.0;
    sound_speed_ = 0.0;
    alpha_ = 0.0;
    beta_ = 0.0;
    
    std::fill(rho_comp_.begin(), rho_comp_.end(), 0.0);
    std::fill(chemical_diffusivity_comp_.begin(), chemical_diffusivity_comp_.end(), 0.0);
    std::fill(internal_energy_comp_.begin(), internal_energy_comp_.end(), 0.0);
  }


  ///this prints properties in fluid_mc_isothermal for debugging purpose
  void print_props_mc_isothermal()
  {
       print_in_file("properties", "rho", rho_);
       print_in_file("properties", "mu", mu_);
       print_in_file("properties", "sound_speed", sound_speed_);
       print_in_file("properties", "beta", beta_);
       print_in_file("properties", "rho_comp", rho_comp_);
       print_in_file("properties", "chemical_diffusivity_comp", chemical_diffusivity_comp_);
       print_in_file("properties", " ");
  }


  ///this prints properties in fluid_mc for debugging purpose
  void print_props_mc()
  {
       print_in_file("properties", "rho", rho_);
       print_in_file("properties", "mu", mu_);
       print_in_file("properties", "kappa", kappa_);
       print_in_file("properties", "sound_speed", sound_speed_);
       print_in_file("properties", "alpha", alpha_);
       print_in_file("properties", "beta", beta_);
       print_in_file("properties", "cv", cv_);
       print_in_file("properties", "cp", cp_);
       print_in_file("properties", "rho_comp", rho_comp_);
       print_in_file("properties", "internal_energy_comp", internal_energy_comp_);
       print_in_file("properties", "chemical_diffusivity_comp", chemical_diffusivity_comp_);
       print_in_file("properties", " ");
  }




  int nb_comp_;

  vec wf_;
  vec vf_;
  
  vec wf_ph_ = vec(2, 0.0);                       // wf_ph_[liquid=0],  wf_ph_[gas=1]
  vec vf_ph_ = vec(2, 0.0);                       // vf_ph_[liquid=0],  vf_ph_[gas=1]

  vec rho_comp_;
  vec internal_energy_comp_;
  vec chemical_diffusivity_comp_;

  std::vector<std::shared_ptr<base_chemical>> chemical_ptrs_;
  std::vector<std::shared_ptr<base_mixture>> mixture_ptrs_;

};


}

#endif
