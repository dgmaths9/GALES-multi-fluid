#ifndef GALES_CHEMICALS_HPP
#define GALES_CHEMICALS_HPP



#include <string>
#include <vector>
#include "base_chemical.hpp"



namespace GALES {





/// this is for custom material
class custom_chemical : public base_chemical
{
   public:
   
   custom_chemical(double rho, double mu, double cp, double kappa, double alpha, double beta)
   {
      rho_ = rho;
      mu_ = mu;
      cp_ = cp;
      alpha_ = alpha;         // this will be zero for fully incompressible fluid and its flow
      beta_ = beta;           // this will be zero for fully incompressible fluid and its flow
      kappa_ = kappa;      
      if(beta_==0.0) sound_speed_ = 1.e20;   // We set very high value for fully incompressible fluid and its flow(in reality it should be infinite)
      else sound_speed_ = 1.0/sqrt(rho_*beta_); 
   } 

  //-------------------------------------------------------------------------------------------------------------------------------------        
  /// Deleting the copy and move constructors - no duplication/transfer in anyway
  custom_chemical(const custom_chemical&) = delete;               //copy constructor
  custom_chemical& operator=(const custom_chemical&) = delete;    //copy assignment operator
  custom_chemical(custom_chemical&&) = delete;                    //move constructor  
  custom_chemical& operator=(custom_chemical&&) = delete;         //move assignment operator 
  //-------------------------------------------------------------------------------------------------------------------------------------

  void properties(double p, double T, double strain_rate) final
  {
    if(beta_==0.0)  cv_ = cp_;	// this is true for fully incompressible fluid and its flow
    else cv_ = cp_ - T*alpha_*alpha_/(rho_*beta_);
  }

};





/// this is for melt to be used in magma mixture
class melt : public base_chemical 
{ 
  using vec = std::vector<double>; 
  
  public:    	

  melt(const vec& oxides_wf, double kappa) :  oxides_wf_(oxides_wf)
  {
    const double sum = std::accumulate(oxides_wf_.begin(), oxides_wf_.end(), 0.0);
    for(auto& a: oxides_wf_) a /= sum;
    
    
    // kg/mol            //SiO2   TiO2    Al2O3    Fe2O3    FeO     MnO     MgO     CaO     Na2O    K2O
    oxides_molar_mass_ = {60.084, 79.876, 101.956, 159.692, 71.846, 70.937, 40.304, 56.078, 61.979, 94.195};
    for(auto& a: oxides_molar_mass_) a *= 1.e-3;


    molar_mass_ = 0.0;
    for (std::size_t i=0; i<10; i++)   molar_mass_ += oxides_wf_[i]/oxides_molar_mass_[i];
    molar_mass_ = 1./molar_mass_;     //kg

    
    for(std::size_t i=0; i<10; i++)   oxides_molar_fraction_[i] =  oxides_wf_[i]*molar_mass_/oxides_molar_mass_[i];

            
    vec oxides_dV_dT = {0.0,   7.24,   0.0,   9.09,   2.92, 0.0, 3.27, 3.74, 7.68, 12.08};     //m^3/mol.K 
    for(auto& a: oxides_dV_dT) a *= 1.e-9;
    dV_dT_ = 0.0;
    for (std::size_t i=0; i<10; i++)    dV_dT_ += oxides_molar_fraction_[i]*oxides_dV_dT[i];


    vec oxides_dV_dp = {-1.89, -2.31, -2.26,  -2.53, -0.45, 0.0, 0.27, 0.34, -2.4, -6.75};     //m^3/mol.Pa
    for(auto& a: oxides_dV_dp) a *= 1.e-15;
    dV_dp_ = 0.0;
    for (std::size_t i=0; i<10; i++)    dV_dp_ += oxides_molar_fraction_[i]*oxides_dV_dp[i];


    vec oxides_molar_volume = {26.86, 23.16, 37.42, 42.13, 13.65, 0.0, 11.69, 16.53, 28.88, 45.07};     //m^3/mol
    for(auto& a: oxides_molar_volume) a *= 1.e-6;
    molar_volume_ = 0.0;
    for (std::size_t i=0; i<10; i++)    molar_volume_ += oxides_molar_fraction_[i]*oxides_molar_volume[i];      


                  //SiO2  TiO2   Al2O3  Fe2O3  FeO   MnO  MgO   CaO   Na2O    K2O
    // oxides_cp = {80.0, 111.8 ,157.6, 229.0, 78.9, 0.0, 99.7, 99.9, 102.3, 97.0};     // In J{gfw}^{-1}K^{-1}  ;  gfw---gram formula weight
    //  1 {gfw}^{-1} = 1/molar_mass g^{-1} = 1000/molar_mass kg^{-1}
                          //SiO2   TiO2   Al2O3   Fe2O3    FeO    MnO   MgO     CaO     Na2O    K2O
    const vec oxides_cp = {1331.0, 1392.0, 1545.0, 1434.0, 1100.0, 0.0, 2424.0, 1781.0, 1651.0, 1030.0};    // In Jkg^{-1}K^{-1}
    cp_ = 0.0;
    for (std::size_t i=0; i<oxides_cp.size() ; i++)  cp_ +=  oxides_molar_fraction_[i]*oxides_cp[i];

    kappa_ = kappa;	
  }


  //-------------------------------------------------------------------------------------------------------------------------------------        
  /// Deleting the copy and move constructors - no duplication/transfer in anyway
  melt(const melt&) = delete;               //copy constructor
  melt& operator=(const melt&) = delete;    //copy assignment operator
  melt(melt&&) = delete;                    //move constructor  
  melt& operator=(melt&&) = delete;         //move assignment operator 
  //-------------------------------------------------------------------------------------------------------------------------------------


  void properties(double p, double T, double strain_rate) final
  {       
    const double molar_vol =  molar_volume_ + dV_dT_*(T-1673.0) + dV_dp_*(p-1.e5);     //Lange Carmichal
    alpha_ = 1./molar_vol*dV_dT_; 
    beta_ = -1./molar_vol*dV_dp_;           
    rho_ = molar_mass_/molar_vol;
    sound_speed_ = 1.0/sqrt(rho_*beta_);
    cv_ = cp_ - T*alpha_*alpha_/(rho_*beta_);	
  } 
  
  vec oxides_wf()final {return oxides_wf_;}
  vec oxides_molar_mass()final {return oxides_molar_mass_;}
  
  private:  
  vec oxides_wf_ = vec(10, 0.0);
  vec oxides_molar_mass_ = vec(10, 0.0); 
  vec oxides_molar_fraction_ = vec(10, 0.0);
  double dV_dT_ = 0.0;
  double dV_dp_ = 0.0;
  double molar_volume_ = 0.0;  
  double molar_mass_ = 0.0;
};






/// this is for dissolved_volatiles to be used in magma mixture
class dissolved_volatiles : public base_chemical 
{
public:  	
  dissolved_volatiles(std::string s, double cp, double kappa): a_(10) 
  {
    if(s=="H2O" || s == "h2o" || s == "H2o")         molar_mass_ = 0.01802;
    else if(s=="CO2" || s=="co2" || s=="Co2")    molar_mass_ = 0.04401;

    cp_ = cp;
    kappa_ = kappa;
    a_ = {9.144e-06, 3.685e-09, 1.168e-11, -1.990e-15, 1.220e-16, -1.945e-17, -1.580e-21, 4.680e-24, 1.144e-26, -3.960e-33};
  }


  //-------------------------------------------------------------------------------------------------------------------------------------        
  /// Deleting the copy and move constructors - no duplication/transfer in anyway
  dissolved_volatiles(const dissolved_volatiles&) = delete;               //copy constructor
  dissolved_volatiles& operator=(const dissolved_volatiles&) = delete;    //copy assignment operator
  dissolved_volatiles(dissolved_volatiles&&) = delete;                    //move constructor  
  dissolved_volatiles& operator=(dissolved_volatiles&&) = delete;         //move assignment operator 
  //-------------------------------------------------------------------------------------------------------------------------------------


  void properties(double p, double T, double strain_rate) final
  {
    double mole_vol = a_[0] + a_[1]*T + a_[2]*T*T + a_[3]*T*T*T + p*(a_[4]+a_[5]*T+a_[6]*T*T) + p*p*(a_[7]+a_[8]*T) + a_[9]*p*p*p;
    rho_ = molar_mass_/mole_vol;

    alpha_ = a_[1] + 2*a_[2]*T + 3*a_[3]*T*T + p*(a_[5]+2*a_[6]*T) + p*p*a_[8];
    alpha_ /= mole_vol;
  
    beta_ = -((a_[4]+a_[5]*T+a_[6]*T*T) + 2*p*(a_[7]+a_[8]*T) + 3*a_[9]*p*p);   
    beta_ /= mole_vol;  
 
    sound_speed_ = 1.0/sqrt(rho_*beta_); 
    cv_ = cp_ - T*alpha_*alpha_/(rho_*beta_);	
  }
  
  private:
  std::vector<double> a_;
  double molar_mass_ = 0.0;
};






/// ideal gas class
class ideal_gas : public base_chemical 
{
  public:  	
  ideal_gas(std::string s)
  { 
    if(s=="H2O" || s == "h2o" || s == "H2o")
    {  
      molar_mass_ = 0.01802;
      gamma_ = 1.4;
    }      
    else if(s=="CO2" || s=="co2" || s=="Co2")
    {  
      molar_mass_ = 0.04401;
      gamma_ = 1.4;
    }  
    else if(s=="AIR" || s=="air" || s=="Air")
    {  
      molar_mass_ = 0.02897;
      gamma_ = 1.4;
    }  
    else if(s=="HE" || s=="he" || s=="He")
    {  
      molar_mass_ = 0.004;
      gamma_ = 1.67;
    }  
    else if(s=="ARGON" || s=="argon" || s=="Argon")
    {  
      molar_mass_ = 0.06;
      gamma_ = 1.67;
    }  
    
    R_ = 8.314462/molar_mass_;
    cv_ = R_/(gamma_-1.0);
    cp_ = gamma_*cv_;
    mu_ = 1.e-5;    
    kappa_ = mu_*cp_/0.72;     // Pr = 0.72
   }


  ideal_gas(double molar_mass, double gamma): base_chemical() 
  { 
    molar_mass_ = molar_mass;
    gamma_ = gamma;
    R_ = 8.314462/molar_mass_;
    cv_ = R_/(gamma_-1.0);
    cp_ = gamma_*cv_;
    mu_ = 1.e-5;    
    kappa_ = mu_*cp_/0.72;      // Pr = 0.72
  }


  //-------------------------------------------------------------------------------------------------------------------------------------        
  /// Deleting the copy and move constructors - no duplication/transfer in anyway
  ideal_gas(const ideal_gas&) = delete;               //copy constructor
  ideal_gas& operator=(const ideal_gas&) = delete;    //copy assignment operator
  ideal_gas(ideal_gas&&) = delete;                    //move constructor  
  ideal_gas& operator=(ideal_gas&&) = delete;         //move assignment operator 
  //-------------------------------------------------------------------------------------------------------------------------------------


  void properties(double p, double T, double strain_rate) final
  { 
    rho_ = p/(R_*T);
    beta_ = 1.0/p;
    alpha_ = 1.0/T;	 
    sound_speed_ = sqrt(gamma_*R_*T);
  }


private:
  double gamma_ = 0.0;
  double R_ = 0.0;
  double molar_mass_ = 0.0;
};



}


#endif
 
