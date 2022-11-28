#ifndef GALES_GIORDANOETALVISCOSITYMODEL_HPP
#define GALES_GIORDANOETALVISCOSITYMODEL_HPP


#include <math.h>
#include <vector>
#include <numeric>


namespace GALES {


  /**
      Girdano et el viscosity model as in
      Giordano, D. et al. Viscosity of magmatic liquids: a model. Earth and Planetary Science Letters 271, 123–134 (2008)
  */


  class GiordanoEtAlViscosityModel 
  {


    public:
    
    
    //-------------------------------------------------------------------------------------------------------------------------------------        
    GiordanoEtAlViscosityModel() = default; //constructor
    /// Deleting the copy and move constructors - no duplication/transfer in anyway
    GiordanoEtAlViscosityModel(const GiordanoEtAlViscosityModel&) = delete;               //copy constructor
    GiordanoEtAlViscosityModel& operator=(const GiordanoEtAlViscosityModel&) = delete;    //copy assignment operator
    GiordanoEtAlViscosityModel(GiordanoEtAlViscosityModel&&) = delete;                    //move constructor  
    GiordanoEtAlViscosityModel& operator=(GiordanoEtAlViscosityModel&&) = delete;         //move assignment operator 
    //-------------------------------------------------------------------------------------------------------------------------------------
    

    void set(const std::vector<double>& oxides_wf, const std::vector<double>& oxides_molar_mass)
    {
      auto ox_molar_mass = oxides_molar_mass;               //in kg/mol
      for(auto& a: ox_molar_mass) a *= 1.e3;              // in g/mol
      auto wf_magma = oxides_wf;                       //SiO2[0]  TiO2[1]  Al2O3[2]  Fe2O3[3]  FeO[4]  MnO[5]  MgO[6]  CaO[7]  Na2O[8]  K2O[9]
      wf_magma[4] = 0.9309*(wf_magma[3]+wf_magma[4]);         //compute the Fe total (Fe2O3 + FeO)
      wf_magma[3] = 0.0;
      const double sum = std::accumulate(wf_magma.begin(), wf_magma.end(), 0.0);          //closeToOne(wf_magma);
      for(auto& a: wf_magma) a /= sum;
      
      for(std::size_t i=0; i<wf_magma.size(); i++)
        molar_mass_ += wf_magma[i]/ox_molar_mass[i];        
      molar_mass_ = 1.0/molar_mass_;           // in g w.r.t. new wf_magma

      std::vector<double> mole_percent(10, 0.0);             
      for(std::size_t i=0; i<wf_magma.size(); i++)
       mole_percent[i] = 100.0*wf_magma[i]*molar_mass_/ox_molar_mass[i];      //mole_percent[i] = mole_fraction[i]*100

      parameters_[0] = -4.55;
                                                                                                                                     
      // parameters are computed for anhydrous mixture (without H2O)      
      parameters_[1] = b_[0]*(mole_percent[0]+mole_percent[1])    // SiO2+TiO2
                     + b_[1]*mole_percent[2]                      // Al2O3
                     + b_[2]*(mole_percent[4]+mole_percent[5])    // FeO+MnO
                     + b_[3]*mole_percent[6]                      // MgO
                     + b_[4]*mole_percent[7]                      // CaO
                     + b_[5]*mole_percent[8];                     // Na2O
      
      parameters_[2] = b1_[0]*(mole_percent[0]+mole_percent[1])*(mole_percent[4]+mole_percent[5]+mole_percent[6])        // (SiO2+TiO2)*(FeO(T)+MnO+MgO) 
                     + b1_[1]*(mole_percent[0]+mole_percent[1]+mole_percent[2])*(mole_percent[8]+mole_percent[9])        // (SiO2+TiO2+Al2O3)⁎(Na2O+K2O)
                     + b1_[2]*mole_percent[2]*(mole_percent[8]+mole_percent[9]);                                         // (Al2O3)⁎(Na2O+K2O)
            
      parameters_[3] = b1_[1]*(mole_percent[0]+mole_percent[1]+mole_percent[2]);            //(SiO2+TiO2+Al2O3) to be multiplied with mole_percent of H2O

      parameters_[4] = c_[0]*mole_percent[0]                                    //SiO2
                     + c_[1]*(mole_percent[1]+mole_percent[2])                  //TiO2+Al2O3
                     + c_[2]*(mole_percent[4]+mole_percent[5]+mole_percent[6])  //FeO(T)+MnO+MgO
                     + c_[3]*mole_percent[7]                                    //CaO
                     + c_[4]*(mole_percent[8]+mole_percent[9]);                 //Na2O+K2O
                         
      parameters_[5] = 0.3*(mole_percent[2]+mole_percent[4]+mole_percent[5]+mole_percent[6]+mole_percent[7])*(mole_percent[8]+mole_percent[9]);    // (Al2O3+FeO(T)+MnO+MgO+CaO)*(Na2O+K2O)

      parameters_[6] = 0.3*(mole_percent[2]+mole_percent[4]+mole_percent[5]+mole_percent[6]+mole_percent[7]);      // (Al2O3+FeO(T)+MnO+MgO+CaO) to be multiplied with mole_percent of H2O      
    }





    double computeViscosity(double T, double wf_h2o_l, double wf_co2_l)
    {
      const double h2o_mola_mass =  18.015 ;          // in grams
      
      if(1.0 - wf_co2_l == 0.0)
      {
         Error("1-wf_co2_l is zero in denominator");                    
      }      
      const double new_wf_h2o = wf_h2o_l/(1-wf_co2_l);
      double mole_fr_h2o = 1 + h2o_mola_mass*(1-new_wf_h2o)/(new_wf_h2o*molar_mass_);
      mole_fr_h2o = 1.0/mole_fr_h2o;
      const double mole_fr_no_h2o = 1-mole_fr_h2o;      
      const double h2o_mole_percent = 100.0*mole_fr_h2o;
      
      const double A = parameters_[0];
      const double B = parameters_[1]*mole_fr_no_h2o + parameters_[2]*mole_fr_no_h2o*mole_fr_no_h2o + parameters_[3]*mole_fr_no_h2o*h2o_mole_percent 
                 + (b_[5]+b_[6])*h2o_mole_percent + b_[6]*log(1+h2o_mole_percent);                
                 
      const double C = parameters_[4]*mole_fr_no_h2o + c_[5]*log(1+h2o_mole_percent) + parameters_[5]*mole_fr_no_h2o*mole_fr_no_h2o + parameters_[6]*mole_fr_no_h2o*h2o_mole_percent;
      
      if(T-C == 0.0)
      {
         Error("T-C is zero in denominator");   
      }      
      const double log_mu = A+B/(T-C);
      double mu = pow(10., log_mu);
      return mu;
    }


    private:

    double molar_mass_ = 0.0;

          // SiO2+TiO2   Al2O3   FeO+MnO    MgO    CaO   Na2O+H2O  H2O+ln(1+H2O)
    std::vector<double> b_ = {159.6,    -173.3,    72.1,    75.7, -39.0, -84.1,    141.5};

       // (SiO2+TiO2)⁎(FeO(T)+MnO+MgO)    (SiO2+TiO2+Al2O3)⁎(Na2O+K2O+H2O)      (Al2O3)⁎(Na2O+K2O)
    std::vector<double> b1_ = {-2.43,                            -0.91,                                  17.6};

        //  SiO2   TiO2+Al2O3   FeO(T)+MnO+MgO    CaO    Na2O+K2O     ln(1+H2O)
    std::vector<double> c_ = {2.75,   15.7,         8.3,            10.2,   -12.3,       -99.5};

    std::vector<double> parameters_ = std::vector<double>(7, 0.0);
    
  };


}
#endif
