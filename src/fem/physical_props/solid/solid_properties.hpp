#ifndef GALES_SOLID_PROPERTIES_HPP
#define GALES_SOLID_PROPERTIES_HPP



#include "solid_properties_reader.hpp"



namespace GALES {


 /**
      This class defined the physical properties of material.
      Propeties are read from "props.txt" by solid_properties_reader class       
 */



class solid_properties
{
  
  public:

  //--------------------------------constructor------------------------------------
  solid_properties()
  {  
      solid_properties_reader reader(*this);
  }
  //-------------------------------------------------------------------------------



  //-------------------------------------------------------------------------------------------------------------------------------------        
  /// Deleting the copy and move constructors - no duplication/transfer in anyway
  solid_properties(const solid_properties&) = delete;               //copy constructor
  solid_properties& operator=(const solid_properties&) = delete;    //copy assignment operator
  solid_properties(solid_properties&&) = delete;                    //move constructor  
  solid_properties& operator=(solid_properties&&) = delete;         //move assignment operator 
  //-------------------------------------------------------------------------------------------------------------------------------------



  
  //-------------------------------------------------------------------------------
  void layer_properties(double el_y_min, double el_y_max)
  {
     for(int i=0; i<num_layers_; i++)
        if(layers_y_min_[i] <= el_y_min && el_y_max <= layers_y_max_[i])
        {
            rho_ = layers_rho_[i];
            E_ = layers_E_[i];
            nu_ = layers_nu_[i];
        }
  } 
  //-------------------------------------------------------------------------------
   
   
  
  
  
   
  //-------------------------------------------------------------------------------  
  void material_type(const std::string& s){material_type_ = s;}   
  auto material_type()const {return material_type_;}   
  
  void condition_type(const std::string& s){condition_type_ = s;}   
  auto condition_type()const {return condition_type_;}   
  
  void rho(double x){rho_ = x;}
  auto rho() const {return rho_;}
    
  void E(double x){E_ = x;}
  auto E() const {return E_;}  

  void nu(double x){nu_ = x;}
  auto nu() const {return nu_;}  

  void a_damping(double x){a_damping_ = x;}
  auto a_damping() const {return a_damping_;}  

  void b_damping(double x){b_damping_ = x;}
  auto b_damping() const {return b_damping_;}  
  
  void nb_Maxwell_el(int x){nb_Maxwell_el_ = x;}
  auto nb_Maxwell_el() const {return nb_Maxwell_el_;}
    
  void Max_el_E(const std::vector<double>& v){Max_el_E_ = v;}
  auto Max_el_E() const {return Max_el_E_;}  
  
  void Max_el_eta(const std::vector<double>& v){Max_el_eta_ = v;}
  auto Max_el_eta() const {return Max_el_eta_;}  

  void heterogeneous(bool f) {heterogeneous_ = f;}  
  auto heterogeneous() const {return heterogeneous_;}  
  
  void num_layers(int i){num_layers_ = i;}
  auto num_layers() const {return num_layers_;}  
  
  void layers_y_min(double d){layers_y_min_.push_back(d);}
  auto layers_y_min() const {return layers_y_min_;}  
  
  void layers_y_max(double d){layers_y_max_.push_back(d);}
  auto layers_y_max() const {return layers_y_max_;}  

  void layers_rho(double d){layers_rho_.push_back(d);}
  auto layers_rho() const {return layers_rho_;}  

  void layers_E(double d){layers_E_.push_back(d);}
  auto layers_E() const {return layers_E_;}  

  void layers_nu(double d){layers_nu_.push_back(d);}
  auto layers_nu() const {return layers_nu_;}  
  //-------------------------------------------------------------------------------
    
   
   

  private:
   
  std::string material_type_ = "Hookes";

  // Plain_Strain: one dimension very large in comparison to other two; e.g. beam  
  // Plain_Stress: one dimesion is very small in comparison to other two; e.g. thin plate
  std::string condition_type_ = "";  
      
  double rho_ = 0.0;
  double E_ = 0.0;
  double nu_ = 0.0;    
  
  double a_damping_ = 0.0;
  double b_damping_ = 0.0;

  int nb_Maxwell_el_ = 1;
  std::vector<double> Max_el_E_;
  std::vector<double> Max_el_eta_;

  bool heterogeneous_ = false;
  int num_layers_ = 0;
  std::vector<double> layers_y_min_;
  std::vector<double> layers_y_max_;
  std::vector<double> layers_rho_;
  std::vector<double> layers_E_;
  std::vector<double> layers_nu_;
    
   
};


}

#endif
