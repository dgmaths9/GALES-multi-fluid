#ifndef GALES_FLUID_PROPERTIES_HPP
#define GALES_FLUID_PROPERTIES_HPP



#include "fluid_properties_reader.hpp"



namespace GALES {


 /**
      This class defined the physical properties of fluid.
      Propeties along with their configuration such as "mix_of_mix" etc. are read from an input "props.txt" by fluid_properties_reader class       
 */



class fluid_properties : public base_mixture
{

  using vec = boost::numeric::ublas::vector<double>;
  
  public:

  fluid_properties()
  {  
      nb_comp_ = 1;
      fluid_properties_reader reader(*this);
  }



  //-------------------------------------------------------------------------------------------------------------------------------------        
  /// Deleting the copy and move constructors - no duplication/transfer in anyway
  fluid_properties(const fluid_properties&) = delete;               //copy constructor
  fluid_properties& operator=(const fluid_properties&) = delete;    //copy assignment operator
  fluid_properties(fluid_properties&&) = delete;                    //move constructor  
  fluid_properties& operator=(fluid_properties&&) = delete;         //move assignment operator 
  //-------------------------------------------------------------------------------------------------------------------------------------
  



  /// This is called from fluid_sc_isothermal
  void properties(double p, double strain_rate)final
  {
     if(mix_of_mix_) 
     {
       mixture_ptrs_[0]->properties(p, isothermal_T_, strain_rate);
       rho_ = mixture_ptrs_[0]->rho_;
       mu_ = mixture_ptrs_[0]->mu_;    
       sound_speed_ = mixture_ptrs_[0]->sound_speed_;
       beta_ = mixture_ptrs_[0]->beta_;       
     }     
     else 
     {
       chemical_ptrs_[0]->properties(p, isothermal_T_, strain_rate);
       rho_ = chemical_ptrs_[0]->rho_;
       mu_ = chemical_ptrs_[0]->mu_;    
       sound_speed_ = chemical_ptrs_[0]->sound_speed_;
       beta_ = chemical_ptrs_[0]->beta_;     
     }
  }






  /// This is called from fluid_sc
  void properties(double p, double T, double strain_rate)final
  {
     if(mix_of_mix_) 
     {
       mixture_ptrs_[0]->properties(p,T,strain_rate);
       rho_ = mixture_ptrs_[0]->rho_;
       mu_ = mixture_ptrs_[0]->mu_;    
       kappa_ = mixture_ptrs_[0]->kappa_;     
       sound_speed_ = mixture_ptrs_[0]->sound_speed_;
       cv_ = mixture_ptrs_[0]->cv_;
       cp_ = mixture_ptrs_[0]->cp_;
       alpha_ = mixture_ptrs_[0]->alpha_;	 
       beta_ = mixture_ptrs_[0]->beta_;       
     }     
     else 
     {
       chemical_ptrs_[0]->properties(p,T,strain_rate);
       rho_ = chemical_ptrs_[0]->rho_;
       mu_ = chemical_ptrs_[0]->mu_;    
       kappa_ = chemical_ptrs_[0]->kappa_;     
       sound_speed_ = chemical_ptrs_[0]->sound_speed_;
       cv_ = chemical_ptrs_[0]->cv_;
       cp_ = chemical_ptrs_[0]->cp_;
       alpha_ = chemical_ptrs_[0]->alpha_;	 
       beta_ = chemical_ptrs_[0]->beta_;     
     }
    
     if(Sutherland_law_)
     {
         mu_ = Sutherland_law_a_*sqrt(T*T*T)/(T+Sutherland_law_b_);
         kappa_ = mu_*cp_/0.72;   //Pr = 0.72           
     }
  }
  


   

  /// This is called from fluid_mc_isothermal
  void properties(double p, const vec &Y, double strain_rate)final
  {
     if(mix_of_mix_)  /// This is for mixture of submixtures
     {
       clear();
       for (int c=0; c<nb_comp_; c++)
       {
         mixture_ptrs_[c]->properties(p, isothermal_T_, strain_rate);
         rho_ += Y[c]/mixture_ptrs_[c]->rho_;
         mu_ += Y[c]*mixture_ptrs_[c]->mu_;
         sound_speed_ += Y[c]/(pow(mixture_ptrs_[c]->rho_,2) * pow(mixture_ptrs_[c]->sound_speed_,2));

         rho_comp_[c] = mixture_ptrs_[c]->rho_;
         chemical_diffusivity_comp_[c] = mixture_ptrs_[c]->chemical_diffusivity_;
       }
       rho_ = 1.0/rho_;
       sound_speed_ = sqrt(1./sound_speed_);
       sound_speed_ /= rho_;       

       for (int c=0; c<nb_comp_; c++)
       {
         beta_ += Y[c]*mixture_ptrs_[c]->beta_*rho_/mixture_ptrs_[c]->rho_;
       }
     }
     
     else  /// This is for mixture of chemical components
     {
       clear();
       for(int c=0; c<nb_comp_; c++)
       {
         chemical_ptrs_[c]->properties(p, isothermal_T_, strain_rate);
         rho_ += Y[c]/chemical_ptrs_[c]->rho_;
         mu_ += Y[c]*chemical_ptrs_[c]->mu_;
         chemical_diffusivity_ += Y[c]*chemical_ptrs_[c]->chemical_diffusivity_;
         sound_speed_ += Y[c]/(pow(chemical_ptrs_[c]->rho_,2) * pow(chemical_ptrs_[c]->sound_speed_,2));

         rho_comp_[c] = chemical_ptrs_[c]->rho_;
         chemical_diffusivity_comp_[c] = chemical_ptrs_[c]->chemical_diffusivity_;
       }
       rho_ = 1./rho_;
       sound_speed_ = sqrt(1.0/sound_speed_);
       sound_speed_ /= rho_;     

       for (int c=0; c<nb_comp_; c++)
       {
         beta_ += Y[c]*chemical_ptrs_[c]->beta_*rho_/chemical_ptrs_[c]->rho_;
       }
     }  
  }







  /// This is called from fluid_mc
  void properties(double p, double T, const vec &Y, double strain_rate)final
  {
     if(mix_of_mix_)  /// This is for mixture of submixtures
     {
       clear();
       for (int c=0; c<nb_comp_; c++)
       {
         mixture_ptrs_[c]->properties(p,T,strain_rate);
         rho_ += Y[c]/mixture_ptrs_[c]->rho_;
         mu_ += Y[c]*mixture_ptrs_[c]->mu_;
         cp_ += Y[c]*mixture_ptrs_[c]->cp_;
         cv_ += Y[c]*mixture_ptrs_[c]->cv_;
         kappa_ += Y[c]*mixture_ptrs_[c]->kappa_;
         sound_speed_ += Y[c]/(pow(mixture_ptrs_[c]->rho_,2) * pow(mixture_ptrs_[c]->sound_speed_,2));

         rho_comp_[c] = mixture_ptrs_[c]->rho_;
         internal_energy_comp_[c] = mixture_ptrs_[c]->cv_*T;
         chemical_diffusivity_comp_[c] = mixture_ptrs_[c]->chemical_diffusivity_;
       }
       rho_ = 1.0/rho_;
       sound_speed_ = sqrt(1./sound_speed_);
       sound_speed_ /= rho_;       

       for (int c=0; c<nb_comp_; c++)
       {
         alpha_ += Y[c]*mixture_ptrs_[c]->alpha_*rho_/mixture_ptrs_[c]->rho_;
         beta_ += Y[c]*mixture_ptrs_[c]->beta_*rho_/mixture_ptrs_[c]->rho_;
       }
     }
     
     else  /// This is for mixture of chemical components
     {
       clear();
       for(int c=0; c<nb_comp_; c++)
       {
         chemical_ptrs_[c]->properties(p,T,strain_rate);
         rho_ += Y[c]/chemical_ptrs_[c]->rho_;
         mu_ += Y[c]*chemical_ptrs_[c]->mu_;
         cp_ += Y[c]*chemical_ptrs_[c]->cp_;
         cv_ += Y[c]*chemical_ptrs_[c]->cv_;
         kappa_ += Y[c]*chemical_ptrs_[c]->kappa_;
         chemical_diffusivity_ += Y[c]*chemical_ptrs_[c]->chemical_diffusivity_;
         sound_speed_ += Y[c]/(pow(chemical_ptrs_[c]->rho_,2) * pow(chemical_ptrs_[c]->sound_speed_,2));
         
         rho_comp_[c] = chemical_ptrs_[c]->rho_;
         internal_energy_comp_[c] = chemical_ptrs_[c]->cv_*T;
         chemical_diffusivity_comp_[c] = chemical_ptrs_[c]->chemical_diffusivity_;
       }
       rho_ = 1./rho_;
       sound_speed_ = sqrt(1.0/sound_speed_);
       sound_speed_ /= rho_;     

       for (int c=0; c<nb_comp_; c++)
       {
         alpha_ += Y[c]*chemical_ptrs_[c]->alpha_*rho_/chemical_ptrs_[c]->rho_;
         beta_ += Y[c]*chemical_ptrs_[c]->beta_*rho_/chemical_ptrs_[c]->rho_;
       }
     }  
  }








   void mixture_of_mix(bool f) {mix_of_mix_ = f;}

   void Sutherland_law(bool f) {Sutherland_law_ = f;}
   void Sutherland_law_a(double a) {Sutherland_law_a_ = a;}
   void Sutherland_law_b(double b) {Sutherland_law_b_ = b;}

   void isothermal_T(double b) 
   {
     isothermal_T_ = b;
     if(isothermal_T_ == 0.0)
     {
       Error("isothermal_T is set to 0.0;  it should not be zero");
     }
   }
   
   int nb_comp() const{return nb_comp_;}
   


   private:
   bool mix_of_mix_ = false;      

   bool Sutherland_law_ = false;
   double Sutherland_law_a_ = 0.0;
   double Sutherland_law_b_ = 0.0;
    
   double isothermal_T_ = 273.0;   
   
};


}

#endif
