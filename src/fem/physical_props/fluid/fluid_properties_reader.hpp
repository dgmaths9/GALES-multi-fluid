#ifndef GALES_FLUID_PROPERTIES_READER_HPP
#define GALES_FLUID_PROPERTIES_READER_HPP



#include "chemicals.hpp"
#include "mixtures.hpp"



namespace GALES {



   /**
        This class reads fluid properties.
        To distinguish between the type of a material we use keywords: mix_of_mix, mix_of_ch, magma_mix, ch_as_mix, ch
        The lines of the file are read by function "read_one_line" defined in algorithm file.
        Whatever is read from the file is printed in outfile file at runtime to check if everything was defined in the right way.
   */






class fluid_properties_reader
{

  using vec = boost::numeric::ublas::vector<double>;
  
  public:

  template<typename T>
  fluid_properties_reader(T& props)
  {  
      std::vector<std::string> s;
      std::ifstream file("props.txt");
      if(!file.is_open())
      {
         Error("unable to open and read  'props.txt' "); 
      }

      while(file)
      {
         read_one_line(file, s);         
         if(s[0] == "fluid")
         {
            read_f_material(file, s, props);
         }
      }
      file.close();                  
  }





  //-------------------------------------------------------------------------------------------------------------------------------------        
  /// Deleting the copy and move constructors - no duplication/transfer in anyway
  fluid_properties_reader(const fluid_properties_reader&) = delete;               //copy constructor
  fluid_properties_reader& operator=(const fluid_properties_reader&) = delete;    //copy assignment operator
  fluid_properties_reader(fluid_properties_reader&&) = delete;                    //move constructor  
  fluid_properties_reader& operator=(fluid_properties_reader&&) = delete;         //move assignment operator 
  //-------------------------------------------------------------------------------------------------------------------------------------




  /**  
       This function reads fluid material properties.  
  */
  template<typename T>
  void read_f_material(std::ifstream& file, std::vector<std::string>& s, T& props)
  {
     while(file)
     {
         read_one_line(file, s); 

         if(s[0] == "{")
         {
              print_only_pid<0>(std::cerr)<<"\n"<<"---------------------------------- fluid material ---------------------------------------"<<"\n";                      
         } 
         else if(s[0] == "Sutherland_law")
         {
            props.Sutherland_law(true);
            props.Sutherland_law_a(stod(s[1])); 
            props.Sutherland_law_b(stod(s[2])); 
         }   
         else if(s[0] == "isothermal_T" || s[0] == "Isothermal_T" || s[0] == "ISOTHERMAL_T")
         {
            props.isothermal_T(stod(s[1])); 
         }   
         else if(s[0] == "mix_of_mix")
         {
            props.mixture_of_mix(true);
            read_mix_of_mix(file, s, props);
         }
         else if(s[0] == "mix_of_ch")
         {
            read_mix_of_ch(file, s, props);            
         }
         else if(s[0] == "ch")
         {
            std::shared_ptr<base_chemical> ch_ptr;
            read_ch(file, s, ch_ptr);            
            props.chemical_ptrs_.push_back(ch_ptr);
         }
         else if(s[0] == "}")
         {
            print_only_pid<0>(std::cerr)<<"----------------------------------------------------------------------------------------"<<"\n";             
            return;
         }
     }    
  }




   
  /**
      This functions reads mixture made of submixtures. submixture is typically either a magma mix, a custom mixture or a chemical considered as a mixture.
  */
  template<typename T>
  void read_mix_of_mix(std::ifstream& file, std::vector<std::string>& s, T& props)
  {           
     while(file)
     {
         read_one_line(file, s);         
         
         if(s[0] == "{")
         {
            print_only_pid<0>(std::cerr)<<"  mix of submix:  "<<"\n\n";                      
         }
         else if(s[0] == "ch_as_mix")
         {
            read_ch_as_mix(file, s, props);                                  
         }
         else if(s[0] == "custom_mix")
         {
            read_custom_mix(file, s, props);                                  
         }
         else if(s[0] == "magma_mix")
         {
            read_magma_mix(file, s, props);
         }
         else if(s[0] == "}")
         {
             const int nb_comp = props.mixture_ptrs_.size();
             props.nb_comp_ = nb_comp;
             props.rho_comp_.resize(nb_comp);
             props.internal_energy_comp_.resize(nb_comp);
             props.chemical_diffusivity_comp_.resize(nb_comp);                  
             return;
         }    
     }            
  }











  /**
         This reads chemical as one mixture component.
         This makes sense if we have to make a whole mixture where one component is already a mixture and another is a chemical.
         In this case it become necessary to project ch as a mixture.
  */
  template<typename T>
  void read_ch_as_mix(std::ifstream& file, std::vector<std::string>& s, T& props)
  {  
     std::vector<std::shared_ptr<base_chemical>> chemical_ptrs;
     while(file)
     {
       read_one_line(file, s);     

       if(s[0] == "{")
       {
          print_only_pid<0>(std::cerr)<<"    ch as mix:  "<<"\n\n";                
       }
       else if(s[0] == "ch")
       {
         std::shared_ptr<base_chemical> ch_ptr;
         read_ch(file, s, ch_ptr);
         chemical_ptrs.push_back(ch_ptr);
       }
       else if(s[0] == "}")
       {
          /**  
             After reading the ch as mixture at "}" we create a shared pointer of mixture type with address of ch_as_mix and push the pointer in a vector "mixture_" defined in mixture.hpp
             Here we are using run time polymorphism in which a pointer to base class can be assigned to the address of the derived class.
             The virtual function will make sure that call goes to the "ch_as_mix" class
          */
          std::shared_ptr<base_mixture> mix_ptr = std::make_shared<ch_as_mix>(chemical_ptrs);
          props.mixture_ptrs_.push_back(mix_ptr);
          return;
       }   
     }
  }













  /**
         This reads a custom mixture as one mixture component.
  */
  template<typename T>
  void read_custom_mix(std::ifstream& file, std::vector<std::string>& s, T& props)
  {       
     double rho(0.0), mu(0.0), cp(0.0), alpha(0.0), beta(0.0), kappa(0.0);
     while(file)
     {
         read_one_line(file, s);          

         if(s[0] == "{")
         {
           print_only_pid<0>(std::cerr)<<"    custom mix:  "<<"\n\n";         
         }
         else if(s[0] == "rho")
         {
           rho = stod(s[1]);
           print_only_pid<0>(std::cerr)<<"            rho: "<<rho<<"\n";
         }               
         else if(s[0] == "mu")
         {
           mu = stod(s[1]);
           print_only_pid<0>(std::cerr)<<"            mu: "<<mu<<"\n";
         }               
         else if(s[0] == "cp")
         {
           cp = stod(s[1]);             
           print_only_pid<0>(std::cerr)<<"            cp: "<<cp<<"\n";
         }  
         else if(s[0] == "kappa")
         {
           kappa = stod(s[1]);   
           print_only_pid<0>(std::cerr)<<"            kappa: "<<kappa<<"\n";
         }  
         else if(s[0] == "alpha")
         {
           alpha = stod(s[1]);             
           print_only_pid<0>(std::cerr)<<"            alpha: "<<alpha<<"\n";
         }  
         else if(s[0] == "beta")
         {
           beta = stod(s[1]);             
           print_only_pid<0>(std::cerr)<<"            beta: "<<beta<<"\n\n";
         }  
         else if(s[0] == "}")
         {    
          /**  
             After reading the custom mixture at "}" we create a shared pointer of mixture type with address of custom_mixture and push the pointer in a vector "mixture_" defined in mixture.hpp
            Here we are using run time polymorphism in which a pointer to base class can be assigned to the address of the derived class.
             The virtual function will make sure that call goes to the "custom_mixture" class
          */
          std::shared_ptr<base_mixture> mix_ptr = std::make_shared<custom_mixture>(rho, mu, cp, kappa, alpha, beta);
          props.mixture_ptrs_.push_back(mix_ptr);
          return;
         }   
     }               
  }













  /**
       This function reads magma mixture
  */
  template<typename T>
  void read_magma_mix(std::ifstream& file, std::vector<std::string>& s, T& props)
  {
     std::string name;
     std::vector<double> oxide_wf(10,0.0);
     double melt_kappa(1.5);
     double h2o_wf(0.0), co2_wf(0.0);
     double h2o_l_cp(2278.0), h2o_l_kappa(0.64);
     double co2_l_cp(2278.0), co2_l_kappa(0.64);
     double crystal_fraction(0.0), crystal_rho(0.0);
     std::string gas_on_mu_model, crystal_on_mu_model;
     
     while(file)
     {
       read_one_line(file, s);     

       if(s[0] == "{")
       {
          print_only_pid<0>(std::cerr)<<"    magma mix:  "<<"\n\n";                    
       }
       else if(s[0] == "name")
       {         
          name = s[1]; 
          print_only_pid<0>(std::cerr)<<"        name: "<<name<<"\n";
       }   
       else if(s[0] == "oxide_wf")
       {  
         for(int i=0; i<10; i++)  
           oxide_wf[i] = stod(s[1+i]);

         std::stringstream ss;
         for(int i=0; i<10; i++)  ss<< oxide_wf[i]<<"  "; 
         print_only_pid<0>(std::cerr)<<"        oxide_wf: "<<ss.str()<<"\n";
       }    
       else if(s[0] == "melt_kappa")   
       {
          melt_kappa = stod(s[1]);
          print_only_pid<0>(std::cerr)<<"        melt_kappa: "<<melt_kappa<<"\n\n";   
       }                 
       else if(s[0] == "h2o_wf")
       {         
          h2o_wf = stod(s[1]); 
          print_only_pid<0>(std::cerr)<<"        h2o_wf: "<<h2o_wf<<"\n";
       }   
       else if(s[0] == "h2o_l_cp")
       {         
          h2o_l_cp = stod(s[1]); 
          print_only_pid<0>(std::cerr)<<"        h2o_l_cp: "<<h2o_l_cp<<"\n";
       }   
       else if(s[0] == "h2o_l_kappa")
       {         
          h2o_l_kappa = stod(s[1]); 
          print_only_pid<0>(std::cerr)<<"        h2o_l_kappa: "<<h2o_l_kappa<<"\n";
       }   
       else if(s[0] == "co2_wf")
       {         
          co2_wf = stod(s[1]); 
          print_only_pid<0>(std::cerr)<<"        co2_wf: "<<co2_wf<<"\n";
       }   
       else if(s[0] == "co2_l_cp")
       {         
          co2_l_cp = stod(s[1]); 
          print_only_pid<0>(std::cerr)<<"        co2_l_cp: "<<co2_l_cp<<"\n";
       }   
       else if(s[0] == "co2_l_kappa")
       {         
          co2_l_kappa = stod(s[1]); 
          print_only_pid<0>(std::cerr)<<"        co2_l_kappa: "<<co2_l_kappa<<"\n";
       }   
       else if(s[0] == "crystal_fraction")
       {         
          crystal_fraction = stod(s[1]); 
          print_only_pid<0>(std::cerr)<<"        crystal_fraction: "<<crystal_fraction<<"\n";
       }   
       else if(s[0] == "crystal_rho")
       { 
          crystal_rho = stod(s[1]); 
          print_only_pid<0>(std::cerr)<<"        crystal_rho: "<<crystal_rho<<"\n";
       }  
       else if(s[0] == "gas_on_mu_model")
       {       
          gas_on_mu_model = s[1];
          print_only_pid<0>(std::cerr)<<"        gas_on_mu_model: "<< gas_on_mu_model<<"\n";         
       }  
       else if(s[0] == "crystal_on_mu_model")   
       {
          crystal_on_mu_model = s[1];
          print_only_pid<0>(std::cerr)<<"        crystal_on_mu_model: "<<crystal_on_mu_model<<"\n\n";   
       }                 
       else if(s[0] == "}")
       {
          /**  
             After reading the magma mixture at "}" we create a shared pointer of mixture type with address of mama mixture and push the pointer in a vector "mixture_" defined in mixture.hpp
             Here we are using run time polymorphism in which a pointer to base class can be assigned to the address of the derived class.
             The virtual function will make sure that call goes to the "magma_mixture" class
          */
          std::shared_ptr<base_mixture> mix_ptr = std::make_shared<magma_mixture>(name, oxide_wf, melt_kappa, h2o_wf, h2o_l_cp, h2o_l_kappa, co2_wf, co2_l_cp, 
                                                                         co2_l_kappa, gas_on_mu_model, crystal_on_mu_model, crystal_fraction, crystal_rho);
          props.mixture_ptrs_.push_back(mix_ptr);
          return;
       }   
     }  
  }













  /**
      This function reads mix of chemicals.
  */
  template<typename T>
  void read_mix_of_ch(std::ifstream& file, std::vector<std::string>& s, T& props)
  {       
     std::vector<std::shared_ptr<base_chemical>> chemical_ptrs;
     while(file)
     {
       read_one_line(file, s);     

       if(s[0] == "{")
       {
         print_only_pid<0>(std::cerr)<<"    mix of ch:  "<<"\n\n";       
       }
       else if(s[0] == "ch")
       {
         std::shared_ptr<base_chemical> ch_ptr;
         read_ch(file, s, ch_ptr);
         chemical_ptrs.push_back(ch_ptr);
       }
       else if(s[0] == "}")
       {
          const int nb_comp = chemical_ptrs.size();
          props.chemical_ptrs_ = chemical_ptrs;
          props.nb_comp_ = nb_comp;
          props.rho_comp_.resize(nb_comp);
          props.internal_energy_comp_.resize(nb_comp);
          props.chemical_diffusivity_comp_.resize(nb_comp);                  
          return;
       }   
     }
  }













  /**
      This function reads ch. After reading from the file a shared pointer of type ch with address of defined type is retuned as output in the function argument "ch_ptr"
      Here again run time polymorphism is used.
  */
  void read_ch(std::ifstream& file, std::vector<std::string>& s, std::shared_ptr<base_chemical>& ch_ptr)
  {       
     while(file)
     {
         read_one_line(file, s); 
         
         if(s[0] == "{")
         {
            // do nothing
         }                          
         else if(s[0] == "ideal_gas" && (s[1] == "air" || s[1] == "Air" || s[1] == "AIR"))                         /// ideal_gas    air
         {
             ch_ptr = std::make_shared<ideal_gas>("AIR");
             print_only_pid<0>(std::cerr)<<"        ch: ideal_gas(AIR)"<<"\n\n";              
         }             
         else if(s[0] == "ideal_gas" && (s[1] == "h2o" || s[1] == "H2o" || s[1] == "H2O"))                   /// ideal_gas    h2o
         {
             ch_ptr = std::make_shared<ideal_gas>("H2O");
             print_only_pid<0>(std::cerr)<<"        ch: ideal_gas(H2O)"<<"\n\n";
         }             
         else if(s[0] == "ideal_gas" && (s[1] == "co2" || s[1] == "Co2" || s[1] == "CO2"))                   /// ideal_gas    co2
         {
             ch_ptr = std::make_shared<ideal_gas>("CO2");
             print_only_pid<0>(std::cerr)<<"        ch: ideal_gas(CO2)"<<"\n\n";
         }    
         else if(s[0] == "ideal_gas" && (s[1] == "custom" || s[1] == "Custom" || s[1] == "CUSTOM"))          /// ideal_gas    custom
         {
             double molar_mass(0.0), gamma(0.0);
             for(int i=0; i<2; i++)
             {
                read_one_line(file, s);
                if(s[0] == "molar_mass")     molar_mass = stod(s[1]);            /// molar_mass    --
                if(s[0] == "gamma")               gamma = stod(s[1]);            ///      gamma    -- 
             }   
             ch_ptr = std::make_shared<ideal_gas>(molar_mass, gamma);
             print_only_pid<0>(std::cerr)<<"        ch: ideal_gas(custom) with molar_mass: "<<molar_mass<<" and gamma: "<<gamma<<"\n\n";
         }    
         else if(s[0] == "custom")                                                                      /// custom material with user defined 
         {                                                                                              /// rho, mu, cp, kappa, alpha, beta
             read_custom_ch(file, s, ch_ptr);                          
             return; 
         }                  
         else if(s[0] == "}")
         {    
             return;
         }  
     }
  }      











  /**
      This function reads custom ch. This is called from "read_ch" function.  
      After reading from the file a shared pointer of type ch with address of defined type is retuned as output in the function argument "ch_ptr"
      Here again run time polymorphism is used.
  */
  void read_custom_ch(std::ifstream& file, std::vector<std::string>& s, std::shared_ptr<base_chemical>& ch_ptr)
  {
     double rho(0.0), mu(0.0), cp(0.0), alpha(0.0), beta(0.0), kappa(0.0);
      
     while(file)
     {
         read_one_line(file, s); 
      
         if(s[0] == "{")
         {
             print_only_pid<0>(std::cerr)<<"        ch: custom"<<"\n";            
         }
         else if(s[0] == "rho")
         {
            rho = stod(s[1]);
            print_only_pid<0>(std::cerr)<<"            rho: "<<rho<<"\n";
         }               
         else if(s[0] == "mu")
         {
            mu = stod(s[1]);
            print_only_pid<0>(std::cerr)<<"            mu: "<<mu<<"\n";
         }               
         else if(s[0] == "cp")
         {
            cp = stod(s[1]);             
            print_only_pid<0>(std::cerr)<<"            cp: "<<cp<<"\n";
         }  
         else if(s[0] == "kappa")
         {
            kappa = stod(s[1]);   
            print_only_pid<0>(std::cerr)<<"            kappa: "<<kappa<<"\n";
         }  
         else if(s[0] == "alpha")
         {
            alpha = stod(s[1]);             
            print_only_pid<0>(std::cerr)<<"            alpha: "<<alpha<<"\n";
         }  
         else if(s[0] == "beta")
         {
            beta = stod(s[1]);             
            print_only_pid<0>(std::cerr)<<"            beta: "<<beta<<"\n\n";
         }  
         else if(s[0] == "}")
         {    
            ch_ptr = std::make_shared<custom_chemical>(rho, mu, cp, kappa, alpha, beta);
            return;
         }   
     }               
  }

     
};


}

#endif
