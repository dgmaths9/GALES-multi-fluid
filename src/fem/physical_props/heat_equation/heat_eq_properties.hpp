#ifndef GALES_HEAT_EQ_PROPERTIES_HPP
#define GALES_HEAT_EQ_PROPERTIES_HPP




namespace GALES {


 /**
      This class defined the physical properties of heat equation.
 */



class heat_eq_properties 
{
  
  public:

  heat_eq_properties()
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
         if(s[0] == "heat_equation")
         {
            read_heat_eq_props(file, s);
         }
      }
      file.close();                  
  }
  
  
  
  
  
  
  
  /**  
       This function reads heat equation properties.  
  */
  void read_heat_eq_props(std::ifstream& file, std::vector<std::string>& s)
  {
     while(file)
     {
         read_one_line(file, s); 
         
         if(s[0] == "{")
         {
              print_only_pid<0>(std::cerr)<<"\n"<<"---------------------------------- heat equation properties ---------------------------------------"<<"\n";                               
         }
         else if(s[0] == "rho")
         {
            rho(stod(s[1]));
            print_only_pid<0>(std::cerr)<<"            rho: "<<rho()<<"\n";
         }                                 
         else if(s[0] == "cp")
         {
            cp(stod(s[1]));
            print_only_pid<0>(std::cerr)<<"            cp: "<<cp()<<"\n";
         }               
         else if(s[0] == "kappa")
         {
            kappa(stod(s[1]));
            print_only_pid<0>(std::cerr)<<"            kappa: "<<kappa()<<"\n";
         }               
         else if(s[0] == "}")
         {             
            print_only_pid<0>(std::cerr)<<"----------------------------------------------------------------------------------------"<<"\n";             
            return;
         }   
     }    
  }
  
  
  
  
  
  //-------------------------------------------------------------------------------------------------------------------------------------        
  /// Deleting the copy and move constructors - no duplication/transfer in anyway
  heat_eq_properties(const heat_eq_properties&) = delete;               //copy constructor
  heat_eq_properties& operator=(const heat_eq_properties&) = delete;    //copy assignment operator
  heat_eq_properties(heat_eq_properties&&) = delete;                    //move constructor  
  heat_eq_properties& operator=(heat_eq_properties&&) = delete;         //move assignment operator 
  //-------------------------------------------------------------------------------------------------------------------------------------



    
  void rho(double x){rho_ = x;}
  double rho() const {return rho_;}

  void cp(double x){cp_ = x;}
  double cp() const {return cp_;}  

  void kappa(double x){kappa_ = x;}
  double kappa() const {return kappa_;}  
  
  private:
      
  double rho_ = 0.0;
  double cp_ = 0.0;
  double kappa_ = 0.0;
   
};


}

#endif
