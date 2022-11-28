#ifndef GALES_ADV_DIFF_PROPERTIES_HPP
#define GALES_ADV_DIFF_PROPERTIES_HPP




namespace GALES {


 /**
      This class defined the physical properties of ADV_DIFF equation.
 */



class adv_diff_properties 
{
  
  public:

  adv_diff_properties()
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
         if(s[0] == "adv_diff")
         {
            read_adv_diff_props(file, s);
         }
      }
      file.close();                  
  }
  
  
  
  
  
  /**  
       This function reads adv_diff properties.  
  */
  void read_adv_diff_props(std::ifstream& file, std::vector<std::string>& s)
  {
     while(file)
     {
         read_one_line(file, s); 
         
         if(s[0] == "{")
         {
              print_only_pid<0>(std::cerr)<<"\n"<<"---------------------------------- adv diff properties ---------------------------------------"<<"\n";                               
         }
         else if(s[0] == "rho")
         {
            rho(stod(s[1]));
            print_only_pid<0>(std::cerr)<<"            rho: "<<rho()<<"\n";
         }                                 
         else if(s[0] == "mu")
         {
            mu(stod(s[1]));
            print_only_pid<0>(std::cerr)<<"            mu: "<<mu()<<"\n";
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
  adv_diff_properties(const adv_diff_properties&) = delete;               //copy constructor
  adv_diff_properties& operator=(const adv_diff_properties&) = delete;    //copy assignment operator
  adv_diff_properties(adv_diff_properties&&) = delete;                    //move constructor  
  adv_diff_properties& operator=(adv_diff_properties&&) = delete;         //move assignment operator 
  //-------------------------------------------------------------------------------------------------------------------------------------

    
  void rho(double x){rho_ = x;}
  double rho() const {return rho_;}

  void mu(double x){mu_ = x;}
  double mu() const {return mu_;}  
  
  private:
      
  double rho_ = 0.0;
  double mu_ = 0.0;
   
};


}

#endif
