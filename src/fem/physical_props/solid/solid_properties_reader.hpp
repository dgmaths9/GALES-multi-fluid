#ifndef GALES_SOLID_PROPERTIES_READER_HPP
#define GALES_SOLID_PROPERTIES_READER_HPP





namespace GALES {



   /**
        This class reads solid properties.
        The lines of the file are read by function "read_one_line" defined in algorithm file.
        Whatever is read from the file is printed in outfile file at runtime to check if everything was defined in the right way.
   */






class solid_properties_reader
{

  public:

  template<typename T>
  solid_properties_reader(T& props)
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
         if(s[0] == "solid")
         {
            read_s_material(file, s, props);
         }
      }
      file.close();                  
  }




  //-------------------------------------------------------------------------------------------------------------------------------------        
  /// Deleting the copy and move constructors - no duplication/transfer in anyway
  solid_properties_reader(const solid_properties_reader&) = delete;               //copy constructor
  solid_properties_reader& operator=(const solid_properties_reader&) = delete;    //copy assignment operator
  solid_properties_reader(solid_properties_reader&&) = delete;                    //move constructor  
  solid_properties_reader& operator=(solid_properties_reader&&) = delete;         //move assignment operator 
  //-------------------------------------------------------------------------------------------------------------------------------------




  /**  
       This function reads solid material properties.  
  */
  template<typename T>
  void read_s_material(std::ifstream& file, std::vector<std::string>& s, T& props)
  {
     while(file)
     {
         read_one_line(file, s); 
         
         if(s[0] == "{")
         {
              print_only_pid<0>(std::cerr)<<"\n"<<"---------------------------------- solid material ---------------------------------------"<<"\n";                               
         }
         else if(s[0] == "material" || s[0] == "Material" || s[0] == "MATERIAL")
         {
           if(s[1]=="hookes" || s[1]=="Hookes" || s[1]=="HOOKES")
           {           
             props.material_type("Hookes");
             print_only_pid<0>(std::cerr)<<"        material: Hookes"<<"\n";
           }  
           else if(s[1]=="svk" || s[1]=="Svk" || s[1]=="SVK")
           {  
             props.material_type("St_V_Kirchhoff");
             print_only_pid<0>(std::cerr)<<"        material: Saint Venant Kirchhoff"<<"\n";
           }  
           else if(s[1]=="neo_hookean" || s[1]=="Neo_Hookean" || s[1]=="NEO_HOOKEAN" || s[1]=="Neo Hookean"|| s[1]=="Neo hookean")
           {  
             props.material_type("Neo_Hookean");
             print_only_pid<0>(std::cerr)<<"        material: Neo Hookean"<<"\n";
           }             
         }              
         else if(s[0] == "plane_strain" || s[0] == "Plane_Strain" || s[0] == "Plane_strain" || s[0] == "PLAIN_STRAIN" || s[0] == "PLAIN STRAIN" || s[0] == "plain strain")
         {
           if(s[1]=="T" || s[1]=="t" || s[1]=="true"|| s[1]=="True" || s[1]=="TRUE")
           {           
             props.condition_type("Plain_Strain");
             print_only_pid<0>(std::cerr)<<"        condition: Plain_Strain"<<"\n";
           } 
         }   
         else if(s[0] == "plane_stress" || s[0] == "Plane_Stress" || s[0] == "Plane_stress" || s[0] == "PLAIN_STRESS" || s[0] == "PLAIN STRESS" || s[0] == "plain stress")
         {
           if(s[1]=="T" || s[1]=="t" || s[1]=="true"|| s[1]=="True" || s[1]=="TRUE")
           {           
             props.condition_type("Plain_Stress");          
             print_only_pid<0>(std::cerr)<<"        condition: Plain_Stress"<<"\n";
           }  
         }              
         else if(s[0] == "axisymmetric" || s[0] == "AxiSymmetric" || s[0] == "Axisymmetric" || s[0] == "AXISYMMETRIC")
         {
           if(s[1]=="T" || s[1]=="t" || s[1]=="true"|| s[1]=="True" || s[1]=="TRUE")
           {           
             props.condition_type("Axisymmetric");          
             print_only_pid<0>(std::cerr)<<"        condition: Axisymmetric"<<"\n";
           }  
         }              
         else if(s[0] == "rho")
         {
            props.rho(stod(s[1]));
            print_only_pid<0>(std::cerr)<<"            rho: "<<props.rho()<<"\n";
         }                                 
         else if(s[0] == "E")
         {
            props.E(stod(s[1]));
            print_only_pid<0>(std::cerr)<<"            E: "<<props.E()<<"\n";
         }               
         else if(s[0] == "nu")
         {
            props.nu(stod(s[1]));
            print_only_pid<0>(std::cerr)<<"            nu: "<<props.nu()<<"\n";
         }               
         else if(s[0] == "a_damping")
         {
            props.a_damping(stod(s[1]));
            print_only_pid<0>(std::cerr)<<"            a_damping: "<<props.a_damping()<<"\n";
         }               
         else if(s[0] == "b_damping")
         {
            props.b_damping(stod(s[1]));
            print_only_pid<0>(std::cerr)<<"            b_damping: "<<props.b_damping()<<"\n";
         }               
         else if(s[0] == "nb_Maxwell_el")
         {
            props.nb_Maxwell_el(stoi(s[1]));
            print_only_pid<0>(std::cerr)<<"            nb_Maxwell_el: "<<props.nb_Maxwell_el()<<"\n";
         }               
         else if(s[0] == "Max_el_E")
         {
            const int nb_el = props.nb_Maxwell_el();
            std::vector<double> v(nb_el,0.0);
            for(int i=0; i<nb_el; i++)  v.push_back(stod(s[1+i]));
            props.Max_el_E(v);
            
            std::stringstream ss;
            for(int i=0; i<nb_el; i++)  ss<< (props.Max_el_E())[i]<<"  "; 
            print_only_pid<0>(std::cerr)<<"         Max_el_E: "<<ss.str()<<"\n";
         }               
         else if(s[0] == "Max_el_eta")
         {
            const int nb_el = props.nb_Maxwell_el();
            std::vector<double> v(nb_el,0.0);
            for(int i=0; i<nb_el; i++)  v.push_back(stod(s[1+i]));
            props.Max_el_eta(v);

            std::stringstream ss;
            for(int i=0; i<nb_el; i++)  ss<< (props.Max_el_eta())[i]<<"  "; 
            print_only_pid<0>(std::cerr)<<"         Max_el_eta: "<<ss.str()<<"\n";
         }               
         else if(s[0] == "heterogeneous")
         {
            props.heterogeneous(true);
            read_layers(file, s, props);
            props.num_layers(props.layers_rho().size());
         }
         else if(s[0] == "}")
         {             
            print_only_pid<0>(std::cerr)<<"----------------------------------------------------------------------------------------"<<"\n";             
            return;
         }   
     }    
  }




  


  template<typename T>
  void read_layers(std::ifstream& file, std::vector<std::string>& s, T& props)
  {
     print_only_pid<0>(std::cerr)<<"\n"<<"         Heterogeneous layers: "<<"\n\n";

     while(file)
     {
       read_one_line(file, s);
       
       if(s[0] == "layer")
       {
          read_layer(file, s, props);
       } 
       else if(s[0] == "}")
       {
          return;
       }   
     }
  }
  
  
  
  
  template<typename T>
  void read_layer(std::ifstream& file, std::vector<std::string>& s, T& props)
  {
     print_only_pid<0>(std::cerr)<<"         layer: "<<"\n";
     while(file)
     {
       read_one_line(file, s); 
           
       if(s[0]=="y_min")
       {
         double y_min = stod(s[1]);
         props.layers_y_min(y_min);       
         print_only_pid<0>(std::cerr)<<"         y_min: "<< y_min<<"\n";
       }
       else if(s[0]=="y_max")
       {
         double y_max = stod(s[1]);
         props.layers_y_max(y_max);       
         print_only_pid<0>(std::cerr)<<"         y_max: "<< y_max<<"\n";
       }
       else if(s[0]=="rho")
       {
         double rho = stod(s[1]);
         props.layers_rho(rho);       
         print_only_pid<0>(std::cerr)<<"         rho: "<< rho<<"\n";
       }
       else if(s[0]=="E")
       {
         double E = stod(s[1]);
         props.layers_E(E);       
         print_only_pid<0>(std::cerr)<<"         E: "<< E<<"\n";
       }
       else if(s[0]=="nu")
       {
         double nu = stod(s[1]);
         props.layers_nu(nu);       
         print_only_pid<0>(std::cerr)<<"         nu: "<< nu<<"\n\n";
       }
       else if(s[0] == "}")
       {    
         return;
       }
     }     
  }
  
  
     
};


}

#endif
