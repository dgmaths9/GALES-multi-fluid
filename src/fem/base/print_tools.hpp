#ifndef GALES_PRINT_TOOLS_HPP
#define GALES_PRINT_TOOLS_HPP


#include <iostream>
#include <iomanip>
#include <vector>
#include <iterator>
#include <sstream>
#include <fstream>
#include <mpi.h>
#include <initializer_list>



namespace GALES{






   ///----------------------------------------------------------------------------------------- 
   /// This is to simply place an empty line.
   void print()
   {
      std::cerr<<std::endl;
   }
   ///----------------------------------------------------------------------------------------- 






   ///----------------------------------------------------------------------------------------- 
   /// This is to print data value.
   template<typename data_type>
   void print(const data_type& data)
   {
      std::cerr<<data<<std::endl;
   }
   ///----------------------------------------------------------------------------------------- 






   ///----------------------------------------------------------------------------------------- 
   /// This is to print data name and its value.
   template<typename data_type>
   void print(std::string varname, const data_type& data)
   {
      std::cerr<<varname<<" = "<< data<<std::endl;
   }
   ///----------------------------------------------------------------------------------------- 





   ///----------------------------------------------------------------------------------------- 
   /// This is to print std vector. 
   template<typename T>
   void print(const std::vector<T>& v)
   {
     std::stringstream ss;
     copy(v.begin(), v.end(), std::ostream_iterator<T>(ss, "  "));
     print(ss.str());
   }
   ///----------------------------------------------------------------------------------------- 






   ///----------------------------------------------------------------------------------------- 
   /// This is to print std vector with std::cout or std::cerr. 
    template<typename E, typename T, typename data_type>
    std::basic_ostream<E, T> &operator << (std::basic_ostream<E,T>& os, const std::vector<data_type>& v) 
    {
        std::size_t size = v.size();        
        std::basic_ostringstream<E, T, std::allocator<E>> s;        
        s.flags(os.flags());
        s.imbue(os.getloc());
        s.precision(os.precision());
        s << "[" << size << "](";
        if(size > 0) s << v[0];
        for(std::size_t i=1; i<size; i++)  s << ", " << v[i];
        s << ")";
        return os << s.str().c_str();
    }
   ///----------------------------------------------------------------------------------------- 





   ///----------------------------------------------------------------------------------------- 
   /// This is to print data value on process 0 only.
   template<typename data_type>
   void print0(const data_type& data)
   {
      const int rank = get_rank();
      if(rank==0)
        print(data);
   }
   ///----------------------------------------------------------------------------------------- 






   ///----------------------------------------------------------------------------------------- 
   /// This is to print data name and value on process 0 only.
   template<typename data_type>
   void print0(std::string varname, const data_type& data)
   {
      const int rank = get_rank();
      if(rank==0)
        print(varname, data);
   }
   ///----------------------------------------------------------------------------------------- 






   ///----------------------------------------------------------------------------------------- 
   /// This is to print std vector on process 0 only. 
   template<typename T>
   void print0(const std::vector<T>& v)
   {
      const int rank = get_rank();
      if(rank==0)
      {
        std::stringstream ss;
        copy(v.begin(), v.end(), std::ostream_iterator<T>(ss, "  "));
        print(ss.str());
      }  
   }
   ///----------------------------------------------------------------------------------------- 

   
   
   
   
   
   
   
   ///----------------------------------------------------------------------------------------- 
   /// This is to for debugging purpose. 
   template<typename data_type>
   void debug(const data_type& d)
   {
      const int rank = get_rank();
      const int size = get_size();
      for(int pid=0; pid<size; pid++)      
        if(rank==pid)
        {      
          std::stringstream ss;
          ss <<"debug_time_"<<time::get().t()<<"_pid_"<<pid<<std::flush;
          std::ofstream os(ss.str(), std::ios::app);
          os<< d <<std::endl;
          os.close();
        }  
   }
   ///----------------------------------------------------------------------------------------- 
 

   
   
   
   
   ///----------------------------------------------------------------------------------------- 
   /// This is to print data value in file on specified process. 
   template<typename data_type>
   void print_in_file(int pid, std::string f, const data_type& d)
   {
      const int rank = get_rank();
        if(rank==pid)
        {      
          std::stringstream ss;
          ss << f <<"_time_"<<time::get().t()<<"_pid_"<<pid<<std::flush;
          std::ofstream os(ss.str(), std::ios::app);
          os<< d <<std::endl;
          os.close();
        }  
   }
   ///----------------------------------------------------------------------------------------- 








   ///----------------------------------------------------------------------------------------- 
   /// This is to print std vector in file on specified process. 
   template<typename T>
   void print_in_file(int pid, std::string f, const std::vector<T>& v)
   {
      const int rank = get_rank();
      if(rank==pid)
      {      
          std::stringstream ss;
          ss << f <<"_time_"<<time::get().t()<<"_pid_"<<pid<<std::flush;
          std::ofstream os(ss.str(), std::ios::app);
          std::stringstream ss1;
          copy(v.begin(), v.end(), std::ostream_iterator<T>(ss1, "  "));
          os<< ss1.str() <<std::endl;
          os.close();
      }  
   }
   ///----------------------------------------------------------------------------------------- 







   ///----------------------------------------------------------------------------------------- 
   /// This is to print data name and value in file on specified process. 
   template<typename data_type>
   void print_in_file(int pid, std::string f, std::string d_name, const data_type& d)
   {
      const int rank = get_rank();
        if(rank==pid)
        {      
          std::stringstream ss;
          ss << f <<"_time_"<<time::get().t()<<"_pid_"<<pid<<std::flush;
          std::ofstream os(ss.str(), std::ios::app);
          os<< d_name <<" = "<< d <<std::endl;
          os.close();
        }  
   }
   ///----------------------------------------------------------------------------------------- 








   ///----------------------------------------------------------------------------------------- 
   /// This is to print data value in file on specified process for non linear iteration(it_num). 
   template<typename data_type>
   void print_in_file(int pid, std::string f, int it_num, const data_type& d)
   {
      const int rank = get_rank();
        if(rank==pid)
        {      
          std::stringstream ss;
          ss << f <<"_time_"<<time::get().t()<<"_pid_"<<pid<<"_it_"<<it_num<<std::flush;
          std::ofstream os(ss.str(), std::ios::app);
          os<< d <<std::endl;
          os.close();
        }  
   }
   ///----------------------------------------------------------------------------------------- 








   ///----------------------------------------------------------------------------------------- 
   /// This is to print data name and value in file on specified process for non linear iteration(it_num). 
   template<typename data_type>
   void print_in_file(int pid, std::string f, int it_num, std::string d_name, const data_type& d)
   {
      const int rank = get_rank();
        if(rank==pid)
        {      
          std::stringstream ss;
          ss << f<<"_time_"<<time::get().t()<<"_pid_"<<pid<<"_it_"<<it_num<<std::flush;
          std::ofstream os(ss.str(), std::ios::app);
          os<< d_name <<" = "<< d <<std::endl;
          os.close();
        }  
   }
   ///----------------------------------------------------------------------------------------- 









   ///----------------------------------------------------------------------------------------- 
   /// This is to print data value in file on all processes.
   template<typename data_type>
   void print_in_file(std::string f, const data_type& d)
   {
      const int size = get_size();
      for(int pid=0; pid<size; pid++)
          print_in_file(pid, f, d);
   }
   ///----------------------------------------------------------------------------------------- 







   ///----------------------------------------------------------------------------------------- 
   /// This is to print data name and value in file on all processes.
   template<typename data_type>
   void print_in_file(std::string f, std::string d_name, const data_type& d)
   {
      const int size = get_size();
      for(int pid=0; pid<size; pid++)
          print_in_file(pid, f, d_name, d);
   }
   ///----------------------------------------------------------------------------------------- 






   ///----------------------------------------------------------------------------------------- 
   /// This is to print std vector in file on all processes. 
   template<typename T>
   void print_in_file(std::string f, const std::vector<T>& v)
   {
      const int size = get_size();
      for(int pid=0; pid<size; pid++)
          print_in_file(pid, f, v);
   }
   ///----------------------------------------------------------------------------------------- 








   ///----------------------------------------------------------------------------------------- 
   /// This is to print data value in file on all processes for non linear iteration(it_num).
   template<typename data_type>
   void print_in_file(std::string f, int it_num, const data_type& d)
   {
      const int size = get_size();
      for(int pid=0; pid<size; pid++)
          print_in_file(pid, f, it_num, d);
   }
   ///----------------------------------------------------------------------------------------- 








   ///----------------------------------------------------------------------------------------- 
   /// This is to print data name and value in file on all processes for non linear iteration(it_num).
   template<typename data_type>
   void print_in_file(std::string f, int it_num, std::string d_name, const data_type& d)
   {
      const int size = get_size();
      for(int pid=0; pid<size; pid++)
          print_in_file(pid, f, it_num, d_name, d);
   }
   ///----------------------------------------------------------------------------------------- 








   ///----------------------------------------------------------------------------------------- 
  /**
    This class uses the << operator on a given object only for the process specified by the template, while the others do nothing.
    It is not meant to be used directly but by the means of ad hoc operator <<.
  */  
  template<int pid_id=0>
  struct print_only_pid
  {
    print_only_pid(std::ostream& o):stream_m(o){}

    template<typename s_type>
    void operator()(const s_type& s)
    {
      const int rank = get_rank();      
      if (rank==pid_id)
      {
        stream_m.flush();
        stream_m<<s;
        stream_m.flush();
      }
    }
    std::ostream& stream_m;
  };



  /// The operator to be actually used by the user. For example: print_only_pid<1>(std::cerr)<<"error!!! \n";
  template<typename d_type, int PID_ID>
  inline print_only_pid<PID_ID> operator<<(print_only_pid<PID_ID> p, d_type const & s)
  {
    p(s);
    return p;
  }
   ///----------------------------------------------------------------------------------------- 



}   
   
   
#endif   
