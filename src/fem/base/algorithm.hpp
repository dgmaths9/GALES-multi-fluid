#ifndef GALES_ALGORITHM_HPP
#define GALES_ALGORITHM_HPP


#include <vector>
#include <algorithm>
#include <iterator>
#include <dirent.h>
#include <fstream>
#include <string>
#include <map>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/algorithm/string.hpp>
#include <type_traits>
#include <typeinfo>
#include <mpi.h>



namespace GALES{





  ///----------------------------------------------------------------------------------------- 
  // Error(" err_msg") can now be used to produce a error message that includes the file and line in which the error occurrred
  void error(const char* file, int line, const std::string& message)
  {
    print_only_pid<0>(std::cerr)<<"\n" << file << ":" << line << ": error: "  << message <<"\n\n";
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  #define Error( message ) error( __FILE__, __LINE__, message );
  ///----------------------------------------------------------------------------------------- 









  ///----------------------------------------------------------------------------------------- 
  /// Sort the vector. 
  template<typename T>
  void sort(std::vector<T>& v)
  {
      sort(v.begin(), v.end());
  }
  ///----------------------------------------------------------------------------------------- 
  







  ///----------------------------------------------------------------------------------------- 
  /// Sort the vector. 
  template<typename T>
  auto sort_vector(const std::vector<T>& v)
  {
      auto v1 = v;
      sort(v1.begin(), v1.end());
      return v1;
  }
  ///----------------------------------------------------------------------------------------- 







  ///----------------------------------------------------------------------------------------- 
  /// unique the vector. 
  template<typename T>
  void unique(std::vector<T>& v)
  {
      auto last = std::unique(v.begin(), v.end()) ;
      v.erase(last, v.end());             
  }
  ///----------------------------------------------------------------------------------------- 








  ///----------------------------------------------------------------------------------------- 
  /// unique the vector. 
  template<typename T>
  auto unique_vector(const std::vector<T>& v)
  {
      auto v1 = v;
      auto last = std::unique(v1.begin(), v1.end()) ;
      v1.erase(last, v1.end());             
      return v1;
  }
  ///----------------------------------------------------------------------------------------- 






  ///----------------------------------------------------------------------------------------- 
  /// Sort and unique the vector. 
  template<typename T>
  void sort_unique(std::vector<T>& v)
  {
      sort(v);
      unique(v);            
  }
  ///----------------------------------------------------------------------------------------- 








  ///----------------------------------------------------------------------------------------- 
  /// Sort and unique the vector. 
  template<typename T>
  auto sort_unique_vector(const std::vector<T>& v)
  {
      auto v1 = v;
      sort(v1);
      unique(v1);            
      return v1;
  }
  ///----------------------------------------------------------------------------------------- 











  ///----------------------------------------------------------------------------------------- 
   /// This is to format data; it returns a string with data alligned according to the input width.
   template<typename data_type>
   auto parse(const data_type& d, int width)
   {
     std::stringstream ss;
     ss<< std::left;
     ss.fill(' ');      
     ss.width(width); 
     ss << d;    
     return ss.str();
   }
  ///----------------------------------------------------------------------------------------- 








  ///----------------------------------------------------------------------------------------- 
  /// check if vecs: v1 and v2 are equal irrespective of the order of the numbers
  template<typename T>
  bool are_equal(std::vector<T>& v1, std::vector<T>& v2)
  {
      sort(v1.begin(), v1.end());
      sort(v2.begin(), v2.end());
      if(v1==v2) return true;
      else return false;
  }
  ///----------------------------------------------------------------------------------------- 








  ///----------------------------------------------------------------------------------------- 
  /// return diagonal matrix of a matrix
  template<typename mat>
  mat diag_mat_of_mat(const mat& m)
  {
     mat m1(m.size1(),m.size1(),0.0);
     for(int i=0; i<m.size1(); i++) 
       m1(i,i) = m(i,i);
     return m1;  
  }
  ///----------------------------------------------------------------------------------------- 







  ///----------------------------------------------------------------------------------------- 
  /// This function reads only one line from the file (input) and return a vector of string as output
  void read_one_line(std::ifstream& file, std::vector<std::string>& split_result)
  {
      std::string s;
      getline(file,s);
      boost::trim(s);
      boost::split(split_result, s, boost::is_any_of(", "), boost::token_compress_on);
  }
  ///----------------------------------------------------------------------------------------- 







  ///----------------------------------------------------------------------------------------- 
  /// This function reads only one line from the file (input) to skip
  void skip_one_line(std::ifstream& file)
  {
      std::string s;
      getline(file,s);
  }
  ///----------------------------------------------------------------------------------------- 







  ///----------------------------------------------------------------------------------------- 
  /// This function returns a vector upto num-1;   e.g.  range(3) = {0,1,2} like python style
  auto range(int num)
  {
    std::vector<int> v(num);
    for(int i=0; i<num; i++)
       v[i] = i;
    return v;
  }
  ///----------------------------------------------------------------------------------------- 








  ///----------------------------------------------------------------------------------------- 
  /// This function generate directories
  void make_dirs(std::string s)
  {
    std::string str = "mkdir -p " + s;
    int p = system(str.c_str());
  }  
  ///----------------------------------------------------------------------------------------- 







  ///----------------------------------------------------------------------------------------- 
  /// This function returns a vector of all files contained in a directory
  std::vector<std::string> dir_files(const char *path) 
  {  
    std::vector<std::string> files; 

    DIR *dir = opendir(path);   
    if (dir == NULL) 
    {
       Error("Directory is empty");
    }
    struct dirent *entry;
    while ((entry = readdir(dir)) != NULL) 
    {
        files.push_back(entry->d_name);
    }
    closedir(dir);


    std::vector<std::string> Files; 
    for(auto s: files)
    {
      if(s == "." || s == ".."){}
      else Files.push_back(s);
    }
    
    return Files;        
  }
  ///----------------------------------------------------------------------------------------- 









  ///----------------------------------------------------------------------------------------- 
  /// This function writes vector "v" in a file in binary format
  template<typename T>
  void write_bin(std::string file, const std::vector<T>& v)
  {
    std::ofstream f(file, std::ios::out | std::ios::binary);
    if(!f.is_open())
    {
         Error("file is not opened");
    }       
    f.write(reinterpret_cast<const char*>(&v[0]), v.size()*sizeof(T));
    f.close();
  }
  ///----------------------------------------------------------------------------------------- 







  ///----------------------------------------------------------------------------------------- 
  /// This function reads vector "v" from a binary file
  template<typename T>
  void read_bin(std::string file, std::vector<T>& v)
  {
    std::ifstream f(file, std::ios::in | std::ios::binary);
    if(!f.is_open())
    {
         Error("file is not opened");
    }        
    
    f.seekg(0, std::ios::end);      // go to end of file
    int length = f.tellg();    // get length in bytes of file
    f.seekg (0, std::ios::beg);     // go to beginning of file
    
    int size_of_v = length/sizeof(T);        
    v.resize(size_of_v);
       
    f.read(reinterpret_cast<char*>(&v[0]), v.size()*sizeof(T));
    f.close();
  }
  ///----------------------------------------------------------------------------------------- 








  ///----------------------------------------------------------------------------------------- 
  /// This function (type bool) computes inverse of a square matrix using LU-decomposition with backsubstitution of unit vectors
  template<typename mat>
  bool inverse(const mat& from, mat& to)
  {
    using namespace boost::numeric::ublas;
    typedef permutation_matrix<std::size_t> pmatrix;
    mat A(from);                          // create a working copy of the input
    pmatrix pm(A.size1());                // create a permutation matrix for the LU-factorization
    const int res = lu_factorize(A,pm);   // perform LU-factorization
    if(res != 0) 
    {
       Error("LU factorization fails in inverse function");   
    }

    to.assign(identity_matrix<double>(A.size1()));        // create identity matrix of "inverse"
    lu_substitute(A, pm, to);                             // backsubstitute to get the inverse
    return true;
  }
  ///----------------------------------------------------------------------------------------- 






  ///----------------------------------------------------------------------------------------- 
  /// This function (type mat) computes inverse of a square matrix using LU-decomposition with backsubstitution of unit vectors
  template<typename mat>
  mat inverse(const mat& from)
  {
    using namespace boost::numeric::ublas;
    typedef permutation_matrix<std::size_t> pmatrix;
    mat A(from);                           // create a working copy of the input
    pmatrix pm(A.size1());                 // create a permutation matrix for the LU-factorization
    const int res = lu_factorize(A, pm);   // perform LU-factorization
    if(res != 0)
    {
       Error("LU factorization fails in inverse function");   
    }        
    mat to(identity_matrix<double>(A.size1()));    // create identity matrix of "inverse"
    lu_substitute(A, pm, to);                      // backsubstitute to get the inverse
    return to;
  }
  ///----------------------------------------------------------------------------------------- 







  ///----------------------------------------------------------------------------------------- 
  /// Denman_Beavers_method to compute sqrt inverse of a matrix
  template<typename mat, typename vec>
  mat mat_sqrt_inv(const mat& A)
  {
    boost::numeric::ublas::identity_matrix<double> I(A.size1()); 
    mat P = A;
    mat Q = I;
   
    const double err(1.0);
    const double tol(1.e-2);
   
    mat Q_inv(A.size1(), A.size1(), 0.0);
    mat P_old_inv(A.size1(), A.size1(),0.0);
      
    int it=1;
    while(err > tol)
    {
        mat P_old = P;
        inverse(Q, Q_inv);
        inverse(P_old, P_old_inv);            
        
        P = 0.5*(P_old + Q_inv);        
        Q = 0.5*(Q + P_old_inv);
        
        mat err_mat = P - P_old;                        
        err = 0.0;
        
        for(std::size_t i=0; i<A.size1(); i++)
        {
           vec r = row(err_mat, i);
           double temp_err = norm_inf(r);
           if(temp_err > err)
              err = temp_err;         
        }
        it++;
    }   
   return Q;    //Q = sqrt_inv(A);    P = sqrt(A);
  }  
  ///----------------------------------------------------------------------------------------- 



















  ///----------------------------------------------------------------------------------------- 
   /// This function reads the matching gids of the nodes at fluid-solid interface from input/fsi_nd.txt and returns the corresponding gids in form of map 
   auto read_fsi_matching_nd_gid(std::string file)
   {
     std::ifstream if_stream(file);
     if(!if_stream.is_open())
     {
       Error("unable to open and read   'input/fsi_nd.txt' ");    
     }     

     int rows;
     if_stream >> rows;
    
     int a,b;
     std::map<int,int> matching_nd_gid;    
     for(int i=0; i<rows; i++)
     {
       if_stream >> a >> b;
       matching_nd_gid[a] = b;
     }       
     if_stream.close();       
     return matching_nd_gid;
   }
  ///----------------------------------------------------------------------------------------- 
  









  ///----------------------------------------------------------------------------------------- 
  /// Here we average the v1 and v2 over the repeated nodes (2D version)
   void avg_over_repeated_nodes(
     const std::vector<int>& repeated_nds, const std::vector<double>& v1, const std::vector<double>& v2,
     std::vector<int>& unique_nds, std::vector<double>& avg_v1, std::vector<double>& avg_v2
     )
   { 
        avg_v1.resize(unique_nds.size());
        avg_v2.resize(unique_nds.size());
   
        for(int i=0; i<unique_nds.size(); i++)
        {
          double d1(0.0), d2(0.0), k(0.0);          
          for(int j=0; j<repeated_nds.size(); j++)
            if(unique_nds[i] == repeated_nds[j])
            {
                d1 += v1[j];
                d2 += v2[j];
                k++;
            }
          
          avg_v1[i] = d1/k; 
          avg_v2[i] = d2/k; 
       }   
  }
  ///----------------------------------------------------------------------------------------- 
  








  ///----------------------------------------------------------------------------------------- 
  /// Here we average v1, v2 and v3 over the repeated nodes (3D version)
   void avg_over_repeated_nodes(
     const std::vector<int>& repeated_nds, const std::vector<double>& v1, const std::vector<double>& v2, const std::vector<double>& v3,
     std::vector<int>& unique_nds, std::vector<double>& avg_v1, std::vector<double>& avg_v2, std::vector<double>& avg_v3
     )
   { 
        avg_v1.resize(unique_nds.size());
        avg_v2.resize(unique_nds.size());
        avg_v3.resize(unique_nds.size());
   
        for(int i=0; i<unique_nds.size(); i++)
        {
          double d1(0.0), d2(0.0), d3(0.0), k(0.0);          
          for(int j=0; j<repeated_nds.size(); j++)
            if(unique_nds[i] == repeated_nds[j])
            {
                d1 += v1[j];
                d2 += v2[j];
                d3 += v3[j];
                k++;
            }
          
          avg_v1[i] = d1/k; 
          avg_v2[i] = d2/k; 
          avg_v3[i] = d3/k; 
       }   
  }
  ///----------------------------------------------------------------------------------------- 








  ///----------------------------------------------------------------------------------------- 
  /// This function computes a map (nd_var_map) from my_nd, my_v1 and my_v2 (2D version)
  /// this is called from   sc_f_tr<2>,   sc_isothermal_f_tr<2>, mc_f_tr<2>,  mc_isothermal_f_tr<2> and mc_f_heat_flux<2>
   void all_pid_avg_over_repeated_nodes(const std::vector<int>& my_nd, const std::vector<double>& my_v1, const std::vector<double>& my_v2, std::map<int, std::vector<double>>& nd_var_map)
   {   
         /**  example demonstration for tractions on two processes
         
                            P0                                                   P1
               4       5  5       6   6       7                    7       8   8       9  
               o-------o  o-------o   o-------o                    o-------o   o-------o   
              nd1    nd2  nd1    nd2  nd1    nd2                   nd1    nd2  nd1    nd2 
        
              my_nd = [                                         my_nd = [
                         4,5,                                             7,8, 
                         5,6,                                             8,9 
                         6,7                                            ]
                      ]                          

              my_v1 = [                                         my_v1 = [
                     gn1_tx_4, gn2_tx_5,                                 gn1_tx_7, gn2_tx_8,
                     gn1_tx_5, gn2_tx_6,                                 gn1_tx_8, gn2_tx_9,
                     gn1_tx_6, gn2_tx_7,                                ]
                    ]
              
              my_v2 = [                                         my_v2 = [
                     gn1_ty_4, gn2_ty_5,                                 gn1_ty_7, gn2_ty_8,
                     gn1_ty_5, gn2_ty_6,                                 gn1_ty_8, gn2_ty_9,
                     gn1_ty_6, gn2_ty_7,                                ]
                    ]
         */


        
        auto my_unique_nd = unique_vector(my_nd);           /// my_nd has repeated enteries; we make it unique
        std::vector<double> my_avg_v1, my_avg_v2;        

        /// Here we average my_v1 and my_v2 over the repeated nodes
        avg_over_repeated_nodes(my_nd, my_v1, my_v2, my_unique_nd, my_avg_v1, my_avg_v2);
        
        /*
                    P0                                                P1
            my_unique_nd = [4,5,6,7]                             my_unique_nd = [7,8,9]
            
            my_avg_v1 = [                                           my_avg_v1 = [
                            gn1_tx_4                                              gn1_tx_7
                           (gn2_tx_5 + gn1_tx_5)/2                               (gn2_tx_8 + gn1_tx_8)/2
                           (gn2_tx_6 + gn1_tx_6)/2                                gn2_tx_9 
                            gn2_tx_7                                             ]
                         ]

            //similarly my_avg_v2

            At this point each process has its nodes and the traction forces computed at nodes.
            Now we need to combine this for all processes
        */




        MPI_Barrier(MPI_COMM_WORLD);
        int rank = get_rank();
        
        const int send_buffer_size = my_unique_nd.size();      
        int sum;
        std::vector<int> num_el_buffer, u;   
        MPI_NumElbuffer_U_Sum(send_buffer_size, num_el_buffer, u, sum);   

        /**            
                        pid =  0   1   
                       data = ----|---
              num_el_buffer =  4   3    
                 position u = 0---4   
                        sum = 7    
        */               
        

        std::vector<int> all_pid_nd(sum);
        MPI_Allgatherv(&my_unique_nd[0], num_el_buffer[rank], MPI_INT, &all_pid_nd[0], &num_el_buffer[0], &u[0], MPI_INT, MPI_COMM_WORLD);

        std::vector<double> all_pid_v1(sum);
        MPI_Allgatherv(&my_avg_v1[0], num_el_buffer[rank], MPI_DOUBLE, &all_pid_v1[0], &num_el_buffer[0], &u[0], MPI_DOUBLE, MPI_COMM_WORLD);

        std::vector<double> all_pid_v2(sum);
        MPI_Allgatherv(&my_avg_v2[0], num_el_buffer[rank], MPI_DOUBLE, &all_pid_v2[0], &num_el_buffer[0], &u[0], MPI_DOUBLE, MPI_COMM_WORLD);


        /** 
                    P0                                                                    P1
            all_pid_nd = [4, 5, 6, 7, 7, 8, 9]                                        same as P0
            all_pid_v1 = [tx_4, tx_5, tx_6, tx_7, tx_7, tx_8, tx_9]
            all_pid_v2 = [ty_4, ty_5, ty_6, ty_7, ty_7, ty_8, ty_9]                                        
        */



        auto all_pid_unique_nds = unique_vector(all_pid_nd);     /// all_pid_nd has repeated enteries; we make it unique
        std::vector<double> all_pid_avg_v1, all_pid_avg_v2;

        /// Here we average the tractions over the repeated nodes
        avg_over_repeated_nodes(all_pid_nd, all_pid_v1, all_pid_v2, all_pid_unique_nds, all_pid_avg_v1, all_pid_avg_v2);
        
       
        for(int i=0; i<all_pid_unique_nds.size(); i++)
          nd_var_map[all_pid_unique_nds[i]] = std::vector<double>{all_pid_avg_v1[i], all_pid_avg_v2[i]};
   }
  ///----------------------------------------------------------------------------------------- 









  ///----------------------------------------------------------------------------------------- 
  /// This function computes a map (nd_var_map) from my_nd, my_v1, my_v2 and my_v3 (3D version)
  /// this is called from   sc_f_tr<3>,   sc_isothermal_f_tr<3>, mc_f_tr<3>,  mc_isothermal_f_tr<3> and mc_f_heat_flux<3>
   void all_pid_avg_over_repeated_nodes(const std::vector<int>& my_nd, const std::vector<double>& my_v1, const std::vector<double>& my_v2, 
                        const std::vector<double>& my_v3, std::map<int, std::vector<double>>& nd_var_map)
   {   
        /**  example demonstration for traction on two processes with quadrangle sides
        
                            P0                                                        P1
            12(n4)  11(n3) 11(n4) 10(n3)  10(n4)   9(n3)                  9(n4)   8(n3)  8(n4)   7(n3) 
               o-------o    o-------o     o-------o                          o-------o    o-------o    
               |       |    |       |     |       |                          |       |    |       |    
               |       |    |       |     |       |                          |       |    |       |    
               o-------o    o-------o     o-------o                          o-------o    o-------o    
            1(n1)   2(n2)  2(n1) 3(n2)   3(n1)    4(n2)                  4(n1)    5(n2)   5(n1)    6(n2) 
        
              my_nd = [                                                my_nd = [
                         1,2,11,12                                               4,5,8,9 
                         2,3,10,11                                               5,6,7,8 
                         3,4,9,10                                              ]
                      ]                          

              my_v1 = [                                                my_v1 = [
                     gn1_tx_1, gn2_tx_2, gn3_tx_11, gn4_tx_12                   gn1_tx_4, gn2_tx_5, gn3_tx_8, gn4_tx_9
                     gn1_tx_2, gn2_tx_3, gn3_tx_10, gn4_tx_11                   gn1_tx_5, gn2_tx_6, gn3_tx_7, gn4_tx_8
                     gn1_tx_3, gn2_tx_4, gn3_tx_9, gn4_tx_10                   ]
                    ]
              
              my_v2 = [                                                my_v2 = [
                     gn1_ty_1, gn2_ty_2, gn3_ty_11, gn4_ty_12                   gn1_ty_4, gn2_ty_5, gn3_ty_8, gn4_ty_9
                     gn1_ty_2, gn2_ty_3, gn3_ty_10, gn4_ty_11                   gn1_ty_5, gn2_ty_6, gn3_ty_7, gn4_ty_8
                     gn1_ty_3, gn2_ty_4, gn3_ty_9, gn4_ty_10                   ]
                    ]
              
              my_v3 = [                                                my_v3 = [
                     gn1_tz_1, gn2_tz_2, gn3_tz_11, gn4_tz_12                   gn1_tz_4, gn2_tz_5, gn3_tz_8, gn4_tz_9
                     gn1_tz_2, gn2_tz_3, gn3_tz_10, gn4_tz_11                   gn1_tz_5, gn2_tz_6, gn3_tz_7, gn4_tz_8
                     gn1_tz_3, gn2_tz_4, gn3_tz_9, gn4_tz_10                   ]
                    ]
         */        


        auto my_unique_nd = unique_vector(my_nd);           /// my_nd has repeated enteries; we make it unique
        std::vector<double> my_avg_v1, my_avg_v2, my_avg_v3;        

        /// Here we average my_v1, my_v2 and my_v3 over the repeated nodes
        avg_over_repeated_nodes(my_nd, my_v1, my_v2, my_v3, my_unique_nd, my_avg_v1, my_avg_v2, my_avg_v3);


        /*
                    P0                                                P1
            my_unique_nd = [1,2,11,12,3,10,4,9]                     my_unique_nd = [4,5,8,9,6,7]
            
            my_avg_v1 = [                                           my_avg_v1 = [
                            gn1_tx_1                                              gn1_tx_4
                           (gn2_tx_2 + gn1_tx_2)/2                               (gn2_tx_5 + gn1_tx_5)/2
                           (gn3_tx_11 + gn4_tx_11)/2                             (gn3_tx_8 + gn4_tx_8)/2
                            gn4_tx_12                                             gn4_tx_9
                           (gn2_tx_3 + gn1_tx_3)/2                                gn2_tx_6
                           (gn3_tx_10 + gn4_tx_10)/2                              gn3_tx_7
                            gn2_tx_4                                            ]
                            gn3_tx_9                             
                       ]

            //similarly for my_avg_v2 and my_avg_v3

            At this point each process has its nodes and the traction forces computed at nodes.
            Now we need to combine this for all processes
        */




        MPI_Barrier(MPI_COMM_WORLD);
        int rank = get_rank();
        
        const int send_buffer_size = my_unique_nd.size();      
        int sum;
        std::vector<int> num_el_buffer, u;   
        MPI_NumElbuffer_U_Sum(send_buffer_size, num_el_buffer, u, sum);   

        /**            
                        pid =  0   1   
                       data = ----|---
              num_el_buffer =  4   3    
                 position u = 0---4   
                        sum = 7    
        */               


        std::vector<int> all_pid_nd(sum);
        MPI_Allgatherv(&my_unique_nd[0], num_el_buffer[rank], MPI_INT, &all_pid_nd[0], &num_el_buffer[0], &u[0], MPI_INT, MPI_COMM_WORLD);

        std::vector<double> all_pid_v1(sum);
        MPI_Allgatherv(&my_avg_v1[0], num_el_buffer[rank], MPI_DOUBLE, &all_pid_v1[0], &num_el_buffer[0], &u[0], MPI_DOUBLE, MPI_COMM_WORLD);

        std::vector<double> all_pid_v2(sum);
        MPI_Allgatherv(&my_avg_v2[0], num_el_buffer[rank], MPI_DOUBLE, &all_pid_v2[0], &num_el_buffer[0], &u[0], MPI_DOUBLE, MPI_COMM_WORLD);

        std::vector<double> all_pid_v3(sum);
        MPI_Allgatherv(&my_avg_v3[0], num_el_buffer[rank], MPI_DOUBLE, &all_pid_v3[0], &num_el_buffer[0], &u[0], MPI_DOUBLE, MPI_COMM_WORLD);

        /** 
                    P0                                                                    P1
            all_pid_nd = [4, 5, 6, 7, 7, 8, 9]                                        same as P0
            all_pid_tx = [tx_4, tx_5, tx_6, tx_7, tx_7, tx_8, tx_9]
            all_pid_ty = [ty_4, ty_5, ty_6, ty_7, ty_7, ty_8, ty_9]                                        
        */

        auto all_pid_unique_nds = unique_vector(all_pid_nd);     /// all_pid_nd has repeated enteries; we make it unique
        std::vector<double> all_pid_avg_v1, all_pid_avg_v2, all_pid_avg_v3;

        /// Here we average the tractions over the repeated nodes
        avg_over_repeated_nodes(all_pid_nd, all_pid_v1, all_pid_v2, all_pid_v3, all_pid_unique_nds, all_pid_avg_v1, all_pid_avg_v2, all_pid_avg_v3);
        
       
        for(int i=0; i<all_pid_unique_nds.size(); i++)
          nd_var_map[all_pid_unique_nds[i]] = std::vector<double>{all_pid_avg_v1[i], all_pid_avg_v2[i], all_pid_avg_v3[i]};
   }   
  ///----------------------------------------------------------------------------------------- 




} //namespace GALES

#endif
