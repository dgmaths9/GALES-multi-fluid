#ifndef _GALES_TIME_HPP_
#define _GALES_TIME_HPP_




namespace GALES {

  
  /**
      This class is defined with a singleton design pattern. This means only one instance of this class is generated throughout the code.
      This is good for special class like time to be defined globally so that time related functions are available in all the files.      
      This is sometimes called the Meyers’ Singleton after Scott Meyers publication “C++ and the Perils of Double-Checked Locking”. 
      The constructor is made private so that objects can not be instantiation outside of the class. 
      Both the copy and move constructors are deleted, so that this class cannot be duplicated/transferred in anyway.
      Since C++11, the initialization of magic statics is guaranteed to be thread safe. 
      Internally, the first thread that calls get() will initialize the pointer, blocking until initialization is complete. All other threads will then use the initialized value. 
      Before C++11, calling get() from multiple threads was classified as undefined behaviour (UB). 
      Most C++ compilers support magic statics, but you have to be careful about the thread safety aspect, as this was only introduced into the C++ standard in C++11! 
      Any C++11 compliant version of GCC will support thread-safe magic statics.
      On the Windows side of things, thread safety was implemented for magic statics in Visual Studio 2015. The earlier versions are not thread safe.
                
      The functions of this class can be called as
      
      time::get().t()  -----> to get current time
      time::get().delta_t() ------> to get time step
      time::get().tick(dt) -----> to update current time t by dt
      
      similarly others
  */


  class time
  {
 
    //----------------------- singleton makeing functions -------------------------------------
    private:                       
    time(){};          /// Disallow instantiation outside of the class.

    public:
    /// Deleting the copy and move constructors - no duplication/transfer in anyway
    time(const time&) = delete;
    time& operator=(const time&) = delete;
    time(time&&) = delete;
    time& operator=(time&&) = delete;

    static auto & get()
    {
        /// Here we create magic static
        static time instance;  
        return instance;
    }
    //-----------------------------------------------------------------------------------------
        
   
    public:
    
    /// This function sets the time step "delta_t_" 
    void delta_t(double dt) {delta_t_ = dt; }

    /// This function returns the time step "delta_t_"
    double delta_t()const {return delta_t_; }

    /// This function sets the time "t_" 
    void t(double t) {t_ = t; }

    /// This function returns the time "t_" 
    double t()const {return t_; }

    /// This function updates the curent time "t_ = t_ + delta_t_" 
    void tick() {t_ += delta_t_;}


    private: 	
    double delta_t_ = 0.0;
    double t_ = 0.0;
  };







}
#endif
