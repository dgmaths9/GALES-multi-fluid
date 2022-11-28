#ifndef POINT_HPP
#define POINT_HPP



namespace GALES{


 
  /**
      Base point class for node  
  */
 
  template<int dim>  class point;
  


  template<>
  class point<2>
  {
    public:
    
    point() : x_(0.0), y_(0.0) {}
    point(double x, double y) : x_(x), y_(y) {}
         
    void set_x(double x) { x_ = x;}                                  /// set x-coord of the point
    double get_x() const {return x_;}                                /// get x-coord of the point

    void set_y(double y) { y_ = y;}                                  /// set y-coord of the point
    double get_y() const {return y_;}                                /// get y-coord of the point

    private:
    double x_ = 0.0;
    double y_ = 0.0;
    static const int dimension = 2;    
  };







  template<>
  class point<3>
  {
    public:

    point() : x_(0.0), y_(0.0), z_(0.0) {}
    point(double x, double y, double z) : x_(x), y_(y), z_(z) {}

    void set_x(double x) { x_ = x;}                                  /// set x-coord of the point
    double get_x() const {return x_;}                                /// get x-coord of the point

    void set_y(double y) { y_ = y;}                                  /// set y-coord of the point
    double get_y() const {return y_;}                                /// get y-coord of the point

    void set_z(double z) { z_ = z;}                                  /// set z-coord of the point
    double get_z() const {return z_;}                                /// get z-coord of the point

    private:
    double x_ = 0.0;
    double y_ = 0.0;
    double z_ = 0.0;
    static const int dimension = 3;    
  };


  typedef point<1> point_1d;
  typedef point<2> point_2d;
  typedef point<3> point_3d;



}

#endif
