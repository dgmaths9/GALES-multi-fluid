#ifndef GALES_NON_NEWTONIAN_HPP
#define GALES_NON_NEWTONIAN_HPP


#include <numeric>
#include <algorithm>
#include <math.h>



namespace GALES {


/**
   This file defines calculation of shear rate and some general laws for non newtonian flows such as power law and Herschel_Bulkley model
*/



using vec = boost::numeric::ublas::vector<double>;


double compute_shear_rate(const vec& dI_dx, const vec& dI_dy)
{
   const double dv1_1(dI_dx[1]), dv2_1(dI_dx[2]);
   const double dv1_2(dI_dy[1]), dv2_2(dI_dy[2]);
   const double shear_rate = sqrt(std::pow(dv1_1,2) + std::pow(dv2_2,2) + 2*std::pow(dv1_2 + dv2_1,2) );
   return shear_rate;
}


double compute_shear_rate(const vec& dI_dx, const vec& dI_dy, const vec& dI_dz)
{
   const double dv1_1(dI_dx[1]), dv2_1(dI_dx[2]), dv3_1(dI_dx[3]);
   const double dv1_2(dI_dy[1]), dv2_2(dI_dy[2]), dv3_2(dI_dy[3]);
   const double dv1_3(dI_dz[1]), dv2_3(dI_dz[2]), dv3_3(dI_dz[3]);
   const double shear_rate = sqrt(std::pow(dv1_1,2) + std::pow(dv2_2,2) + std::pow(dv3_3,2) + 2*std::pow(dv1_2 + dv2_1,2) + 2*std::pow(dv1_3 + dv3_1,2) + 2*std::pow(dv2_3 + dv3_2,2));
   return shear_rate;
}


double Power_Law(double consistency_index, double power_law_index, double shear_rate)
{
  if(shear_rate < std::numeric_limits<double>::epsilon())  return 0.0;
  else  return consistency_index*pow(shear_rate, power_law_index-1);
}



double Herschel_Bulkley(double consistency_index, double power_law_index, double yield_stress, double critical_shear_rate, double shear_rate)
{
  if(shear_rate < std::numeric_limits<double>::epsilon())  return 0.0;
  else if(shear_rate < critical_shear_rate)  return consistency_index*pow(critical_shear_rate, power_law_index-1) + yield_stress/critical_shear_rate;
  else  return consistency_index*pow(shear_rate, power_law_index-1) + yield_stress/shear_rate;
}



}
#endif
