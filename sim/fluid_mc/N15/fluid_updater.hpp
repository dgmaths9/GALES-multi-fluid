#ifndef _FLUID_UPDATER_HPP_
#define _FLUID_UPDATER_HPP_

#include "../../../src/fem/fem.hpp"


namespace GALES {


  struct fluid_updater 
  {

    using vec = boost::numeric::ublas::vector<double>;

    fluid_updater(dof_state& state): state_(state) {}


    void predictor()
    {
      state_.dofs(1) = state_.dofs(0); //P = last corrected
    }



    void corrector(const vec& correction)
    {
      state_.dofs(0) += correction;             //I   = I + correction
    }

    private:
    dof_state& state_;


 };
}// namespace GALES
#endif

