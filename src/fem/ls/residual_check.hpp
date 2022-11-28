#ifndef RESIDUAL_CHECK_HPP
#define RESIDUAL_CHECK_HPP


#include <vector>
#include "Epetra_FEVector.h"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_StandardCatchMacros.hpp"
#include "Teuchos_RCPDecl.hpp"
#include "BelosConfigDefs.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosEpetraAdapter.hpp"



namespace GALES{


  /**
      This class computes the non linear residual which is the norm of the rhs vector.
      At the first iteration of each time step we compute absolute norm of the rhs vector
      and for the next iterations we compute the relative norm of the rhs vector.
      We return true if "non_linear_res_it" is small means the solution converges.  
  */


  class residual_check
  {
     public:
      
      //------------------------------------------------------- ------------------------------------------------------------------------------        
      residual_check() = default;                       // default constructor
      /// Deleting the copy and move constructors - no duplication/transfer in anyway
      residual_check(const residual_check&) = delete;               //copy constructor
      residual_check& operator=(const residual_check&) = delete;    //copy assignment operator
      residual_check(residual_check&&) = delete;                    //move constructor  
      residual_check& operator=(residual_check&&) = delete;         //move assignment operator 
      //-------------------------------------------------------------------------------------------------------------------------------------
      
      
               
      bool execute(Epetra_FEVector& b, int it_count, double& non_linear_res)
      {
        typedef double                            ST;
        typedef Epetra_MultiVector                MV;
        typedef Belos::MultiVecTraits<ST,MV>     MVT;
  
        const Teuchos::RCP<Epetra_FEVector> b_rcp(new Epetra_FEVector(b));
        Teuchos::RCP<Epetra_MultiVector> B = Teuchos::rcp_implicit_cast<Epetra_MultiVector>(b_rcp); 
        std::vector<double> rhs_norm(1);        
        MVT::MvNorm( *B, rhs_norm );

 
        if(it_count==0)
        {  
            non_linear_res_it0_ = rhs_norm[0];
            non_linear_res = rhs_norm[0];  
            return false;
        }  
        else
        {
            const double non_linear_res_it = rhs_norm[0];            
            non_linear_res = non_linear_res_it/non_linear_res_it0_;            
            
            /**
                                M dx = -R
                 criterion:
                 
                 ||R||  < 1.e-5                 absolute err
             or   
                 ||R||
               ----------  <  1.e-8             relative err
                 ||R||_0
            */
            
            if(non_linear_res_it <= std::max(1.e-8*non_linear_res_it0_, 1.e-5))
              return true;
        }
        return false;
     }
     
        
     
    private:    
      double non_linear_res_it0_ = 0.0;  

  };



}

#endif
