#ifndef AMESOS_LINEAR_SOLVER_HPP
#define AMESOS_LINEAR_SOLVER_HPP




#include "Amesos_ConfigDefs.h"
#include "Amesos.h"


#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Amesos.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Import.h"
#include "Epetra_Export.h"
#include "Epetra_Map.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "Epetra_LinearProblem.h"
#include "Teuchos_ParameterList.hpp"





  /**
       This class is a wrapper to amesos direct linear solver factory
       
       Specifies the solver. String "SolverType" can assume one of the following values:
                     _
       - Lapack       |
       - Klu          |<------- These are prebuilt as a part of trilinos   
       - Paraklete   _|
                     
                     _
       - Umfpack      |
       - Pardiso      |
       - Taucs        |
       - Superlu      |<------- These are third party libraries to trilinos 
       - Superludist  |
       - Mumps        |
       - Dscpack      |
       - Scalapack   _|                     
  */


 
namespace GALES {





  class amesos_linear_solver
  {

    public:


    explicit amesos_linear_solver(read_setup& setup)
    : 
      solver_name_(setup.ls_solver()), 
      X_(nullptr) 
    {}



      //-------------------------------------------------------------------------------------------------------------------------------------        
      /// Deleting the copy and move constructors - no duplication/transfer in anyway
      amesos_linear_solver(const amesos_linear_solver&) = delete;               //copy constructor
      amesos_linear_solver& operator=(const amesos_linear_solver&) = delete;    //copy assignment operator
      amesos_linear_solver(amesos_linear_solver&&) = delete;                    //move constructor  
      amesos_linear_solver& operator=(amesos_linear_solver&&) = delete;         //move assignment operator 
      //-------------------------------------------------------------------------------------------------------------------------------------



      int execute(Epetra_FECrsMatrix& A, Epetra_FEVector& b)
      {
        X_ = std::make_shared<Epetra_Vector>(A.RowMap());        
        X_->PutScalar(0.0);
      
        Epetra_LinearProblem Problem(&A, &*X_, &b);
        Amesos_BaseSolver* Solver = 0;
        Amesos Factory;
        std::string SolverType = solver_name_;
        Solver = Factory.Create(SolverType, Problem);
        if (Solver == 0) 
        {
          Error("Specified solver is not available");
        }
  
  
        Teuchos::ParameterList List;
        List.set("PrintTiming", true);
        List.set("PrintStatus", true);  
        Solver->SetParameters(List);


        Solver->SymbolicFactorization();
        Solver->NumericFactorization();
        Solver->Solve();

   
        Epetra_MultiVector Ax(*X_);
        A.Multiply(false, *X_, Ax);
        Ax.Update(1.0, b, -1.0);
        Ax.Norm2(&abs_res_err_);
        double rhs_norm;
        b.Norm2(&rhs_norm);
        rel_res_err_ = abs_res_err_/rhs_norm;

  
        delete Solver;          
        return 0;
      }




      auto solution()
      {
        return X_;
      }



      auto RRE()  /// relative residual error
      {
        return rel_res_err_;
      }


      auto ARE()  /// absolute residual error
      {
        return abs_res_err_;
      }



      auto num_it()
      {
        return 1;
      }





    private:
      std::string solver_name_;
      double abs_res_err_;
      double rel_res_err_;
      std::shared_ptr<Epetra_Vector> X_;

    }; 


}

  #endif
