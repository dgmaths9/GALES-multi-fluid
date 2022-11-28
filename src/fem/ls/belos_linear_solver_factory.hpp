#ifndef BELOS_LINEAR_SOLVER_HPP
#define BELOS_LINEAR_SOLVER_HPP



#include "BelosSolverFactory.hpp"

#include "BelosConfigDefs.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosEpetraAdapter.hpp"

#include "EpetraExt_readEpetraLinearSystem.h"
#include "Epetra_Map.h"
#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_CrsMatrix.h"

#include "Ifpack.h"


#include "BelosGCRODRSolMgr.hpp"
#include "BelosPseudoBlockGmresSolMgr.hpp"
#include "BelosPseudoBlockCGSolMgr.hpp"
#include "BelosBlockGmresSolMgr.hpp"
#include "BelosBlockCGSolMgr.hpp"
#include "BelosFixedPointSolMgr.hpp"
#include "BelosLSQRSolMgr.hpp"
#include "BelosPCPGSolMgr.hpp"
#include "BelosRCGSolMgr.hpp"
#include "BelosBiCGStabSolMgr.hpp"

#include "ml_include.h" 
#include "ml_MultiLevelPreconditioner.h"

#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_StandardCatchMacros.hpp"






/**  This is a generalised class for belos solvers. 
     User can select the solver and preconditioner by setting solver_name and precond
     Options are:
  
    solver_name
    1.  GMRES
    2.  Flexible GMRES
    3.  Block CG 
    4.  PseudoBlockCG
    5.  Stochastic CG
    6.  Recycling GMRES 
    7.  Recycling CG
    8.  MINRES 
    9.  LSQR 
    10. TFQMR
    11. Pseudoblock
    12. GmresPoly
    13. CGPoly 
*/

namespace GALES {



  class belos_linear_solver
  {

    public:

    belos_linear_solver(read_setup& setup, linear_system& lp)
    :
      rel_res_tol_(setup.ls_rel_res_tol()),  
      maxsubspace_(setup.ls_maxsubspace()),  
      maxrestarts_(setup.ls_maxrestarts()),
      maxiters_(setup.ls_maxiters()),
      fill_(setup.ls_fill()), 
      solver_name_(setup.ls_solver()), 
      precond_(setup.ls_precond()),
      details_(setup.ls_details()),
      lp_(lp)
    {}




    //-------------------------------------------------------------------------------------------------------------------------------------        
    /// Deleting the copy and move constructors - no duplication/transfer in anyway
    belos_linear_solver(const belos_linear_solver&) = delete;               //copy constructor
    belos_linear_solver& operator=(const belos_linear_solver&) = delete;    //copy assignment operator
    belos_linear_solver(belos_linear_solver&&) = delete;                    //move constructor  
    belos_linear_solver& operator=(belos_linear_solver&&) = delete;         //move assignment operator 
    //-------------------------------------------------------------------------------------------------------------------------------------




      double execute()
      {
        double start(MPI_Wtime());

        typedef Epetra_MultiVector                MV;
        typedef Epetra_Operator                   OP;
        typedef Belos::MultiVecTraits<double,MV>     MVT;
        typedef Belos::OperatorTraits<double,MV,OP>  OPT;

        using Teuchos::ParameterList;
        using Teuchos::parameterList;
        using Teuchos::RCP;
        using Teuchos::rcp; 
        using Teuchos::rcp_implicit_cast;



        // ---------------------Get the problem--------------------------------------------------------
        const RCP<Epetra_FECrsMatrix> A(new Epetra_FECrsMatrix(*(lp_.matrix())));

        
        const RCP<Epetra_FEVector> b_rcp(new Epetra_FEVector(*(lp_.rhs())));
        RCP<Epetra_MultiVector> B = rcp_implicit_cast<Epetra_MultiVector>(b_rcp);
                
        
        const RCP<Epetra_FEVector> x_rcp(new Epetra_FEVector(A->RowMap()));
        x_rcp->PutScalar(0.0);
        X_ = rcp_implicit_cast<Epetra_MultiVector>(x_rcp);
        //----------------------------------------------------------------------------------------------



        if(maxiters_==-1)  maxiters_ = B->GlobalLength()-1;
        
 



        //-------------------------Construct preconditioner------------------------------
        RCP<Belos::EpetraPrecOp> M;
        if(precond_ == "ILU" || precond_ == "ICT")
        {
            ParameterList ifpackList;
            Ifpack Factory;
            std::string PrecType = precond_;
            int OverlapLevel = 1; // must be >= 0.
            RCP<Ifpack_Preconditioner> Prec = rcp( Factory.Create(PrecType, &*A, OverlapLevel) );
            assert(Prec != Teuchos::null);
            ifpackList.set("fact: level-of-fill", fill_);
            ifpackList.set("schwarz: combine mode", "Add");
            IFPACK_CHK_ERR(Prec->SetParameters(ifpackList));
            IFPACK_CHK_ERR(Prec->Initialize());
            IFPACK_CHK_ERR(Prec->Compute());
            M = rcp( new Belos::EpetraPrecOp( Prec ) );
        }
        else if(precond_ == "ML")
        {
            ParameterList MLList; // Set MLList for Smoothed Aggregation
            ML_Epetra::SetDefaults("SA", MLList); // reset parameters ML User's Guide
            MLList.set("smoother: type","Chebyshev"); // Chebyshev smoother  ... aztec??
            MLList.set("smoother: sweeps",3);
            MLList.set("smoother: pre or post", "both"); // both pre- and post-smoothing
            MLList.set("coarse: type","Amesos-KLU"); // solve with serial direct solver KLU
            RCP<Epetra_Operator> Prec = rcp(  new ML_Epetra::MultiLevelPreconditioner(*A, MLList) );
            assert(Prec != Teuchos::null);
            M = rcp( new Belos::EpetraPrecOp( Prec ) );
        }    
        //-----------------------------------------------------------------------------------


 


        //------------------------set parameters----------------------------------------------
        Belos::SolverFactory<double, MV, OP> factory;
        RCP<ParameterList> solverParams = parameterList();
        solverParams->set ("Num Blocks", maxsubspace_);
        solverParams->set ("Block Size", 1);
        solverParams->set ("Maximum Iterations", maxiters_);
        solverParams->set ("Maximum Restarts", maxrestarts_);
        solverParams->set ("Convergence Tolerance", rel_res_tol_);
        if (details_)
          solverParams->set ("Verbosity", Belos::Errors + Belos::Warnings + Belos::TimingDetails + Belos::StatusTestDetails );
        //---------------------------------------------------------------------------









        // Solve the linear problem. A, X, B, and M are passed by (smart) pointer, not copied.
        RCP<Belos::SolverManager<double, MV, OP> > solver = factory.create (solver_name_, solverParams);        
        RCP<Belos::LinearProblem<double, MV, OP> > problem = rcp (new Belos::LinearProblem<double, MV, OP> (A, X_, B));        
        problem->setRightPrec(M);
        problem->setProblem();
        solver->setProblem (problem);        
        solver->solve();
        num_iterations_ = solver->getNumIters();


        std::vector<double> rhs_norm(1);        
        MVT::MvNorm( *B, rhs_norm );
        Epetra_MultiVector resid(A->RowMap(), 1);
        OPT::Apply( *A, *X_, resid );
        MVT::MvAddMv( -1.0, resid, 1.0, *B, resid );
        std::vector<double> res_error(1);        
        MVT::MvNorm( resid, res_error);
        abs_res_err_ = res_error[0];
        rel_res_err_ = res_error[0]/rhs_norm[0];


        return MPI_Wtime()-start;            
      }





      auto solution()
      {
        Epetra_Vector* x = (*X_)(0);
        return x;
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
        return num_iterations_;
      }





    private:
      Teuchos::RCP<Epetra_MultiVector> X_;
      double rel_res_err_;
      double abs_res_err_;
      double rel_res_tol_;
      int maxsubspace_;
      int maxrestarts_;
      int maxiters_;
      int fill_;
      std::string solver_name_;
      std::string precond_;
      int num_iterations_;
      bool details_ = false;
      linear_system& lp_;
    }; 



}

  #endif
