#include <mpi.h>
#include <iostream>
 

#include "../../../src/fem/fem.hpp"
#include "fluid_ic_bc.hpp"
#include "fluid_assembly.hpp"
#include "fluid_updater.hpp"
#include "fluid_dofs_ic_bc.hpp"
#include "fluid_2d.hpp"






using namespace GALES;


int main(int argc, char* argv[])
{

  MPI_Init(&argc, &argv);




  //---------- creating directories ----------------------------------------------------------------------
  make_dirs("results/fluid_dofs");







  const int dim = 2;
  point<dim> dummy;



  //---------------------- ic bc -------------------------------------------------------------------------------------------------------
  typedef fluid_ic_bc<dim> fluid_ic_bc_type;
  
  //----------------------- dofs ic bc---------------------------------------------------------------------------------------------------------
  typedef fluid_dofs_ic_bc<fluid_ic_bc_type, dim> fluid_dofs_ic_bc_type;

  //--------------------- scatters --------------------------------------------------------------------------------------------------------
  typedef scatter<fluid_dofs_ic_bc_type>  fluid_scatter_type; 

  //------------------ pressure profile -------------------------------------------------------------------------------------------------
  typedef pressure_profile_mc<fluid_ic_bc_type, dim> pressure_profile_type;

  //-------------------- integrals -------------------------------------------------------------------------------------------------------
  typedef fluid<fluid_ic_bc_type, dim> fluid_integral_type;

  //-------------------- assembly ------------------------------------------------------------------------------------------------------
  typedef fluid_assembly<fluid_integral_type, fluid_scatter_type, dim> fluid_assembly_type;













  //---------------- read setup file -------------------------------
  read_setup setup;

  
  //-----------------fluid props------------------------------------------
  fluid_properties fluid_props;



  const int nb_comp = fluid_props.nb_comp();
  const int nb_dofs_f = 2 + dim + nb_comp - 1;  


  

  //------------------ mesh building --------------------------------------------------------------------
  print_only_pid<0>(std::cerr)<< "FLUID MESH\n";
  double mesh_read_start = MPI_Wtime();
  Mesh<dim> fluid_mesh(setup.fluid_mesh_file(), nb_dofs_f, setup);
  print_only_pid<0>(std::cerr)<<"Mesh reading took: "<<std::setprecision(setup.precision()) << MPI_Wtime()-mesh_read_start<<" s\n\n";



  
  //-----------------------------------------maps------------------------------------------------------------
  epetra_maps fluid_maps(fluid_mesh);
  
  //----------------------------------------dofs state-----------------------------------------------------  
  auto fluid_dof_state = std::make_unique<dof_state>(fluid_maps.state_map()->NumMyElements(), 2);
  
  //---------------------------------------model---------------------------------------------------
  model<dim> fluid_model(fluid_mesh, *fluid_dof_state, fluid_maps, setup);  

  //----------------------------------------ic bc---------------------------------------------------------
  fluid_ic_bc_type fluid_ic_bc;

  //----------------------------------------dofs ic bc---------------------------------------------------------
  fluid_dofs_ic_bc_type fluid_dofs_ic_bc(fluid_model, fluid_ic_bc, nb_comp);
  
  //--------------pressure profile-----------------------------------------------------
  pressure_profile_type pressure_profile(fluid_ic_bc, fluid_model, fluid_props);


  for(const auto& nd : fluid_mesh.nodes())
    if(nd->get_y() <= -3420.0)
    {
      auto p = fluid_dof_state->get_dof(nd->first_dof_lid());
      fluid_dof_state->set_dof(nd->first_dof_lid(), p+15.e6);
    }

  //----------------------scatter------------------------------------------------------
  fluid_scatter_type fluid_scatter(fluid_dofs_ic_bc);   

  
  //---------------------updater------------------------------------------------------
  fluid_updater fluid_updater(*fluid_dof_state);

      
  //--------------------linear system-------------------------------------------------
  linear_system fluid_lp(fluid_model);
  
        
  //--------------------linear solver----------------------------------------------------
  belos_linear_solver fluid_ls_solver(setup, fluid_lp); 


  //---------------------non linear residual check-------------------------------------------
  residual_check fluid_res_check;
  
  //------------------------integral---------------------------------------------------------  
  fluid_integral_type fluid_integral(fluid_ic_bc, fluid_props, fluid_model); 

  //-------------------------assembly--------------------------------------------------------
  fluid_assembly_type fluid_assembly(fluid_integral, fluid_scatter, fluid_model, fluid_lp);


  // ------------------parallel I/O---------------------------------------------------  
  IO fluid_io(*fluid_maps.state_map(), *fluid_maps.dof_map());

  //----------adaptive time step-------------------------------------------------
  adaptive_time_criteria<dim> criteria;































  
  
  time::get().t(0.0);
  time::get().delta_t(setup.delta_t());
    
  double restart_or_t_zero = MPI_Wtime();  
  if(!setup.restart())
  {
    //--------- setting up dofs for time = 0 ----------------
    fluid_dofs_ic_bc.dirichlet_bc();
    fluid_io.write("fluid_dofs/", *fluid_dof_state);
    print_only_pid<0>(std::cerr)<<"Dofs writing for time 0 took: "<<std::setprecision(setup.precision()) << MPI_Wtime()-restart_or_t_zero<<" s\n\n";       
  }
  else
  {
    //--------------read restart-----------------------------
    time::get().t(setup.restart_time());
    fluid_io.read("fluid_dofs/", *fluid_dof_state);
    print_only_pid<0>(std::cerr)<<"Dofs reading at restart time took: "<<std::setprecision(setup.precision()) << MPI_Wtime()-restart_or_t_zero<<" s\n\n";       
  }
  time::get().tick();



  // --------------------------------TIME LOOP-----------------------------------------------------------
  double time_loop_start = MPI_Wtime();
  double non_linear_res(0.0);
  int tot_time_iterations(0);

  for(; time::get().t()<=setup.end_time(); time::get().tick())
  {
    double t_iteration(MPI_Wtime());

    print_only_pid<0>(std::cerr)<<"simulation time:  "<<std::setprecision(setup.precision()) << time::get().t()<<" s\n";

    fluid_updater.predictor();



      for(int it_count=0; it_count < setup.n_max_it(); it_count++)   //fluid iterations loop
      {
       //------------------------ SOLUTION For fluid ------------------------
        if(setup.t_dependent_dirichlet_bc())
           fluid_dofs_ic_bc.dirichlet_bc();
	
	double fill_time = fluid_assembly.execute();   


        bool residual_converged = fluid_res_check.execute(*fluid_lp.rhs(), it_count, non_linear_res);      
        if(residual_converged)
        {
          print_only_pid<0>(std::cerr) << "FLUID  it: " << parse(it_count,5) << "Assembly: " << parse(fill_time,15) << "NonLinearRes: " << parse(non_linear_res,15) << "\n";        
          break;
        } 	              
        double solver_time = fluid_ls_solver.execute();


	fluid_updater.corrector(fluid_io.state_vec(*fluid_ls_solver.solution()));
	fluid_model.dof_trimmer(4, 0, 1);


        print_only_pid<0>(std::cerr) 
        << "FLUID  it: " << parse(it_count,5) << "Assembly: " << parse(fill_time,15) << "NonLinearRes: " << parse(non_linear_res,15)<< "Solver: " << parse(solver_time,15)
        << "AbsResErr: " << parse(fluid_ls_solver.ARE(),15) << "RelResErr: "<< parse(fluid_ls_solver.RRE(),15)
        << "Num_it: " << parse(fluid_ls_solver.num_it(),8) << "dt: "<< time::get().delta_t()
        << "\n";

        // -----------------------end of fluid solution---------------------------------------------------------
       }
       



    if((tot_time_iterations+1)%(setup.print_freq()) == 0)
    {
       fluid_io.write("fluid_dofs/", *fluid_dof_state);
    }
       


     //-------- adaptive time step-----------------
     if(setup.adaptive_time_step() == true && tot_time_iterations > setup.N())
     {
        criteria.f_mc(fluid_model, fluid_props, setup, nb_dofs_f);
     }



    print_only_pid<0>(std::cerr)<<"Time step took: "<<MPI_Wtime()-t_iteration<<" s\n\n";
    tot_time_iterations++;    	
  }

  print_only_pid<0>(std::cerr)<<"End time reached!!!! "<<"\n\n";
  print_only_pid<0>(std::cerr)<<"Total run time: "<<MPI_Wtime()-time_loop_start<<" s\n\n\n";

  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
  return 0;
}
