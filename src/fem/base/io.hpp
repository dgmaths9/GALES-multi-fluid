#ifndef _GALES_IO_HPP_
#define _GALES_IO_HPP_


#include <string>
#include <sstream>
#include <fstream>

#include "mpi.h"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_MpiComm.h"




namespace GALES{
  
  



  /**
    This class is for parallel I/O.
    We assume I/O to be performed using a unique and simple data map, in which a contiguous slice of the distributed vector is handled by each processor. 
  */    
  
  class IO
  {
    public:
    
    //-------------------------------------------------------------------------------------------------------------------------------------
    /// This constructor use state_map.MaxAllGID()+1 to construct a map and an epetra vector "tmp_" to be used in read and write functions.
    /// This also builds the dofMap_to_stateMap_importer_ to be used in state_vec.
    /// state_map and dof_map are implementaed in map class.
    IO(const Epetra_Map& state_map, const Epetra_Map& dof_map): state_map_(&state_map)
    {
      global_dim_ = state_map.MaxAllGID()+1;                                      //global_dim = mesh.tot_nodes()*nb_dofs
      Epetra_Map map(global_dim_, 0, Epetra_MpiComm(MPI_COMM_WORLD));             //this creates Epetra_Map with num of enteries = global_dim  
      tmp_ = std::make_unique<Epetra_Vector>(map);                                //a poiner to a temporary distributed Epetra_vector with map defined in line above
      local_dim_ = map.NumMyElements();                                           //local number of elements of Epetra_vector 
      local_start_ = map.MinMyGID();                                             //lowest gid of elements in Epetra_vector on each pid     
      stateMap_to_tmpMap_importer_ = std::make_unique<Epetra_Import>(tmp_->Map(), state_map);        // state_map ----> tmp_->map    importer
      tmpMap_to_stateMap_importer_ = std::make_unique<Epetra_Import>(state_map, tmp_->Map());    // tmp_->map ----> state_map    importer
      dofMap_to_stateMap_importer_ = std::make_unique<Epetra_Import>(state_map, dof_map);      // dof_map  ----> state_map    importer
    }
    //-------------------------------------------------------------------------------------------------------------------------------------
    




    //-------------------------------------------------------------------------------------------------------------------------------------    
    /// Here map is shared_node_map (used in mesh_writer and secondary dofs).
    explicit IO(const Epetra_Map& shared_node_map): shared_node_map_(&shared_node_map)
    {
        global_dim_ = shared_node_map.MaxAllGID()+1;                                //global_dim = mesh.tot_nodes()
        Epetra_Map map(global_dim_, 0, Epetra_MpiComm(MPI_COMM_WORLD));             //this creates Epetra_Map with num of enteries = global_dim  
        tmp_ = std::make_unique<Epetra_Vector>(map);                                //a poiner to a temporary distributed Epetra_vector with map defined in line above
        local_dim_ = map.NumMyElements();                                           //local number of elements of Epetra_vector 
        local_start_ = map.MinMyGID();                                             //lowest gid of elements in Epetra_vector on each pid     
     	sharedNodeMap_to_tmpMap_importer_ = std::make_unique<Epetra_Import>(tmp_->Map(), shared_node_map);  // shared_node_map ----> tmp_->map    importer
    }
    //-------------------------------------------------------------------------------------------------------------------------------------
    
    
    
    
    
    //-------------------------------------------------------------------------------------------------------------------------------------        
    /// Deleting the copy and move constructors - no duplication/transfer in anyway
    IO(const IO&) = delete;               //copy constructor
    IO& operator=(const IO&) = delete;    //copy assignment operator
    IO(IO&&) = delete;                    //move constructor  
    IO& operator=(IO&&) = delete;         //move assignment operator 
    //-------------------------------------------------------------------------------------------------------------------------------------
        



    //-------------------------------------------------------------------------------------------------------------------------------------
    /// This function open the file with name "fname" as a mpi file to read write in parallel
    /// MPI_Type_create_subarray creates an MPI datatype representing a subset of an array.
    void set_p_file(std::string fname, MPI_File& handle, MPI_Datatype& mysubarray, unsigned int flags=0)const
    { 
      if(MPI_File_open(MPI_COMM_WORLD, &fname[0], flags, MPI_INFO_NULL, &handle) != MPI_SUCCESS)    //opening the MPI_FILE
      {
         Error("MPI process " + std::to_string(get_rank()) + "Failure in opening the file.");
      }                   
      //it creates an MPI type which describes the memory layout of a subarray given: a larger array of some given type; a set of subsizes; and a "corner" from which to start.
      MPI_Type_create_subarray(1, &global_dim_, &local_dim_, &local_start_, MPI_ORDER_C, MPI_DOUBLE, &mysubarray);  

      MPI_Type_commit(&mysubarray);                                                                                      //It is necessary to commit mysubarray before using it      
      /*
           The MPI_File_set_view routine changes the process’s view of the data in the file 
           The beginning of the data accessible in the file through that view is set to disp (we set it 0 by calling MPI_Offset(0)); 
           The type of data is set to etype (MPI_DOUBLE); 
           Tthe distribution of data to processes is set to filetype (mysubarray). 
           In addition, MPI_File_set_view resets the independent file pointers and the shared file pointer to zero. MPI_File_set_view is collective across the handle; 
           The data types passed in etype (MPI_DOUBLE is standard  and already commited) and filetype must be committed.           
           The disp (MPI_Offset(0)) displacement argument specifies the position (absolute offset in bytes from the beginning of the file) where the view begins.           
      */
      if(MPI_File_set_view(handle, MPI_Offset(0), MPI_DOUBLE, mysubarray, &"native"[0], MPI_INFO_NULL) != MPI_SUCCESS)   //building the MPI view for the file
      {
         Error("MPI process " + std::to_string(get_rank()) + "Failure in setting the MPI_File_set_view.");      
      }        
    }    
    //-------------------------------------------------------------------------------------------------------------------------------------
    
    



    //-------------------------------------------------------------------------------------------------------------------------------------
    /// This function reads in parallel the data from the file to tmp_ Epetra_Vector 
    /// and then with that defines the state vector. 
    template<typename state_type>
    void read(std::string dir, state_type& state, const int state_level = 0)const
    {
      std::stringstream ss;
      ss<<"results/"<< dir<<time::get().t();
      std::ifstream ifile(ss.str());    	
      if(!(bool)ifile)
      {
         Error("The file   '" + ss.str() + "'  to read does not exist!!!!!");
      }  
      MPI_File handle;                                                  //MPI_File represents a file handle in MPI, it is used in MPI IO routines such as MPI_File_open for instance.
      MPI_Datatype mysubarray;                                          //MPI type to describe the memory layout to read data on local processes
      set_p_file(ss.str(), handle, mysubarray, MPI_MODE_RDONLY);        //opening and setting file view in parallel
      MPI_Status status;
            
      if(MPI_File_read_all(handle, &(*tmp_)[0], tmp_->MyLength(), MPI_DOUBLE, &status) != MPI_SUCCESS)   // here tmp_ is filled with data  
      {
         Error("MPI process " + std::to_string(get_rank()) + "Failure in MPI_File_read_all");      
      }        

      //------- copy data:    unique dofs --->shared dofs      tmp_ ---> state     
      Epetra_Vector app(*state_map_);
      app.Import(*tmp_, *tmpMap_to_stateMap_importer_, Insert);
      app.ExtractCopy(&(state.dofs(state_level))[0]);
      //----------------------------------------------------------

      MPI_Type_free(&mysubarray);
      MPI_File_close(&handle);
    }    
    //-------------------------------------------------------------------------------------------------------------------------------------
      





    //-------------------------------------------------------------------------------------------------------------------------------------
    /// This function first opens a file in parallel using set_p_file() function
    /// Then data is copied from state vector to tmp_ vector 
    /// and then write that data in parallel in file
    template<typename state_type>
    void write(std::string dir, const state_type& state, const int state_level = 0)const
    {
      std::stringstream ss;
      ss<<"results/"<< dir<<time::get().t();
      MPI_File handle;                                        //MPI_File represents a file handle in MPI, it is used in MPI IO routines such as MPI_File_open for instance.
      MPI_Datatype mysubarray;                                //MPI type to describe the memory layout to read data on local processes
      set_p_file(ss.str(), handle, mysubarray, MPI_MODE_WRONLY| MPI_MODE_CREATE);              //opening and setting file view in parallel

      //-------     copy data:  shared dofs---->unique dofs     state ----> tmp_   
      Epetra_Vector app(Copy, *state_map_, const_cast<double*>(&(state.dofs(state_level))[0]));	  	  
      tmp_->Import(app, *stateMap_to_tmpMap_importer_, Insert);
      //---------------------------------------------------------

      MPI_Status status;  
      // I noticed with krafla_sim (144 cores) and andesite_dacite (216 cores) MPI_File_write_all is deadlocking while MPI_File_write is working fine
      // Need to understand why????
      if( MPI_File_write(handle, &(*tmp_)[0], tmp_->MyLength(), MPI_DOUBLE, &status) != MPI_SUCCESS)   //write tmp_ in parallel  
      {
         Error("MPI process " + std::to_string(get_rank()) + "Failure in MPI_File_write_all");      
      }        
      
      MPI_Type_free(&mysubarray);
      MPI_File_close(&handle);
    }
    //-------------------------------------------------------------------------------------------------------------------------------------





    //-------------------------------------------------------------------------------------------------------------------------------------
    /// This has the same functionality as the previous function. It writes vector "v" in place of state vector.
    /// This is called to write the mesh and secondary dofs in parallel.
    /// Here map m is shared_node_map.
    void write2(std::string s, const std::vector<double>& v)const
    {
      MPI_File handle;                                       //MPI_File represents a file handle in MPI, it is used in MPI IO routines such as MPI_File_open for instance.
      MPI_Datatype mysubarray;                               //MPI type to describe the memory layout to read data on local processes
      set_p_file(s, handle, mysubarray, MPI_MODE_WRONLY| MPI_MODE_CREATE);        //opening and setting file view in parallel

      //-------    copy data:  shared dofs---->unique dofs     v  ----> tmp_
      Epetra_Vector app(Copy, *shared_node_map_, const_cast<double*>(&v[0]));	  	  
      tmp_->Import(app, *sharedNodeMap_to_tmpMap_importer_, Insert);
      //--------------------------------------------------------

      MPI_Status status;      
      if(MPI_File_write_all(handle, &(*tmp_)[0], tmp_->MyLength(), MPI_DOUBLE, &status) != MPI_SUCCESS)   //write tmp_ in parallel   
      {
         Error("MPI process " + std::to_string(get_rank()) + "Failure in MPI_File_write_all");      
      }        
      
      MPI_Type_free(&mysubarray);
      MPI_File_close(&handle);
    }        
    //-------------------------------------------------------------------------------------------------------------------------------------
    




    //-------------------------------------------------------------------------------------------------------------------------------------
    /// This is used in solvers to get a boost vector (such as dofs vector with state_map) from linear solver vector(unique_dof_map)
    auto state_vec(const Epetra_Vector& from)const
    {
      Epetra_Vector app(*state_map_, true);
      app.Import(from, *dofMap_to_stateMap_importer_, Insert);
      boost::numeric::ublas::vector<double> v(state_map_->NumMyElements());
      app.ExtractCopy(&v[0]);
      return v;	  
    }
    //-------------------------------------------------------------------------------------------------------------------------------------
    
    


    private:    

    std::unique_ptr<Epetra_Vector> tmp_;
    int global_dim_, local_dim_, local_start_;
    const Epetra_Map* state_map_ = nullptr;  
    std::unique_ptr<Epetra_Import> dofMap_to_stateMap_importer_ = nullptr;
    std::unique_ptr<Epetra_Import> stateMap_to_tmpMap_importer_ = nullptr;
    std::unique_ptr<Epetra_Import> tmpMap_to_stateMap_importer_ = nullptr;
    std::unique_ptr<Epetra_Import> sharedNodeMap_to_tmpMap_importer_ = nullptr;
    const Epetra_Map* shared_node_map_ = nullptr;    
  };





    
  
}
#endif

