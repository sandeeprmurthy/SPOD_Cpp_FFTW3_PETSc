#include <mpi.h>
#include <fftw3-mpi.h>
#include <hdf5.h>
#include <iostream>
#include <map> 
#include <vector>
#include <string>
#include <sstream>
#include <cassert>
#include <unistd.h>
#include <fstream>
#include <iomanip>

#include <numeric>
#include <algorithm>
#include <tuple>

// Include PETSc and SLEPc headers
#include <petscmat.h>
#include <slepceps.h>

// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void CreatePetscVecFromStdVector(const std::vector<double>& std_vec, Vec& petsc_vec, MPI_Comm comm = PETSC_COMM_WORLD) {
    PetscInt size = static_cast<PetscInt>(std_vec.size());
    
    // Create a PETSc vector
    VecCreate(comm, &petsc_vec);
    VecSetSizes(petsc_vec, PETSC_DECIDE, size);
    VecSetFromOptions(petsc_vec);

    // Set values in the PETSc vector
    PetscInt start, end;
    VecGetOwnershipRange(petsc_vec, &start, &end);

    for (PetscInt i = start; i < end; ++i) {
        PetscScalar value = std_vec[i];
        VecSetValue(petsc_vec, i, value, INSERT_VALUES);
    }

    VecAssemblyBegin(petsc_vec);
    VecAssemblyEnd(petsc_vec);
}

// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

std::tuple<std::vector<double>, std::vector<double>, double, int, int, int, int, std::vector<double>>
processHDF5File(const std::string &file_name) {
    try {
        // Open the HDF5 file
        hid_t file_id = H5Fopen(file_name.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
        if (file_id < 0) {
            throw std::runtime_error("Failed to open HDF5 file.");
        }

        // Open 'data' dataset
        hid_t data_id = H5Dopen(file_id, "data", H5P_DEFAULT);
        if (data_id < 0) {
            throw std::runtime_error("Failed to open 'data' dataset.");
        }

        // Get dimensions of the 'data' dataset
        hid_t data_space = H5Dget_space(data_id);
        hsize_t data_dims[2];
        H5Sget_simple_extent_dims(data_space, data_dims, nullptr);
        int nt = data_dims[0]; // Number of snapshots
        int nx = data_dims[1]; // Number of grid points * number of variables
        hsize_t data_size = nt * nx;

        // Read 'data' into a vector
        std::vector<double> data(data_size);
        if (H5Dread(data_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data.data()) < 0) {
            throw std::runtime_error("Failed to read 'data' dataset.");
        }

        // Close 'data' resources
        H5Sclose(data_space);
        H5Dclose(data_id);

        // Open 'grid' dataset
        hid_t grid_id = H5Dopen(file_id, "grid", H5P_DEFAULT);
        if (grid_id < 0) {
            throw std::runtime_error("Failed to open 'grid' dataset.");
        }

        // Get dimensions of the 'grid' dataset
        hid_t grid_space = H5Dget_space(grid_id);
        hsize_t grid_dims[2];
        H5Sget_simple_extent_dims(grid_space, grid_dims, nullptr);
        int ng = grid_dims[0]; // Number of grid points
        int grid_cols = grid_dims[1];
        hsize_t grid_size = ng * grid_cols;

        // Read 'grid' into a vector
        std::vector<double> grid(grid_size);
        if (H5Dread(grid_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, grid.data()) < 0) {
            throw std::runtime_error("Failed to read 'grid' dataset.");
        }

        // Close 'grid' resources
        H5Sclose(grid_space);
        H5Dclose(grid_id);

        // Open 'dt' dataset
        hid_t dt_id = H5Dopen(file_id, "dt", H5P_DEFAULT);
        if (dt_id < 0) {
            throw std::runtime_error("Failed to open 'dt' dataset.");
        }

        double dt;
        if (H5Dread(dt_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &dt) < 0) {
            throw std::runtime_error("Failed to read 'dt' dataset.");
        }

        // Close 'dt' resources
        H5Dclose(dt_id);

        // Calculate derived values
        int nvar = nx / ng;

        // Extract weight_grid (third column of the grid dataset)
        std::vector<double> weight_grid(ng);
        for (int i = 0; i < ng; ++i) {
            weight_grid[i] = grid[i * grid_cols + 2]; // Accessing the 3rd column
        }

        // Normalize weight_grid by its mean
        double mean_weight_grid = std::accumulate(weight_grid.begin(), weight_grid.end(), 0.0) / ng;
        for (auto &value : weight_grid) {
            value /= mean_weight_grid;
        }

        // Create weight_phy (all ones)
        std::vector<double> weight_phy(ng, 1.0);

        // Calculate the final weight (element-wise multiplication)
        std::vector<double> weight(ng);
        std::transform(weight_grid.begin(), weight_grid.end(), weight_phy.begin(), weight.begin(),
                       [](double wg, double wp) { return wg * wp; });

        // Close file resources
        H5Fclose(file_id);

        // Return results as a tuple
        return {data, grid, dt, ng, nt, nx, nvar, weight};

    } catch (const std::exception &e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return {};
    }
}

// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

// Deallocate allocated memory on error
void deallocate_memory(std::vector<double*>& buffers) {
    for (auto& ptr : buffers) {
        if (ptr) {
            delete[] ptr;
            ptr = nullptr;
        }
    }
}

// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

// Function to compute the magnitude of a vector
PetscReal computeMagnitude(const std::vector<PetscScalar>& vec) {
    PetscReal magnitude = 0.0;
    for (const auto& val : vec) {
        magnitude += PetscRealPart(val) * PetscRealPart(val) + PetscImaginaryPart(val) * PetscImaginaryPart(val);
    }
    return std::sqrt(magnitude);
}

// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++}

// Function to estimate required memory in bytes
size_t estimate_memory_requirement(int num_freq_bins, int num_segments, int NUM_VARIABLES, size_t total_local_size) {
    return static_cast<size_t>(num_freq_bins) *
           static_cast<size_t>(num_segments) *
           static_cast<size_t>(NUM_VARIABLES) *
           static_cast<size_t>(total_local_size) * sizeof(PetscScalar);
}

// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

// Function to get available system memory in bytes
size_t get_available_memory() {
#if defined(__linux__)
    // On Linux, read available memory from /proc/meminfo
    long pages = sysconf(_SC_AVPHYS_PAGES);  // Available pages
    long page_size = sysconf(_SC_PAGE_SIZE); // Page size in bytes
    return static_cast<size_t>(pages) * static_cast<size_t>(page_size);
#else
    // For other platforms, assume a safe lower bound (e.g., 150GB)
    return 150L * 1024 * 1024 * 1024; // 150 GB
#endif
}

// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

// Helper function to concatenate arguments
template <typename T>
void concat_args(std::ostringstream& oss, T arg) {
    oss << arg;
}

// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//template <typename T, typename... Args>
//void concat_args(std::ostringstream& oss, T arg, Args... args) {
//    oss << arg;
//    concat_args(oss, args...);
//}

// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//// Stream-like function for rank 0 output
//template <typename... Args>
//void rank_0_cout(int rank, Args... args) {
//    if (rank == 0) {
//        std::ostringstream oss;
//        concat_args(oss, args...); // Concatenate all arguments
//        std::cout << oss.str() << std::endl;
//    }
//}

// Helper function to append arguments to the output stream
template <typename T>
void append_to_stream(std::ostringstream& oss, const T& arg) {
    if constexpr (std::is_same_v<T, double>) {
        oss << std::setprecision(16) << std::fixed << arg;
    } else {
        oss << arg;
    }
}

template <typename T, typename... Args>
void concat_args(std::ostringstream& oss, const T& first, const Args&... rest) {
    append_to_stream(oss, first);
    if constexpr (sizeof...(rest) > 0) {
        oss << " ";  // Add space between arguments
        concat_args(oss, rest...);
    }
}

template <typename... Args>
void rank_0_cout(int rank, Args... args) {
    if (rank == 0) {
        std::ostringstream oss;
        concat_args(oss, args...); // Concatenate all arguments
        std::cout << oss.str() << std::endl;
    }
}

// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

// Utility function to check if file exists
bool file_exists(const std::string& filename) {
    return (access(filename.c_str(), F_OK) != -1);
}

// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

// Function to read data from HDF5 files in parallel
void read_hdf5_parallel(const std::string& filename,
                        const std::vector<std::string>& grid_names,
                        const std::vector<std::string>& variables,
                        std::vector<std::vector<double*>>& data_buffers,
                        const std::vector<ptrdiff_t>& local_n0s,
                        const std::vector<ptrdiff_t>& local_0_starts,
                        const std::vector<ptrdiff_t>& N1s, const std::vector<ptrdiff_t>& N2s,
                        MPI_Comm comm) {
    int rank;
    MPI_Comm_rank(comm, &rank);

    int file_exists_flag = file_exists(filename) ? 1 : 0;
    MPI_Bcast(&file_exists_flag, 1, MPI_INT, 0, comm);
    if (!file_exists_flag) {
        if (rank == 0) {
            std::cerr << "Error: File " << filename << " does not exist." << std::endl;
        }
        MPI_Abort(comm, 1);
    }

    // Set up file access property list with MPI-IO driver
    hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_all_coll_metadata_ops(plist_id, 1);
    H5Pset_coll_metadata_write(plist_id, 1);    
    MPI_Info info;
    MPI_Info_create(&info);
    MPI_Info_set(info, "cb_buffer_size", "16777216"); // Set collective buffer size to 16MB
    MPI_Info_set(info, "cb_nodes", "4");              // Number of aggregator nodes
    // Add more hints as needed based on your system
    H5Pset_fapl_mpio(plist_id, comm, info);
    MPI_Info_free(&info);
    // H5Pset_fapl_mpio(plist_id, comm, MPI_INFO_NULL);
    hid_t file_id = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, plist_id);
    H5Pclose(plist_id);

    if (file_id < 0) {
        std::cerr << "Error: Unable to open file " << filename << " by rank " << rank << std::endl;
        MPI_Abort(comm, 1);
    }

    for (size_t grid_idx = 0; grid_idx < grid_names.size(); ++grid_idx) {
        ptrdiff_t local_n0 = local_n0s[grid_idx];
        ptrdiff_t local_0_start = local_0_starts[grid_idx];
        ptrdiff_t N1 = N1s[grid_idx];
        ptrdiff_t N2 = N2s[grid_idx];        
        for (size_t var_idx = 0; var_idx < variables.size(); ++var_idx) {
            std::string dataset_path = "/PlasCom2/Simulation/Single_Jet_3D/" + grid_names[grid_idx] + "/" + variables[var_idx];
            hid_t dataset_id = H5Dopen(file_id, dataset_path.c_str(), H5P_DEFAULT);

            if (dataset_id < 0) {
                std::cerr << "Error: Unable to open dataset " << dataset_path << " by rank " << rank << std::endl;
                MPI_Abort(comm, 1);
            }

            // Select hyperslab in the file
            hid_t filespace = H5Dget_space(dataset_id);
            hsize_t offset[3] = {static_cast<hsize_t>(local_0_start), 0, 0};
            hsize_t count[3] = {static_cast<hsize_t>(local_n0), static_cast<hsize_t>(N1), static_cast<hsize_t>(N2)};
            H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, count, NULL);

            // Define the memory space
            hsize_t dimsm[3] = {count[0], count[1], count[2]};
            hid_t memspace = H5Screate_simple(3, dimsm, NULL);

            // Set up the collective read property list
            plist_id = H5Pcreate(H5P_DATASET_XFER);
            H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

            // Perform the read operation
            herr_t status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, memspace, filespace, plist_id, data_buffers[grid_idx][var_idx]);

            if (status < 0) {
                std::cerr << "Error: Reading dataset " << dataset_path << " by rank " << rank << std::endl;
                MPI_Abort(comm, 1);
            }

            // Close resources
            H5Pclose(plist_id);
            H5Sclose(memspace);
            H5Sclose(filespace);
            H5Dclose(dataset_id);
        }
    }
    H5Fclose(file_id);
}


// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

// Main function
int main(int argc, char** argv) {
    
    // ####################################################################################
    
    PetscErrorCode ierr = SlepcInitialize(&argc, &argv, NULL, NULL);
    if (ierr) return ierr;

    MPI_Comm comm = PETSC_COMM_WORLD;
    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    rank_0_cout(rank, "MPI Initialized ...");

	// ####################################################################################

    fftw_mpi_init();
    rank_0_cout(rank, "FFTW MPI initialized ...");
    
    // ####################################################################################

    // Define the grid dimensions (adjust as per your data)
    const int NUM_VARIABLES = 1;
    const int NUM_GRIDS = 5;

  // ####################################################################################
    
		int oad = 0; // offset_axail_downstream
	
    std::vector<ptrdiff_t> N0s = {300,  300,     300,      129,      300}; // Adjust these dimensions
    std::vector<ptrdiff_t> N1s = {231,  593, 895-oad, 1098-oad, 1460-oad};
    std::vector<ptrdiff_t> N2s = {70,    73,     126,      129,      196};    
    
    ptrdiff_t total_global_size = 0;
    std::vector<ptrdiff_t> global_sizes(NUM_GRIDS);
    std::vector<ptrdiff_t> start_global_indices(NUM_GRIDS);
    
    for (int grid_index = 0; grid_index < NUM_GRIDS; ++grid_index) {
        ptrdiff_t N0 = N0s[grid_index];
        ptrdiff_t N1 = N1s[grid_index];
        ptrdiff_t N2 = N2s[grid_index];

        ptrdiff_t global_size = N0 * N1 * N2;
        global_sizes[grid_index] = global_size;

        if (grid_index == 0) {
            start_global_indices[grid_index] = 0;
        } else {
            start_global_indices[grid_index] = start_global_indices[grid_index - 1] + global_sizes[grid_index - 1];
        }

        total_global_size += global_size;
    }    

    std::vector<std::string> grid_names = {"grid1", "grid2", "grid3", "grid4", "grid5"};
    //std::vector<std::string> variables = {"rho", "rhoV-1", "rhoV-2", "rhoV-3", "rhoE"};	
    std::vector<std::string> variables = {"pressure"};

	// ####################################################################################

    // Process command line arguments
    int total_time_steps = 1000;  // Default value
    int segment_length = 128;     // Length of each segment (N)
    double overlap = 0.5;         // Overlap percentage between segments
    double dt = 1e-3;
    int ITERATION_START = 373000;
    double DELTA_ITERATION = 200;
    double DELTA_T = DELTA_ITERATION * dt; // Time step size
// TESTING
    const double f0 = 2.0;       // Frequency of the sine wave in Hz

	// ####################################################################################

    for (int i = 1; i < argc; ++i) {
        if (std::string(argv[i]) == "-total_time_steps") {
            total_time_steps = std::atoi(argv[++i]);
        } else if (std::string(argv[i]) == "-segment_length") {
            segment_length = std::atoi(argv[++i]);
        } else if (std::string(argv[i]) == "-overlap") {
            overlap = std::atof(argv[++i]);
        } else if (std::string(argv[i]) == "-delta_t") {
            DELTA_T = std::atof(argv[++i]);
        }
    }



// ####################################################################################
// Setup for the FFT to begin
// ####################################################################################



    int window_shift = static_cast<int>(segment_length * (1 - overlap));
    int num_segments = static_cast<int>(std::floor((total_time_steps - segment_length) / window_shift)) + 1;

    rank_0_cout(rank, "Total number of segments: num_segments = ", num_segments);

	// ####################################################################################

    int required_time_steps = segment_length + (num_segments - 1) * window_shift;
    if (required_time_steps > total_time_steps) {
        if (rank == 0) {
            std::cerr << "Not enough time steps for the given parameters. Adjusting total_time_steps to "
                      << required_time_steps << std::endl;
        }
        total_time_steps = required_time_steps;
    }

	// ####################################################################################

// TESTING
		// List of time step filenames
		std::vector<std::string> time_step_files;
		//    // Adjust the filename pattern as needed
		//    for (int t = 0; t < total_time_steps; ++t) {
		//        std::stringstream ss;
		//        ss << "/scratch1/03119/srmurth2/new_box/Sandeep_PhD_Aeroacoustics/LES_PlasCom2_Simulation_Data/Full_Domain_Single_Jet_TTR_1_With_Downstream_Sponge/PlasCom2_" << t * DELTA_ITERATION + ITERATION_START << ".h5";
		//        time_step_files.push_back(ss.str());
		//    }
		time_step_files.push_back("./PlasCom2_000890400.h5");  
		
	// ####################################################################################    		
	
    // Setup mean data file name   
		std::vector<std::string> mean_file_name;
		mean_file_name.push_back("./mean_flow_TTR_1.h5");

	// ####################################################################################

    // Frequency array (in Hz) based on DELTA_T
    rank_0_cout(rank, "Frequency array (in Hz) based on DELTA_T ...");
    int num_freq_bins = segment_length / 2 + 1;
    std::vector<double> frequencies(num_freq_bins);
    for (int f = 0; f < num_freq_bins; ++f) {
        frequencies[f] = f / (DELTA_T * segment_length);
    }
    
	// ####################################################################################

    // Compute local sizes and total size
    rank_0_cout(rank, "Compute local sizes and total size ...");
    std::vector<ptrdiff_t> local_n0s(NUM_GRIDS);
    std::vector<ptrdiff_t> local_0_starts(NUM_GRIDS);
    std::vector<ptrdiff_t> local_sizes(NUM_GRIDS);
    std::vector<ptrdiff_t> start_indices(NUM_GRIDS);
    ptrdiff_t total_local_size = 0;

    for (int grid_index = 0; grid_index < NUM_GRIDS; ++grid_index) {
        ptrdiff_t N0 = N0s[grid_index];
        ptrdiff_t N1 = N1s[grid_index];
        ptrdiff_t N2 = N2s[grid_index];

        ptrdiff_t local_n0, local_0_start;
        ptrdiff_t alloc_local = fftw_mpi_local_size_3d(N0, N1, N2, comm, &local_n0, &local_0_start);

        ptrdiff_t local_size = local_n0 * N1 * N2;
        local_n0s[grid_index] = local_n0;
        local_0_starts[grid_index] = local_0_start;
        local_sizes[grid_index] = local_size;

        if (grid_index == 0) {
            start_indices[grid_index] = 0;
        } else {
            start_indices[grid_index] = start_indices[grid_index - 1] + local_sizes[grid_index - 1];
        }

        total_local_size += local_size;
    }

	// ####################################################################################

    // Get the available memory
    rank_0_cout(rank, "Get the available memory ...");
  	size_t available_memory = get_available_memory();
  	
	// ####################################################################################  	
    
    // Check if we can allocate memory for time_data
    rank_0_cout(rank, "Check if we can allocate memory for time_data ...");
    size_t required_memory = NUM_VARIABLES * total_local_size * segment_length * sizeof(double);
    if (required_memory > available_memory) {
        std::cerr << "Error: Insufficient memory available for time_data allocation.\n";
        std::cerr << "Required: " << required_memory / (1024 * 1024) << " MB, Available: " << available_memory / (1024 * 1024) << " MB\n";
        return EXIT_FAILURE;
    }

	// Allocate memory for time_data
	rank_0_cout(rank, "Allocate memory for time_data ...");
	std::vector<double*> time_data(NUM_VARIABLES);
    for (int var = 0; var < NUM_VARIABLES; ++var) {
        time_data[var] = new double[total_local_size * segment_length];
        if (!time_data[var]) {
            std::cerr << "Error: Memory allocation failed for time_data[" << var << "].\n";
            deallocate_memory(time_data);
            return EXIT_FAILURE;
        }
    }
    available_memory -= required_memory;    

	// ####################################################################################

	// Allocate data_buffers for each grid and variable
	rank_0_cout(rank, "Allocate data_buffers for each grid and variable...");
	std::vector<std::vector<double*>> data_buffers(NUM_GRIDS, std::vector<double*>(NUM_VARIABLES));

	required_memory = sizeof(double) * total_local_size * NUM_VARIABLES;

	if (required_memory > available_memory) {
		std::cerr << "Error: Insufficient memory available for data_buffers allocation.\n";
		std::cerr << "Required: " << required_memory / (1024 * 1024) << " MB, Available: "
				    << available_memory / (1024 * 1024) << " MB\n";
		deallocate_memory(time_data);
		return EXIT_FAILURE;
	}

	for (int grid_index = 0; grid_index < NUM_GRIDS; ++grid_index) {
		for (int var = 0; var < NUM_VARIABLES; ++var) {
			data_buffers[grid_index][var] = new double[local_sizes[grid_index]];
			if (!data_buffers[grid_index][var]) {
				  std::cerr << "Error: Memory allocation failed for data_buffers[" << grid_index
				            << "][" << var << "].\n";
				  deallocate_memory(time_data);
				  // Deallocate previously allocated data_buffers
				  for (int g = 0; g <= grid_index; ++g) {
				      for (int v = 0; v < NUM_VARIABLES; ++v) {
				          delete[] data_buffers[g][v];
				      }
				  }
				  return EXIT_FAILURE;
			}
		}
	}
	available_memory -= required_memory;
	
	// ####################################################################################	

	// Allocate mean_data_buffers for mean flow grid and variable
	rank_0_cout(rank, "Allocate mean_data_buffers for each grid and variable...");
	std::vector<std::vector<double*>> mean_data_buffers(NUM_GRIDS, std::vector<double*>(NUM_VARIABLES));

	required_memory = sizeof(double) * total_local_size * NUM_VARIABLES;

	if (required_memory > available_memory) {
		std::cerr << "Error: Insufficient memory available for mean_data_buffers allocation.\n";
		std::cerr << "Required: " << required_memory / (1024 * 1024) << " MB, Available: "
				    << available_memory / (1024 * 1024) << " MB\n";
		deallocate_memory(time_data);
		return EXIT_FAILURE;
	}

	for (int grid_index = 0; grid_index < NUM_GRIDS; ++grid_index) {
		for (int var = 0; var < NUM_VARIABLES; ++var) {
			mean_data_buffers[grid_index][var] = new double[local_sizes[grid_index]];
			if (!mean_data_buffers[grid_index][var]) {
				  std::cerr << "Error: Memory allocation failed for mean_data_buffers[" << grid_index
				            << "][" << var << "].\n";
				  deallocate_memory(time_data);
				  // Deallocate previously allocated mean_data_buffers
				  for (int g = 0; g <= grid_index; ++g) {
				      for (int v = 0; v < NUM_VARIABLES; ++v) {
				          delete[] mean_data_buffers[g][v];
				      }
				  }
				  return EXIT_FAILURE;
			}
		}
	}
	available_memory -= required_memory;

	// ####################################################################################

    // Estimate required memory
    rank_0_cout(rank, "Estimate required memory ...");
    required_memory = estimate_memory_requirement(num_freq_bins, num_segments, NUM_VARIABLES, total_local_size);

    if (required_memory > available_memory) {
        std::cerr << "Error: Insufficient memory available to allocate Q_frequencies.\n";
        std::cerr << "Required memory: " << required_memory / (1024 * 1024) << " MB\n";
        std::cerr << "Available memory: " << available_memory / (1024 * 1024) << " MB\n";
        return EXIT_FAILURE;
    }

    // Allocate Q_frequencies as a 1D vector
    // Corrected size calculation to account for grids
    rank_0_cout(rank, "Allocate Q_frequencies as a 1D vector ...");
    std::vector<PetscScalar> Q_frequencies(num_freq_bins * num_segments * NUM_VARIABLES * total_local_size);
    rank_0_cout(rank, "Q_frequencies allocated successfully!");



// ####################################################################################
// Begin the FFT section
// ####################################################################################	


    
    // Prepare window function
    rank_0_cout(rank, "Prepare window function ...");
    std::vector<double> window(segment_length);
    double mean_window = 0.0;
    for (int t = 0; t < segment_length; ++t) {
        window[t] = 0.54 - 0.46 * cos(2 * M_PI * t / (segment_length - 1)); // Hamming window
        mean_window += window[t]/segment_length;
    }

	// ####################################################################################

	// Load the mean_data_buffers 
  	read_hdf5_parallel(mean_file_name[0], grid_names, variables, mean_data_buffers,
                     local_n0s, local_0_starts, N1s, N2s, comm);           
    
	// ####################################################################################

	// Main loop over segments and time steps
    rank_0_cout(rank, "Main loop over segments and time steps ...");
    for (int seg = 0; seg < num_segments; ++seg) {
    
        for (int t = 0; t < segment_length; ++t) {
        	rank_0_cout(rank, "seg = ",seg,": t = ",t," of ",segment_length);
        	
            int time_idx = seg * window_shift + t;
            if (time_idx >= total_time_steps) break;
            
            if(t==0){
            	read_hdf5_parallel(time_step_files[0], grid_names, variables, data_buffers,
                               local_n0s, local_0_starts, N1s, N2s, comm);   
            }

            for (int var = 0; var < NUM_VARIABLES; ++var) {

                for (int grid_index = 0; grid_index < NUM_GRIDS; ++grid_index) {
                    ptrdiff_t local_size = local_sizes[grid_index];
                    ptrdiff_t start_idx = start_indices[grid_index];
                    
                    for (ptrdiff_t i = 0; i < local_size; ++i) {
// TESTING
                        time_data[var][(start_idx + i) * segment_length + t] = (data_buffers[grid_index][var][i] - mean_data_buffers[grid_index][var][i]) * window[t] * sin(2 * M_PI * f0 * time_idx * DELTA_T);
                    
//                        time_data[var][(start_idx + i) * segment_length + t] = (data_buffers[grid_index][var][i] - mean_data_buffers[grid_index][var][i]) * window[t];
     
                    }
                }
            }
        }

        // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        // Perform FFT over time for each spatial point and variable
        rank_0_cout(rank, "Perform FFT over time for each spatial point and variable ...");
        for (int var = 0; var < NUM_VARIABLES; ++var) {
        	rank_0_cout(rank, "seg = ",seg,": var = ",var," of ",NUM_VARIABLES);
        	
            for (ptrdiff_t i = 0; i < total_local_size; ++i) {

								// Allocate input and output arrays
								fftw_complex* out = fftw_alloc_complex(num_freq_bins);  
								double* in = &time_data[var][i * segment_length];
								
								// Setup the plan
								fftw_plan plan_forward = fftw_plan_dft_r2c_1d(segment_length, in, out, FFTW_ESTIMATE);
								assert(plan_forward != nullptr);

								// Execute the FFT
								fftw_execute(plan_forward);

                // Store FFT output for each frequency
                for (int freq = 0; freq < num_freq_bins; ++freq) {
                    PetscScalar value = out[freq][0] + PETSC_i * out[freq][1];
                    size_t index = ((size_t)freq) * ((size_t)num_segments) * ((size_t)total_local_size) * NUM_VARIABLES +
                                   ((size_t)seg) * ((size_t)total_local_size) * NUM_VARIABLES +
                                   ((size_t)i) * NUM_VARIABLES + var;
                    Q_frequencies[index] = (1.0/(mean_window*segment_length))*value;
                }
                
								// Cleanup
								fftw_destroy_plan(plan_forward);
								fftw_free(out);

            }
        }
    }
    
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    // Loop through the data structure
    rank_0_cout(rank, "Loop through the data structure to find peak frequency ...");
    PetscReal max_magnitude = 0.0;
    PetscReal avg_magnitude = 0.0;
    double peak_frequency = 0.0;    
    int var = 0;

    rank_0_cout(rank, "Loop through the data structure to find peak frequency...");

    for (int freq = 0; freq < num_freq_bins; ++freq) {
        double sum_magnitude_local = 0.0;

        // Compute sum_magnitude on each processor
        for (int seg = 0; seg < num_segments; ++seg) {
            for (int i = 0; i < total_local_size; ++i) {
                size_t index = static_cast<size_t>(freq) * num_segments * total_local_size * NUM_VARIABLES +
                               static_cast<size_t>(seg) * total_local_size * NUM_VARIABLES +
                               static_cast<size_t>(i) * NUM_VARIABLES + var;

                std::vector<PetscScalar> vec = { Q_frequencies[index] };
                sum_magnitude_local += computeMagnitude(vec);
            }
        }

        // Reduce sum_magnitude across all processes
        double sum_magnitude_global = 0.0;
        MPI_Reduce(&sum_magnitude_local, &sum_magnitude_global, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

        double avg_magnitude = 0.0;
        if (rank == 0) {
            avg_magnitude = sum_magnitude_global / num_segments;

            // Update the maximum magnitude and peak frequency
            if (avg_magnitude > max_magnitude) {
                max_magnitude = avg_magnitude;
                peak_frequency = frequencies[freq];
            }

            // Print the results
            std::cout << "At frequencies[" << freq << "] = " << frequencies[freq]
                      << ": avg_magnitude = " << avg_magnitude << std::endl;
        }
    }

    if (rank == 0) {
        std::cout << "Peak frequency: " << peak_frequency
                  << " with max magnitude: " << max_magnitude << std::endl;
    }	

//    // Verify that peak frequency matches f0
//    rank_0_cout(rank, "Verify that peak frequency matches f0 ...");
//    double frequency_tolerance = 1e-2;  // Acceptable error in Hz
//    bool success = std::abs(peak_frequency - f0) < frequency_tolerance;

//    if (rank == 0) {
//        if (success) {
//            std::cout << "FFT unit test passed. Peak frequency: " << peak_frequency << " Hz" << std::endl;
//        } else {
//            std::cout << "FFT unit test failed. Expected frequency: " << f0 << " Hz, Detected frequency: "
//                      << peak_frequency << " Hz" << std::endl;
//        }
//    }

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

    for (int var = 0; var < NUM_VARIABLES; ++var) {
        delete[] time_data[var];
    }

    for (int grid_index = 0; grid_index < NUM_GRIDS; ++grid_index) {
        for (int var = 0; var < NUM_VARIABLES; ++var) {
            delete[] data_buffers[grid_index][var];
            delete[] mean_data_buffers[grid_index][var];
        }
    }



// ####################################################################################
// Begin the SPOD calculation
// ####################################################################################	

	// Read and load the weight matrix
	rank_0_cout(rank, "Read and load the weight matrix ...");
	Mat W;
	PetscViewer viewer;
  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"./W.dat",FILE_MODE_READ,&viewer);CHKERRQ(ierr);
  ierr = MatCreate(PETSC_COMM_WORLD,&W);CHKERRQ(ierr);
  ierr = MatSetFromOptions(W);CHKERRQ(ierr);
  ierr = MatLoad(W,viewer);CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);	

	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	// Loop over frequencies to compute SPOD
	rank_0_cout(rank, "Loop over frequencies to compute SPOD ...");
	for (int freq = 0; freq < num_freq_bins; ++freq) {
		
		rank_0_cout(rank, "freq = ",freq);
		
		// Assemble data matrix Q(f) of size ( total_global_size * NUM_VARIABLES) x (num_segments)
		Mat Q;
		MatCreateDense(comm, PETSC_DECIDE, PETSC_DECIDE, total_global_size * NUM_VARIABLES, num_segments, NULL, &Q);
		MatSetUp(Q);
		
		// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

		// Fill Q matrix with complex values
		rank_0_cout(rank, "Fill Q matrix with complex values ...");
		for (int seg = 0; seg < num_segments; ++seg) {

			rank_0_cout(rank, "		seg = ",seg);
		
			for (int grid_index = 0; grid_index < NUM_GRIDS; ++grid_index) {
			
				for(int i_grid=0; i_grid < local_sizes[grid_index]; ++i_grid){

					int i_local = start_indices[grid_index] + i_grid;

					int i_global = start_global_indices[grid_index] + local_0_starts[grid_index]*N1s[grid_index]*N2s[grid_index] + i_grid;

					for (int var = 0; var < NUM_VARIABLES; ++var) {

						PetscInt row = i_global * NUM_VARIABLES + var;

						// Calculate the index into the 1D Q_frequencies vector
						size_t index = ((size_t)freq) * ((size_t)num_segments) * ((size_t)total_local_size) * NUM_VARIABLES +
						((size_t)seg) * ((size_t)total_local_size) * NUM_VARIABLES +
						((size_t)i_local) * NUM_VARIABLES + var;

						PetscScalar value = Q_frequencies[index];
			
						MatSetValue(Q, row, seg, value, INSERT_VALUES);
					}
				}
			}
		}
		rank_0_cout(rank, "Begin Assembly ...");
		MatAssemblyBegin(Q, MAT_FINAL_ASSEMBLY);
		MatAssemblyEnd(Q, MAT_FINAL_ASSEMBLY);
		rank_0_cout(rank, "Finalize assembly ...");
		
//		if(freq==0){
//		  // Write the Q matrix to a binary file
//		  PetscViewer viewer;
//		  PetscViewerBinaryOpen(comm, "Q_matrix.dat", FILE_MODE_WRITE, &viewer);
//		  MatView(Q, viewer);
//		  PetscViewerDestroy(&viewer);
//		}
		
		// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

		// Compute Hermitian transpose of Q
		rank_0_cout(rank, "Compute Hermitian transpose of Q ...");
		Mat Q_H;
		MatHermitianTranspose(Q, MAT_INITIAL_MATRIX, &Q_H);

		// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

		// Compute correlation matrix C(f) = Q^H * W * Q
		rank_0_cout(rank, "Compute correlation matrix C(f) = Q^H * W * Q ...");
		Mat C;
		MatMatMatMult(Q_H, W, Q, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &C);
		
		// Scale the matrix C
		MatScale(C, 1.0/num_segments);

		// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

		// Clean up Q_H as it's no longer needed
		rank_0_cout(rank, "Clean up Q_H as it's no longer needed ...");
		MatDestroy(&Q_H);
		
		// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

		// Perform eigenvalue decomposition on C(f)
		rank_0_cout(rank, "Perform eigenvalue decomposition on C(f) ...");
		EPS eps;
		EPSCreate(comm, &eps);
		EPSSetOperators(eps, C, NULL);
		EPSSetProblemType(eps, EPS_HEP);
		EPSSetFromOptions(eps);
		EPSSolve(eps);
		
		// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

		// Get number of converged eigenvalues
		PetscInt nconv;
		EPSGetConverged(eps, &nconv);
		
		// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

		// Get eigenvalues and eigenvectors
		rank_0_cout(rank, "Get eigenvalues and eigenvectors ...");
		std::vector<PetscScalar> eigenvalues(nconv);
		std::vector<Vec> eigenvectors(nconv);
		for (PetscInt i = 0; i < nconv; ++i) {
		    EPSGetEigenpair(eps, i, &eigenvalues[i], NULL, NULL, NULL);
		    VecCreate(comm, &eigenvectors[i]);
		    VecSetSizes(eigenvectors[i], PETSC_DECIDE, num_segments);
		    VecSetFromOptions(eigenvectors[i]);
		    EPSGetEigenvector(eps, i, eigenvectors[i], NULL);
		}
		
		// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

		// Compute SPOD modes: Phi(f) = Q(f) * V(f) * Lambda^{-1/2}(f)
		// Since Q is large, we will compute modes one by one
		rank_0_cout(rank, "Compute SPOD modes: Phi(f) = Q(f) * V(f) * Lambda^{-1/2}(f) ...");
		for (PetscInt mode = 0; mode < nconv; ++mode) {
		    PetscScalar lambda_sqrt_inv = 1.0 / PetscSqrtScalar(eigenvalues[mode]);
		    Vec phi;
		    VecCreate(comm, &phi);
		    VecSetSizes(phi, PETSC_DECIDE, total_global_size * NUM_VARIABLES);
		    VecSetFromOptions(phi);

		    // Multiply Q * v
		    MatMult(Q, eigenvectors[mode], phi);
		    // Scale by lambda^{-1/2}
		    VecScale(phi, lambda_sqrt_inv);

		    // Save mode to file
		    std::stringstream ss_mode;
		    ss_mode << "./results/spod_modes/spod_mode_freq_" << freq << "_mode_" << mode << ".dat";
		    PetscViewer viewer;
		    PetscViewerBinaryOpen(comm, ss_mode.str().c_str(), FILE_MODE_WRITE, &viewer);
		    VecView(phi, viewer);
		    PetscViewerDestroy(&viewer);

		    VecDestroy(&phi);
		    VecDestroy(&eigenvectors[mode]);
		}
		
		// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

		// Save eigenvalues
		rank_0_cout(rank, "Save eigenvalues ...");
		if (rank == 0) {
		    std::stringstream ss_eigval;
		    ss_eigval << "./results/eigenvalues/eigenvalues_freq_" << freq << ".txt";
		    std::ofstream eigenvalues_file(ss_eigval.str());
		    for (PetscInt i = 0; i < nconv; ++i) {
		        eigenvalues_file << eigenvalues[i] << std::endl;
		    }
		    eigenvalues_file.close();
		}
		
		// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

		// Clean up
		MatDestroy(&Q);
		MatDestroy(&C);
		EPSDestroy(&eps);
	}
	
	// ####################################################################################	
 
    fftw_mpi_cleanup();
    SlepcFinalize();
    return 0;

	// ####################################################################################
	// 									END	
	// ####################################################################################	
	
	rank_0_cout(rank, "Program Done!");    
    
}

