#include <hdf5.h>
#include <iostream>
#include <vector>
#include <numeric>
#include <algorithm>
#include <string>

int main() {
    try {
        std::string data_path = "./"; // Update with your actual data path
        std::string file_name = data_path + "/jetLES.h5";

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

        // Output results
        std::cout << "Time step (dt): " << dt << " seconds" << std::endl;
        std::cout << "Number of grid points (ng): " << ng << std::endl;
        std::cout << "Number of snapshots (nt): " << nt << std::endl;
        std::cout << "Number of variables (nvar): " << nvar << std::endl;
        std::cout << "First 5 weights: ";
        for (int i = 0; i < std::min(5, static_cast<int>(weight.size())); ++i) {
            std::cout << weight[i] << " ";
        }
        std::cout << std::endl;

        // Close file resources
        H5Fclose(file_id);

    } catch (const std::exception &e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return -1;
    }

    return 0;
}

