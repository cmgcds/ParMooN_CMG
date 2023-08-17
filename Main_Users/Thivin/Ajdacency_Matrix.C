#include<stdio.h>
#include<cstring>
#include<vector>
#include<cmath>
#include<fstream>
#include<sstream>
#include<iostream>
#include<queue>
#include<string>
#include<algorithm>

using namespace std;

int main()
{

    std:string filename = "/home/thivin/ITC_Output/Siminhale_Particle_PolyDisperse_15LPM/adjacency_matrix_for_siminhale_reflevel_1.txt";

// / read the adjacency matrix file
    // First line contains the Size of rowPtr and colInd
    std::ifstream file(filename);
    std::vector<int> row_ptr;
    std::vector<int> col_ind;

    // Check if object is valid
    if(!file)
    {
        std::cerr << "Cannot open the File : "<<filename<<std::endl;
        std::cout << "Error on file: " << __FILE__ << " on line: " << __LINE__ << std::endl;
        exit(0);
    }

    std::string str;
    int line_count = 0;

    int num_cells = 276577;

    // Read the first line
    while (std::getline(file, str) ) {
        line_count++;
        // std::cout  << str << '\n';
        // Check if the line contains the word "Tetrahedra"
        if (line_count == 1) {
            // convert the string to integer
            std::istringstream iss(str);
            
            int row_ptr_size;
            int col_ind_size;
            iss >> row_ptr_size >> col_ind_size;  // Add the number of cells to the variable num_cells

            std::cout << "row_ptr size: " << row_ptr_size << '\n';
            std::cout << "col_ind size: " << col_ind_size << '\n';
            

            // Resize the row_ptr vector
            row_ptr.resize(row_ptr_size);

            // Resize the col_ind vector
            col_ind.resize(col_ind_size);
        }
        // Read the row pointer Array
        if(line_count == 2)
        {
            std::istringstream iss(str);
            for(int i = 0; i < row_ptr.size(); ++i)
            {
                iss >> row_ptr[i];
            }
        }

        // Read the col_ind Array
        if(line_count == 3)
        {
            std::istringstream iss(str);
            for(int i = 0; i < col_ind.size(); ++i)
            {
                iss >> col_ind[i];
            }
        }
    }

    // Sanity check
        if (num_cells != row_ptr.size() - 1)
        {
            std::cout << "num cells coded : " << num_cells << " row_ptr size: " << row_ptr.size() << '\n';
            std::cout << "Error on file: " << __FILE__ << " on line: " << __LINE__ << std::endl;
            std::cout << "Number of cells and row_ptr size mismatch" << std::endl;
            exit(0);
        }

        // Sanity check, Check the last element of row_ptr with size of col_ind
        if (row_ptr[row_ptr.size() - 1] != col_ind.size())
        {
            std::cout << "Error on file: " << __FILE__ << " on line: " << __LINE__ << std::endl;
            std::cout << "Last element of row_ptr " << row_ptr[row_ptr.size() - 1] << "  and col_ind " << col_ind.size()   <<" size mismatch" << std::endl;
            exit(0);
        }
    
     //-----  read a VTK file 
    std::string filename_vtk = "/home/thivin/ITC_Output/Siminhale_Fluid_Solutions_refined/VTK/unitcube.00001.vtk";

    
    // Read all the lines, untill it reaches the line "CELL_TYPES"
    std::ifstream reference_vtk_file(filename_vtk);
    if (!reference_vtk_file.is_open()) {
        std::cerr << "Error: could not open file" << std::endl;
        return 1;
    }

    std::vector<std::string> lines;
    std::string line;
    while (std::getline(reference_vtk_file, line)) {
        if (line.find("POINT_DATA") != std::string::npos) {
            break;
        }
        lines.push_back(line);
    }

    reference_vtk_file.close();



    //open a file stream to write the output
    std::ofstream output_file;
    output_file.open("/home/thivin/ITC_Output/Siminhale_Particle_PolyDisperse_15LPM/level_neibhours.txt");
    
    int search_depth = 4;
    // Now. lets compute the level neibhours for each cell 
    for (int i = 0 ; i < num_cells ; i++)
    {
        int n_neibhours[3] = {0,0,0};
        
        if(i %200 !=0)
        {
            continue;
        }
        std::queue<std::pair<int, int>> cell_queue;
        std::vector<int> visited_cells;

        std::vector<int> visited_levels;

        cell_queue.push({i,0});


        while(!cell_queue.empty())
        {
            int current_cell = cell_queue.front().first;
            int current_depth = cell_queue.front().second;
            cell_queue.pop();

            // Get the current levels neighbours of the cell 
            int start_index = row_ptr[current_cell];
            int end_index = row_ptr[current_cell+1];


            // For cells in the current level
            for(int index = start_index; index < end_index; index++)
            {
                int neighbourCellNo = col_ind[index];
                
                if(current_depth + 1 < search_depth)
                {
                    // check if the neighbour cell is already in the queue
                    if(std::find(visited_cells.begin(), visited_cells.end(), neighbourCellNo) == visited_cells.end())
                    {
                        cell_queue.push({neighbourCellNo, current_depth + 1});
                        visited_cells.push_back(neighbourCellNo);
                        visited_levels.push_back(current_depth+1);
                        n_neibhours[current_depth]++;
                    }

                }
                    
            }
        }

        // output_file << i << " " << n_neibhours[0] << " " << n_neibhours[1] << " " << n_neibhours[2] << std::endl;
        cout << i << " " << n_neibhours[0] << " " << n_neibhours[1] << " " << n_neibhours[2] << std::endl;
        std::string filename_vtk = "/home/thivin/ITC_Output/Siminhale_Particle_PolyDisperse_15LPM/VTK/zonal." + std::to_string(i) + ".vtk";

        std::ofstream output_file_vtk;
        output_file_vtk.open(filename_vtk);

        //check if the file is open
        if(!output_file_vtk)
        {
            std::cerr << "Cannot open the File : "<<filename_vtk<<std::endl;
            std::cout << "Error on file: " << __FILE__ << " on line: " << __LINE__ << std::endl;
            exit(0);
        }

        // Write the stored lines to the file
        for (const auto& line : lines) {
            output_file_vtk << line << std::endl;
        }

        // Write the cell data
        output_file_vtk << "CELL_DATA " << 276577 << std::endl; // HARDCODED THIVIN
        output_file_vtk << "SCALARS cell_level int 1" << std::endl;
        output_file_vtk << "LOOKUP_TABLE default" << std::endl;

        cout << "[INFO] Writing cell level data" << endl;


        // Write the cell level data
        // Logic : For a cell , write the level as 0 and if the cell exists in the visited_cells vector, then write the level as visited_levels
        for(int j = 0; j < num_cells; j++)
        {
            if(j == i)
            {
                output_file_vtk << 1 << std::endl;
            }
            else
            {
                if(std::find(visited_cells.begin(), visited_cells.end(), j) != visited_cells.end())
                {
                    int index = std::find(visited_cells.begin(), visited_cells.end(), j) - visited_cells.begin();
                    output_file_vtk << visited_levels[index]+1 << std::endl;
                }
                else
                {
                    output_file_vtk << 0 << std::endl;
                }
            }
        }


        output_file_vtk.close();
        cout << "[INFO] Writing cell level data done" << endl;
    }
    
    output_file.close();    

    // Here if 9 is the cell, and 4,5,6,7 are its immediate neibh, then 4,5,6,7 are level 1 neibhours of 9 as per queue 
    // and level 0 neibhours of 9 as per the vector visited_cells and visited_levels

   

    
    


    return 0;
}