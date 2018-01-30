/*
#############################################################################
# If you use PhysiCell in your project, please cite PhysiCell and the ver-  #
# sion number, such as below:                                               #
#                                                                           #
# We implemented and solved the model using PhysiCell (Version 1.2.1) [1].  #
#                                                                           #
# [1] A Ghaffarizadeh, SH Friedman, SM Mumenthaler, and P Macklin,          #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for            #
#     Multicellular Systems, PLoS Comput. Biol. 2017 (in revision).         #
#     preprint DOI: 10.1101/088773                                          #
#                                                                           #
# Because PhysiCell extensively uses BioFVM, we suggest you also cite       #
#     BioFVM as below:                                                      #
#                                                                           #
# We implemented and solved the model using PhysiCell (Version 1.2.1) [1],  #
# with BioFVM [2] to solve the transport equations.                         #
#                                                                           #
# [1] A Ghaffarizadeh, SH Friedman, SM Mumenthaler, and P Macklin,          #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for            #
#     Multicellular Systems, PLoS Comput. Biol. 2017 (in revision).         #
#     preprint DOI: 10.1101/088773                                          #
#                                                                           #
# [2] A Ghaffarizadeh, SH Friedman, and P Macklin, BioFVM: an efficient     #
#    parallelized diffusive transport solver for 3-D biological simulations,#
#    Bioinformatics 32(8): 1256-8, 2016. DOI: 10.1093/bioinformatics/btv730 #
#                                                                           #
#############################################################################
#                                                                           #
# BSD 3-Clause License (see https://opensource.org/licenses/BSD-3-Clause)   #
#                                                                           #
# Copyright (c) 2015-2017, Paul Macklin and the PhysiCell Project           #
# All rights reserved.                                                      #
#                                                                           #
# Redistribution and use in source and binary forms, with or without        #
# modification, are permitted provided that the following conditions are    #
# met:                                                                      #
#                                                                           #
# 1. Redistributions of source code must retain the above copyright notice, #
# this list of conditions and the following disclaimer.                     #
#                                                                           #
# 2. Redistributions in binary form must reproduce the above copyright      #
# notice, this list of conditions and the following disclaimer in the       #
# documentation and/or other materials provided with the distribution.      #
#                                                                           #
# 3. Neither the name of the copyright holder nor the names of its          #
# contributors may be used to endorse or promote products derived from this #
# software without specific prior written permission.                       #
#                                                                           #
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS       #
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED #
# TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A           #
# PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER #
# OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,  #
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,       #
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR        #
# PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF    #
# LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING      #
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS        #
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.              #
#                                                                           #
#############################################################################
*/

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <ctime>
#include <cmath>
#include <omp.h>
#include <fstream>

#include "./core/PhysiCell.h"
#include "./modules/PhysiCell_standard_modules.h" 

// custom user modules 

#include "./custom_modules/cancer_immune_3D.h" 
	
using namespace BioFVM;
using namespace PhysiCell;

// set number of threads for OpenMP (parallel computing)
int omp_num_threads = 8; // set this to # of CPU cores x 2 (for hyperthreading)

int main( int argc, char* argv[] )
{
	/* Process command line options here */
	
	bool detailed_output = false; 
	
	// time setup 
	std::string time_units = "min"; 
	double t = 0.0; // current simulation time 

	double immune_activation_time = 60 * 24 * 14; // activate immune response at 14 days 
	
	double t_output_interval = 60; // 5; // output once per hour 
	double t_max = immune_activation_time + 60*24*7;  // activation time + one week
	double t_next_output_time = t; 
	int output_index = 0; // used for creating unique output filenames 	
	
	
	if( argc != 10 )
	{
		std::cout << "use: " << std::endl  
		<< "cancer-immune-EMEWS [attachment_rate] [attachment_lifetime] ... " << std::endl  
		
		<< "\t[oncoprotein_detection_threshold] [oncoprotein_saturation]" << std::endl << std::endl  
		
		<< "\t[migration_bias] [kill_rate]" << std::endl << std::endl  
		
		<< "\t[detailed_output (0 or 1)] [omp_num_threads] [random seed value]" << std::endl << std::endl  
		
		<< "Good defaults: " << std::endl 
		<< "cancer-immune-EMEWS  0.2 60.0  0.5 2.0  0.5 0.067  0 8 0" << std::endl << std::endl; 

		return -1; 
	}
	cancer_immune_options.attachment_rate = strtod( argv[1] , NULL ); 			
	cancer_immune_options.attachment_lifetime = strtod( argv[2] , NULL ); 			

	cancer_immune_options.oncoprotein_detection_threshold = strtod( argv[3] , NULL ); 			
	cancer_immune_options.oncoprotein_saturation = strtod( argv[4] , NULL ); 		

	cancer_immune_options.migration_bias = strtod( argv[5] , NULL ); 			
	cancer_immune_options.kill_rate = strtod( argv[6] , NULL ); 			
	
	detailed_output = (bool) strtol( argv[7] , NULL , 10 ); 			
	omp_num_threads = strtol( argv[8] , NULL , 10 ); 
	int rand_seed = strtol( argv[9] , NULL , 10 ); 
	
	// OpenMP setup
	omp_set_num_threads(omp_num_threads);

	// PNRG setup 
	SeedRandom(rand_seed); 

	// int initial_cell_count = strtol( argv[1] , NULL , 10 ); 
	// double initial_oxygen = strtod( argv[2] , NULL ); 	
	
	/* end of processing command line options */ 

	/* Microenvironment setup */ 

	setup_microenvironment(); 

	/* PhysiCell setup */ 
 	
	// set mechanics voxel size, and match the data structure to BioFVM

	double mechanics_voxel_size = 30; 
	Cell_Container* cell_container = create_cell_container_for_microenvironment( microenvironment, mechanics_voxel_size );

	create_cell_types();

	setup_tissue();

	/* Users typically start modifying here. START USERMODS */ 


	/* Users typically stop modifying here. END USERMODS */ 
	
	// set MultiCellDS save options 

	set_save_biofvm_mesh_as_matlab( true ); 
	set_save_biofvm_data_as_matlab( true ); 
	set_save_biofvm_cell_data( true ); 
	set_save_biofvm_cell_data_as_custom_matlab( true );

	// save a simulation snapshot 

	save_PhysiCell_to_MultiCellDS_xml_pugi( "initial" , microenvironment , t ); 

	// save a quick SVG cross section through z = 0, after setting its 
	// length bar to 200 microns 

	PhysiCell_SVG_options.length_bar = 200; 

	// for simplicity, set a pathology coloring function 
	
	std::vector<std::string> (*cell_coloring_function)(Cell*) = cancer_immune_coloring_function;
	
	SVG_plot( "initial.svg" , microenvironment, 0.0 , t, cell_coloring_function );
	
	// set the performance timers 

	BioFVM::RUNTIME_TIC();
	BioFVM::TIC();
	
	std::ofstream report_file ("simulation_report.txt"); 	// create the data log file 
	report_file<<"simulated time\tnum cells\tnum division\tnum death\twall time"<<std::endl;
	
	// main loop 
	
	try 
	{	
		while( t < t_max + 0.1*diffusion_dt )
		{
			// save data if it's time. 
			if(  fabs( t - t_next_output_time ) < 0.01 * diffusion_dt )
			{
				log_output(t, output_index, microenvironment, report_file);
				
				char filename[1024]; 
				if( detailed_output == true )
				{
					sprintf( filename , "output%08u" , output_index ); 
					save_PhysiCell_to_MultiCellDS_xml_pugi( filename , microenvironment , t ); 
				}
				 
				sprintf( filename , "snapshot%08u.svg" , output_index ); 
				SVG_plot( filename , microenvironment, 0.0 , t, cell_coloring_function );
				
				static bool immune_cells_introduced = false; 
				if( t > immune_activation_time - 0.01*diffusion_dt && immune_cells_introduced == false )
				{
					immune_cells_introduced = true; 
					t_output_interval = 3.0; 
					introduce_immune_cells();
				} 
		
				output_index++; 
				t_next_output_time += t_output_interval;
			}
			// update the microenvironment
			microenvironment.simulate_diffusion_decay( diffusion_dt );
			if( default_microenvironment_options.calculate_gradients )
			{ microenvironment.compute_all_gradient_vectors(); }
			
			// run PhysiCell 
			((Cell_Container *)microenvironment.agent_container)->update_all_cells(t);
			
			t += diffusion_dt; 
		}
		log_output(t, output_index, microenvironment, report_file);
		report_file.close();
	}
	catch( const std::exception& e )
	{ // reference to the base of a polymorphic object
		std::cout << e.what(); // information from length_error printed
	}
	
	// save a final simulation snapshot 
	
	save_PhysiCell_to_MultiCellDS_xml_pugi( "final" , microenvironment , t ); 
	SVG_plot( "final.svg" , microenvironment, 0.0 , t, cell_coloring_function );
	
	// timer 
	
	std::cout << std::endl << "Total simulation runtime: " << std::endl; 
	BioFVM::display_stopwatch_value( std::cout , BioFVM::runtime_stopwatch_value() ); 

	return 0; 
}
