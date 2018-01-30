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

#ifndef __PhysiCell_phenotype_h__
#define __PhysiCell_phenotype_h__

#include <vector>
#include <string>
#include <unordered_map>
#include <map> 

#include "../BioFVM/BioFVM.h" 

using namespace BioFVM; 

namespace PhysiCell{

class Cell;
class Cycle_Model; 
class Phenotype; 

/*
// future use?
class BM_Point
{
 private:
 public:
	std::vector<double> position; 
	std::vector<double> unit_normal; 
		
	BM_Point(); 
	
	double distance_from_point( std::vector<double>& point ); 
	void displacement_from_point( std::vector<double>& point, std::vector<double>& displacement ); 
};
*/

class Phase
{
 public:
	int index; // an internal index for the cycle model
	int code; // a global identifier code 
	std::string name; 
	
	bool division_at_phase_exit; // does this phase trigger division? 
	bool removal_at_phase_exit; // does this phase trigger removal? 
	
	void (*entry_function)( Cell* pCell, Phenotype& phenotype, double dt ); 
	
	Phase(); // done
};

class Phase_Link
{
 public:
	int start_phase_index;
	int end_phase_index; 
	
	bool fixed_duration; 
	
	bool (*arrest_function)( Cell* pCell, Phenotype& phenotype, double dt ); 
		// return true if arrested, false if not 
		
	void (*exit_function)( Cell* pCell, Phenotype& phenotype, double dt );
		// function to be excecuted when completing the phase transition 
	
	Phase_Link(); // done
};

class Cycle_Data
{
 private:
 
	// this maps the end_phase_index to the link index in each 
	// phase_links[i]
	// So, index_inverse_map[i][j] = k, corresponds to 
	// phases[i], phase_links[i][k] (which links from phase i to phase j)
	// transition_rates[i][k] (the transition rate from phase i to phase j) 
	std::vector< std::unordered_map<int,int> > inverse_index_maps; 
	
 public:
	Cycle_Model* pCycle_Model; 

	std::string time_units; 
	
	std::vector< std::vector<double> > transition_rates; 
	
	int current_phase_index; 
	double elapsed_time_in_phase; 
	
	Cycle_Data(); // done 
	
	// return current phase (by reference)
	Phase& current_phase( void ); // done 
	
	// make the data structures consistent with the corresponding cell cycle model
	void sync_to_cycle_model( void ); // done 

	// access the transition rate from phase i to phase j (by reference)
	double& transition_rate(int start_phase_index, int end_phase_index ); // done

	double& exit_rate(int phase_index ); // This returns the first transition rate out of 
		// phase # phase_index. It is only relevant if the phase has only one phase link 
		// (true for many cycle models). 
};

class Cycle_Model
{
 private:
 
	// this maps the end_phase_index to the link index in each 
	// phase_links[i]
	// So, index_inverse_map[i][j] = k, corresponds to 
	// phases[i], phase_links[i][k] (which links from phase i to phase j)
	// transition_rates[i][k] (the transition rate from phase i to phase j)
	std::vector< std::unordered_map<int,int> > inverse_index_maps; 
 
 public:
	std::string name; 
	int code; 

	std::vector<Phase> phases; 
	std::vector< std::vector<Phase_Link> > phase_links; 
	
	int default_phase_index; 
	
	Cycle_Data data; // this will be copied to individual cell agents 

	Cycle_Model(); 
	
	void advance_model( Cell* pCell, Phenotype& phenotype, double dt ); // done 
	
	int add_phase( int code, std::string name ); // done 
	
	int add_phase_link( int start_index, int end_index , 
		bool (*arrest_function)( Cell* pCell, Phenotype& phenotype, double dt ) ); // done 
	int add_phase_link( int start_index, int end_index , double rate , 
		bool (*arrest_function)( Cell* pCell, Phenotype& phenotype, double dte ) ); // done 
	
	int find_phase_index( int code ); // done 
	int find_phase_index( std::string name ); // done 
	
	double& transition_rate( int start_index , int end_index ); // done 
	Phase_Link& phase_link(int start_index,int end_index ); // done 
	
	std::ostream& display( std::ostream& os ); // done 
};

class Cycle
{
 private:
 
 public:
	Cycle_Model* pCycle_Model; 
	Cycle_Data data; 
	
	Cycle();// done 
	
	void advance_cycle( Cell* pCell, Phenotype& phenotype, double dt ); // done 
	
	Cycle_Model& model( void ); // done 
	Phase& current_phase( void ); // done 
	int& current_phase_index( void ); // done 
	
	void sync_to_cycle_model( Cycle_Model& cm ); // done 
};

class Death_Parameters
{
 public:
	std::string time_units; 
 
	double unlysed_fluid_change_rate;
	double lysed_fluid_change_rate; 
	
	double cytoplasmic_biomass_change_rate;
	double nuclear_biomass_change_rate; 
	
	double calcification_rate; 
	
	double relative_rupture_volume; 
	
	Death_Parameters(); // done 
};

class Death
{
 private:
 public:
	std::vector<double> rates; 
	std::vector<Cycle_Model*> models; 
	std::vector<Death_Parameters> parameters; 
	
	bool dead; 
	int current_death_model_index;
	
	Death(); // done 
	
	int add_death_model( double rate, Cycle_Model* pModel );  // done
	int add_death_model( double rate, Cycle_Model* pModel, Death_Parameters& death_parameters); // done 
	
	int find_death_model_index( int code ); // done 
	int find_death_model_index( std::string name ); // done 
	
	bool check_for_death( double dt ); // done
	void trigger_death( int death_model_index ); // done 
	
	Cycle_Model& current_model( void ); // done
	Death_Parameters& current_parameters( void ); // done 
};

class Volume
{
 public:
	//
	// state variables 
	//
	double total;
	double solid;
	double fluid;
	double fluid_fraction; 
	
	double nuclear;
	double nuclear_fluid;
	double nuclear_solid; 

	double cytoplasmic;
	double cytoplasmic_fluid; 
	double cytoplasmic_solid; 
	
	double calcified_fraction;
	
	double cytoplasmic_to_nuclear_ratio;
	
	double rupture_volume; // in volume units 
	
	//
	// a function that can be set by the user. 
	//
	// void (*volume_update_function)( Cell* pCell, Phenotype& phenotype, double dt ); 
	
	//
	// parameters that can be set by users 
	//
	double cytoplasmic_biomass_change_rate; 
	double nuclear_biomass_change_rate; 
	double fluid_change_rate;

	double calcification_rate; 
	
	double target_solid_cytoplasmic;
	double target_solid_nuclear;
	double target_fluid_fraction;
	
	double target_cytoplasmic_to_nuclear_ratio;

	double relative_rupture_volume; 
	// the volume ratio (compared to initial volume at time of death) 
	// at which a cell ruptures / lyses / bursts. 

	//
	// functions 
	//
	Volume(); // done 
	
	void divide( void ); // done 
	void multiply_by_ratio(double); // done 
	
	void update( Cell* pCell, Phenotype& phenotype, double dt ); // done 
};

class Geometry
{
 public:
	double radius; 
	double nuclear_radius; 
	double surface_area; 
	
	double polarity; 
	
	Geometry(); // done 
	
	void update_radius( Cell* pCell, Phenotype& phenotype, double dt ); // done 
	void update_nuclear_radius( Cell* pCell, Phenotype& phenotype, double dt ); // done 
	void update_surface_area( Cell* pCell, Phenotype& phenotype, double dt ); // done 
	
	void update( Cell* pCell, Phenotype& phenotype, double dt ); // done 
};

class Mechanics
{
 public:
	double cell_cell_adhesion_strength; 
	double cell_BM_adhesion_strength;
	double cell_cell_repulsion_strength;
	double cell_BM_repulsion_strength; 
	
	// this is a multiple of the cell (equivalent) radius
	double relative_maximum_adhesion_distance; 
	double maximum_adhesion_distance; // needed? 
	
	Mechanics(); // done 
};

class Motility
{
 public:
	bool is_motile; 
 
	double persistence_time; // mean time to keep going in one direction 
		// before resampling for a new direction. 
	double migration_speed; // migration speed along chosen direction, 
		// in absence of all other adhesive / repulsive forces 
	
	std::vector<double> migration_bias_direction; // a unit vector
		// random motility is biased in this direction (e.g., chemotaxis)
	double migration_bias; // how biased is motility
		// if 0, completely random. if 1, deterministic along the bias vector 
		
	bool restrict_to_2D; 
		// if true, set random motility to 2D only. 
		
	std::vector<double> motility_vector; 
		
	Motility(); // done 
};

class Secretion
{
 private:
 public:
	Microenvironment* pMicroenvironment; 
	
	std::vector<double> secretion_rates; 
	std::vector<double> uptake_rates; 
	std::vector<double> saturation_densities;
	
	// in the default constructor, we'll size to the default microenvironment, if 
	// specified. (This ties to BioFVM.) 
	Secretion(); // done 

	// use this to properly size the secretion parameters to the microenvironment in 
	// pMicroenvironment
	void sync_to_current_microenvironment( void ); // done 
	
	void advance( Basic_Agent* pCell, Phenotype& phenotype , double dt ); 
	
	// use this to properly size the secretion parameters to the microenvironment 
	void sync_to_microenvironment( Microenvironment* pNew_Microenvironment ); // done 
	
	void set_all_secretion_to_zero( void ); // NEW
	void set_all_uptake_to_zero( void ); // NEW
	void scale_all_secretion_by_factor( double factor ); // NEW
	void scale_all_uptake_by_factor( double factor ); // NEW
};

class Cell_Functions
{
 private:
 public:
	Cycle_Model cycle_model; 

	void (*volume_update_function)( Cell* pCell, Phenotype& phenotype , double dt ); // used in cell 
	void (*update_migration_bias)( Cell* pCell, Phenotype& phenotype, double dt ); 
	
	void (*custom_cell_rule)( Cell* pCell, Phenotype& phenotype, double dt ); 
	void (*update_phenotype)( Cell* pCell, Phenotype& phenotype, double dt ); // used in celll
	
	void (*update_velocity)( Cell* pCell, Phenotype& phenotype, double dt ); 
	
	void (*add_cell_basement_membrane_interactions)(Cell* pCell, Phenotype& phenotype, double dt );
	double (*calculate_distance_to_membrane)( Cell* pCell, Phenotype& phenotype, double dt );
	
	void (*set_orientation)(Cell* pCell, Phenotype& phenotype, double dt );
	
	Cell_Functions(); // done 
};

class Phenotype
{
 private:
 public:
	bool flagged_for_division; 
	bool flagged_for_removal; 
 
	Cycle cycle; 
	Death death; 
	Volume volume; 
	Geometry geometry; 
	Mechanics mechanics; 
	Motility motility; 
	Secretion secretion; 
	
	Phenotype(); // done 
	
	void sync_to_functions( Cell_Functions& functions ); // done 
	
	// make sure cycle, death, etc. are synced to the defaults. 
	void sync_to_default_functions( void ); // done 
};

/*
class Microenvironment_Options
{
 private:
 
 public: 
	Microenvironment* pMicroenvironment;
	std::string name; 
 
	std::string time_units; 
	std::string spatial_units; 
	double dx;
	double dy; 
	double dz; 
	
	bool outer_Dirichlet_conditions; 
	std::vector<double> Dirichlet_condition_vector; 
	
	bool simulate_2D; 
	std::vector<double> X_range; 
	std::vector<double> Y_range; 
	std::vector<double> Z_range; 
	
	Microenvironment_Options(); 
	
	bool calculate_gradients; 
};

extern Microenvironment_Options default_microenvironment_options; 
extern Microenvironment microenvironment;

void initialize_microenvironment( void ); 
*/

};

#endif
