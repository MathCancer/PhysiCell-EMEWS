#include "../core/PhysiCell.h"
#include "../modules/PhysiCell_standard_modules.h" 

using namespace BioFVM; 
using namespace PhysiCell;

class Cancer_Immune_Options
{
 private: 
	
 public:
	// member data: 
	
	// immune cell options 
	double elastic_coefficient; 
	double attachment_rate; 
	double kill_rate; // kill_rate_index
	double attachment_lifetime; 
	double immune_cell_apoptosis_rate; 
	double relative_immune_cell_oxygen_uptake_rate; 
	double relative_immune_cell_adhesion; 
	double relative_immune_cell_repulsion; 
	
	// cancer cell options
	double oncoprotein_mean;
	double oncoprotein_standard_deviation; 
	double immunostimulatory_factor_secretion_rate_live;
	double immunostimulatory_factor_secretion_rate_dead; 
	
	// overall problem options 
	double initial_tumor_radius; 
	double oxygen_conditions; 
	int number_of_immune_cells; 
	double domain_size; 
	
	// immunostimulatory factor 
	double immunostimulatory_diffusion_coefficient; 
	double immunostimulatory_decay_rate; 
	
	// for immune_cell_attempt_attachment
	double oncoprotein_saturation; 
	double oncoprotein_detection_threshold; // oncoprotein_threshold
	double max_attachment_distance; 
	double min_attachment_distance; 
	
	double persistence_time; 
	double migration_speed;
	double migration_bias; 
	
	double break_adhesion_distance; // need to implement -- break the bond if separation is too big
 
	// default constructor with default values; 
	Cancer_Immune_Options();
	
/*

	immune_cell.phenotype.motility.is_motile = true; 
	immune_cell.phenotype.motility.persistence_time = 10.0; 
	immune_cell.phenotype.motility.migration_speed = 1;  
	immune_cell.phenotype.motility.migration_bias = 0.5;
	
	immune_cell.phenotype.mechanics.cell_cell_adhesion_strength *= 0.0;
	immune_cell.phenotype.mechanics.cell_cell_repulsion_strength *= 5.0;
	
	// set functions 
	
	immune_cell.functions.update_phenotype = NULL; 
	immune_cell.functions.custom_cell_rule = immune_cell_rule; 
	immune_cell.functions.update_migration_bias = immune_cell_motility;	
	
	// set custom data values 
	
	immune_cell.custom_data[ "oncoprotein" ] = 0.0; 
	immune_cell.custom_data[ "kill rate" ] = 1.0/15.0; // how often it tries to kill
	immune_cell.custom_data[ "attachment lifetime" ] = 60.00; // how long it can stay attached 
	immune_cell.custom_data[ "attachment rate" ] = 1.0/5.0; // how long it wants to wander before attaching	
*/	
	
	
};

extern Cancer_Immune_Options cancer_immune_options; 

// custom cell phenotype function to scale immunostimulatory factor with hypoxia 
void tumor_cell_phenotype_with_and_immune_stimulation( Cell* pCell, Phenotype& phenotype, double dt ); 

extern Cell_Definition immune_cell; 

void create_immune_cell_type( void ); 

// set the tumor cell properties, then call the function 
// to set up the tumor cells 
void create_cell_types( void );

void setup_tissue(); 

void introduce_immune_cells( void ); 

// set up the microenvironment to include the immunostimulatory factor 
void setup_microenvironment( void );   

std::vector<std::string> cancer_immune_coloring_function( Cell* );

// cell rules for extra elastic adhesion

void attach_cells( Cell* pCell_1, Cell* pCell_2 );
void detach_cells( Cell* pCell_1 , Cell* pCell_2 );

void add_elastic_velocity( Cell* pActingOn, Cell* pAttachedTo , double elastic_constant ); 
void extra_elastic_attachment_mechanics( Cell* pCell, Phenotype& phenotype, double dt );

// immune cell functions for attacking a cell 
Cell* immune_cell_check_neighbors_for_attachment( Cell* pAttacker , double dt ); 
bool immune_cell_attempt_attachment( Cell* pAttacker, Cell* pTarget , double dt ); // only attack if oncoprotein 
bool immune_cell_attempt_apoptosis( Cell* pAttacker, Cell* pTarget, double dt ); 
bool immune_cell_trigger_apoptosis( Cell* pAttacker, Cell* pTarget ); 

void immune_cell_rule( Cell* pCell, Phenotype& phenotype, double dt ); 



void immune_cell_attach( Cell* pAttacker, Cell* pTarget ); // use attach_cells?? 
void immune_cell_detach( Cell* pAttacker, Cell* pTarget ); // use detach_cells ?? 


// immune cell functions for motility 

void immune_cell_motility( Cell* pCell, Phenotype& phenotype, double dt ); 


