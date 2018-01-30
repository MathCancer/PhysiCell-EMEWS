#include "./cancer_immune_3D.h"

Cancer_Immune_Options::Cancer_Immune_Options()
{
	// immune cell options 
	elastic_coefficient = 0.01; // 1/min 
	attachment_rate = 1.0 / 5.0; // wander 5 minutes before attempting to attach 
	kill_rate = 1.0 / 15.0; // 15 minutes to kill a cell in normal circumstances (o.p. = 1)
	attachment_lifetime = 60.0; // stay adhered 60 minutes before detaching and trying again 
	immune_cell_apoptosis_rate = 6.94e-5; // 1.0 / (10.0 * 24.0 * 60.0 ); 
	relative_immune_cell_oxygen_uptake_rate = 0.1; 
	relative_immune_cell_adhesion = 0.0; 
	relative_immune_cell_repulsion = 5.0; 

	// cancer cell options
	oncoprotein_mean = 1.0;
	oncoprotein_standard_deviation = 0.25; 
	immunostimulatory_factor_secretion_rate_live = 10.0;
	immunostimulatory_factor_secretion_rate_dead = 10.0; 

	// overall problem options 
	initial_tumor_radius = 250.0; 
	oxygen_conditions = 38.0; // physioxic conditions 
	number_of_immune_cells = 7500; 
	domain_size = 750.0; 

	// immunostimulatory factor 
	immunostimulatory_diffusion_coefficient = 1e3; 
	immunostimulatory_decay_rate = 0.016; 

	// for immune_cell_attempt_attachment
	oncoprotein_saturation = 2.0; // expression over this level does not help recognition 
	oncoprotein_detection_threshold = 0.5; // expression below this level are not recognized 
	max_attachment_distance = 18.0; // 
	min_attachment_distance = 14.0; // 

	persistence_time = 10.0; 
	migration_speed = 1.0; 
	migration_bias = 0.5; 

	break_adhesion_distance = 20.0; // need to implement -- break the bond if separation is too big

	return; 
}

Cancer_Immune_Options cancer_immune_options; 

Cell_Definition immune_cell; 

void create_immune_cell_type( void )
{
	immune_cell = cell_defaults; 
	
	immune_cell.name = "immune cell";
	immune_cell.type = 1; 

	// turn off proliferation; 
	
	int cycle_start_index = live.find_phase_index( PhysiCell_constants::live ); 
	int cycle_end_index = live.find_phase_index( PhysiCell_constants::live ); 
	
	immune_cell.phenotype.cycle.data.transition_rate(cycle_start_index,cycle_end_index) = 0.0; 	
	
	int apoptosis_index = cell_defaults.phenotype.death.find_death_model_index( PhysiCell_constants::apoptosis_death_model ); 
	int oxygen_i = microenvironment.find_density_index( "oxygen" ); 
	
	// reduce o2 uptake 
	
	immune_cell.phenotype.secretion.uptake_rates[oxygen_i] *= cancer_immune_options.relative_immune_cell_oxygen_uptake_rate; 
		// 0.1; 
	
	// set apoptosis to survive 10 days (on average) 
	
	immune_cell.phenotype.death.rates[apoptosis_index] = cancer_immune_options.immune_cell_apoptosis_rate; 
		//	1.0 / (10.0 * 24.0 * 60.0 ); 
	
	// turn on motility; 
	immune_cell.phenotype.motility.is_motile = true; 
	immune_cell.phenotype.motility.persistence_time = cancer_immune_options.persistence_time; 
		// 10.0; 
	immune_cell.phenotype.motility.migration_speed = cancer_immune_options.migration_speed; 
		// 1;  
	immune_cell.phenotype.motility.migration_bias = cancer_immune_options.migration_bias; 
		// 0.5;
	
	immune_cell.phenotype.mechanics.cell_cell_adhesion_strength *= cancer_immune_options.relative_immune_cell_adhesion; 
		// 0.0;
	immune_cell.phenotype.mechanics.cell_cell_repulsion_strength *= cancer_immune_options.relative_immune_cell_repulsion; 
		// 5.0;
	
	// set functions 
	
	immune_cell.functions.update_phenotype = NULL; 
	immune_cell.functions.custom_cell_rule = immune_cell_rule; 
	immune_cell.functions.update_migration_bias = immune_cell_motility;	
	
	// set custom data values 
	
	immune_cell.custom_data[ "oncoprotein" ] = 0.0; 
	immune_cell.custom_data[ "kill rate" ] = cancer_immune_options.kill_rate; 
		// 1.0/15.0; // how often it tries to kill
	immune_cell.custom_data[ "attachment lifetime" ] = cancer_immune_options.attachment_lifetime; 
		// 60.00; // how long it can stay attached 
	immune_cell.custom_data[ "attachment rate" ] = cancer_immune_options.attachment_rate; 
		// 1.0/5.0; // how long it wants to wander before attaching	
	
	return; 
}

void create_cell_types( void )
{
	// use the same random seed so that future experiments have the 
	// same initial histogram of oncoprotein, even if threading means 
	// that future division and other events are still not identical 
	// for all runs 
	SeedRandom(0); 
	
	// housekeeping 
	
	initialize_default_cell_definition();
	cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment ); 
	
	// turn the default cycle model to live, 
	// so it's easier to turn off proliferation
	
	cell_defaults.phenotype.cycle.sync_to_cycle_model( live ); 
	
	// Make sure we're ready for 2D
	
	cell_defaults.functions.set_orientation = up_orientation; 
	// cell_defaults.phenotype.geometry.polarity = 1.0; 
	cell_defaults.phenotype.motility.restrict_to_2D = false; // true; 
	
	// set to no motility for cancer cells 
	cell_defaults.phenotype.motility.is_motile = false; 
	
	// use default proliferation and death 
	
	int cycle_start_index = live.find_phase_index( PhysiCell_constants::live ); 
	int cycle_end_index = live.find_phase_index( PhysiCell_constants::live ); 
	int apoptosis_index = cell_defaults.phenotype.death.find_death_model_index( PhysiCell_constants::apoptosis_death_model ); 
	
	cell_defaults.parameters.o2_proliferation_saturation = cancer_immune_options.oxygen_conditions; 
		// 38.0;  
	cell_defaults.parameters.o2_reference = cancer_immune_options.oxygen_conditions; 
		// 38.0; 
	
	// set default uptake and secretion 
	
	int oxygen_i = microenvironment.find_density_index( "oxygen" ); 
	
	// oxygen 
	cell_defaults.phenotype.secretion.secretion_rates[oxygen_i] = 0; 
	cell_defaults.phenotype.secretion.uptake_rates[oxygen_i] = 10; 
	cell_defaults.phenotype.secretion.saturation_densities[oxygen_i] = cancer_immune_options.oxygen_conditions; // 38; 

	// immunostimulatory 
	
	int immune_factor_i = microenvironment.find_density_index( "immunostimulatory factor" ); 
	
	cell_defaults.phenotype.secretion.saturation_densities[immune_factor_i] = 1; 

	// set the default cell type to o2-based proliferation with the effect of the 
	// on oncoprotein, and secretion of the immunostimulatory factor 
	
	cell_defaults.functions.update_phenotype = tumor_cell_phenotype_with_and_immune_stimulation; 
	
	// add the extra bit of "attachment" mechanics 
	cell_defaults.functions.custom_cell_rule = extra_elastic_attachment_mechanics; 
	
	cell_defaults.name = "cancer cell"; 
	cell_defaults.type = 0; 
	
	// add custom data 
	
	cell_defaults.custom_data.add_variable( "oncoprotein" , "dimensionless", cancer_immune_options.oncoprotein_mean ); 
//	cell_defaults.custom_data.add_variable( "oncoprotein" , "dimensionless", 1.0 ); 
	
	cell_defaults.custom_data.add_variable( "elastic coefficient" , "1/min" , cancer_immune_options.elastic_coefficient ); 
//	cell_defaults.custom_data.add_variable( "elastic coefficient" , "1/min" , 0.01 ); 
	cell_defaults.custom_data.add_variable( "kill rate" , "1/min" , 0 ); // how often it tries to kill
	cell_defaults.custom_data.add_variable( "attachment lifetime" , "min" , 0 ); // how long it can stay attached 
	cell_defaults.custom_data.add_variable( "attachment rate" , "1/min" ,0 ); // how long it wants to wander before attaching
	
	// create the immune cell type 
	create_immune_cell_type(); 
	
	return; 
}

void setup_microenvironment( void )
{
	// set domain parameters
/*
	default_microenvironment_options.X_range = {-1000, 1000}; 
	default_microenvironment_options.Y_range = {-1000, 1000}; 
	default_microenvironment_options.Z_range = {-1000, 1000}; 
*/	
	double domain_size = cancer_immune_options.domain_size; 
	default_microenvironment_options.X_range = {-domain_size, domain_size}; 
	default_microenvironment_options.Y_range = {-domain_size, domain_size}; 
	default_microenvironment_options.Z_range = {-domain_size, domain_size}; 


	default_microenvironment_options.simulate_2D = false; 

	// gradients are needed for this example 

	default_microenvironment_options.calculate_gradients = true; 

	// add the immunostimulatory factor 

	microenvironment.add_density( "immunostimulatory factor", "dimensionless" ); 

	int oxygen_i = 0; // microenvironment.find_density_index( "oxygen" ); 
	int immune_factor_i = microenvironment.find_density_index( "immunostimulatory factor" ); 

	microenvironment.diffusion_coefficients[immune_factor_i] = cancer_immune_options.immunostimulatory_diffusion_coefficient; 
		// 1e3; 
	microenvironment.decay_rates[immune_factor_i] = cancer_immune_options.immunostimulatory_decay_rate; 
		// .016; 
	
	// let BioFVM use oxygen as the default 

	default_microenvironment_options.use_oxygen_as_first_field = true; 

	// set Dirichlet conditions 

	default_microenvironment_options.outer_Dirichlet_conditions = true;
	default_microenvironment_options.Dirichlet_condition_vector[oxygen_i] = cancer_immune_options.oxygen_conditions; // 38; // physioxic conditions 
	default_microenvironment_options.Dirichlet_condition_vector[immune_factor_i] = 0; 

	default_microenvironment_options.Dirichlet_activation_vector[oxygen_i] = true;  
	default_microenvironment_options.Dirichlet_activation_vector[immune_factor_i] = false;  // no Dirichlet for the immunostimulatory factor 

	initialize_microenvironment(); 	

	return; 
}	

void introduce_immune_cells( void )
{
	double tumor_radius = -9e9;   
	double temp_radius = 0.0; 
	
	// for the loop, deal with the (faster) norm squared 
	for( int i=0; i < (*all_cells).size() ; i++ )
	{
		temp_radius = norm_squared( (*all_cells)[i]->position ); 
		if( temp_radius > tumor_radius )
		{ tumor_radius = temp_radius; }
	}
	// now square root to get to radius 
	tumor_radius = sqrt( tumor_radius ); 
	
	// if this goes wackadoodle, choose cancer_immune_options.initial_tumor_radius 
	if( tumor_radius < cancer_immune_options.initial_tumor_radius )
	{ tumor_radius = cancer_immune_options.initial_tumor_radius; }
	
	std::cout << "current tumor radius: " << tumor_radius << std::endl; 
	
	// now seed immune cells 
	
	int number_of_immune_cells = cancer_immune_options.number_of_immune_cells; // 7500; // 100; // 40; 
	double radius_inner = tumor_radius + 30.0; // 75 // 50; 
	double radius_outer = radius_inner + 75.0; // 100; // 1000 - 50.0; 
	
	double mean_radius = 0.5*(radius_inner + radius_outer); 
	double std_radius = 0.33*( radius_outer-radius_inner)/2.0; 
	
	for( int i=0 ;i < number_of_immune_cells ; i++ )
	{
		double theta = UniformRandom() * 6.283185307179586476925286766559; 
		double phi = acos( 2.0*UniformRandom() - 1.0 );  
		
		double radius = NormalRandom( mean_radius, std_radius ); 
		
		Cell* pCell = create_cell( immune_cell ); 
		pCell->assign_position( radius*cos(theta)*sin(phi), radius*sin(theta)*sin(phi), radius*cos(phi) ); 
	}
	
	return; 
}


std::vector<std::vector<double>> create_cell_sphere_positions(double cell_radius, double sphere_radius)
{
	std::vector<std::vector<double>> cells;
	int xc=0,yc=0,zc=0;
	double x_spacing= cell_radius*sqrt(3);
	double y_spacing= cell_radius*2;
	double z_spacing= cell_radius*sqrt(3);
	
	std::vector<double> tempPoint(3,0.0);
	// std::vector<double> cylinder_center(3,0.0);
	
	for(double z=-sphere_radius;z<sphere_radius;z+=z_spacing, zc++)
	{
		for(double x=-sphere_radius;x<sphere_radius;x+=x_spacing, xc++)
		{
			for(double y=-sphere_radius;y<sphere_radius;y+=y_spacing, yc++)
			{
				tempPoint[0]=x + (zc%2) * 0.5 * cell_radius;
				tempPoint[1]=y + (xc%2) * cell_radius;
				tempPoint[2]=z;
				
				if(sqrt(norm_squared(tempPoint))< sphere_radius)
				{ cells.push_back(tempPoint); }
			}
			
		}
	}
	return cells;
	
}

void setup_tissue( void )
{
	// static int oncoprotein_i = pCell->custom_data.find_variable_index( "oncoprotein" ); 
	static int oncoprotein_i = cell_defaults.custom_data.find_variable_index( "oncoprotein" ); 
	
	// place a cluster of tumor cells at the center 
	
	double cell_radius = cell_defaults.phenotype.geometry.radius; 
	double cell_spacing = 0.95 * 2.0 * cell_radius; 
	
	double tumor_radius = cancer_immune_options.initial_tumor_radius; // 250.0; 
	
	Cell* pCell = NULL; 
	
	std::vector<std::vector<double>> positions = create_cell_sphere_positions(cell_radius,tumor_radius); 
	std::cout << "creating " << positions.size() << " closely-packed tumor cells ... " << std::endl; 
	
	for( int i=0; i < positions.size(); i++ )
	{
		pCell = create_cell(); // tumor cell 
		pCell->assign_position( positions[i] );
		pCell->custom_data[oncoprotein_i] = NormalRandom( cancer_immune_options.oncoprotein_mean, 
			cancer_immune_options.oncoprotein_standard_deviation );
//		pCell->custom_data[oncoprotein_i] = NormalRandom( 1.0, 0.25 );
	}
	
	double sum = 0.0; 
	double min = 9e9; 
	double max = -9e9; 
	for( int i=0; i < all_cells->size() ; i++ )
	{
		double r = (*all_cells)[i]->custom_data[oncoprotein_i]; 
		sum += r;
		if( r < min )
		{ min = r; } 
		if( r > max )
		{ max = r; }
	}
	double mean = sum / ( all_cells->size() + 1e-15 ); 
	// compute standard deviation 
	sum = 0.0; 
	for( int i=0; i < all_cells->size(); i++ )
	{
		sum +=  ( (*all_cells)[i]->custom_data[oncoprotein_i] - mean )*( (*all_cells)[i]->custom_data[oncoprotein_i] - mean ); 
	}
	double standard_deviation = sqrt( sum / ( all_cells->size() - 1.0 + 1e-15 ) ); 
	
	std::cout << std::endl << "Oncoprotein summary: " << std::endl
			  << "===================" << std::endl; 
	std::cout << "mean: " << mean << std::endl; 
	std::cout << "standard deviation: " << standard_deviation << std::endl; 
	std::cout << "[min max]: [" << min << " " << max << "]" << std::endl << std::endl; 
	
	return; 
}

// custom cell phenotype function to scale immunostimulatory factor with hypoxia 
void tumor_cell_phenotype_with_and_immune_stimulation( Cell* pCell, Phenotype& phenotype, double dt )
{
	static int cycle_start_index = live.find_phase_index( PhysiCell_constants::live ); 
	static int cycle_end_index = live.find_phase_index( PhysiCell_constants::live ); 
	
	// static int oncoprotein_i = pCell->custom_data.find_variable_index( "oncoprotein" ); 
	static int oncoprotein_i = cell_defaults.custom_data.find_variable_index( "oncoprotein" ); 
	
	// update secretion rates based on hypoxia 
	
	static int o2_index = microenvironment.find_density_index( "oxygen" ); 
	static int immune_factor_index = microenvironment.find_density_index( "immunostimulatory factor" ); 
	double o2 = pCell->nearest_density_vector()[o2_index];	

/*	
	if( o2 > pCell->parameters.o2_hypoxic_response )
	{
		phenotype.secretion.secretion_rates[immune_factor_index] = 0.0; 
	}
	else
	{
		double hypoxia = ( pCell->parameters.o2_hypoxic_response - o2 ) / ( pCell->parameters.o2_hypoxic_response + 1e-13 ); 
		phenotype.secretion.secretion_rates[ immune_factor_index ] = 10.0 * hypoxia; 	
	}
*/
	// new 
	phenotype.secretion.secretion_rates[immune_factor_index] = cancer_immune_options.immunostimulatory_factor_secretion_rate_live; 
		// 10.0; 
	
	update_cell_and_death_parameters_O2_based(pCell,phenotype,dt);
	
	// if cell is dead, don't bother with future phenotype changes. 
	// set it to secrete the immunostimulatory factor 
	if( phenotype.death.dead == true )
	{
		phenotype.secretion.secretion_rates[immune_factor_index] = cancer_immune_options.immunostimulatory_factor_secretion_rate_dead; 
			// 10; 
		pCell->functions.update_phenotype = NULL; 		
		return; 
	}

	// multiply proliferation rate by the oncoprotein 
	phenotype.cycle.data.transition_rate( cycle_start_index ,cycle_end_index ) *= pCell->custom_data[oncoprotein_i] ; 
	
	return; 
}

std::vector<std::string> cancer_immune_coloring_function( Cell* pCell )
{
	static int oncoprotein_i = pCell->custom_data.find_variable_index( "oncoprotein" ); 
	
	// immune are black
	std::vector< std::string > output( 4, "black" ); 
	
	if( pCell->type == 1 )
	{ 
		output[0] = "lime";
		output[1] = "lime";
		output[2] = "green"; 
		return output;
	} 

	// if I'm under attack, color me 
	if( pCell->state.neighbors.size() > 0 )
	{
		output[0] = "darkcyan"; // orangered // "purple"; // 128,0,128
		output[1] = "black"; // "magenta"; 
		output[2] = "cyan"; // "magenta"; //255,0,255
		return output; 
	}
	
	// live cells are green, but shaded by oncoprotein value 
	if( pCell->phenotype.death.dead == false )
	{
		int oncoprotein = (int) round( 0.5 * pCell->custom_data[oncoprotein_i] * 255.0 ); 
		char szTempString [128];
		sprintf( szTempString , "rgb(%u,%u,%u)", oncoprotein, oncoprotein, 255-oncoprotein );
		output[0].assign( szTempString );
		output[1].assign( szTempString );

		sprintf( szTempString , "rgb(%u,%u,%u)", (int)round(output[0][0]/2.0) , (int)round(output[0][1]/2.0) , (int)round(output[0][2]/2.0) );
		output[2].assign( szTempString );
		
		return output; 
	}

	// if not, dead colors 
	
	if (pCell->phenotype.cycle.current_phase().code == PhysiCell_constants::apoptotic )  // Apoptotic - Red
	{
		output[0] = "rgb(255,0,0)";
		output[2] = "rgb(125,0,0)";
	}
	
	// Necrotic - Brown
	if( pCell->phenotype.cycle.current_phase().code == PhysiCell_constants::necrotic_swelling || 
		pCell->phenotype.cycle.current_phase().code == PhysiCell_constants::necrotic_lysed || 
		pCell->phenotype.cycle.current_phase().code == PhysiCell_constants::necrotic )
	{
		output[0] = "rgb(250,138,38)";
		output[2] = "rgb(139,69,19)";
	}	
	
	return output; 
}

void add_elastic_velocity( Cell* pActingOn, Cell* pAttachedTo , double elastic_constant )
{
	std::vector<double> displacement = pAttachedTo->position - pActingOn->position; 
	axpy( &(pActingOn->velocity) , elastic_constant , displacement ); 
	
	return; 
}

void extra_elastic_attachment_mechanics( Cell* pCell, Phenotype& phenotype, double dt )
{
	for( int i=0; i < pCell->state.neighbors.size() ; i++ )
	{
		add_elastic_velocity( pCell, pCell->state.neighbors[i], pCell->custom_data["elastic coefficient"] ); 
	}

	return; 
}	

void attach_cells( Cell* pCell_1, Cell* pCell_2 )
{
	#pragma omp critical
	{
		
	bool already_attached = false; 
	for( int i=0 ; i < pCell_1->state.neighbors.size() ; i++ )
	{
		if( pCell_1->state.neighbors[i] == pCell_2 )
		{ already_attached = true; }
	}
	if( already_attached == false )
	{ pCell_1->state.neighbors.push_back( pCell_2 ); }
	
	already_attached = false; 
	for( int i=0 ; i < pCell_2->state.neighbors.size() ; i++ )
	{
		if( pCell_2->state.neighbors[i] == pCell_1 )
		{ already_attached = true; }
	}
	if( already_attached == false )
	{ pCell_2->state.neighbors.push_back( pCell_1 ); }

	}

	return; 
}

void detach_cells( Cell* pCell_1 , Cell* pCell_2 )
{
	#pragma omp critical
	{
		bool found = false; 
		int i = 0; 
		while( !found && i < pCell_1->state.neighbors.size() )
		{
			// if cell 2 is in cell 1's list, remove it
			if( pCell_1->state.neighbors[i] == pCell_2 )
			{
				int n = pCell_1->state.neighbors.size(); 
				// copy last entry to current position 
				pCell_1->state.neighbors[i] = pCell_1->state.neighbors[n-1]; 
				// shrink by one 
				pCell_1->state.neighbors.pop_back(); 
				found = true; 
			}
			i++; 
		}
	
		found = false; 
		i = 0; 
		while( !found && i < pCell_2->state.neighbors.size() )
		{
			// if cell 1 is in cell 2's list, remove it
			if( pCell_2->state.neighbors[i] == pCell_1 )
			{
				int n = pCell_2->state.neighbors.size(); 
				// copy last entry to current position 
				pCell_2->state.neighbors[i] = pCell_2->state.neighbors[n-1]; 
				// shrink by one 
				pCell_2->state.neighbors.pop_back(); 
				found = true; 
			}
			i++; 
		}

	}
	
	return; 
}

void immune_cell_motility( Cell* pCell, Phenotype& phenotype, double dt )
{
	// if attached, biased motility towards director chemoattractant 
	// otherwise, biased motility towards cargo chemoattractant 
	
	static int immune_factor_index = microenvironment.find_density_index( "immunostimulatory factor" ); 

	// if not docked, attempt biased chemotaxis 
	if( pCell->state.neighbors.size() == 0 )
	{
		// phenotype.motility.migration_bias = 0.25; 
		phenotype.motility.is_motile = true; 
		
		phenotype.motility.migration_bias_direction = pCell->nearest_gradient(immune_factor_index);	
		normalize( &( phenotype.motility.migration_bias_direction ) );			
	}
	else
	{
		phenotype.motility.is_motile = false; 
	}
	
	return; 
}

Cell* immune_cell_check_neighbors_for_attachment( Cell* pAttacker , double dt )
{
	std::vector<Cell*> nearby = pAttacker->cells_in_my_container(); 
	int i = 0; 
	while( i < nearby.size() )
	{
		// don't try to kill yourself 
		if( nearby[i] != pAttacker )
		{
			if( immune_cell_attempt_attachment( pAttacker, nearby[i] , dt ) )
			{ return nearby[i]; }
		}
		i++; 
	}
	
	return NULL; 
}

bool immune_cell_attempt_attachment( Cell* pAttacker, Cell* pTarget , double dt )
{
	static int oncoprotein_i = pTarget->custom_data.find_variable_index( "oncoprotein" ); 
	static int attach_rate_i = pAttacker->custom_data.find_variable_index( "attachment rate" ); 

	static double oncoprotein_saturation = cancer_immune_options.oncoprotein_saturation; 
		// 2.0; 
	static double oncoprotein_threshold = cancer_immune_options.oncoprotein_detection_threshold; 
		// 0.5; // 0.1; 
	static double oncoprotein_difference = oncoprotein_saturation - oncoprotein_threshold;
	
	static double max_attachment_distance = cancer_immune_options.max_attachment_distance; 
		// 18.0; 
	static double min_attachment_distance = cancer_immune_options.min_attachment_distance; 
		// 14.0; 
	static double attachment_difference = max_attachment_distance - min_attachment_distance; 
	
	if( pTarget->custom_data[oncoprotein_i] > oncoprotein_threshold && pTarget->phenotype.death.dead == false )
	{
		std::vector<double> displacement = pTarget->position - pAttacker->position;
		double distance_scale = norm( displacement ); 
		if( distance_scale > max_attachment_distance )
		{ return false; } 
	
		double scale = pTarget->custom_data[oncoprotein_i];
		scale -= oncoprotein_threshold; 
		scale /= oncoprotein_difference;
		if( scale > 1.0 )
		{ scale = 1.0; } 
		
		distance_scale *= -1.0; 
		distance_scale += max_attachment_distance; 
		distance_scale /= attachment_difference; 
		if( distance_scale > 1.0 )
		{ distance_scale = 1.0; } 
		
		if( UniformRandom() < pAttacker->custom_data[attach_rate_i] * scale * dt * distance_scale )
		{
			// std::cout << "\t attach!" << " " << pTarget->custom_data[oncoprotein_i] << std::endl; 
			attach_cells( pAttacker, pTarget ); 
		}
		
		return true; 
	}
	
	return false; 
}

bool immune_cell_attempt_apoptosis( Cell* pAttacker, Cell* pTarget, double dt )
{
	static int oncoprotein_i = pTarget->custom_data.find_variable_index( "oncoprotein" ); 
	static int apoptosis_model_index = pTarget->phenotype.death.find_death_model_index( "apoptosis" );	
	static int kill_rate_index = pAttacker->custom_data.find_variable_index( "kill rate" ); 
	
	static double oncoprotein_saturation = cancer_immune_options.oncoprotein_saturation;
		// 2.0; 
	static double oncoprotein_threshold = cancer_immune_options.oncoprotein_detection_threshold;
		// 0.5; // 0.1; 
	static double oncoprotein_difference = oncoprotein_saturation - oncoprotein_threshold;
	
	// new 
	if( pTarget->custom_data[oncoprotein_i] < oncoprotein_threshold )
	{ return false; }
	
	// new 
	double scale = pTarget->custom_data[oncoprotein_i];
	scale -= oncoprotein_threshold; 
	scale /= oncoprotein_difference;
	if( scale > 1.0 )
	{ scale = 1.0; } 
	
	
//	if( UniformRandom() < pAttacker->custom_data[kill_rate_index] * pTarget->custom_data[oncoprotein_i] * dt )
	if( UniformRandom() < pAttacker->custom_data[kill_rate_index] * scale * dt )
	{ 
		// std::cout << "\t\t kill!" << " " << pTarget->custom_data[oncoprotein_i] << std::endl; 
		return true; 
	}
	return false; 
}

bool immune_cell_trigger_apoptosis( Cell* pAttacker, Cell* pTarget )
{
	static int apoptosis_model_index = pTarget->phenotype.death.find_death_model_index( "apoptosis" );	
	
	// if the Target cell is already dead, don't bother!
	if( pTarget->phenotype.death.dead == true )
	{ return false; }

	pTarget->start_death( apoptosis_model_index );
	return true; 
}

void immune_cell_rule( Cell* pCell, Phenotype& phenotype, double dt )
{
	static int attach_lifetime_i = pCell->custom_data.find_variable_index( "attachment lifetime" ); 
	
	if( phenotype.death.dead == true )
	{
		// the cell death functions don't automatically turn off custom functions, 
		// since those are part of mechanics. 
		
		// Let's just fully disable now. 
		pCell->functions.custom_cell_rule = NULL; 
		return; 
	}
	
	// if I'm docked
	if( pCell->state.neighbors.size() > 0 )
	{
		extra_elastic_attachment_mechanics( pCell, phenotype, dt );
		
		// attempt to kill my attached cell
		
		bool detach_me = false; 
		
		if( immune_cell_attempt_apoptosis( pCell, pCell->state.neighbors[0], dt ) )
		{
			immune_cell_trigger_apoptosis( pCell, pCell->state.neighbors[0] ); 
			detach_me = true; 
		}
		
		// decide whether to detach 
		
		if( UniformRandom() < dt / ( pCell->custom_data[attach_lifetime_i] + 1e-15 ) )
		{ detach_me = true; }
	
		// detach if pulled too far apart 
		
		std::vector<double> displacement = pCell->position - pCell->state.neighbors[0]->position;
		static double detachment_distance_squared = cancer_immune_options.break_adhesion_distance * cancer_immune_options.break_adhesion_distance; 
		if( norm_squared(displacement) > detachment_distance_squared )
		{ detach_me = true; }
		
		// if I detach, resume motile behavior 
		
		if( detach_me )
		{
			detach_cells( pCell, pCell->state.neighbors[0] ); 
			phenotype.motility.is_motile = true; 
		}
		return; 
	}
	
	// I'm not docked, look for cells nearby and try to docked
	
	// if this returns non-NULL, we're now attached to a cell 
	if( immune_cell_check_neighbors_for_attachment( pCell , dt) )
	{
		// set motility off 
		phenotype.motility.is_motile = false; 
		return; 
	}
	phenotype.motility.is_motile = true; 
	
	return; 
}

