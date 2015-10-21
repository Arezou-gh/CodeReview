#include "StokesProblem.hh"
#include <cmath>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <set>
#include <string>
#include <tuple>
#include "SchurComplement.hh"

/*......................................................................................*/

template <int dim>
// Typo: exemple -> example
StokesProblem<dim>::StokesProblem(int exemple): dof_handler(triangulation),  max_degree (10), Tolerance (1e-16)

{
	if (exemple==1)
		exact_solution = new ExactSolutionEx1<dim>();
	else if (exemple==2)
		exact_solution = new ExactSolutionEx2<dim>();
	// I would use else if here, and then make an else that throws an error or simply return from the program
	// Currently, you start Example 3, even is someone chose to call StokesProblem(5)
	else 
		exact_solution = new ExactSolutionEx3<dim>();
	
	for (unsigned int degree=1; degree<=max_degree; ++degree)
	{
		fe_collection.push_back (FESystem<dim>(FE_Q<dim> (degree + 1), dim,
				FE_Q<dim> (degree), 1));

		quadrature_collection.push_back(QGauss<dim> (degree+2));
		face_quadrature_collection.push_back (QGauss<dim-1> (degree+2));

		quadrature_collection_Err.push_back(QGauss<dim> (degree+3));
		face_quadrature_collection_Err.push_back (QGauss<dim-1> (degree+3));
	}
	fe_collection.push_back (FESystem<dim>(FE_Nothing<dim>(), dim,
			FE_Nothing<dim>(), 1));
	quadrature_collection.push_back(QGauss<dim>(1));
	face_quadrature_collection.push_back (QGauss<dim-1>(1));
}

/*.....................................................................................*/
template <int dim>
StokesProblem <dim>::~StokesProblem()
{
 dof_handler.clear();
 delete exact_solution;
}

/*......................................................................................*/
template <int dim>
bool StokesProblem <dim>::decreasing (const std::pair<double,typename hp::DoFHandler<dim>::active_cell_iterator > &i,
                                      const std::pair<double,typename hp::DoFHandler<dim>::active_cell_iterator > &j)
{
 return ((i.first) > (j.first));
}

/*......................................................................................*/
// Generate mesh


template <int dim>
void StokesProblem <dim>::generate_mesh(){

  // If this is for somebody else than you, then better also mention
  // year and title, otherwise I will not find the paper
// example 3.1, paper Morin-Nocheto- Uzawa
 GridGenerator::hyper_cube (triangulation, -1, 1);
 triangulation.refine_global (2);

 std::ofstream out ("grid-hyper_cube.eps");
 GridOut grid_out;
 grid_out.write_eps (triangulation, out);


 // Is this commented code any longer necessary?
/*
 std::vector<Point<dim> > vertices (8);

 vertices [0]=Point<dim> (-1,-1);
 vertices [1]=Point<dim> (0,-1);
 vertices [2]=Point<dim> (-1,0);
 vertices [3]=Point<dim> (0,0);
 vertices [4]=Point<dim> (1,0);
 vertices [5]=Point<dim> (-1,1);
 vertices [6]=Point<dim> (0,1);
 vertices [7]=Point<dim> (1,1);

 const unsigned int n_cells=3;
 std::vector<CellData<dim> > cell(n_cells);
 cell[0].vertices[0]=0;
 cell[0].vertices[1]=1;
 cell[0].vertices[2]=2;
 cell[0].vertices[3]=3;

 cell[1].vertices[0]=2;
 cell[1].vertices[1]=3;
 cell[1].vertices[2]=5;
 cell[1].vertices[3]=6;

 cell[2].vertices[0]=3;
 cell[2].vertices[1]=4;
 cell[2].vertices[2]=6;
 cell[2].vertices[3]=7;

 triangulation.create_triangulation(vertices,cell,SubCellData());
 triangulation.refine_global (1);
 std::ofstream out ("grid-L-Shape.eps");
 GridOut grid_out;
 grid_out.write_eps (triangulation, out);
*/
 
}
//.....................................................................................
// Set the finite element order to the minimum on every cell.
template <int dim>
void StokesProblem <dim>::set_global_active_fe_indices (hp::DoFHandler<dim> &dof_handler)
{
 typename hp::DoFHandler<dim>::active_cell_iterator cell= dof_handler.begin_active(), end_cell = dof_handler.end();
 for (; cell!=end_cell; ++cell)
   cell->set_active_fe_index (0);
}
/*......................................................................................*/
// setup system()

template <int dim>
void StokesProblem <dim>::setup_system(){

	system_matrix.clear();

	dof_handler.distribute_dofs (fe_collection);

	DoFRenumbering::Cuthill_McKee (dof_handler);

	std::vector<unsigned int> block_component (dim+1, 0);
	block_component[dim]=1;
	DoFRenumbering::component_wise(dof_handler, block_component);

	{
		constraints.clear ();
		FEValuesExtractors::Vector velocities(0);
		DoFTools::make_hanging_node_constraints (dof_handler, constraints);
		// Make the spacing uniform in a way you like (either space between every parameter or no space at all
		VectorTools::interpolate_boundary_values (dof_handler,0,*exact_solution,constraints,fe_collection.component_mask(velocities));


		// Since with Dirichlet velocity Bdry condition, pressure will be defined up to a constant, in order to make the
		// solution to be unique, we need to add an additional constraint for Pressure
		// We choose for example the first cell of triangulation and do as follow:
		//
		typename hp::DoFHandler<dim>::active_cell_iterator
		first_cell = dof_handler.begin_active();
		// Is this comment still needed?
		//  std::cout<< "vertex _0  of the first_cell   " << first_cell->vertex(0) << std::endl;
		std::vector<types::global_dof_index> local_dof_indices (first_cell->get_fe().dofs_per_cell);
		first_cell->get_dof_indices(local_dof_indices);
		Point<dim> Pnt_in_ref = first_cell->get_fe().unit_support_point(first_cell->get_fe().component_to_system_index(dim,0));
		MappingQ1<dim> mapping;
		Point<dim> Pnt_in_real = mapping.transform_unit_to_real_cell(first_cell,Pnt_in_ref);

		types::global_dof_index first_pressure_dof = local_dof_indices[first_cell->get_fe().component_to_system_index(dim,0)];
		
		// component_to_system_index: "Compute the shape function for the given vector component and index."
		constraints.add_line (first_pressure_dof);
		Vector<double> values(3);
		 exact_solution->vector_value(Pnt_in_real,values);
		constraints.set_inhomogeneity (first_pressure_dof,values[dim]);


		//****************/

	}
	constraints.close();
	// What is this comment line for? If you want to separate it logically just add an empty line or
	// describe in a comment what comes next
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
	std::vector<types::global_dof_index> dofs_per_block (2);
	DoFTools::count_dofs_per_block (dof_handler, dofs_per_block, block_component);
	const unsigned int n_u = dofs_per_block[0],
			n_p = dofs_per_block[1];

	std::cout   << "   Number of degrees of freedom: "
			<< dof_handler.n_dofs()
			<< " (" << n_u << '+' << n_p << ')'
			<< std::endl;
	{
		BlockCompressedSimpleSparsityPattern csp (2,2);
		csp.block(0,0).reinit (n_u, n_u);
		csp.block(1,0).reinit (n_p, n_u);
		csp.block(0,1).reinit (n_u, n_p);
		csp.block(1,1).reinit (n_p, n_p);
		csp.collect_sizes();
		DoFTools::make_sparsity_pattern (dof_handler, csp, constraints, false);
		sparsity_pattern.copy_from (csp);
	}
	system_matrix.reinit (sparsity_pattern);
	solution.reinit (2);
	solution.block(0).reinit (n_u);
	solution.block(1).reinit (n_p);
	solution.collect_sizes ();
	system_rhs.reinit (2);
	system_rhs.block(0).reinit (n_u);
	system_rhs.block(1).reinit (n_p);
	system_rhs.collect_sizes ();
}
// No need for the following line, and the "assemble system" comment can be removed as well,
// it is the function name after all ;-)
/*......................................................................................*/
// assemble system

template <int dim>
void StokesProblem <dim>::assemble_system () {
	hp::FEValues<dim> hp_fe_values (fe_collection, quadrature_collection, update_values|update_quadrature_points|update_JxW_values|update_gradients);

	FullMatrix<double> local_matrix;
	Vector<double> local_rhs;
	std::vector<types::global_dof_index> local_dof_indices;

	std::vector<Vector<double> >  rhs_values;

	const FEValuesExtractors::Vector velocities (0);
	const FEValuesExtractors::Scalar pressure (dim);

	std::vector<SymmetricTensor<2,dim> > symgrad_phi_u;
	// if grad_phi_u is not needed remove it
	//std::vector<Tensor<2,dim> > grad_phi_u;
	std::vector<double> div_phi_u;
	std::vector<Tensor<1,dim> > phi_u;
	std::vector<double> phi_p;

	typename hp::DoFHandler<dim>::active_cell_iterator
	cell = dof_handler.begin_active(),
	endc = dof_handler.end();
	for (; cell!=endc; ++cell)
	{
		const unsigned int   dofs_per_cell = cell->get_fe().dofs_per_cell;
		local_matrix.reinit (dofs_per_cell, dofs_per_cell);
		local_rhs.reinit (dofs_per_cell);
		local_matrix=0;
		local_rhs=0;

		hp_fe_values.reinit (cell);
                // I am confused by the & signs in the next 2 lines. The functions return stuff by reference,
                // but you assign it to a real variable, so it should be ok to drop them, right?
		const FEValues<dim> &fe_values = hp_fe_values.get_present_fe_values ();
		const std::vector<double>& JxW_values = fe_values.get_JxW_values ();
		const unsigned int n_q_points = fe_values.n_quadrature_points;

		rhs_values.resize(n_q_points, Vector<double>(dim+1));
		rhs_function.vector_value_list (fe_values.get_quadrature_points(), rhs_values);

		symgrad_phi_u.resize(dofs_per_cell);
		// Same here, can you remove the comment?
		//grad_phi_u.resize(dofs_per_cell);
		div_phi_u.resize(dofs_per_cell);
		phi_u.resize (dofs_per_cell);
		phi_p.resize(dofs_per_cell);

		for (unsigned int q=0; q<n_q_points; ++q)
		{
			for (unsigned int k=0; k<dofs_per_cell; ++k)
			{
				symgrad_phi_u[k] = fe_values[velocities].symmetric_gradient (k, q);
				// remove following line?
				//grad_phi_u[k] = fe_values[velocities].gradient (k, q);
				div_phi_u[k] = fe_values[velocities].divergence (k, q);
				phi_u[k] = fe_values[velocities].value (k, q);
				phi_p[k] = fe_values[pressure].value (k, q);
			}

			for (unsigned int i=0; i<dofs_per_cell; ++i)
			{
				for (unsigned int j=0; j<dofs_per_cell; ++j)
				{
                                        local_matrix(i,j) += (2 * (symgrad_phi_u[i] * symgrad_phi_u[j]) - div_phi_u[i] * phi_p[j]
                                         - phi_p[i] * div_phi_u[j])*JxW_values[q];
				} 

				local_rhs(i) += (phi_u[i][0] * rhs_values[q](0) + phi_u[i][1] * rhs_values [q](1)) * JxW_values[q];
			} 
		} 

		local_dof_indices.resize (dofs_per_cell);
		cell->get_dof_indices (local_dof_indices);
		constraints.distribute_local_to_global (local_matrix, local_rhs, local_dof_indices, system_matrix, system_rhs);
	} 
}

// The following two lines can be dropped, and the comment below should be rewritten to full sentences and
// moved to the header file.
/*......................................................................................*/
// Solve

// pressure mass matrix to precondition the Schur complement
//  Later, when solving, we then precondition the Schur complement with M^{1}B^{T} where the matrix
//A is related to the Laplace operator
// Thus, solving with A is a lot more complicated: the matrix is badly conditioned and we know that we
//need many iterations unless we have a very good preconditioner

//in 2d, we use the ultimate preconditioner, namely a direct sparse LU decomposition of the matrix.
//This is implemented using the SparseDirectUMFPACK class that uses the UMFPACK direct solver to compute the decomposition.

template <int dim>
void StokesProblem <dim>::solve ()
{
	SparseDirectUMFPACK A_inverse;
	A_inverse.initialize (system_matrix,
			SparseDirectUMFPACK::AdditionalData());
	A_inverse.vmult (solution, system_rhs);


	constraints.distribute (solution);
	// Are the following lines needed?
/*
	solution.block (1).add (-1.0 * pressure_mean_value ());
	constraints.distribute (solution);	
	*/
}

// I guess this is something you are planning to implement later?
// then it is ok to keep it here for now. Depending on your intention what to
// do with this function, you could also rename it to iterative_solve and simply
// not call it, if it is not working yet.
// Iterative Solver
/*
template <int dim>
void StokesProblem <dim>::solve ()
{
 SparseDirectUMFPACK A_inverse;
 A_inverse.initialize (system_matrix.block(0,0),
     SparseDirectUMFPACK::AdditionalData());
 Vector<double> tmp (solution.block(0).size());
 {
   Vector<double> schur_rhs (solution.block(1).size());
   A_inverse.vmult (tmp, system_rhs.block(0));
   system_matrix.block(1,0).vmult (schur_rhs, tmp);
   schur_rhs -= system_rhs.block(1);

   SchurComplement schur_complement (system_matrix, A_inverse);
   SolverControl solver_control (solution.block(1).size(),
   1e-6*schur_rhs.l2_norm());

  // SolverCG<>    cg (solver_control);
   SolverGMRES<>     gmres (solver_control);
   SparseDirectUMFPACK preconditioner;
   preconditioner.initialize (system_matrix.block(1,1),
	SparseDirectUMFPACK::AdditionalData());
   gmres.solve(schur_complement, solution.block(1), schur_rhs,
			preconditioner);
  // cg.solve (schur_complement, solution.block(1), schur_rhs, preconditioner);
   
   //cout<<" residuals of each step " << solver_control.enable_history_data() << endl;
   constraints.distribute (solution);
   std::cout << "   "
                << solver_control.last_step()
                << " GMRES iterations for Stokes subsystem."
                << std::endl; 
   //  std::cout << "  "
   //<< solver_control.last_step()
   //<< " outer CG Schur complement iterations for pressure"
   //<< std::endl;
 }
 system_matrix.block(0,1).vmult (tmp, solution.block(1));
 tmp *= -1.0;
 tmp += system_rhs.block(0);
 A_inverse.vmult (solution.block(0), tmp);
 constraints.distribute (solution);
 
 // Normalized the solution with keeping pressure's mean value=0.
 solution.block (1).add (-1.0 * pressure_mean_value ());
 constraints.distribute (solution);
}
*/

// This function is currently nowhere called. If you think you need it, keep it around
// otherwise delete it.
/*......................................................................................*/
template <int dim>
double StokesProblem <dim>::pressure_mean_value () const
{

	// get pressure such that satisfies mean value property:
	hp::FEValues<dim> hp_fe_values (fe_collection, quadrature_collection, update_values|update_JxW_values);
	const FEValuesExtractors::Scalar pressure (dim);

	std::vector<double> values;
	double domain_mean_val_p=0;
	double measure_domain=0;
	typename hp::DoFHandler<dim>::active_cell_iterator
	cell = dof_handler.begin_active(),
	endc = dof_handler.end();
	for (; cell!=endc; ++cell)
	{
		hp_fe_values.reinit (cell);

		const FEValues<dim> &fe_values = hp_fe_values.get_present_fe_values ();
		const std::vector<double>& JxW_values = fe_values.get_JxW_values ();
		const unsigned int n_q_points = fe_values.n_quadrature_points;
		values.resize(n_q_points);
		fe_values[pressure].get_function_values(solution, values);
		for (unsigned int q=0; q<n_q_points; ++q)
		{
			domain_mean_val_p += values[q]*JxW_values[q];
			measure_domain += JxW_values[q];
		}
	}

	return domain_mean_val_p / measure_domain;
}

/*......................................................................................*/

template <int dim>
double StokesProblem <dim>::exact_pressure_mean_value () const
{

	hp::FEValues<dim> hp_fe_values (fe_collection, quadrature_collection, update_values | update_quadrature_points|update_JxW_values);
	const FEValuesExtractors::Scalar pressure (dim);

	std::vector<Vector<double> > values;
	double domain_mean_val_p=0;
	double measure_domain=0;
	typename hp::DoFHandler<dim>::active_cell_iterator
	cell = dof_handler.begin_active(),
	endc = dof_handler.end();
	for (; cell!=endc; ++cell)
	{
		hp_fe_values.reinit (cell);

		const FEValues<dim> &fe_values = hp_fe_values.get_present_fe_values ();
		const std::vector<double>& JxW_values = fe_values.get_JxW_values ();
		const unsigned int n_q_points = fe_values.n_quadrature_points;
		values.resize(n_q_points,Vector<double>(dim+1));
		exact_solution->vector_value_list(fe_values.get_quadrature_points(), values);
		for (unsigned int q=0; q<n_q_points; ++q)
		{
			domain_mean_val_p += values[q][dim]*JxW_values[q];
			measure_domain += JxW_values[q];
		}
	}
	return domain_mean_val_p / measure_domain;
}

// Comments that only repeat the function name do not help much ;-).
/*......................................................................................*/
// compute_error

template <int dim>
void StokesProblem <dim>::compute_error (Vector<double> &error_per_cell, Vector<double> &Vect_Pressure_Err, Vector<double> &Vect_grad_Velocity_Err, Vector<double> & Vec_Velocity_Err)
{
	hp::FEValues<dim> hp_fe_values (fe_collection, quadrature_collection_Err, update_values|update_quadrature_points|update_JxW_values|update_gradients|update_hessians);
	const FEValuesExtractors::Vector velocities (0);
	const FEValuesExtractors::Scalar pressure (dim);

	std::vector<double> values;
	std::vector<Tensor<2,dim> > gradients;
	std::vector<Tensor<1,dim> > velocity_values;
	
	std::vector<std::vector<Tensor<1,dim> > > exact_solution_gradients;
	std::vector<Vector<double> > exact_solution_values;

	typename hp::DoFHandler<dim>::active_cell_iterator
	cell = dof_handler.begin_active(),
	endc = dof_handler.end();

	const double mean_exact_pressure = exact_pressure_mean_value ();
	// Remove the following line?
	//std::cout << "*** " << mean_exact_pressure << std::endl;

	unsigned int cell_index=0;
	for (; cell!=endc; ++cell,++cell_index)
	{
		double subtract_p=0.;
		double grad_u_vals=0.;
		double u_vals=0.;
		hp_fe_values.reinit (cell);
		const FEValues<dim> &fe_values = hp_fe_values.get_present_fe_values ();
		const std::vector<double>& JxW_values = fe_values.get_JxW_values ();
		const std::vector<Point<dim> >& quadrature_points = fe_values.get_quadrature_points();
		const unsigned int n_q_points =fe_values.n_quadrature_points;

		velocity_values.resize(n_q_points);
		gradients.resize(n_q_points);
		values.resize(n_q_points);
		exact_solution_gradients.resize(n_q_points , std::vector<Tensor<1,dim> > (dim+1));
		exact_solution_values.resize(n_q_points, Vector<double> (dim+1));
 
		fe_values[velocities].get_function_values(solution, velocity_values);
		fe_values[velocities].get_function_gradients(solution, gradients);
		fe_values[pressure].get_function_values(solution, values);

		exact_solution->vector_gradient_list(quadrature_points, exact_solution_gradients);
		exact_solution->vector_value_list(quadrature_points, exact_solution_values);
 
		double diff_laplace_u_grad_p=0;

		for (unsigned int q=0; q<n_q_points; ++q)
		{
			values[q] -= exact_solution_values[q](dim);
			subtract_p +=values[q]*values[q]* JxW_values[q];

			for (unsigned int i=0; i<dim; ++i)
			{
				velocity_values[q][i]-=exact_solution_values[q](i);
				gradients[q][i]-=exact_solution_gradients[q][i];
			}

			grad_u_vals += gradients[q].norm_square() * JxW_values[q];
			u_vals += velocity_values[q].norm_square() * JxW_values[q];

		} 

		error_per_cell(cell_index) = sqrt (subtract_p + grad_u_vals);

		Vect_Pressure_Err(cell_index)=sqrt(subtract_p);
		Vect_grad_Velocity_Err(cell_index)=sqrt(grad_u_vals);
		Vec_Velocity_Err(cell_index)=sqrt(u_vals);

	}// cell
	
	std::cout<< std::endl;
	double L2_norm_grad_velocity_Err= Vect_grad_Velocity_Err.l2_norm();
	std::cout<< "L2_norm_grad_velocity_Err : "<< L2_norm_grad_velocity_Err << std::endl;
	std::cout<< std::endl;
	double L2_norm_velocity_Err= Vec_Velocity_Err.l2_norm();
	std::cout<< "L2_norm_velocity_Err : "<< L2_norm_velocity_Err << std::endl;
	std::cout<< std::endl;
	std::cout<< std::endl;
	double L2_norm_pressure_Err=Vect_Pressure_Err.l2_norm();
	std::cout<< "L2_norm_pressure_Err : "<< L2_norm_pressure_Err << std::endl;
	std::cout<< std::endl;
	std::cout<< std::endl;
	double L2_norm_total_Err= sqrt (std::pow (L2_norm_grad_velocity_Err,2)+ std::pow (L2_norm_pressure_Err,2));
	std::cout<< "L2_norm of Tottal_ERROR is : "<< L2_norm_total_Err << std::endl;
	std::cout<< std::endl;

}
/*......................................................................................*/
// compute_estimator

// Can you specify in the function name what it is that you are estimating?
// Maybe also rename the argument to estimate_per_cell, I was thinking what est
// might mean, until I saw the function name
template <int dim>
void StokesProblem <dim>::estimate (Vector<double> &est_per_cell)
{

	hp::FEValues<dim> hp_fe_values (fe_collection, quadrature_collection,
        update_values|update_quadrature_points|update_JxW_values|update_gradients|update_hessians);
	hp::FEFaceValues<dim> hp_fe_face_values(fe_collection, face_quadrature_collection, 
        update_JxW_values|update_gradients|update_normal_vectors);
	hp::FEFaceValues<dim> hp_neighbor_face_values(fe_collection, 
        face_quadrature_collection, update_gradients);
	hp::FESubfaceValues<dim> hp_subface_values(fe_collection, face_quadrature_collection, 
        update_JxW_values|update_gradients|update_normal_vectors);
	hp::FESubfaceValues<dim> hp_neighbor_subface_values(fe_collection, face_quadrature_collection, update_gradients);

	std::vector<Tensor<1,dim> > gradients_p;
	std::vector<double> divergences;
	std::vector<Tensor<1,dim> > laplacians;

	std::vector<Tensor<2,dim> > gradients;
	std::vector<Tensor<2,dim> > neighbor_gradients;

	const FEValuesExtractors::Vector velocities (0);
	const FEValuesExtractors::Scalar pressure (dim);

	std::vector<Vector<double> >  rhs_values;

	// Maybe rename this variables to residual_estimate_per_cell
	// and jump_estimate_per_cell, it is simpler to read, and variables
	// usually do not start with capital letters
	Vector<double> res_est_per_cell(triangulation.n_active_cells());
	Vector<double> Jump_est_per_cell(triangulation.n_active_cells());

	typename hp::DoFHandler<dim>::active_cell_iterator
	cell = dof_handler.begin_active(),
	endc = dof_handler.end();
	unsigned int cell_index=0;

	for (; cell!=endc; ++cell,++cell_index)
	{
		hp_fe_values.reinit (cell);
		// same here with the references &
		const FEValues<dim> &fe_values = hp_fe_values.get_present_fe_values ();
		const std::vector<double>& JxW_values = fe_values.get_JxW_values ();
		const unsigned int n_q_points = fe_values.n_quadrature_points;

		rhs_values.resize(n_q_points, Vector<double>(dim+1));
		rhs_function.vector_value_list (fe_values.get_quadrature_points(), rhs_values);

		divergences.resize(n_q_points);
		gradients_p.resize(n_q_points);
		laplacians.resize(n_q_points);

		fe_values[pressure].get_function_gradients(solution, gradients_p);
		fe_values[velocities].get_function_divergences(solution, divergences);
		fe_values[velocities].get_function_laplacians(solution, laplacians);

		// can you give these variables more meaningful names?
                double term1=0; // for the residual term in estimator definition
		double term2=0; // divergence term in estimator definition
		
		for (unsigned int q=0; q<n_q_points; ++q){
			term2 += divergences[q] * divergences[q] * JxW_values[q];

			for (unsigned int i=0; i<2; ++i)
				gradients_p[q][i] -= rhs_values[q](i) + laplacians[q][i];

			term1 += contract(gradients_p[q],gradients_p[q])*JxW_values[q];
		}// q
		// What does q mean?
		res_est_per_cell(cell_index)= pow((cell->diameter())/(cell->get_fe().degree), 2.0 ) * (term1) + term2;

		// No need for so many dots ;-)
		// compute jump_est_per_cell
		// self explaining name possible?
		double term3=0;//jumpped part of the estimator
		for (unsigned int face_number=0; face_number<GeometryInfo<2>::faces_per_cell; ++face_number)

			if ((cell->face(face_number)->at_boundary()==false) && (cell->face(face_number)->has_children() == false)
				&& (cell->face(face_number)->level() == cell->level()))
			{
				const unsigned int q_index = std::max (cell->active_fe_index(),
						cell->neighbor(face_number)->active_fe_index());

				hp_fe_face_values.reinit       (cell,                        face_number,                             q_index);
				hp_neighbor_face_values.reinit (cell->neighbor(face_number), cell->neighbor_of_neighbor(face_number), q_index);

				const FEFaceValues<2> &neighbor_face_values = hp_neighbor_face_values.get_present_fe_values ();
				const FEFaceValues<2> &fe_face_values = hp_fe_face_values.get_present_fe_values ();

				const std::vector<double>& JxW_values = fe_face_values.get_JxW_values ();

				const unsigned int n_face_q_points = fe_face_values.n_quadrature_points;

				gradients.resize(n_face_q_points);
				neighbor_gradients.resize(n_face_q_points);

				neighbor_face_values[velocities].get_function_gradients(solution, neighbor_gradients);
				fe_face_values[velocities].get_function_gradients(solution, gradients);

				std::vector<Tensor<1,dim> > jump_per_face;
				jump_per_face.resize(n_face_q_points);

				// please rename to jump_value, these abbreviations are not worth the saved characters
				double jump_val=0;

				for (unsigned int q=0; q<n_face_q_points; ++q)
				{
					for (unsigned int i=0; i<2; ++i){
						for (unsigned int j=0; j<2; ++j){
						    // Please try to uniformly use spaces in calculations and remove unnecessary parentheses
							jump_per_face[q][i] = (gradients[q][i][j] - neighbor_gradients[q][i][j]) * fe_face_values.normal_vector(q)[j];
						}
					}
					jump_val += contract(jump_per_face[q],jump_per_face[q]) * JxW_values[q];
				}
				term3 += cell->face(face_number)->diameter())/(2.0 * cell->get_fe().degree)*jump_val;
			} 

			// else if the neighbor has children
		// Same here, try to use spaces in the same uniform way throughout your code. You use a space after the opening parentheses
		// but not before the closing one. This makes the code harder to read.
			else if ( (cell->face(face_number)->at_boundary()==false) && (cell->face(face_number)->has_children() == true))
			{
				for (unsigned int subface=0;subface<cell->face(face_number)->n_children();++subface)
				{
					const unsigned int q_index = std::max(cell->neighbor_child_on_subface (face_number, subface)->active_fe_index(), cell->active_fe_index());

					hp_neighbor_face_values.reinit (cell->neighbor_child_on_subface (face_number, subface), cell->neighbor_of_neighbor(face_number), q_index);
					hp_subface_values.reinit (cell,face_number, subface, q_index);

					const FEFaceValues<2> &neighbor_face_values  = hp_neighbor_face_values.get_present_fe_values ();
					const FESubfaceValues<2> &fe_subface_values = hp_subface_values.get_present_fe_values ();


					const std::vector<double>& JxW_values = fe_subface_values.get_JxW_values ();

					const unsigned int n_subface_q_points = fe_subface_values.n_quadrature_points;

					gradients.resize(n_subface_q_points);
					neighbor_gradients.resize(n_subface_q_points);

					neighbor_face_values[velocities].get_function_gradients(solution, neighbor_gradients);
					fe_subface_values[velocities].get_function_gradients(solution, gradients);

					std::vector<Tensor<1,dim> > jump_per_subface;
					jump_per_subface.resize(n_subface_q_points);

							double jump_val=0;
							for (unsigned int q=0; q<n_subface_q_points; ++q)
							{
								for (unsigned int i=0; i<2; ++i){
									for (unsigned int j=0; j<2; ++j){
										jump_per_subface[q][j] += (gradients[q][i][j]- neighbor_gradients[q][i][j])*(fe_subface_values.normal_vector(q)[j]);
									}
								}
								jump_val += contract(jump_per_subface[q],jump_per_subface[q])*(JxW_values[q]);
							}
							term3 +=(cell->face(face_number)->child(subface)->diameter())/(2.0 * cell->get_fe().degree)*jump_val;
				}
			}

		// if the neighbor is coarser

			else if ( (cell->face(face_number)->at_boundary()==false) && (cell->neighbor_is_coarser(face_number)))

			{
				const unsigned int q_index = std::max(cell->active_fe_index(),cell->neighbor(face_number)->active_fe_index());
				hp_fe_face_values.reinit(cell, face_number,q_index);
				hp_neighbor_subface_values.reinit(cell->neighbor(face_number),cell->neighbor_of_coarser_neighbor(face_number).first, cell->neighbor_of_coarser_neighbor(face_number).second,q_index);

				const FEFaceValues<dim> &fe_face_values  = hp_fe_face_values.get_present_fe_values ();
				const FESubfaceValues<dim> &neighbor_subface_values = hp_neighbor_subface_values.get_present_fe_values ();


				const std::vector<double>& JxW_values = fe_face_values.get_JxW_values ();

				const unsigned int n_face_q_points = fe_face_values.n_quadrature_points;

				gradients.resize(n_face_q_points);
				neighbor_gradients.resize(n_face_q_points);

				neighbor_subface_values[velocities].get_function_gradients(solution, neighbor_gradients);
				fe_face_values[velocities].get_function_gradients(solution, gradients);

				std::vector<Tensor<1,dim> > jump_per_face;
				jump_per_face.resize(n_face_q_points);
				double jump_val=0;
				for (unsigned int q=0;
						q<n_face_q_points; ++q)
				{
					for (unsigned int i=0; i<2; ++i){
						for (unsigned int j=0; j<2; ++j){
							jump_per_face[q][i] += (gradients[q][i][j]- neighbor_gradients[q][i][j])*(fe_face_values.normal_vector(q)[j]);
						}
					}
					jump_val += contract(jump_per_face[q],jump_per_face[q])*JxW_values[q];
				}
				term3 +=(cell->face(face_number)->diameter())/(2.0 * cell->get_fe().degree)*jump_val;

			} // else if coarse neighbor

		Jump_est_per_cell(cell_index) = term3;
		est_per_cell(cell_index)=sqrt(Jump_est_per_cell(cell_index)+res_est_per_cell(cell_index));
	}
}

/*......................................................................................*/
// get_layers_of_patch_around_cell

template <int dim>
std::vector<typename hp::DoFHandler<dim>::active_cell_iterator> StokesProblem <dim>::get_patch_around_cell(const typename hp::DoFHandler<dim>::active_cell_iterator &cell)
{
 std::vector<typename hp::DoFHandler<dim>::active_cell_iterator> patch;
 std::set<typename hp::DoFHandler<dim>::active_cell_iterator> cells_done;

 // Uniform use of spaces
 patch.push_back (cell);
 cells_done.insert(cell);
 // then call it layer instead of i. ;-)
 // and add a line: const unsigned int max_layers = 1 and let the loop run to max_layers. That makes it more readable.
 //  i counter for the number of patch layers ... n_layers=1 here (1 level of patch around cell)

 // These loops use a very different indentation scheme than all other functions.
 // We can talk tomorrow about the use of the astyle program that can format your
 // code in a uniform way.
 for (unsigned int i=0; i<1; ++i)
   {
     const unsigned int patch_size = patch.size();
     for (unsigned int j=0; j<patch_size; ++j)
       {
         for (unsigned int face_number=0; face_number< GeometryInfo<dim>::faces_per_cell; ++face_number)
           {
             if (patch[j]->face(face_number)->at_boundary()==false)
               {
                 if (patch[j]->face(face_number)->has_children() == false)
                   {
                     typename hp::DoFHandler<dim>::active_cell_iterator celll = patch[j]->neighbor(face_number);
                     if (cells_done.count(celll)==0)
                       {
                         patch.push_back(celll);
                         cells_done.insert(celll);
                       }
                   }
                 else
                   for (unsigned int subface=0; subface< patch[j]->face(face_number)->n_children(); ++subface)
                     {
                       typename hp::DoFHandler<dim>::active_cell_iterator child_cell = patch[j]->neighbor_child_on_subface (face_number, subface);
                       if (cells_done.count(child_cell)==0)
                         {
                           patch.push_back(child_cell);
                           cells_done.insert(child_cell);
                         }
                     }
               } 
           }
       } 
   } 
 return patch;
}

// consider renaming patch to patches or something else that indicates that it is a vector of several things
template <int dim>
// in ASPECT we usually separate the return type from the name, to keep the lines shorter. Like this:
std::vector<typename hp::DoFHandler<dim>::cell_iterator>
StokesProblem <dim>::get_cells_at_coarsest_common_level (const std::vector<typename hp::DoFHandler<dim>::active_cell_iterator>  &patch)
{
 Assert (patch.size() > 0, ExcMessage("vector containing patch cells should not be an empty vector!"));
 // Remove the following line
 //Assert (patch.size() > 0, ExcInternalError());
 unsigned int min_level = static_cast<unsigned int> (patch[0]->level());
 unsigned int max_level = static_cast<unsigned int> (patch[0]->level());
 for (unsigned int i=0; i<patch.size();++i)
   {
     min_level = std::min (min_level, static_cast<unsigned int> (patch[i]->level()) );
     max_level = std::max (max_level, static_cast<unsigned int> (patch[i]->level()) );
   }

 std::set<typename hp::DoFHandler<dim>::cell_iterator>  uniform_cells;

 typename std::vector<typename hp::DoFHandler<dim>::active_cell_iterator>::const_iterator  patch_c;

 for (patch_c=patch.begin(); patch_c!=patch.end () ; ++patch_c){
	  if (static_cast<unsigned int>((*patch_c)->level()) == min_level)
		  uniform_cells.insert (*patch_c);
	  else
	  {
		  typename hp::DoFHandler<dim>::cell_iterator parent = *patch_c;

		  while (static_cast<unsigned int> (parent->level()) > min_level)
			  parent = parent-> parent();
		  uniform_cells.insert (parent);
	  }
 }

 return std::vector<typename hp::DoFHandler<dim>::cell_iterator> (uniform_cells.begin(), uniform_cells.end());

}	

template <int dim>
// The following argument list is much more readable than the original one
// Also what is DoFHandler_active_cell_iterator? There is a type DoFHandler<dim>::active_cell_iterator, why not use that one?
// Does this compile?
void StokesProblem<dim>::build_triangulation_from_patch(const std::vector<DoFHandler_active_cell_iterator> &patch,
		                                        Triangulation<dim> &local_triangulation,
		                                        unsigned int &level_h_refine,
                                                        unsigned int &level_p_refine,
                                                        std::map<Triangulation_active_cell_iterator,DoFHandler_active_cell_iterator> &patch_to_global_tria_map)
{
	std::vector<DoFHandler_cell_iterator> uniform_cells = 
    get_cells_at_coarsest_common_level (patch);

	level_h_refine=static_cast<unsigned int> (patch[0]->level());
	level_p_refine=static_cast<unsigned int> (patch[0]->active_fe_index());

	local_triangulation.clear();
	std::vector<Point<dim> > vertices;
	const unsigned int n_uniform_cells=uniform_cells.size();
	std::vector<CellData<dim> > cells(n_uniform_cells);
	// Please name counting variables in meaningful ways, like cell_index, vertex_index
	// nobody but you remembers what i,j,k,m was supposed to mean inside the loop below ;-)
	unsigned int k=0;// for enumerating cells
	unsigned int i=0;// for enumerating vertices

	typename std::vector<DoFHandler_cell_iterator>::const_iterator uniform_c;

	for (uniform_c=uniform_cells.begin(); uniform_c!=uniform_cells.end(); ++uniform_c)
	{
		bool repeat_vertex;
		for (unsigned int j=0;  j< GeometryInfo<dim>::vertices_per_cell; ++j)
		{
			Point<dim> position=(*uniform_c)->vertex (j);
			repeat_vertex=false;

			for (unsigned int m=0; m<i; ++m)
			{
				if (position == vertices[m]){
					repeat_vertex=true;
					cells[k].vertices[j]=m ;
					break;
				}
			}
			if (repeat_vertex==false)
			{
				vertices.push_back(position);
				cells[k].vertices[j]=i;
				i=i+1;
			}

		}//for vertices_per_cell
		k++;
	}

	local_triangulation.create_triangulation(vertices,cells,SubCellData());
	Assert (local_triangulation.n_active_cells() == uniform_cells.size(), ExcInternalError());

	local_triangulation.clear_user_flags ();
	unsigned int index=0;
	std::map<Triangulation_cell_iterator, DoFHandler_cell_iterator> patch_to_global_tria_map_tmp;
	// What means coarse_cell_t?
	for (Triangulation_cell_iterator coarse_cell_t = local_triangulation.begin(); 
      coarse_cell_t != local_triangulation.end(); ++coarse_cell_t, ++index)
	{
		patch_to_global_tria_map_tmp.insert (std::make_pair(coarse_cell_t, uniform_cells[index]));
		AssertThrow ((std::fabs(coarse_cell_t->center()(0) - uniform_cells[index]->center()(0))<1e-16 && 
          std::fabs(coarse_cell_t->center()(1) - uniform_cells[index]->center()(1)) <1e-16), 
        ExcInternalError());
	}
	
	bool refinement_necessary;
	do
	{
		refinement_necessary = false;
		for (Triangulation_active_cell_iterator cell_tt = local_triangulation.begin_active(); 
        cell_tt != local_triangulation.end(); ++cell_tt)
		{
			if (patch_to_global_tria_map_tmp[cell_tt]->has_children())
			{
				cell_tt -> set_refine_flag();
				refinement_necessary = true;
			}
			else for (unsigned int i=0; i<patch.size(); ++i)
      {
				if (patch_to_global_tria_map_tmp[cell_tt]==patch[i])
        {
					cell_tt->set_user_flag();
					break;
        }
			}
		}

		if (refinement_necessary)
		{
			local_triangulation.execute_coarsening_and_refinement ();

			for (Triangulation_cell_iterator cell_ttt = local_triangulation.begin(); 
          cell_ttt != local_triangulation.end(); ++cell_ttt)
			{

        if(patch_to_global_tria_map_tmp.find(cell_ttt)!=patch_to_global_tria_map_tmp.end())
        {
          if (cell_ttt-> has_children())
          {
            // Note: Since the cell got children, then it should not be in the map anymore...
            // children may be added into the map, instead

            // these children may not yet be in the map
            for (unsigned int c=0; c< cell_ttt ->n_children(); ++c)
            {
              if (patch_to_global_tria_map_tmp.find(cell_ttt->child(c)) == 
                  patch_to_global_tria_map_tmp.end())
              {
                patch_to_global_tria_map_tmp.insert (std::make_pair(cell_ttt ->child(c), 
                      patch_to_global_tria_map_tmp[cell_ttt]->child(c)));

                AssertThrow( (std::fabs (cell_ttt ->child(c)->center()(0) - 
                        patch_to_global_tria_map_tmp[cell_ttt]->child(c)->center()(0)) < 1e-16 && 
                      std::fabs (cell_ttt ->child(c)->center()(1) - 
                        patch_to_global_tria_map_tmp[cell_ttt]->child(c)->center()(1)) < 1e-16), 
                    ExcInternalError());
              }
            }
            patch_to_global_tria_map_tmp.erase(cell_ttt);
          }
        }
			}
		}

	}
	while (refinement_necessary);

	typename std::map<Triangulation_cell_iterator,DoFHandler_cell_iterator>::iterator map_tmp_it =
    patch_to_global_tria_map_tmp.begin(),map_tmp_end = patch_to_global_tria_map_tmp.end();

	for (; map_tmp_it!=map_tmp_end; ++map_tmp_it)
		patch_to_global_tria_map[map_tmp_it->first] = map_tmp_it->second;
}

/*......................................................................................................................................*/
//set_active_fe_indices for cells on and out of each patch


// mark cells in the "local_triangulation" that exist in the patch: go
// through all cells in the "local_triangulation" and see whether the
// corresponding cell in the global triangulation is part of
// the 'patch' list of cells
template <int dim>
void StokesProblem<dim>::set_active_fe_indices(hp::DoFHandler<dim> &local_dof_handler,
    std::map<Triangulation_active_cell_iterator, DoFHandler_active_cell_iterator> 
    &patch_to_global_tria_map)
{
 DoFHandler_active_cell_iterator patch_cell = local_dof_handler.begin_active(), 
          end_patch_cell = local_dof_handler.end();
 for (; patch_cell!=end_patch_cell; ++patch_cell)
 {
   if (patch_cell->user_flag_set()==true)
   {
     DoFHandler_active_cell_iterator global_cell = patch_to_global_tria_map[patch_cell];

     patch_cell->set_active_fe_index(global_cell->active_fe_index());
   }
   else if (patch_cell->user_flag_set()==false)
   {
     // which assigns FE_Nothing for the cells out of patch
     patch_cell->set_active_fe_index (max_degree);
   }
   else
     Assert (false, ExcNotImplemented());
 }
}

//.................................................................................................................................
template <int dim>
void StokesProblem <dim>:: patch_output (unsigned int patch_number, const unsigned int cycle, hp::DoFHandler<dim> &local_dof_handler, BlockVector<double> &local_solu)
{

 std::vector<std::string> solution_names (dim, "patch_velocity");
 solution_names.push_back ("patch_pressure");

 std::vector<DataComponentInterpretation::DataComponentInterpretation> data_component_interpretation
 (dim, DataComponentInterpretation::component_is_part_of_vector);
 data_component_interpretation.push_back (DataComponentInterpretation::component_is_scalar);

 DataOut<dim,hp::DoFHandler<dim> > patch_data_out;

 patch_data_out.attach_dof_handler (local_dof_handler);


 patch_data_out.add_data_vector (local_solu, solution_names, DataOut<dim,hp::DoFHandler<dim> >::type_dof_data, data_component_interpretation);
 patch_data_out.build_patches ();

 std::string filename = "patch_solution-" +
     Utilities::int_to_string (cycle, 2) +
     +"-"+Utilities::int_to_string (patch_number, 2) +".vtu";
 std::ofstream output (filename.c_str());
 patch_data_out.write_vtu (output);
}
/*......................................................................................................................................*/
//Compute h_convergence_estimator   &   h_workload_number  for each patch around cell

template <int dim>
void StokesProblem <dim>::h_patch_conv_load_no ( const unsigned int cycle , double &h_convergence_est_per_cell,
		unsigned int &h_workload_num, const typename hp::DoFHandler<dim>::active_cell_iterator &cell,
		 unsigned int & patch_number)
{

 Triangulation<dim> local_triangulation;
 unsigned int level_h_refine;
 unsigned int level_p_refine;
 std::map<typename Triangulation<dim>::active_cell_iterator, typename hp::DoFHandler<dim>::active_cell_iterator> patch_to_global_tria_map;
 
 std::vector<typename hp::DoFHandler<dim>::active_cell_iterator> patch = get_patch_around_cell(cell);
 build_triangulation_from_patch (patch, local_triangulation, level_h_refine, level_p_refine, patch_to_global_tria_map);
 hp::DoFHandler<dim> local_dof_handler(local_triangulation);

 set_active_fe_indices (local_dof_handler,patch_to_global_tria_map );
 local_dof_handler.distribute_dofs (fe_collection);

 DoFRenumbering::Cuthill_McKee (local_dof_handler);

 h_convergence_est_per_cell=0.;
 double h_solu_norm_per_patch=0.;

 ConstraintMatrix constraints_patch;
 BlockSparsityPattern sparsity_pattern_patch;

 std::vector<unsigned int> block_component_patch (dim+1, 0);
 block_component_patch[dim]=1;
 DoFRenumbering::component_wise (local_dof_handler, block_component_patch);

 std::vector<types::global_dof_index> dofs_per_block_patch (2);
 DoFTools::count_dofs_per_block (local_dof_handler, dofs_per_block_patch, block_component_patch);

// local_solution?
 BlockVector<double> local_solu (dofs_per_block_patch);

 // Here we are trying to project the values of the global vector "solution" into vector "local_solu" which is
 // solution over patch cells corresponding to cell "K".
 // We actually need these projected solutions in order to construct the residual terms and let it to be the
 // right-hand side of our local variational problems (The solution of these local problems are in fact the Ritz-Representation 
 // of the residuals on each patch) 

 // what is cl?
 typename hp::DoFHandler<dim>::active_cell_iterator patch_cl= local_dof_handler.begin_active(), end_patch_cl = local_dof_handler.end();
 for (; patch_cl !=end_patch_cl; ++patch_cl)
 {
	  const unsigned int   dofs_per_cl = patch_cl->get_fe().dofs_per_cell;
	  // we check if the corresponding finite element for this cell is not 'FE_Nothing!' and it takes usual finite element.
	  if (dofs_per_cl!=0)
	  {
		  Vector<double> local_solution_values(dofs_per_cl);
		  typename hp::DoFHandler<dim>::active_cell_iterator global_cell = patch_to_global_tria_map[patch_cl];

		  global_cell-> get_dof_values (solution,local_solution_values);
		  patch_cl->set_dof_values (local_solution_values, local_solu);
	  }
 }

//......................  Solution Transfer .......................  h_refinement of patch cells ................................. //

 bool need_to_refine = false;

 // cc?
 typename hp::DoFHandler<dim>::active_cell_iterator patch_cc= local_dof_handler.begin_active(), end_patch_cc = local_dof_handler.end();
 for (; patch_cc!=end_patch_cc; ++patch_cc)
 {
	  if (static_cast<unsigned int> (patch_cc->level()) <  (level_h_refine+1) )
	  {
		  need_to_refine = true;
		  patch_cc->set_refine_flag();
	  }
 }

 if (need_to_refine == true)
 {
	  // user flags will be overwritten by
	  // execute_coarsening_and_refinement. save their values into
	  // the material_id, since that one not only survives
	  // refinement, but is inherited by the children
	  for (typename hp::DoFHandler<dim>::cell_iterator cell_ttt = local_dof_handler.begin(); cell_ttt != local_dof_handler.end(); ++cell_ttt)
	    if (cell_ttt->user_flag_set())
	      cell_ttt->set_material_id (1);
	    else
	      cell_ttt->set_material_id (0);

	  local_triangulation.prepare_coarsening_and_refinement ();
	  SolutionTransfer<dim,BlockVector<double>, hp::DoFHandler<dim>> solution_transfer(local_dof_handler);
	  solution_transfer.prepare_for_pure_refinement();


	  local_triangulation.execute_coarsening_and_refinement ();

	  // get user flags back out of the material_id field
	  for (typename hp::DoFHandler<dim>::cell_iterator cell_ttt = local_dof_handler.begin(); cell_ttt != local_dof_handler.end(); ++cell_ttt)
	    if (cell_ttt->material_id() == 1)
	      cell_ttt->set_user_flag();
	    else
	      cell_ttt->clear_user_flag();
       
	  local_dof_handler.distribute_dofs (fe_collection);


    std::vector<unsigned int> block_component_patch (dim+1, 0);
    block_component_patch[dim]=1;
    DoFRenumbering::component_wise (local_dof_handler, block_component_patch);

    std::vector<types::global_dof_index> dofs_per_block_patch (2);
        DoFTools::count_dofs_per_block (local_dof_handler, dofs_per_block_patch, block_component_patch);
    // resize the vector temp to the correct size
    BlockVector<double> temp (dofs_per_block_patch);
    solution_transfer.refine_interpolate(local_solu , temp);
    local_solu = temp;	
	  
 }

// Use normal comments instead of these lines. They are unnecessary bold and
 // prevents me from reading through your function in one piece.
// setup_h_patch_system and patch_rhs

 // Why are there spaces between local_dof_handler and n_dofs?
 unsigned int local_system_size = local_dof_handler. n_dofs();
 h_workload_num = local_dof_handler. n_dofs();

 {
   constraints_patch.clear ();
   FEValuesExtractors::Vector velocities(0);
   DoFTools::make_hanging_node_constraints(local_dof_handler, constraints_patch);

   // And please write sentences in the comments, no need to use underscore or abbreviations
   // Zero boundary condition on patch
   {
     // what means patch_cl? is there a way to rename the variable? Otherwise mention it here in a comment
   	typename hp::DoFHandler<dim>::active_cell_iterator patch_cl= local_dof_handler.begin_active(), end_patch_cl = local_dof_handler.end();
   	for (; patch_cl !=end_patch_cl; ++patch_cl)
   	{
   		std::vector<types::global_dof_index> local_face_dof_indices ((patch_cl->get_fe()).dofs_per_face);
   		if (patch_cl->user_flag_set() == true)
   		{
   			for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
   			{
   				bool face_is_on_patch_Bdry = false;
   				if ((patch_cl->face(f)->at_boundary()) ||  ((patch_cl->neighbor(f)->has_children() == false)
   						&&
   						(patch_cl->neighbor(f)->user_flag_set() == false)))
   					face_is_on_patch_Bdry = true;
   				else if ((patch_cl->face(f)->at_boundary()) || (patch_cl->neighbor(f)->has_children() == true))
   				{
   					for (unsigned int sf=0; sf< patch_cl->face(f)->n_children(); ++sf)
   						if (patch_cl->neighbor_child_on_subface (f, sf) -> user_flag_set() == false)
   						{
   							face_is_on_patch_Bdry = true;
   							break;
   						}
   				}

   				if (face_is_on_patch_Bdry)
   				{
   					patch_cl->face(f)->get_dof_indices (local_face_dof_indices, patch_cl->active_fe_index());
   					for (unsigned int i=0; i<local_face_dof_indices.size(); ++i)
   						
   						if ((patch_cl->get_fe()).face_system_to_component_index(i).first < dim)
						  {
   							constraints_patch.add_line (local_face_dof_indices[i]);
						  }
   				}
   			}
   		}

   	}
   }

 }
 constraints_patch.close();

 DoFTools::count_dofs_per_block (local_dof_handler, dofs_per_block_patch, block_component_patch);
 const unsigned int n_u=dofs_per_block_patch[0], n_p=dofs_per_block_patch[1];

 {
   BlockCompressedSetSparsityPattern csp (dofs_per_block_patch, dofs_per_block_patch);
   DoFTools::make_sparsity_pattern (local_dof_handler, csp, constraints_patch, false);
   sparsity_pattern_patch.copy_from(csp);
 }

 BlockSparseMatrix<double> patch_system (sparsity_pattern_patch);
 BlockVector<double> patch_solution (dofs_per_block_patch);
 BlockVector<double> patch_rhs (dofs_per_block_patch);

 // assemble patch_system and patch_rhs

 hp::FEValues<dim> hp_fe_values (fe_collection, quadrature_collection, update_values|update_quadrature_points|update_JxW_values|update_gradients|update_hessians);

 FullMatrix<double> local_matrix_patch;
 Vector<double> local_rhs_patch;
 Vector<double> local_rhs1;
 Vector<double> local_rhs2;
 std::vector<types::global_dof_index> local_dof_indices;

 std::vector<Vector<double> >  rhs_values;
 const RightHandSide<dim> rhs_function;

 const FEValuesExtractors::Vector velocities (0);
 const FEValuesExtractors::Scalar pressure (dim);

 std::vector<Tensor<2,dim> > grad_phi_u;
 std::vector<double> div_phi_u;
 std::vector<Tensor<1,dim> > phi_u;
 std::vector<double> phi_p;

 std::vector<Tensor<1,dim> > gradients_p;
 std::vector<double> divergences;
 std::vector<Tensor<1,dim> > laplacians;

 std::vector<double> values;
 std::vector<Tensor<2,dim> > gradients;


 typename hp::DoFHandler<dim>::active_cell_iterator patch_cll= local_dof_handler.begin_active(), end_patch_cll = local_dof_handler.end();
 for (; patch_cll!=end_patch_cll; ++patch_cll)
{

     const unsigned int   dofs_per_cell = patch_cll->get_fe().dofs_per_cell;
     if (dofs_per_cell!=0) {
         local_matrix_patch.reinit (dofs_per_cell, dofs_per_cell);
         local_rhs_patch.reinit (dofs_per_cell);
         local_rhs1.reinit (dofs_per_cell);
         local_rhs2.reinit (dofs_per_cell);

         hp_fe_values.reinit (patch_cll);
         const FEValues<dim> &fe_values = hp_fe_values.get_present_fe_values ();
         const std::vector<double>& JxW_values = fe_values.get_JxW_values ();
         const unsigned int n_q_points = fe_values.n_quadrature_points;

         rhs_values.resize(n_q_points, Vector<double>(dim+1));
         rhs_function.vector_value_list (fe_values.get_quadrature_points(), rhs_values);

         grad_phi_u.resize(dofs_per_cell);
         div_phi_u.resize(dofs_per_cell);
         phi_u.resize (dofs_per_cell);
         phi_p.resize(dofs_per_cell);

         divergences.resize(n_q_points);
         gradients_p.resize(n_q_points) ;
         laplacians.resize(n_q_points);

         fe_values[pressure].get_function_gradients(local_solu, gradients_p);
         fe_values[velocities].get_function_divergences(local_solu, divergences);
         fe_values[velocities].get_function_laplacians(local_solu, laplacians);


         for (unsigned int q=0; q<n_q_points; ++q)
           {
             for (unsigned int k=0; k<dofs_per_cell; ++k)
               {
                 grad_phi_u[k] = fe_values[velocities].gradient (k, q);
                 phi_u[k] = fe_values[velocities].value (k, q);
                 phi_p[k] = fe_values[pressure].value (k, q);
               }


             for (unsigned int i=0; i<dofs_per_cell; ++i)
             {
           	  for (unsigned int j=0; j<dofs_per_cell; ++j)
           	  {

           	  local_matrix_patch(i,j) += (double_contract (grad_phi_u[i], grad_phi_u[j])  + (phi_p[i] * phi_p[j]))* JxW_values[q];

           	  } 


                 local_rhs1(i)+= ((rhs_values[q](0)+ laplacians[q][0]-gradients_p[q][0])*(phi_u[i][0])+
					  (rhs_values[q](1)+ laplacians[q][1]-gradients_p[q][1])*(phi_u[i][1]))* JxW_values[q];
                 local_rhs2(i)+= (phi_p[i]*divergences[q])*JxW_values[q];
                 local_rhs_patch(i)= local_rhs1(i)+local_rhs2(i);
               }
           }

        

         local_dof_indices.resize (dofs_per_cell);
         patch_cll->get_dof_indices (local_dof_indices);
         constraints_patch.distribute_local_to_global (local_matrix_patch, local_rhs_patch, local_dof_indices, patch_system, patch_rhs);
     }
 }


 SolverControl           solver_control_stiffness (patch_rhs.block(0).size(),1e-8*patch_rhs.block(0).l2_norm());
 SolverCG<>              cg_stiff (solver_control_stiffness);

 PreconditionSSOR<> preconditioner_stiffness;
 preconditioner_stiffness.initialize(patch_system.block(0,0), 1.2);

 cg_stiff.solve (patch_system.block(0,0), patch_solution.block(0), patch_rhs.block(0),
		  preconditioner_stiffness);



 SolverControl           solver_control_mass (patch_rhs.block(1).size(),1e-8*patch_rhs.block(1).l2_norm());
 SolverCG<>              cg_mass (solver_control_mass);

 PreconditionSSOR<> preconditioner_mass;
 preconditioner_mass.initialize(patch_system.block(1,1), 1.2);

 cg_mass.solve (patch_system.block(1,1), patch_solution.block(1), patch_rhs.block(1),
		  preconditioner_mass);
 constraints_patch.distribute (patch_solution);

//.......................................  get the L2 norm of the gradient of velocity solution and pressure value  .....................//

 double pressure_val=0;
 double grad_u_val=0;
 typename hp::DoFHandler<dim>::active_cell_iterator patch_cel= local_dof_handler.begin_active(), end_patch_cel = local_dof_handler.end();
 for (; patch_cel!=end_patch_cel; ++patch_cel)
 {
	  const unsigned int   dofs_per_cel = patch_cel->get_fe().dofs_per_cell;
	  if (dofs_per_cel!=0)
	  {
		  hp_fe_values.reinit (patch_cel);
		  const FEValues<dim> &fe_values = hp_fe_values.get_present_fe_values ();
		  const std::vector<double>& JxW_values = fe_values.get_JxW_values ();
		  const unsigned int n_q_points = fe_values.n_quadrature_points;

		  gradients.resize(n_q_points);
		  values.resize(n_q_points);

		  fe_values[velocities].get_function_gradients(patch_solution, gradients);
		  fe_values[pressure].get_function_values(patch_solution, values);

		  for (unsigned int q=0; q<n_q_points; ++q)
		  {
			  pressure_val +=values[q]*values[q]* JxW_values[q];

			  for (unsigned int i=0; i<dim; ++i)

				  grad_u_val += contract(gradients[q][i],gradients[q][i])* JxW_values[q];
			
		  } 
		  h_solu_norm_per_patch +=pressure_val + grad_u_val;
	  }

 }
 h_convergence_est_per_cell =sqrt(h_solu_norm_per_patch);

		}  

/*......................................................................................................................................*/

//Compute p_convergence_estimator   &   p_workload_number  for each patch around cell

template <int dim>
void StokesProblem <dim>::p_patch_conv_load_no ( const unsigned int cycle , double &p_convergence_est_per_cell,
		unsigned int &p_workload_num, const typename hp::DoFHandler<dim>::active_cell_iterator &cell,
		unsigned int & patch_number)
		{

	Triangulation<dim> local_triangulation;
	unsigned int level_h_refine;
	unsigned int level_p_refine;
	std::map<typename Triangulation<dim>::active_cell_iterator, typename hp::DoFHandler<dim>::active_cell_iterator> patch_to_global_tria_map;

	std::vector<typename hp::DoFHandler<dim>::active_cell_iterator> patch = get_patch_around_cell(cell);

	build_triangulation_from_patch (patch, local_triangulation, level_h_refine, level_p_refine, patch_to_global_tria_map);
	hp::DoFHandler<dim> local_dof_handler(local_triangulation);

	set_active_fe_indices (local_dof_handler,patch_to_global_tria_map);
	local_dof_handler.distribute_dofs (fe_collection);

	DoFRenumbering::Cuthill_McKee (local_dof_handler);
	
	p_convergence_est_per_cell=0.;
	double p_solu_norm_per_patch=0.;

	ConstraintMatrix constraints_patch;
	BlockSparsityPattern sparsity_pattern_patch;

	std::vector<unsigned int> block_component_patch (dim+1, 0);
	block_component_patch[dim]=1;
	DoFRenumbering::component_wise (local_dof_handler, block_component_patch);

	std::vector<types::global_dof_index> dofs_per_block_patch (2);
	DoFTools::count_dofs_per_block (local_dof_handler, dofs_per_block_patch, block_component_patch);

	BlockVector<double> local_solu (dofs_per_block_patch);

	// Here we are trying to project the values of the global vector "solution" into vector "local_solu" which is
	//solution over patch cells corresponding to cell "K".

	typename hp::DoFHandler<dim>::active_cell_iterator patch_cl= local_dof_handler.begin_active(), end_patch_cl = local_dof_handler.end();
	for (; patch_cl !=end_patch_cl; ++patch_cl)
	{

		const unsigned int   dofs_per_cl = patch_cl->get_fe().dofs_per_cell;

		if (dofs_per_cl!=0)
		{

			Vector<double> local_solution_values(dofs_per_cl);
			typename hp::DoFHandler<dim>::active_cell_iterator global_cl = patch_to_global_tria_map[patch_cl];

			global_cl-> get_dof_values (solution,local_solution_values);
			patch_cl->set_dof_values (local_solution_values, local_solu);
		}
	}

	//......................  Solution Transfer .......................  p_refinement of patch cells ................................. //
	bool need_to_p_refine = false;
	typename hp::DoFHandler<dim>::active_cell_iterator patch_cc= local_dof_handler.begin_active(), end_patch_cc = local_dof_handler.end();
	for (; patch_cc!=end_patch_cc; ++patch_cc)
	{

	typename hp::DoFHandler<dim>::active_cell_iterator global_cell = patch_to_global_tria_map[patch_cc];


		if ( ( global_cell->active_fe_index() +1) < (fe_collection.size()-1) && global_cell->active_fe_index() < (level_p_refine+1))
			need_to_p_refine = true;

	}
	if (need_to_p_refine == true)
	{
		
		SolutionTransfer<dim,BlockVector<double>, hp::DoFHandler<dim>> solution_transfer(local_dof_handler);
		solution_transfer.prepare_for_pure_refinement();
		

		typename hp::DoFHandler<dim>::active_cell_iterator patch_cc= local_dof_handler.begin_active(), end_patch_cc = local_dof_handler.end();
		for (; patch_cc!=end_patch_cc; ++patch_cc)
		{

		       typename hp::DoFHandler<dim>::active_cell_iterator global_cc = patch_to_global_tria_map[patch_cc];

                       patch_cc->set_active_fe_index (global_cc->active_fe_index() +1 );

		}

    local_dof_handler.distribute_dofs (fe_collection);

    std::vector<unsigned int> block_component_patch (dim+1, 0);
    block_component_patch[dim]=1;
    DoFRenumbering::component_wise (local_dof_handler, block_component_patch);

    std::vector<types::global_dof_index> dofs_per_block_patch (2);
    DoFTools::count_dofs_per_block (local_dof_handler, dofs_per_block_patch, block_component_patch);
    BlockVector<double> temp (dofs_per_block_patch);
    solution_transfer.refine_interpolate(local_solu , temp);
    local_solu = temp;
	}

	//......................  setup_p_patch_system and  patch_rhs .............. //

	unsigned int local_system_size = local_dof_handler. n_dofs();
	p_workload_num = local_dof_handler. n_dofs();

	{
		constraints_patch.clear ();

		FEValuesExtractors::Vector velocities(0);
		DoFTools::make_hanging_node_constraints(local_dof_handler, constraints_patch);

		//......................... Zero_Bdry_Condition_on_Patch .......................................//
		{
			typename hp::DoFHandler<dim>::active_cell_iterator patch_cl= local_dof_handler.begin_active(), end_patch_cl = local_dof_handler.end();
			for (; patch_cl !=end_patch_cl; ++patch_cl)
			{
				std::vector<types::global_dof_index> local_face_dof_indices ((patch_cl->get_fe()).dofs_per_face);
				if (patch_cl->user_flag_set() == true)
				{
					for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
					{
						bool face_is_on_patch_Bdry = false;
						if ((patch_cl->face(f)->at_boundary()) ||  ((patch_cl->neighbor(f)->has_children() == false)
								&&
								(patch_cl->neighbor(f)->user_flag_set() == false)))
							face_is_on_patch_Bdry = true;
						else if ((patch_cl->face(f)->at_boundary()) || (patch_cl->neighbor(f)->has_children() == true))
						{
							for (unsigned int sf=0; sf< patch_cl->face(f)->n_children(); ++sf)
								if (patch_cl->neighbor_child_on_subface (f, sf) -> user_flag_set() == false)
								{
									face_is_on_patch_Bdry = true;
									break;
								}
						}

						if (face_is_on_patch_Bdry)
						{
							patch_cl->face(f)->get_dof_indices (local_face_dof_indices, patch_cl->active_fe_index());
							for (unsigned int i=0; i<local_face_dof_indices.size(); ++i)
								// system_to_component_index: "Compute vector component and index of this shape function within the shape functions
								// corresponding to this component from the index of a shape function within this finite element"
								if ((patch_cl->get_fe()).face_system_to_component_index(i).first < dim)
								{
									constraints_patch.add_line (local_face_dof_indices[i]);
								}
						}
					}
				}

			}
		}


	}
	constraints_patch.close();
	DoFTools::count_dofs_per_block (local_dof_handler, dofs_per_block_patch, block_component_patch);
	const unsigned int n_u=dofs_per_block_patch[0], n_p=dofs_per_block_patch[1];

	{
		BlockCompressedSetSparsityPattern csp (dofs_per_block_patch, dofs_per_block_patch);
		DoFTools::make_sparsity_pattern (local_dof_handler, csp, constraints_patch, false);
		sparsity_pattern_patch.copy_from(csp);
	}

	BlockSparseMatrix<double> patch_system (sparsity_pattern_patch);
	BlockVector<double> patch_solution (dofs_per_block_patch);
	BlockVector<double> patch_rhs (dofs_per_block_patch);

	// .........................................  assemble  patch_system  and patch_rhs .............................. //

	hp::FEValues<dim> hp_fe_values (fe_collection, quadrature_collection, update_values|update_quadrature_points|update_JxW_values|update_gradients|update_hessians);

	FullMatrix<double> local_matrix_patch;
	Vector<double> local_rhs_patch;
	Vector<double> local_rhs1;
	Vector<double> local_rhs2;
	std::vector<types::global_dof_index> local_dof_indices;

	std::vector<Vector<double> >  rhs_values;
	const RightHandSide<dim> rhs_function;

	const FEValuesExtractors::Vector velocities (0);
	const FEValuesExtractors::Scalar pressure (dim);

	std::vector<Tensor<2,dim> > grad_phi_u;
	std::vector<double> div_phi_u;
	std::vector<Tensor<1,dim> > phi_u;
	std::vector<double> phi_p;

	std::vector<Tensor<1,dim> > gradients_p;
	std::vector<double> divergences;
	std::vector<Tensor<1,dim> > laplacians;

	std::vector<double> values;
	std::vector<Tensor<2,dim> > gradients;


	typename hp::DoFHandler<dim>::active_cell_iterator patch_cll= local_dof_handler.begin_active(), end_patch_cll = local_dof_handler.end();
	for (; patch_cll!=end_patch_cll; ++patch_cll){

		const unsigned int   dofs_per_cell = patch_cll->get_fe().dofs_per_cell;
		if (dofs_per_cell!=0) {
			local_matrix_patch.reinit (dofs_per_cell, dofs_per_cell);
			local_rhs_patch.reinit (dofs_per_cell);
			local_rhs1.reinit (dofs_per_cell);
			local_rhs2.reinit (dofs_per_cell);

			hp_fe_values.reinit (patch_cll);
			const FEValues<dim> &fe_values = hp_fe_values.get_present_fe_values ();
			const std::vector<double>& JxW_values = fe_values.get_JxW_values ();
			const unsigned int n_q_points = fe_values.n_quadrature_points;

			rhs_values.resize(n_q_points, Vector<double>(dim+1));
			rhs_function.vector_value_list (fe_values.get_quadrature_points(), rhs_values);

			grad_phi_u.resize(dofs_per_cell);
			div_phi_u.resize(dofs_per_cell);
			phi_u.resize (dofs_per_cell);
			phi_p.resize(dofs_per_cell);

			divergences.resize(n_q_points);
			gradients_p.resize(n_q_points) ;
			laplacians.resize(n_q_points);

			fe_values[pressure].get_function_gradients(local_solu, gradients_p);
			fe_values[velocities].get_function_divergences(local_solu, divergences);
			fe_values[velocities].get_function_laplacians(local_solu, laplacians);


			for (unsigned int q=0; q<n_q_points; ++q)
			{
				for (unsigned int k=0; k<dofs_per_cell; ++k)
				{
					grad_phi_u[k] = fe_values[velocities].gradient (k, q);
					phi_u[k] = fe_values[velocities].value (k, q);
					phi_p[k] = fe_values[pressure].value (k, q);
				}


				for (unsigned int i=0; i<dofs_per_cell; ++i)
				{
					for (unsigned int j=0; j<dofs_per_cell; ++j)
					{
						local_matrix_patch(i,j) += (double_contract (grad_phi_u[i], grad_phi_u[j])  + (phi_p[i] * phi_p[j]))* JxW_values[q];

					} 
					local_rhs1(i)+= ((rhs_values[q](0)+laplacians[q][0]-gradients_p[q][0])*(phi_u[i][0])+
							(rhs_values[q](1)+laplacians[q][1]-gradients_p[q][1])*(phi_u[i][1]))* JxW_values[q];
					local_rhs2(i)+= (phi_p[i]*divergences[q])*JxW_values[q];
					local_rhs_patch(i)= local_rhs1(i)+local_rhs2(i);
				}
			}

			local_dof_indices.resize (dofs_per_cell);
			patch_cll->get_dof_indices (local_dof_indices);
			constraints_patch.distribute_local_to_global (local_matrix_patch, local_rhs_patch, local_dof_indices, patch_system, patch_rhs);
		}
	}

	

	/*
	// direct solver
	 SparseDirectUMFPACK A_inverse_stiffness;
	 A_inverse_stiffness.initialize (patch_system.block(0,0),
			 SparseDirectUMFPACK::AdditionalData());
	 A_inverse_stiffness.vmult (patch_solution.block(0), patch_rhs.block(0));
	 //constraints_patch.distribute (patch_solution.block(0));

	 SparseDirectUMFPACK A_inverse_mass;
	 A_inverse_mass.initialize (patch_system.block(1,1),
			 SparseDirectUMFPACK::AdditionalData());
	 A_inverse_mass.vmult (patch_solution.block(1), patch_rhs.block(1));
	 //constraints_patch.distribute (patch_solution.block(1));

	 constraints_patch.distribute (patch_solution);
        */

	// iterative solver
	SolverControl           solver_control_stiffness (patch_rhs.block(0).size(),1e-8*patch_rhs.block(0).l2_norm());
	SolverCG<>              cg_stiff (solver_control_stiffness);

	PreconditionSSOR<> preconditioner_stiffness;
	preconditioner_stiffness.initialize(patch_system.block(0,0), 1.2);
	cg_stiff.solve (patch_system.block(0,0), patch_solution.block(0), patch_rhs.block(0),
			preconditioner_stiffness);

	SolverControl           solver_control_mass (patch_rhs.block(1).size(),1e-8*patch_rhs.block(1).l2_norm());
	SolverCG<>              cg_mass (solver_control_mass);

	PreconditionSSOR<> preconditioner_mass;
	preconditioner_mass.initialize(patch_system.block(1,1), 1.2);

	cg_mass.solve (patch_system.block(1,1), patch_solution.block(1), patch_rhs.block(1),
			preconditioner_mass);

	constraints_patch.distribute (patch_solution);

	//.......................................  get the L2 norm of the gradient of velocity solution and pressure value  .....................//

	double pressure_val=0;
	double grad_u_val=0;
	typename hp::DoFHandler<dim>::active_cell_iterator patch_cel= local_dof_handler.begin_active(), end_patch_cel = local_dof_handler.end();
	for (; patch_cel!=end_patch_cel; ++patch_cel)
	{
		const unsigned int   dofs_per_cel = patch_cel->get_fe().dofs_per_cell;
		if (dofs_per_cel!=0)
		{
			hp_fe_values.reinit (patch_cel);
			const FEValues<dim> &fe_values = hp_fe_values.get_present_fe_values ();
			const std::vector<double>& JxW_values = fe_values.get_JxW_values ();
			const unsigned int n_q_points = fe_values.n_quadrature_points;

			gradients.resize(n_q_points);
			values.resize(n_q_points);

			fe_values[velocities].get_function_gradients(patch_solution, gradients);
			fe_values[pressure].get_function_values(patch_solution, values);


			for (unsigned int q=0; q<n_q_points; ++q)
			{
				pressure_val +=values[q]*values[q]* JxW_values[q];

				for (unsigned int i=0; i<dim; ++i)

					grad_u_val += contract(gradients[q][i],gradients[q][i])* JxW_values[q];
				
			} 
			p_solu_norm_per_patch +=pressure_val + grad_u_val;
		}

	}
	p_convergence_est_per_cell =sqrt(p_solu_norm_per_patch);
	
		}

/*..............................................   marking_cells   .....................................................*/
template <int dim>
void StokesProblem <dim>:: marking_cells (const unsigned int cycle, Vector<float> & marked_cells, std::vector<typename hp::DoFHandler<dim>::active_cell_iterator> &candidate_cell_set,
		std::map<typename hp::DoFHandler<dim>::active_cell_iterator, bool > &p_ref_map, Vector<double> &hp_Conv_Est)
		{
	
	std::vector<std::pair<double, typename hp::DoFHandler<dim>::active_cell_iterator> > to_be_sorted;
	std::vector<std::pair<typename hp::DoFHandler<dim>::active_cell_iterator,bool > > to_be_sorted_with_refine_info;

	Vector<double> est_per_cell (triangulation.n_active_cells());
	estimate(est_per_cell);

	Vector<double> convergence_est_per_cell (triangulation.n_active_cells());
	

	hp_Conv_Est.reinit(triangulation.n_active_cells());

	unsigned int cell_index=0;
	unsigned int patch_number=0;

	typename hp::DoFHandler<dim>::active_cell_iterator
	cell = dof_handler.begin_active(),
	endc = dof_handler.end();
	for (; cell!=endc; ++cell , ++cell_index, ++patch_number)
	{
		
		double indicator_per_cell =0.0;

		double h_convergence_est_per_cell;
		unsigned int h_workload_num;
		h_patch_conv_load_no (cycle ,h_convergence_est_per_cell,h_workload_num, cell, patch_number);

		h_convergence_est_per_cell = h_convergence_est_per_cell /est_per_cell(cell_index);

		double p_convergence_est_per_cell;
		unsigned int p_workload_num;
		p_patch_conv_load_no (cycle ,p_convergence_est_per_cell,p_workload_num, cell, patch_number);

		p_convergence_est_per_cell = p_convergence_est_per_cell /est_per_cell(cell_index);

		double h_ratio= h_convergence_est_per_cell /  h_workload_num ;

		double p_ratio= p_convergence_est_per_cell /  p_workload_num ;

		if (h_ratio > p_ratio)
		{
			convergence_est_per_cell(cell_index)=h_convergence_est_per_cell;
			indicator_per_cell= convergence_est_per_cell(cell_index)*est_per_cell(cell_index);
			hp_Conv_Est(cell_index)=indicator_per_cell;
			p_ref_map[cell] = false;

			std::cout<< "H-refinement_marking ...  =>  p_ref_map[cell] = " << p_ref_map[cell] << std::endl;
			p_ref_map.insert (std::make_pair(cell, false));

		}

		else
		{

			convergence_est_per_cell(cell_index)=p_convergence_est_per_cell;
			indicator_per_cell=convergence_est_per_cell(cell_index)*est_per_cell(cell_index);
			hp_Conv_Est(cell_index)=indicator_per_cell;
			p_ref_map[cell] = true;

			std::cout<< "P-refinement_marking ...  =>  p_ref_map[cell] = "  << p_ref_map[cell] << std::endl;
			p_ref_map.insert (std::make_pair(cell, true));

		}  

		to_be_sorted.push_back(std::make_pair(indicator_per_cell,cell));
	}
	

	unsigned int index=0;

	typename hp::DoFHandler<dim>::active_cell_iterator
	celll = dof_handler.begin_active(),
	endcl = dof_handler.end();
	for (; celll!=endcl; ++celll , ++index)
	{
		

		std::cout<<std::endl;
		std::cout<< "convergence_est_per_cell [" << index << " ] ="<<  convergence_est_per_cell(index) << std::endl;
		std::cout<< "est_per_cell [" << index << " ] ="<< est_per_cell(index) << std::endl;
		std::cout<< "  hp_Conv_Est(index) [" << index << " ] ="<<   hp_Conv_Est(index) << " ?= " <<  " hp 2 _Conv_Est(index)[" << index << " ] ="<<   hp_Conv_Est2(index) << std::endl;
		std::cout<< "refinement_marking ...  =>  p-ref==1,  h-ref==0 : "  << p_ref_map[celll] << std::endl;
		std::cout<<std::endl;
	}
	std::cout<<std::endl;
   
	double theta= 0.5;

	std::sort (to_be_sorted.begin(), to_be_sorted.end(), std_cxx1x::bind(&StokesProblem<dim>::decreasing,this,std_cxx1x::_1,std_cxx1x::_2));

	double L2_norm=est_per_cell.l2_norm();
	double sum=0;
	for (unsigned int i=0; i< to_be_sorted. size(); ++i)
	{
		to_be_sorted_with_refine_info.push_back(std::make_pair(to_be_sorted[i].second, p_ref_map[to_be_sorted[i].second]));
		typename hp::DoFHandler<dim>::active_cell_iterator  cell_sort=to_be_sorted[i].second;
		sum+= (to_be_sorted[i].first)*(to_be_sorted[i].first);
		
		candidate_cell_set.push_back (cell_sort);
		
		if (sum >= (theta*(L2_norm))*(theta*(L2_norm)))

			break;
	}
	unsigned int n= candidate_cell_set.size();
	std::cout<< std::endl;
	std::cout<< "number of candidate_cell_set "<< n << std::endl;
	std::cout<< std::endl;
	//..............................................................................................................................

	marked_cells =0.;
	unsigned int i=0;
	typename hp::DoFHandler<dim>::active_cell_iterator  cel = dof_handler.begin_active(),
			endcel = dof_handler.end();
	for (; cel!=endcel; ++cel,++i)
	{
		typename std::vector<typename hp::DoFHandler<dim>::active_cell_iterator>::iterator  mark_candidate;
		for (mark_candidate=candidate_cell_set.begin(); mark_candidate!=candidate_cell_set.end(); ++ mark_candidate)
		{
			if (cel == (*mark_candidate))
				marked_cells(i)=1;

		}
	}

		}

//.................................................................................................................................
//Output_result
template <int dim>

void StokesProblem <dim>::output_results (const unsigned int cycle , Vector<float> & marked_cells , Vector<double> &est_per_cell , Vector<double> &error_per_cell, Vector<double> &Vect_Pressure_Err, Vector<double> &Vect_grad_Velocity_Err , Vector<double> &hp_Conv_Est )

{

 Vector<float> fe_degrees (triangulation.n_active_cells());
 {
   typename hp::DoFHandler<dim>::active_cell_iterator
   cell = dof_handler.begin_active(),
   endc = dof_handler.end();
   for (unsigned int index=0; cell!=endc; ++cell, ++index)
     fe_degrees(index)= fe_collection[cell->active_fe_index()].degree;
 }

 std::vector<std::string> solution_names (dim, "velocity");

 //std::vector<std::string> solution_names;
 solution_names.push_back ("x_velocity");
 solution_names.push_back ("y_velocity");
 solution_names.push_back ("pressure");

 //std::vector<DataComponentInterpretation::DataComponentInterpretation> data_component_interpretation
 //(dim, DataComponentInterpretation::component_is_part_of_vector);

 std::vector<DataComponentInterpretation::DataComponentInterpretation> data_component_interpretation
 (dim, DataComponentInterpretation::component_is_part_of_vector);

 data_component_interpretation.push_back (DataComponentInterpretation::component_is_scalar);
 data_component_interpretation.push_back (DataComponentInterpretation::component_is_scalar);
 data_component_interpretation.push_back (DataComponentInterpretation::component_is_scalar);


 DataOut<dim,hp::DoFHandler<dim> > data_out;

 data_out.attach_dof_handler (dof_handler);
 

 data_out.add_data_vector (solution, solution_names,DataOut<dim,hp::DoFHandler<dim> >::type_dof_data,data_component_interpretation);


 //data_out.add_data_vector (marked_cells, "marked_cells");
 data_out.add_data_vector (fe_degrees, "fe_degree");
 data_out.add_data_vector (est_per_cell, "Error_Estimator");
 data_out.add_data_vector (error_per_cell, "Error");
 data_out.add_data_vector (Vect_Pressure_Err, "Pressure_Error");
 data_out.add_data_vector (Vect_grad_Velocity_Err, "Grad_Velocity_Error");
 data_out.add_data_vector (Vec_Velocity_Err, "Vec_Velocity_Err");
 
/*
 data_out.add_data_vector (h_Conv_Est, "h_refine_Conv_Est");
 data_out.add_data_vector (p_Conv_Est, "p_refine_Conv_Est");
 data_out.add_data_vector (hp_Conv_Est, "hp_refine_Conv_Est");
*/

 data_out.build_patches ();
 std::string filename = "solution-" +
		  Utilities::int_to_string (cycle, 2) +".vtu";
 std::ofstream output (filename.c_str());
 data_out.write_vtu (output);
}

/*..............................................................................................................*/
//refine_in_h_p()

template <int dim>
void StokesProblem <dim>:: refine_in_h_p (const unsigned int cycle, std::vector<typename hp::DoFHandler<dim>::active_cell_iterator> &candidate_cell_set,
		std::map<typename hp::DoFHandler<dim>::active_cell_iterator, bool > &p_ref_map )
		{
	bool need_to_h_refine=false;
	unsigned int candidate_INDEX=0;

	typename std::vector<typename hp::DoFHandler<dim>::active_cell_iterator>::iterator  cell_candidate;
	for (cell_candidate=candidate_cell_set.begin(); cell_candidate!=candidate_cell_set.end(); ++ cell_candidate, ++candidate_INDEX)
	{
		std::vector<typename hp::DoFHandler<dim>::active_cell_iterator> patch_cells = get_patch_around_cell (* cell_candidate);
		unsigned int level_p_refine = static_cast<unsigned int>((*cell_candidate)->active_fe_index());
		unsigned int level_h_refine = static_cast<unsigned int>((*cell_candidate)->level());

		std::cout<< "  p_ref_map[*cell_candidate] :  " << p_ref_map[*cell_candidate] << std::endl;

		if (p_ref_map[*cell_candidate]==false)
		{
			need_to_h_refine = true;
			(*cell_candidate)->set_refine_flag();
			std::cout<< "h-refinement flags" <<std::endl;

		}//if

		else if(p_ref_map[*cell_candidate]==true)
		{

			if ( ((*cell_candidate)-> active_fe_index()+ 1) <  (fe_collection.size()-1)  )
			{
				(*cell_candidate)->set_active_fe_index ((*cell_candidate)->active_fe_index() + 1);

				std::cout<< "p-refinement flags" <<std::endl;
			}

		}

	}
	if ( need_to_h_refine==true)
		triangulation.execute_coarsening_and_refinement();

bool cell_changed=false;
do
{
cell_changed=false;
unsigned int count_cell=0;
typename hp::DoFHandler<dim>::active_cell_iterator
 cell = dof_handler.begin_active(),
 endc = dof_handler.end();
 for (; cell!=endc; ++cell,++count_cell)
   {
	
        std::vector<typename hp::DoFHandler<dim>::active_cell_iterator> patch_cells = get_patch_around_cell (cell);
      
     for (unsigned int i=1; i<patch_cells.size(); ++i)
       {

         if (patch_cells[i]->active_fe_index()+1 < (patch_cells[0]->active_fe_index()))
          {
             patch_cells[i]->set_active_fe_index(patch_cells[0]->active_fe_index()-1);
             cell_changed=true;
          }
          else if (patch_cells[i]->active_fe_index() > (patch_cells[0]->active_fe_index()+1))
                {
                  patch_cells[0]-> set_active_fe_index (patch_cells[i]->active_fe_index()-1);
                  cell_changed=true;

                }
          }
}

}
while (cell_changed==true);
		}
/*......................................................................................................................................................*/
template <int dim>
void StokesProblem <dim>::run()
{
	for (unsigned int cycle=0; cycle<150; ++cycle)
	{

		std::cout<< std::endl;

		std::cout<< "-----------------------------------------------------------" << std::endl;
		std::cout<< std::endl;
		std::cout << "Cycle " << cycle << ':' << std::endl;
		if (cycle == 0)
		{
			generate_mesh();
			set_global_active_fe_indices(dof_handler);

		}
		std::cout<<"Number of active cells: "<< triangulation.n_active_cells() << std::endl;

		std::cout<<"Total number of cells: " << triangulation.n_cells() << std::endl ;
		setup_system ();
		assemble_system();
		solve ();

		Vector<double> error_per_cell (triangulation.n_active_cells());
		Vector<double> Vect_Pressure_Err(triangulation.n_active_cells());
		Vector<double> Vect_grad_Velocity_Err(triangulation.n_active_cells());
		Vector<double> Vect_Velocity_Err(triangulation.n_active_cells());
		
		compute_error  (error_per_cell, Vect_Pressure_Err, Vect_grad_Velocity_Err, Vect_Velocity_Err);
		Vector<double> est_per_cell (triangulation.n_active_cells());
		estimate(est_per_cell);

		double L2_norm_est= est_per_cell.l2_norm();
		std::cout<< "L2_norm of ERROR Estimate is: "<< L2_norm_est << std::endl;

		//triangulation.refine_global (1);
		
		double L2_norm_est= est_per_cell.l2_norm();
		std::cout<< "L2_norm of ERROR Estimate is: "<< L2_norm_est << std::endl;
		std::cout<<std::endl;
		Vector<float> marked_cells(triangulation.n_active_cells());
		std::vector<typename hp::DoFHandler<dim>::active_cell_iterator> candidate_cell_set;
		std::map<typename hp::DoFHandler<dim>::active_cell_iterator, bool > p_ref_map;

		Vector<double> hp_Conv_Est;

		marking_cells(cycle,  marked_cells, candidate_cell_set, p_ref_map, hp_Conv_Est);
		output_results(cycle, marked_cells, est_per_cell, error_per_cell, Vect_Pressure_Err, Vect_grad_Velocity_Err, hp_Conv_Est);

		refine_in_h_p(cycle,  candidate_cell_set, p_ref_map);

		if (L2_norm_est < Tolerance)
			break;
		 


	}// cycle

}//run

//Explicit initialization

template class StokesProblem<2>;

