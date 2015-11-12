#ifndef _STOKESPROBLEM_HH_
#define _STOKESPROBLEM_HH_

#include <vector>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/function.h>
#include <deal.II/base/function_lib.h>
#include <deal.II/base/utilities.h>

#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/full_matrix.h>

#include <deal.II/lac/block_sparse_matrix.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/lac/precondition_block.h>

#include <deal.II/lac/sparse_ilu.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/compressed_sparsity_pattern.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/grid_out.h>

#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_nothing.h>

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/solution_transfer.h>

#include <deal.II/lac/compressed_set_sparsity_pattern.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/hp/dof_handler.h>
#include <deal.II/hp/fe_values.h>

#include <deal.II/base/std_cxx1x/bind.h>

#include "ExactSolutionEx1.hh"
#include "ExactSolutionEx2.hh"
#include "ExactSolutionEx3.hh"

#include "RightHandSide.hh"
//#include "RightHandSideLocal.hh"


using namespace dealii;

template <int dim>
class StokesProblem
{
  public:
    StokesProblem (int exemple);
    ~StokesProblem();
    void run();

    std::vector<typename hp::DoFHandler<dim>::active_cell_iterator> get_patch_around_cell(const typename hp::DoFHandler<dim>::active_cell_iterator &cell);

    void build_triangulation_from_patch (const std::vector<typename hp::DoFHandler<dim>::active_cell_iterator>  &patch, Triangulation<dim> &local_triangulation, unsigned int &level_h_refine, unsigned int &level_p_refine, std::map<typename Triangulation<dim>::active_cell_iterator, typename hp::DoFHandler<dim>::active_cell_iterator> & patch_to_global_tria_map);

  private:
    
    const RightHandSide<dim> rhs_function;
    //const ExactSolution<dim> exact_solution;
    Function<dim>* exact_solution;
    //const RightHandSideLocal<dim> rhs_function_local;

    void generate_mesh ();
    void set_global_active_fe_indices (hp::DoFHandler<dim> &dof_handler);
    void setup_system ();
    void assemble_system ();
    void solve ();
    double pressure_mean_value () const;
    double exact_pressure_mean_value () const;
   
    void compute_error (Vector<double> &error_per_cell, Vector<double> &Vect_Pressure_Err, Vector<double> &Vect_grad_Velocity_Err, Vector<double> &Vect_Velocity_Err);
    void estimate (Vector<double> &est_per_cell);
    
    void set_active_fe_indices (hp::DoFHandler<dim> &local_dof_handler, std::map<typename Triangulation<dim>::active_cell_iterator, typename hp::DoFHandler<dim>::active_cell_iterator> & patch_to_global_tria_map);
      
  
    void patch_output (unsigned int patch_number ,const unsigned int cycle, hp::DoFHandler<dim> &local_dof_handler, BlockVector<double> &local_solu);

    void solution_on_patch_system_1 (const unsigned int cycle, double &h_gradient_velocity_solution,
    		const typename hp::DoFHandler<dim>::active_cell_iterator &cell, hp::DoFHandler<dim> &local_dof_handler, BlockVector<double> &local_solu);
    void solution_on_patch_system_2 (const unsigned int cycle, double &h_pressure_solution,
        		const typename hp::DoFHandler<dim>::active_cell_iterator &cell, hp::DoFHandler<dim> &local_dof_handler, BlockVector<double> &local_solu);

    void h_patch_conv_load_no (const unsigned int cycle ,
    		double &h_convergence_est_per_cell, unsigned int &h_workload_num,
    		const typename hp::DoFHandler<dim>::active_cell_iterator &cell, unsigned int & patch_number);
    
    void p_patch_conv_load_no (const unsigned int cycle ,
        		double &p_convergence_est_per_cell, unsigned int &p_workload_num,
        		const typename hp::DoFHandler<dim>::active_cell_iterator &cell, unsigned int & patch_number);


    std::vector<typename hp::DoFHandler<dim>::cell_iterator> get_cells_at_coarsest_common_level ( const std::vector<typename hp::DoFHandler<dim>::active_cell_iterator>  &patch);

    bool decreasing (const std::pair<double,typename hp::DoFHandler<dim>::active_cell_iterator> &i, const std::pair<double,typename hp::DoFHandler<dim>::active_cell_iterator > &j);


    void marking_cells (const unsigned int cycle,  Vector<float> & marked_cells, std::vector<typename hp::DoFHandler<dim>::active_cell_iterator> &candidate_cell_set, 
    		std::map<typename hp::DoFHandler<dim>::active_cell_iterator, bool > &p_ref_map, Vector<double> &hp_Conv_Est);

    void output_results (const unsigned int cycle , Vector<float> & marked_cells , Vector<double> &est_per_cell , Vector<double> &error_per_cell, Vector<double> &Vect_Pressure_Err, Vector<double> &Vect_grad_Velocity_Err , Vector<double> &hp_Conv_Est );
    

    void refine_in_h_p (const unsigned int cycle, std::vector<typename hp::DoFHandler<dim>::active_cell_iterator> &candidate_cell_set, 
    		std::map<typename hp::DoFHandler<dim>::active_cell_iterator, bool > &p_ref_map );

    Triangulation<dim> triangulation;
    hp::DoFHandler<dim> dof_handler;

    hp::FECollection<dim> fe_collection;
    hp::QCollection<dim> quadrature_collection;
    hp::QCollection<dim-1> face_quadrature_collection;

    hp::QCollection<dim> quadrature_collection_Err;
    hp::QCollection<dim-1> face_quadrature_collection_Err;

    ConstraintMatrix constraints;
    BlockSparsityPattern sparsity_pattern;
    BlockSparseMatrix<double> system_matrix;

    BlockVector<double> solution;
    BlockVector<double> system_rhs;

   // std::vector<typename hp::DoFHandler<dim>::active_cell_iterator> global_cells; 

    const unsigned int max_degree;
    const double Tolerance;
};

#endif
