#include "StokesProblem.hh"
#include "Parameters.hh"

#include <fstream>
#include <vector>

using namespace dealii;


void output_convergence(EXAMPLE example, REFINEMENT refinement, double theta,
                        std::vector<unsigned int> const &n_dofs_per_cycle,
                        std::vector<double> const &error_per_cycle)
{
  std::ofstream output_stream;
  output_stream.open("convergence.txt");
  switch (example)
    {
    case (example_1) :
    {
      output_stream<<"example_1 ";

      break;
    }
    case (example_2) :
    {
      output_stream<<"example_2 ";

      break;
    }
    case (example_3) :
    {
      output_stream<<"example_3 ";

      break;
    }
    case (example_4) :
    {
      output_stream<<"example_4 ";

      break;
    }
    default :
    {
      AssertThrow(false, ExcMessage("Unknow Example"));
    }
    }

  switch (refinement)
    {
    case (h_refine) :
    {
      output_stream<<"h_refinement ";

      break;
    }
    case (p_refine) :
    {
      output_stream<<"p_refinement ";

      break;
    }
    default :
    {
      output_stream<<"hp_refinement ";
    }
    }

  output_stream<<theta<<std::endl;

  for (unsigned int i=0; i<n_dofs_per_cycle.size(); ++i)
    output_stream<<n_dofs_per_cycle[i]<<" "<<error_per_cycle[i]<<std::endl;

  output_stream.close();
}


template <int dim>
void run(Parameters const &parameters, StokesProblem<dim> &stokes_problem)
{
  std::vector<unsigned int> n_dofs_per_cycles;
  std::vector<double> error_per_cycles;
  for (unsigned int cycle=0; cycle<parameters.get_max_n_cycles(); ++cycle)
    {
      std::cout<< "-----------------------------------------------------------" << std::endl;
      std::cout << "Cycle " << cycle << ':' << std::endl;
      if (cycle == 0)
        stokes_problem.generate_mesh();

      std::cout<<"Number of active cells: "<< stokes_problem.n_active_cells() << std::endl;

      std::cout<<"Set up"<<std::endl;
      stokes_problem.setup_system();

      std::cout<<"Number of degrees of freedom: "<<stokes_problem.n_dofs()<<std::endl;
      n_dofs_per_cycles.push_back(stokes_problem.n_dofs());

      std::cout<<"Assemble system"<<std::endl;
      stokes_problem.assemble_system();

      std::cout<<"Solve"<<std::endl;
      stokes_problem.solve();

      std::cout<<"Compute the error"<<std::endl;
      stokes_problem.compute_error();
      std::cout<< "L2_norm of ERROR is: "<< stokes_problem.error_l2_norm() << std::endl;
      error_per_cycles.push_back(stokes_problem.error_l2_norm());

      std::cout<<"Compute the error estimate"<<std::endl;
      stokes_problem.estimate_error();
      std::cout<< "L2_norm of ERROR Estimate is: "<< stokes_problem.error_estimate_l2_norm() << std::endl;

      std::cout<<"Marking Cell"<<std::endl;
      stokes_problem.mark_cells(cycle, parameters.get_theta());

      std::cout<<"Output results"<<std::endl;
      stokes_problem.output_results(cycle);

      std::cout<<"Refinement"<<std::endl;
      stokes_problem.refine_in_h_p();

      if (stokes_problem.error_estimate_l2_norm()<parameters.get_tolerance())
        break;
    }

  output_convergence(parameters.get_example(), parameters.get_refinement(),
                     parameters.get_theta(), n_dofs_per_cycles, error_per_cycles);
}

int main(int argc, char *argv[])
{
  // Silence a warning
  (void) argc;
  deallog.depth_console (0);
  Parameters parameters(argv[1]);

  if (parameters.get_dim()==2)
    {
      StokesProblem<2> stokes_problem(parameters.get_verbose(),
                                      parameters.get_example(), parameters.get_quadrature(),
                                      parameters.get_refinement(), parameters.get_max_degree());
      run(parameters, stokes_problem);
    }
  else
    {
      StokesProblem<3> stokes_problem(parameters.get_verbose(),
                                      parameters.get_example(), parameters.get_quadrature(),
                                      parameters.get_refinement(), parameters.get_max_degree());
      run(parameters, stokes_problem);
    }


  return 0;
}
