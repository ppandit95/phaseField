//beta prime precipitate evolution implementation
//general headers
#include "../../include/dealIIheaders.h"

//precipitate problem headers
#include "parameters.h"
#include "../../src/coupled_CH_AC_Mechanics.h"
 
//initial condition functions
//concentration initial conditions
template <int dim>
double InitialConditionC<dim>::value (const Point<dim> &p, const unsigned int /* component */) const
{
  //set result equal to the concentration initial condition 
  //return 0.03 + 1.0e-3*(2*(0.5 - (double)(std::rand() % 100 )/100.0));
  double dx=spanX/std::pow(2.0,refineFactor);
  double r=0.0;
#if problemDIM==1
  r=p[0];
  return 0.005+0.5*(0.125-0.005)*(1-std::tanh((r-spanX/2.0)/(3*dx)));
#elif problemDIM==2
  r=p.distance(Point<dim>(spanX/2.0,spanY/2.0));
  return 0.005+0.5*(0.125-0.005)*(1-std::tanh((r-spanX/8.0)/(3*dx)));
#elif problemDIM==3
  r=p.distance(Point<dim>(spanX/2.0,spanY/2.0,spanZ/2.0));
  return 0.005+0.5*(0.125-0.005)*(1-std::tanh((r-spanX/8.0)/(3*dx)));
#endif
}

//structural order parameter initial conditions
template <int dim>
double InitialConditionN<dim>::value (const Point<dim> &p, const unsigned int /* component */) const
{
  //set result equal to the structural order paramter initial condition
  double dx=spanX/std::pow(2.0,refineFactor);
  double r=0.0;
#if problemDIM==1
  r=p[0];
  return 0.5*(1.0-std::tanh((r-spanX/2.0)/(6.2*dx)));
#elif problemDIM==2
  //  double r1=p.distance(Point<dim>(spanX/4.0,spanY/4.0));
  //  double r2=p.distance(Point<dim>(3*spanX/4.0,spanY/4.0));
  //  double r3=p.distance(Point<dim>(spanX/4.0,3*spanY/4.0));
  //  double r4=p.distance(Point<dim>(3*spanX/4.0,3*spanY/4.0));
  //  r=std::min(std::min(r1,r2),std::min(r3,r4));
  r=p.distance(Point<dim>(spanX/2.0,spanY/2.0));
  return 0.5*(1.0-std::tanh((r-spanX/8.0)/(3*dx)));
#elif problemDIM==3
  r=p.distance(Point<dim>(spanX/2.0,spanY/2.0,spanZ/2.0));
  return 0.5*(1.0-std::tanh((r-spanX/8.0)/(3*dx)));
#endif
}

//main
int main (int argc, char **argv)
{
  Utilities::System::MPI_InitFinalize mpi_initialization(argc, argv,numbers::invalid_unsigned_int);
  try
    {
      deallog.depth_console(0);
      PrecipitateProblem<problemDIM> precipitateProblem;
      precipitateProblem.run ();
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }
  
  return 0;
}