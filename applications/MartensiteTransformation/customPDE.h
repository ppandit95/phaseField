#include "../../include/matrixFreePDE.h"

template <int dim, int degree>
class customPDE: public MatrixFreePDE<dim,degree>
{
public:
	customPDE(userInputParameters<dim> _userInputs): MatrixFreePDE<dim,degree>(_userInputs) , userInputs(_userInputs) {
		c_dependent_misfit = false;
		for (unsigned int i=0; i<dim; i++){
			for (unsigned int j=0; j<dim; j++){
				if (std::abs(sfts_linear1[i][j])>1.0e-12){
					c_dependent_misfit = true;
				}
			}
		}
	};
    // Function to set the initial conditions (in ICs_and_BCs.h)
    void setInitialCondition(const dealii::Point<dim> &p, const unsigned int index, double & scalar_IC, dealii::Vector<double> & vector_IC);

    // Function to set the non-uniform Dirichlet boundary conditions (in ICs_and_BCs.h)
    void setNonUniformDirichletBCs(const dealii::Point<dim> &p, const unsigned int index, const unsigned int direction, const double time, double & scalar_BC, dealii::Vector<double> & vector_BC);

private:
    #include "../../include/typeDefs.h"

    const userInputParameters<dim> userInputs;

    // Function to set the RHS of the governing equations for explicit time dependent equations (in equations.h)
    void explicitEquationRHS(variableContainer<dim,degree,dealii::VectorizedArray<double> > & variable_list,
                     dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const;

    // Function to set the RHS of the governing equations for all other equations (in equations.h)
    void nonExplicitEquationRHS(variableContainer<dim,degree,dealii::VectorizedArray<double> > & variable_list,
                     dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const;

    // Function to set the LHS of the governing equations (in equations.h)
    void equationLHS(variableContainer<dim,degree,dealii::VectorizedArray<double> > & variable_list,
                     dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const;

    // Function to set postprocessing expressions (in postprocess.h)
    #ifdef POSTPROCESS_FILE_EXISTS
    void postProcessedFields(const variableContainer<dim,degree,dealii::VectorizedArray<double> > & variable_list,
                    variableContainer<dim,degree,dealii::VectorizedArray<double> > & pp_variable_list,
                    const dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const;
    #endif

    // Function to set the nucleation probability (in nucleation.h)
    #ifdef NUCLEATION_FILE_EXISTS
    double getNucleationProbability(variableValueContainer variable_value, double dV) const;
    #endif

	// ================================================================
	// Methods specific to this subclass
	// ================================================================


	// ================================================================
	// Model constants specific to this subclass
	// ================================================================
    double A = userInputs.get_model_constant_double("A");
	double Mn1V = userInputs.get_model_constant_double("Mn1V");

	double ks = userInputs.get_model_constant_double("ks");
	double kg = userInputs.get_model_constant_double("kg");

	double G = userInputs.get_model_constant_double("G");
	double L = userInputs.get_model_constant_double("L");


	const static unsigned int CIJ_tensor_size =2*dim-1+dim/3;
	dealii::Tensor<2,CIJ_tensor_size> CIJ_A = userInputs.get_model_constant_elasticity_tensor("CIJ_A");

	dealii::Tensor<2,dim> sfts_linear = userInputs.get_model_constant_rank_2_tensor("sfts_linear");
	bool n_dependent_stiffness = userInputs.get_model_constant_bool("n_dependent_stiffness");





	// ================================================================

};
