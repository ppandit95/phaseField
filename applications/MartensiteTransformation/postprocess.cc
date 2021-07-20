// =============================================================================================
// loadPostProcessorVariableAttributes: Set the attributes of the postprocessing variables
// =============================================================================================
// This function is analogous to 'loadVariableAttributes' in 'equations.h', but for
// the postprocessing expressions. It sets the attributes for each postprocessing
// expression, including its name, whether it is a vector or scalar (only scalars are
// supported at present), its dependencies on other variables and their derivatives,
// and whether to calculate an integral of the postprocessed quantity over the entire
// domain. Note: this function is not a member of customPDE.

void variableAttributeLoader::loadPostProcessorVariableAttributes(){
	// Variable 0
	set_variable_name				(0,"psi_tot");
	set_variable_type				(0,SCALAR);

    set_dependencies_value_term_RHS(0, "n1, grad(n1), grad(u)");
    set_dependencies_gradient_term_RHS(0, "");

    set_output_integral         	(0,true);

	// Variable 1
	set_variable_name				(1,"von_mises_stress");
	set_variable_type				(1,SCALAR);

    set_dependencies_value_term_RHS(1, "n1, grad(u)");
    set_dependencies_gradient_term_RHS(1, "");

	set_output_integral         	(1,false);

}

// =============================================================================================
// postProcessedFields: Set the postprocessing expressions
// =============================================================================================
// This function is analogous to 'explicitEquationRHS' and 'nonExplicitEquationRHS' in
// equations.h. It takes in "variable_list" and "q_point_loc" as inputs and outputs two terms in
// the expression for the postprocessing variable -- one proportional to the test
// function and one proportional to the gradient of the test function. The index for
// each variable in this list corresponds to the index given at the top of this file (for
// submitting the terms) and the index in 'equations.h' for assigning the values/derivatives of
// the primary variables.

template <int dim,int degree>
void customPDE<dim,degree>::postProcessedFields(const variableContainer<dim,degree,dealii::VectorizedArray<double> > & variable_list,
				variableContainer<dim,degree,dealii::VectorizedArray<double> > & pp_variable_list,
												const dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const {

    // --- Getting the values and derivatives of the model variables ---

		// The first order parameter and its derivatives (names here should match those in the macros above)
		scalarvalueType n1 = variable_list.get_scalar_value(0);
		scalargradType n1x = variable_list.get_scalar_gradient(0);

		// The derivative of the displacement vector (names here should match those in the macros above)
		vectorgradType ux = variable_list.get_vector_gradient(1);

        // --- Setting the expressions for the terms in the postprocessing expressions ---

		double B = 3.0*A + 12.0;
		double C = 2.0*A + 12.0;
		scalarvaluetype fV = 1.0 + constV(A/2)*n1*n1 - constV(B/3)*n1*n1*n1 + constV(C/4)*n1*n1*n1*n1;

        // Start calculating components of the energy density
        scalarvalueType total_energy_density = constV(0.0);

		scalarvalueType psi_sep = constV(ks*(G/L))*fV;

		scalarvalueType psi_grad = constV(0.0);

		for (int i=0; i<dim; i++){
			for (int j=0; j<dim; j++){
				psi_grad += constV(0.5*kg*G*L)*n1x[i]*n1x[j];
			}
		}
		// Calculate the stress-free transformation strain and its derivatives at the quadrature point
		dealii::Tensor<2, dim, dealii::VectorizedArray<double> > sfts1;
		for (unsigned int i=0; i<dim; i++){
		for (unsigned int j=0; j<dim; j++){
			sfts1[i][j] = constV(sfts_linear1[i][j])*n1;
		}
		}

		//compute E2=(E-E0)
		dealii::VectorizedArray<double> E2[dim][dim], S[dim][dim];

		for (unsigned int i=0; i<dim; i++){
			for (unsigned int j=0; j<dim; j++){
				E2[i][j]= constV(0.5)*(ux[i][j]+ux[j][i]) - sfts1[i][j];

			}
		}

		//compute stress
		//S=C*(E-E0)
		dealii::VectorizedArray<double> CIJ_combined[2*dim-1+dim/3][2*dim-1+dim/3];
		dealii::Tensor<2,CIJ_tensor_size> CIJ_M = constV(1.1)*CIJ_A;
		if (n_dependent_stiffness == true){
			for (unsigned int i=0; i<2*dim-1+dim/3; i++){
				for (unsigned int j=0; j<2*dim-1+dim/3; j++){
					CIJ_combined[i][j] = CIJ_A[i][j] + n1*(CIJ_M[i][j] - CIJ_A[i][j]);
				}
			}
			computeStress<dim>(CIJ_combined, E2, S);
		}
		else{
			computeStress<dim>(CIJ_A, E2, S);
		}

		scalarvalueType psi_el = constV(0.0);

		for (unsigned int i=0; i<dim; i++){
			for (unsigned int j=0; j<dim; j++){
				psi_el += constV(0.5) * S[i][j]*E2[i][j];
			}
		}

		total_energy_density = psi_sep + psi_grad + psi_el;

// The Von Mises Stress
dealii::VectorizedArray<double> vm_stress;
if (dim == 3){
    vm_stress = (S[0][0]-S[1][1])*(S[0][0]-S[1][1]) + (S[1][1]-S[2][2])*(S[1][1]-S[2][2]) + (S[2][2]-S[0][0])*(S[2][2]-S[0][0]);
    vm_stress += constV(6.0)*(S[0][1]*S[0][1] + S[1][2]*S[1][2] + S[2][0]*S[2][0]);
    vm_stress *= constV(0.5);
    vm_stress = std::sqrt(vm_stress);
}
else {
    vm_stress = S[0][0]*S[0][0] - S[0][0]*S[1][1] + S[1][1]*S[1][1] + constV(3.0)*S[0][1]*S[0][1];
    vm_stress = std::sqrt(vm_stress);
}


// --- Submitting the terms for the postprocessing expressions ---

pp_variable_list.set_scalar_value_term_RHS(0, total_energy_density);
pp_variable_list.set_scalar_value_term_RHS(1, vm_stress);

}
