// =================================================================================
// Set the attributes of the primary field variables
// =================================================================================
// This function sets attributes for each variable/equation in the app. The
// attributes are set via standardized function calls. The first parameter for each
// function call is the variable index (starting at zero). The first set of
// variable/equation attributes are the variable name (any string), the variable
// type (SCALAR/VECTOR), and the equation type (EXPLICIT_TIME_DEPENDENT/
// TIME_INDEPENDENT/AUXILIARY). The next set of attributes describe the
// dependencies for the governing equation on the values and derivatives of the
// other variables for the value term and gradient term of the RHS and the LHS.
// The final pair of attributes determine whether a variable represents a field
// that can nucleate and whether the value of the field is needed for nucleation
// rate calculations.

void variableAttributeLoader::loadVariableAttributes(){

	// Variable 0
	set_variable_name				(0,"n1");
	set_variable_type				(0,SCALAR);
	set_variable_equation_type		(0,EXPLICIT_TIME_DEPENDENT);

    set_dependencies_value_term_RHS(0, " n1, grad(u)");
    set_dependencies_gradient_term_RHS(0, "grad(n1)");

	// Variable 1
	set_variable_name				(1,"u");
	set_variable_type				(1,VECTOR);
	set_variable_equation_type		(1,TIME_INDEPENDENT);

    set_dependencies_value_term_RHS(1, "");
    set_dependencies_gradient_term_RHS(1, "n1, grad(u)");
    set_dependencies_value_term_LHS(1, "");
    set_dependencies_gradient_term_LHS(1, "n1, grad(change(u))");

}

// =============================================================================================
// explicitEquationRHS (needed only if one or more equation is explict time dependent)
// =============================================================================================
// This function calculates the right-hand-side of the explicit time-dependent
// equations for each variable. It takes "variable_list" as an input, which is a list
// of the value and derivatives of each of the variables at a specific quadrature
// point. The (x,y,z) location of that quadrature point is given by "q_point_loc".
// The function outputs two terms to variable_list -- one proportional to the test
// function and one proportional to the gradient of the test function. The index for
// each variable in this list corresponds to the index given at the top of this file.

template <int dim, int degree>
void customPDE<dim,degree>::explicitEquationRHS(variableContainer<dim,degree,dealii::VectorizedArray<double> > & variable_list,
				 dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const {

// --- Getting the values and derivatives of the model variables ---

// The first order parameter and its derivatives (names here should match those in the macros above)
scalarvalueType n1 = variable_list.get_scalar_value(0);
scalargradType n1x = variable_list.get_scalar_gradient(0);

// The derivative of the displacement vector (names here should match those in the macros above)
vectorgradType ux = variable_list.get_vector_gradient(1);

// --- Setting the expressions for the terms in the governing equations ---
double B = 3.0*A + 12.0;
double C = 2.0*A + 12.0;
scalarvaluetype fV = 1.0 + constV(A/2)*n1*n1 - constV(B/3)*n1*n1*n1 + constV(C/4)*n1*n1*n1*n1;
scalarvaluetype fn1V = constV(A)*n1 - constV(B)*n1*n1 + constV(C)*n1*n1*n1;

//Calculating stiffness matrix for martensite phase
dealii::Tensor<2,CIJ_tensor_size> CIJ_M = constV(1.1)*CIJ_A;


// Calculate the stress-free transformation strain and its derivatives at the quadrature point
dealii::Tensor<2, dim, dealii::VectorizedArray<double> > sfts,sfts1n;

for (unsigned int i=0; i<dim; i++){
for (unsigned int j=0; j<dim; j++){
	// Polynomial fits for the stress-free transformation strains
	sfts[i][j] = constV(sfts_linear[i][j])*n1;
	sfts1n[i][j] = constV(sfts_linear[i][j]);

}
}

//compute E2=(E-E0)
dealii::VectorizedArray<double> E2[dim][dim], S[dim][dim];

for (unsigned int i=0; i<dim; i++){
for (unsigned int j=0; j<dim; j++){
	  E2[i][j]= constV(0.5)*(ux[i][j]+ux[j][i])- sfts1[i][j];
}
}

//compute stress
//S=C*(E-E0)
// Compute stress tensor (which is equal to the residual, Rux)
dealii::VectorizedArray<double> CIJ_combined[CIJ_tensor_size][CIJ_tensor_size];

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

// Compute one of the stress terms in the order parameter chemical potential, nDependentMisfitACp = -C*(E-E0)*(E0_n)
dealii::VectorizedArray<double> nDependentMisfitAC1=constV(0.0);

for (unsigned int i=0; i<dim; i++){
for (unsigned int j=0; j<dim; j++){
	  nDependentMisfitAC1+=-S[i][j]*sfts1n[i][j];
}
}


// Compute the other stress term in the order parameter chemical potential, heterMechACp = 0.5*Hn*(C_beta-C_alpha)*(E-E0)*(E-E0)
dealii::VectorizedArray<double> heterMechAC1=constV(0.0);
dealii::VectorizedArray<double> S2[dim][dim];

if (n_dependent_stiffness == true){
	computeStress<dim>(CIJ_M-CIJ_A, E2, S2);

	for (unsigned int i=0; i<dim; i++){
		for (unsigned int j=0; j<dim; j++){
			heterMechAC1 += constV(0.5)*S2[i][j]*E2[i][j];
		}
	}
}
// The terms in the governing equations
scalarvalueType eq_n1  = (n1-constV(userInputs.dtValue*Mn1V)*( constV(G*ks/L)*fn1V - nDependentMisfitAC1 + heterMechAC1));
scalargradType eqx_n1 = (constV(-userInputs.dtValue*Mn1V*kg*G*L)*n1x);

// --- Submitting the terms for the governing equations ---

variable_list.set_scalar_value_term_RHS(2,eq_n1);
variable_list.set_scalar_gradient_term_RHS(2,eqx_n1);

}

// =============================================================================================
// nonExplicitEquationRHS (needed only if one or more equation is time independent or auxiliary)
// =============================================================================================
// This function calculates the right-hand-side of all of the equations that are not
// explicit time-dependent equations. It takes "variable_list" as an input, which is
// a list of the value and derivatives of each of the variables at a specific
// quadrature point. The (x,y,z) location of that quadrature point is given by
// "q_point_loc". The function outputs two terms to variable_list -- one proportional
// to the test function and one proportional to the gradient of the test function. The
// index for each variable in this list corresponds to the index given at the top of
// this file.

template <int dim, int degree>
void customPDE<dim,degree>::nonExplicitEquationRHS(variableContainer<dim,degree,dealii::VectorizedArray<double> > & variable_list,
				 dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const {

 // --- Getting the values and derivatives of the model variables ---

 // The first order parameter and its derivatives (names here should match those in the macros above)
 scalarvalueType n1 = variable_list.get_scalar_value(0);

 // The derivative of the displacement vector (names here should match those in the macros above)
 vectorgradType ux = variable_list.get_vector_gradient(1);

 // --- Setting the expressions for the terms in the governing equations ---
 double B = 3.0*A + 12.0;
 double C = 2.0*A + 12.0;
 scalarvaluetype fV = 1.0 + constV(A/2)*n1*n1 - constV(B/3)*n1*n1*n1 + constV(C/4)*n1*n1*n1*n1;
 scalarvaluetype fn1V = constV(A)*n1 - constV(B)*n1*n1 + constV(C)*n1*n1*n1;

 //Calculating stiffness matrix for martensite phase
 dealii::Tensor<2,CIJ_tensor_size> CIJ_M = constV(1.1)*CIJ_A;


 // Calculate the stress-free transformation strain and its derivatives at the quadrature point
 dealii::Tensor<2, dim, dealii::VectorizedArray<double> > sfts,sfts1n;

 for (unsigned int i=0; i<dim; i++){
 for (unsigned int j=0; j<dim; j++){
 	// Polynomial fits for the stress-free transformation strains,
 	sfts[i][j] = constV(sfts_linear[i][j])*n1;
 	sfts1n[i][j] = constV(sfts_linear[i][j]);
 }
 }

 //compute E2=(E-E0)
 dealii::VectorizedArray<double> E2[dim][dim], S[dim][dim];

 for (unsigned int i=0; i<dim; i++){
 for (unsigned int j=0; j<dim; j++){
 	  E2[i][j]= constV(0.5)*(ux[i][j]+ux[j][i])- sfts[i][j];
 }
 }

 //compute stress
 //S=C*(E-E0)
 // Compute stress tensor (which is equal to the residual, Rux)
 dealii::VectorizedArray<double> CIJ_combined[CIJ_tensor_size][CIJ_tensor_size];

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


 vectorgradType eqx_u;

 // Fill residual corresponding to mechanics
 // R=-C*(E-E0)

 for (unsigned int i=0; i<dim; i++){
 for (unsigned int j=0; j<dim; j++){
 	  eqx_u[i][j] = - S[i][j];
 }
 }

 // --- Submitting the terms for the governing equations ---

 variable_list.set_vector_gradient_term_RHS(1,eqx_u);

}

// =============================================================================================
// equationLHS (needed only if at least one equation is time independent)
// =============================================================================================
// This function calculates the left-hand-side of time-independent equations. It
// takes "variable_list" as an input, which is a list of the value and derivatives of
// each of the variables at a specific quadrature point. The (x,y,z) location of that
// quadrature point is given by "q_point_loc". The function outputs two terms to
// variable_list -- one proportional to the test function and one proportional to the
// gradient of the test function -- for the left-hand-side of the equation. The index
// for each variable in this list corresponds to the index given at the top of this
// file. If there are multiple elliptic equations, conditional statements should be
// sed to ensure that the correct residual is being submitted. The index of the field
// being solved can be accessed by "this->currentFieldIndex".

template <int dim, int degree>
void customPDE<dim,degree>::equationLHS(variableContainer<dim,degree,dealii::VectorizedArray<double> > & variable_list,
		dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const {

// --- Getting the values and derivatives of the model variables ---

	//n1
	scalarvalueType n1 = variable_list.get_scalar_value(0);

    // --- Setting the expressions for the terms in the governing equations ---

	vectorgradType eqx_Du;

	// Take advantage of E being simply 0.5*(ux + transpose(ux)) and use the dealii "symmetrize" function
	dealii::Tensor<2, dim, dealii::VectorizedArray<double> > E;
	E = symmetrize(variable_list.get_change_in_vector_gradient(1));

	//Calculating stiffness matrix for martensite phase
	 dealii::Tensor<2,CIJ_tensor_size> CIJ_M = constV(1.1)*CIJ_A;

	// Compute stress tensor (which is equal to the residual, Rux)
	if (n_dependent_stiffness == true){
		dealii::Tensor<2, CIJ_tensor_size, dealii::VectorizedArray<double> > CIJ_combined;
		CIJ_combined = CIJ_A + n1*(CIJ_M - CIJ_A);


		computeStress<dim>(CIJ_combined, E, eqx_Du);
	}
	else{
		computeStress<dim>(CIJ_A, E, eqx_Du);
	}

    // --- Submitting the terms for the governing equations ---

	variable_list.set_vector_gradient_term_LHS(1,eqx_Du);

}
