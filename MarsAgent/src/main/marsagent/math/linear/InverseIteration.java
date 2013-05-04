/**
 * 
 */
package marsagent.math.linear;

import java.util.ArrayList;

import marsagent.math.exception.DimensionMissmatchException;
import marsagent.math.exception.NonSquareMatrixException;

import org.jblas.DoubleMatrix;
import org.jblas.Eigen;
import org.jblas.Solve;

/**
 * Solve the lowest eigenvalue <math xmlns="http://www.w3.org/1998/Math/MathML" display="inline"><mi>&lambda;<sub><mn>1</mn></sub><sup><mn>*</mn></sup></mi></math>
 * and it's corresponding eigenvector
 * <math xmlns="http://www.w3.org/1998/Math/MathML" display="inline"><mi><b>&varphi;</b><sub><mn>1</mn></sub><sup><mn>*</mn></sup></mi></math>
 * using the inverse iteration method with Rayleigh quotient.
 * 
 * <p>The <code>InverseIteration</code> solves the following generalized eigenvalue problem:</p>
 * <math xmlns="http://www.w3.org/1998/Math/MathML" display="block">
 *   <mo>(</mo>
 *     <mn><b>A</b></mn>
 *     <mo>-</mo>
 *     <mn>&lambda;<sub><mn>1</mn></sub><sup>*</sup></mn><mn><b>B</b></mn>
 *   <mo>)</mo><mi><b>&varphi;</b><sub><mn>1</mn></sub><sup>*</sup></mi>
 *   <mo>=</mo>
 *   <mn><b>0</b></mn>
 * </math>
 * 
 * <p>where it is assumed that the sought eigenvalue <math xmlns="http://www.w3.org/1998/Math/MathML" display="inline">
 * <mi>&lambda;<sub><mn>1</mn></sub></mi></math> is distinct and that the eigenvalues are ordered according to</p>
 * 
 * <math xmlns="http://www.w3.org/1998/Math/MathML" display="block">
 *   <mi>&lambda;<sub><mn>1</mn></sub></mi><mo>&lt;</mo>
 *   <mi>&lambda;<sub><mn>2</mn></sub></mi><mo>&leq;</mo>
 *   <mi>&lambda;<sub><mn>3</mn></sub></mi><mo>&leq;</mo>
 *   <mo>...</mo><mo>&leq;</mo>
 *   <mi>&lambda;<sub><mi>N</mi></sub></mi>
 * </math>
 * 
 * <p>The eigenvalue <math xmlns="http://www.w3.org/1998/Math/MathML" display="inline"><mi>&lambda;<sub><mn>1</mn></sub><sup><mn>*</mn></sup></mi></math>
 * is approximated such that
 * <math xmlns="http://www.w3.org/1998/Math/MathML" display="inline"><mi>&lambda;<sub><mn>1</mn></sub></mn><mo>&leq;</mo><mi>&lambda;<sub><mn>1</mn></sub><sup><mn>*</mn></sup></mi></math>
 * where <math xmlns="http://www.w3.org/1998/Math/MathML" display="inline"><mi>&lambda;<sub><mn>1</mn></sub></mn></math> is the true lowest eigenvalue.</p>
 * 
 * <p>It is possible to solve for other approximations of the eigenvalue other
 * than the lowest eigenvalue by shifting the eigenvalue spectrum using
 * <math xmlns="http://www.w3.org/1998/Math/MathML" display="inline"><mi>&mu;</mi><mo>&gt;</mo><mi>&lambda;<sub><mi>n</mi></sub></mi><mo>&gt;</mo><mi>&lambda;<sub><mn>1</mn></sub></mi></math>
 * where <math xmlns="http://www.w3.org/1998/Math/MathML" display="inline"><mi>&lambda;<sub><mi>n</mi></sub></mn></math> is the sought eigenvalue.
 * The following problem is then solved:</p>
 * 
 * <math xmlns="http://www.w3.org/1998/Math/MathML" display="block">
 *   <mo>[</mo>
 *     <mo>(</mo><mn><b>A</b></mn><mo>-</mo><mn>&mu;</mn><mn><b>B</b></mn><mo>)</mo>
 *     <mo>-</mo>
 *     <mo>(</mo><mn>&lambda;<sub><mi>n</mi></sub><sup>*</sup></mn><mo>-</mo><mn>&mu;</mn><mo>)</mo>
 *     <mn><b>B</b></mn>
 *   <mo>]</mo><mi><b>&varphi;</b><sub><mi>n</mi></sub><sup>*</sup></mi>
 *   <mo>=</mo>
 *   <mn><b>0</b></mn>
 * </math>
 * 
 * @author 	AlexanderSehlstrom
 * @since	1.0
 * @version	1.0.2
 *
 */
public class InverseIteration {

	/**
	 * Structure holding the computed results.
	 * @author 	AlexanderSehlstrom
	 * @since	1.0
	 * @version	1.0
	 *
	 */
	public class ResultStructure {
		/** Eigenvalue */
		public double lambda;

		/** Eigenvector */
		public DoubleMatrix phi;
	}

	/** Default maximum number of iterations */
	private static final int MAX_ITER = 10;

	/** Default change tolerance */
	private static final double TOL = 1e-9;

	public static void main(String[] args) {
		InverseIteration ii = new InverseIteration();

		DoubleMatrix a = new DoubleMatrix(new double[][]{
				{ 2, -1,  0,  0},
				{-1,  2, -1,  0},
				{ 0, -1,  2, -1},
				{ 0,  0, -1,  2}});

		DoubleMatrix b = DoubleMatrix.eye(a.rows);

		System.out.println("GENERALIZED SYSTEM");
		ResultStructure rs1 = ii.compute(a, b);

		System.out.println("Inverse iteration results");
		System.out.println("lambda = " + rs1.lambda);
		System.out.println("phi(1) = " + rs1.phi.get(0));
		System.out.println("phi(2) = " + rs1.phi.get(1));

		DoubleMatrix rs2lambda = Eigen.symmetricGeneralizedEigenvalues(a, b);
		DoubleMatrix[] rs2phi = Eigen.symmetricGeneralizedEigenvectors(a, b);

		System.out.println("\njBlas solver results");
		DoubleMatrix rs2lambdapositive = rs2lambda.ge(0).dup();
		rs2lambda = rs2lambda.get(rs2lambdapositive);
		System.out.println("lambda = " + rs2lambda.min() + " (the " + rs2lambda.argmin() + "th eigenvalue)");
		System.out.println("phi(1) = " + rs2phi[rs2lambda.argmin()].get(0));
		System.out.println("phi(2) = " + rs2phi[rs2lambda.argmin()].get(1));

		System.out.println("\nRatio: " + (rs1.lambda/rs2lambda.min()));
	}

	/** Maximum number of iterations */
	private int maxiter;

	/** Change tolerance */
	private double tol;

	/**
	 * Simple constructor with default convergency tolerance.
	 */
	public InverseIteration() {
		tol = TOL;
		maxiter = MAX_ITER;
	}

	/**
	 * Constructor giving the option to tweek the convergency criterias.
	 * @param tol change tolerance
	 * @param maxiter maximum number of iterations
	 */
	public InverseIteration(double tol, int maxiter) {
		this.tol = tol;
		this.maxiter = maxiter;
	}

	/**
	 * Solves the lowest eigenvalue and it's corresponding eigenvector for the
	 * eigenvalue problem.
	 * @param a		square matrix <code>(m x m)</code>
	 * @return result
	 * @throws NonSquareMatrixException if <code>a</code> or <code>b</code> is not square
	 * @throws DimensionMissmatchException if any input has missmatching dimensions
	 */
	public ResultStructure compute(DoubleMatrix a) {
		return compute(a, DoubleMatrix.eye(a.rows), 0, DoubleMatrix.ones(a.rows,1));
	}

	/**
	 * Solves the lowest eigenvalue and it's corresponding eigenvector for the
	 * eigenvalue problem with spectrum shift.
	 * @param a		square matrix <code>(m x m)</code>
	 * @param mu	spectrum shift
	 * @return result
	 * @throws NonSquareMatrixException if <code>a</code> or <code>b</code> is not square
	 * @throws DimensionMissmatchException if any input has missmatching dimensions
	 */
	public ResultStructure compute(DoubleMatrix a, double mu) {
		return compute(a, DoubleMatrix.eye(a.rows), mu, DoubleMatrix.ones(a.rows,1));
	}

	/**
	 * Solves the lowest eigenvalue and it's corresponding eigenvector for the
	 * eigenvalue problem with spectrum shift and initial guess mode.
	 * @param a		square matrix <code>(m x m)</code>
	 * @param mu	spectrum shift
	 * @param phi0	initial guess mode <code>(m x 1)</code>
	 * @return result
	 * @throws NonSquareMatrixException if <code>a</code> or <code>b</code> is not square
	 * @throws DimensionMissmatchException if any input has missmatching dimensions
	 */
	public ResultStructure compute(DoubleMatrix a, double mu, DoubleMatrix phi0) {
		return compute(a, DoubleMatrix.eye(a.rows), mu, phi0);
	}

	/**
	 * Solves the lowest eigenvalue and it's corresponding eigenvector for the
	 * generalized eigenvalue problem.
	 * @param a		square matrix <code>(m x m)</code>
	 * @param b		sqare matrix <code>(m x m)</code>
	 * @return result
	 * @throws NonSquareMatrixException if <code>a</code> or <code>b</code> is not square
	 * @throws DimensionMissmatchException if any input has missmatching dimensions
	 */
	public ResultStructure compute(DoubleMatrix a, DoubleMatrix b) {
		return compute(a, b, 0, DoubleMatrix.ones(a.rows,1));
	}

	/**
	 * Solves the first eigenvalue that is smaller than <code>mu</code> and
	 * it's corresponding eigenvector for the generalized eigenvalue problem
	 * with spectrum shift.
	 * @param a		square matrix <code>(m x m)</code>
	 * @param b		sqare matrix <code>(m x m)</code>
	 * @param mu	spectrum shift
	 * @return result
	 * @throws NonSquareMatrixException if <code>a</code> or <code>b</code> is not square
	 * @throws DimensionMissmatchException if any input has missmatching dimensions
	 */
	public ResultStructure compute(DoubleMatrix a, DoubleMatrix b, double mu) {
		return compute(a, b, mu, DoubleMatrix.ones(a.rows,1));
	}

	/**
	 * Solves the first eigenvalue that is smaller than <code>mu</code> and
	 * it's corresponding eigenvector for the generalized eigenvalue problem
	 * with spectrum shift.
	 * @param a		square matrix <code>(m x m)</code>
	 * @param b		sqare matrix <code>(m x m)</code>
	 * @param mu	spectrum shift
	 * @param phi0	initial guess mode <code>(m x 1)</code>
	 * @return result
	 * @throws NonSquareMatrixException if <code>a</code> or <code>b</code> is not square
	 * @throws DimensionMissmatchException if any input has missmatching dimensions
	 */
	public ResultStructure compute(DoubleMatrix a, DoubleMatrix b, double mu, DoubleMatrix phi0) {
		// Check input
		if (!a.isSquare()) {
			throw new NonSquareMatrixException("The a-matrix has to be square");
		}
		if (!b.isSquare()) {
			throw new NonSquareMatrixException("The b-matrix has to be square");
		}
		if (!a.sameSize(b)) {
			throw new DimensionMissmatchException("The a-matrix and the b-matrix has to have the same size");
		}
		if (!phi0.isColumnVector()) {
			throw new DimensionMissmatchException("The phi0-vector has to be a column vector");
		}
		if (phi0.rows != a.rows) {
			throw new DimensionMissmatchException("The phi0-vector has to have the same number of rows as the a-matrix");
		}

		// Setup iteration arrays with initial guesses.
		ArrayList<DoubleMatrix> phi = new ArrayList<DoubleMatrix>(maxiter);
		phi.add(0, phi0.dup());

		ArrayList<Double> lambda = new ArrayList<Double>(maxiter);
		lambda.add(0, -1e3);

		// Shift the a-matix
		// MATLAB: a_hat = a - mu * b;
		DoubleMatrix a_hat = a.sub(b.mul(mu));

		// Give name to variables and initiate some of them
		DoubleMatrix phi_hat, phi_s;
		ResultStructure result = new ResultStructure();
		boolean run = true;
		int s = 0;

		// Iterate
		while (run) {
			// Guess mode of step s
			// MATLAB: phi_s = phi(:,s);
			phi_s = phi.get(s);

			// Shifted guess mode of step s
			// MATLAB: phi_hat = a_hat \ (b * phi_s);
			phi_hat = Solve.solve(a_hat, b.mmul(phi_s));

			// Eigenvalue of step s+1
			// MATLAB: lambda(s+1) = (phi_hat' * b * phi_s) / (phi_hat' * b * phi_hat);
			lambda.add(s+1, phi_hat.transpose().mmul(b).mmul(phi_s).div(phi_hat.transpose().mmul(b).mmul(phi_hat)).get(0));

			// Eigenvector of step s+1
			// MATLAB: phi(:,s+1) = phi_hat / (phi_hat' * b * phi_hat)^(1/2);
			phi.add(s+1, phi_hat.div(Math.sqrt(phi_hat.transpose().mmul(b).mmul(phi_hat).get(0))));

			// Check stop criterion
			run = ((Math.abs(lambda.get(s+1)-lambda.get(s)) < tol) || (s > maxiter) ? false : true );

			// Count
			s ++;
		}

		// Shift the spectrum back again and save results.
		result.lambda = lambda.get(s) + mu;
		result.phi = phi.get(s);

		// Return results
		return result;
	}

	/**
	 * Solves the lowest eigenvalue and it's corresponding eigenvector for the
	 * generalized eigenvalue problem.
	 * @param a		square matrix <code>(m x m)</code>
	 * @param b		sqare matrix <code>(m x m)</code>
	 * @param phi0	initial guess mode <code>(m x 1)</code>
	 * @return result
	 * @throws NonSquareMatrixException if <code>a</code> or <code>b</code> is not square
	 * @throws DimensionMissmatchException if any input has missmatching dimensions
	 */
	public ResultStructure compute(DoubleMatrix a, DoubleMatrix b, DoubleMatrix phi0) {
		return compute(a, b, 0, phi0);
	}

	/**
	 * Get the maximum allowed number of iterations
	 * @return maximum number of iterations
	 */
	public int getMaxiter() {
		return maxiter;
	}

	/**
	 * Get the tolerance
	 * @return change tolerance
	 */
	public double getTol() {
		return tol;
	}

	/**
	 * Set the maximum allowed number of iterations
	 * @param maxiter maximum number of iterations
	 */
	public void setMaxiter(int maxiter) {
		this.maxiter = maxiter;
	}

	/**
	 * Set the tolerance
	 * @param tol change tolerance
	 */
	public void setTol(double tol) {
		this.tol = tol;
	}
}