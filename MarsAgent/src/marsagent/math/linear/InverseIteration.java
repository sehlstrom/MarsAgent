/**
 * 
 */
package marsagent.math.linear;

import org.jblas.DoubleMatrix;

/**
 * Solve the lowest eigenvalue <math xmlns="http://www.w3.org/1998/Math/MathML" display="inline"><mi>&lambda;<sub><mn>1</mn></sub><sup><mn>*</mn></sup></mi></math>
 * and it's corresponding eigenvector
 * <math xmlns="http://www.w3.org/1998/Math/MathML" display="inline"><mi><b>&varphi;</b><sub><mn>1</mn></sub><sup><mn>*</mn></sup></mi></math>
 * using the inverse iteration method with Rayleigh quotient.
 * 
 * <p>The <code>InverseIteration</code> solves the following generalized eigenvalue problem:</p>
 * 
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
 * <p>The eigenvalue <math xmlns="http://www.w3.org/1998/Math/MathML" display="inline"><mi>&lambda;<sub><mn>1</mn></sub><sup><mn>*</mn></sup></mi></math>
 * is approximated such that
 * <math xmlns="http://www.w3.org/1998/Math/MathML" display="inline"><mi>&lambda;<sub><mn>1</mn></sub><sup><mn>*</mn></sup></mi><mo>&gt;</mo><mi>&lambda;<sub><mn>1</mn></sub></mn></math>
 * where <math xmlns="http://www.w3.org/1998/Math/MathML" display="inline"><mi>&lambda;<sub><mn>1</mn></sub></mn></math> is the true lowest eigenvalue</p>
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
 * @since	$
 * @version	$
 *
 */
public class InverseIteration {

	/**
	 * Structure holding the computed results.
	 * @author 	AlexanderSehlstrom
	 * @since	$
	 * @version	$
	 *
	 */
	public class ResultStructure {
		/** Sought eigenvalue */
		public double lambda;

		/** Sought eigenvector */
		public DoubleMatrix phi;
	}

	/** Default maximum number of iterations */
	private static final int MAX_ITER = 10;

	/** Default change tolerance */
	private static final double TOL = 1e-3;

	/** Maximum number of iterations */
	private double maxiter;

	/** Change tolerance */
	private double tol;

	/**
	 * Simple constructor.
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
	 * generalized eigenvalue problem.
	 * @param a		matrix
	 * @param b		matrix
	 * @return result
	 */
	public ResultStructure compute(DoubleMatrix a, DoubleMatrix b) {
		return compute(a, b, 0, DoubleMatrix.ones(a.rows));
	}

	/**
	 * Solves the first eigenvalue that is smaller than <code>mu</code> and
	 * it's corresponding eigenvector for the generalized eigenvalue problem
	 * with spectrum shift.
	 * @param a		matrix
	 * @param b		matrix
	 * @param mu	spectrum shift
	 * @return result
	 */
	public ResultStructure compute(DoubleMatrix a, DoubleMatrix b, double mu) {
		return compute(a, b, mu, DoubleMatrix.ones(a.rows));
	}

	/**
	 * Solves the lowest eigenvalue and it's corresponding eigenvector for the
	 * generalized eigenvalue problem.
	 * @param a		matrix
	 * @param b		matrix
	 * @param phi0	initial guess mode.
	 * @return result
	 */
	public ResultStructure compute(DoubleMatrix a, DoubleMatrix b, DoubleMatrix phi0) {
		return compute(a, b, 0, phi0);
	}

	/**
	 * Solves the first eigenvalue that is smaller than <code>mu</code> and
	 * it's corresponding eigenvector for the generalized eigenvalue problem
	 * with spectrum shift.
	 * @param a		matrix
	 * @param b		matrix
	 * @param mu	spectrum shift
	 * @param phi0	initial guess mode.
	 * @return result
	 */
	private ResultStructure compute(DoubleMatrix a, DoubleMatrix b, double mu, DoubleMatrix phi0) {

		// TODO: Add check of matrix and vector sizes.

		ResultStructure result = new ResultStructure();
		DoubleMatrix phi = phi0.dup();
		DoubleMatrix phi_s;
		DoubleMatrix phi_hat;

		DoubleMatrix K_hat = a.sub(b.mul(mu));

		boolean run = true;
		int s = 0;
		double lambda = 0;
		double lambda_s1 = 0;

		while (run) {
			phi_s = phi;

			phi_hat = K_hat.div(b.mmul(phi_s));

			lambda_s1 = phi_hat.transpose().mmul(b).mmul(phi_s).div(phi_hat.transpose().mmul(b).mmul(phi_hat)).get(0);

			phi = phi_hat.div(Math.pow(phi_hat.transpose().mmul(b).mmul(phi_hat).get(0), 1/2));

			// Check stop criterion
			run = ((Math.abs(lambda_s1-lambda) < tol) || (s > maxiter) ? false : true );

			// Count
			s ++;
		}

		result.lambda = lambda_s1 + mu;
		result.phi = phi;

		return result;
	}
}