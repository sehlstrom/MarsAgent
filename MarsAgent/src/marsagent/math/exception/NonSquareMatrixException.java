/* Copyright 2013 AlexanderSehlstrom
 */
package marsagent.math.exception;

/**
 * Exception to throw when matrix is not square.
 * @author 	AlexanderSehlstrom
 * @since	1.0
 * @version	1.0
 */
public class NonSquareMatrixException extends DimensionMissmatchException {

	/**
	 * Exception to throw when matrix is not square.
	 * @param msg
	 */
	public NonSquareMatrixException(String msg) {
		super(msg);
	}
}
