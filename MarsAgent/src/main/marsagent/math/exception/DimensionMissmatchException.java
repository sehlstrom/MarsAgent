/* Copyright 2013 AlexanderSehlstrom
 */
package marsagent.math.exception;

/**
 * Exception to throw when input has unexpected dimensions.
 * @author 	AlexanderSehlstrom
 * @since	$
 * @version	$
 *
 */
public class DimensionMissmatchException extends IllegalArgumentException {

	/**
	 * Exception to throw when input has unexpected dimensions.
	 * @param msg
	 */
	public DimensionMissmatchException(String msg) {
		super(msg);
	}
}
