package marsagent.math.geometry;

import java.util.Vector;

/**
 * A point is a primitive notation and is zero-dimensional.
 * 
 * @author 	AlexanderSehlstrom
 * @since	$
 * @version	$
 *
 */
public class Point {

	public static Point[] convexHull(final Point[] s) {
		Point pointOnHull;

		// Find left most point in S
		pointOnHull = s[0];
		boolean notfound = true;
		int j = 0;
		for (final Point element : s) {
			while (notfound && (j < s.length)) {
				if (s[j].getX() < pointOnHull.getX()) {
					pointOnHull = s[j];
					notfound = false;
				}
				j++;
			}

			j = 0;
			notfound = true;
		}

		final Vector<Point> c = new Vector();

		// Find convex hull

		//position = sign( (Bx-Ax)*(Y-Ay) - (By-Ay)*(X-Ax) )

		//Math.signum( ()*() - ()*() );

		Point endpoint;
		int i = 0;
		double position;
		do {
			c.add(i, pointOnHull);
			endpoint = s[0];
			for (j = 0; j < (s.length-1); j++) {
				position = Math.signum( ((endpoint.getX()-c.get(i).getX())*(s[j].getY()-c.get(i).getY())) - ((endpoint.getY()-c.get(i).getY())*(s[j].getX()-c.get(i).getX())) );
				if ((endpoint == pointOnHull) || (position == -1)) {
					endpoint = s[j];
				}
			}
			i++;
			pointOnHull = endpoint;
		} while (endpoint != c.get(0));

		final Point[] carray = new Point[c.size()];
		return c.toArray(carray);
	}

	/**
	 * Tests if {@link Point} c is on the the
	 * line going from a to b. 2D space is considered,
	 * i.e. only the x- and y-coordinates will be used.
	 * 
	 * @param a Start point of line
	 * @param b End point of line
	 * @param c Point to test
	 * @return <code>true</code> if the point c is on the line between a and b, <code>false</code> otherwise
	 */
	public static boolean isInLine(final Point a, final Point b, final Point c)
	{
		return (((b.x - a.x)*(c.y - a.y)) - ((b.y - a.y)*(c.x - a.x))) == 0;
	}

	/**
	 * Tests if {@link Point} c is to the left of the
	 * line going from a to b. 2D space is considered,
	 * i.e. only the x- and y-coordinates will be used.
	 * 
	 * <p>Examples:
	 * <pre>
	 * +-----------------+-------------------+
	 * | false:    b     | true:      b      |
	 * |           o     |            o      |
	 * |          /      |           /       |
	 * |         /    o  |    o     /        |
	 * |        /     c  |    c    /         |
	 * |       o         |        o          |
	 * |       a         |        a          |
	 * +-----------------+-------------------+
	 * | true:     a     | false:     a      |
	 * |           o     |            o      |
	 * |          /      |           /       |
	 * |         /    o  |    o     /        |
	 * |        /     c  |    c    /         |
	 * |       o         |        o          |
	 * |       b         |        b          |
	 * +-----------------+-------------------+
	 * </pre>
	 * </p>
	 * @param a Start point of line
	 * @param b End point of line
	 * @param c Point to test
	 * @return <code>true</code> if the point c is to the left of the line from a to b, <code>false</code> otherwise
	 */
	public static boolean isToTheLeft(final Point a, final Point b, final Point c)
	{
		return (((b.x - a.x)*(c.y - a.y)) - ((b.y - a.y)*(c.x - a.x))) < 0;
	}

	/**
	 * Tests if {@link Point} c is to the right of the
	 * line going from a to b. 2D space is considered,
	 * i.e. only the x- and y-coordinates will be used.
	 * 
	 * <p>Examples:
	 * <pre>
	 * +-----------------+-------------------+
	 * | true:     b     | false:     b      |
	 * |           o     |            o      |
	 * |          /      |           /       |
	 * |         /    o  |    o     /        |
	 * |        /     c  |    c    /         |
	 * |       o         |        o          |
	 * |       a         |        a          |
	 * +-----------------+-------------------+
	 * | false:    a     | true:      a      |
	 * |           o     |            o      |
	 * |          /      |           /       |
	 * |         /    o  |    o     /        |
	 * |        /     c  |    c    /         |
	 * |       o         |        o          |
	 * |       b         |        b          |
	 * +-----------------+-------------------+
	 * </pre>
	 * </p>
	 * @param a Start point of line
	 * @param b End point of line
	 * @param c Point to test
	 * @return <code>true</code> if the point is to the right of the line from a to b, <code>false</code> otherwise
	 */
	public static boolean isToTheRight(final Point a, final Point b, final Point c)
	{
		return (((b.x - a.x)*(c.y - a.y)) - ((b.y - a.y)*(c.x - a.x))) > 0;
	}

	/** x-coordinate component */
	double x;

	/** y-coordinate component */
	double y;

	/** z-coordinate component */
	double z;

	/**
	 * Simple constructor.
	 * 
	 * @param x	x-coordinate
	 * @param y y-coordinate
	 * @param z z-coordinate
	 */
	public Point(final double x, final double y, final double z) {
		this.x = x;
		this.y = y;
		this.z = z;
	}

	/**
	 * @return the x
	 */
	public double getX() {
		return x;
	}

	/**
	 * @return the y
	 */
	public double getY() {
		return y;
	}

	/**
	 * @return the z
	 */
	public double getZ() {
		return z;
	}

	/**
	 * @param x the x to set
	 */
	public void setX(final double x) {
		this.x = x;
	}

	/**
	 * @param y the y to set
	 */
	public void setY(final double y) {
		this.y = y;
	}

	/**
	 * @param z the z to set
	 */
	public void setZ(final double z) {
		this.z = z;
	}

	@Override
	public String toString() {
		return "p=(" +x+ ", " + y + ", " + z + ")";
	}
}
