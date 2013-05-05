package marsagent.math.geometry;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;

/**
 * ConvexHull2d finds the convex hull of a set of points in 2D space.
 * 
 * @author 	AlexanderSehlstrom
 *
 */
public class ConvexHull2D {

	private static class XCompare implements Comparator<Point>
	{
		@Override
		public int compare(final Point o1, final Point o2)
		{
			return (new Double(o1.x)).compareTo(new Double(o2.x));
		}
	}

	/**
	 * Find the convex hull of the given points.
	 * 
	 * @param points	ArrayList with {@link Point}s to find the convex hull of
	 * @return An ordered ArrayList of points defining the convex hull.
	 */
	public static ArrayList<Point> find(final ArrayList<Point> points) {
		final ArrayList<Point> xSorted = (ArrayList<Point>) points.clone();
		Collections.sort(xSorted, new XCompare());

		final int n = xSorted.size();

		final Point[] lUpper = new Point[n];

		lUpper[0] = xSorted.get(0);
		lUpper[1] = xSorted.get(1);

		int lUpperSize = 2;

		for (int i = 2; i < n; i++)
		{
			lUpper[lUpperSize] = xSorted.get(i);
			lUpperSize++;

			while ((lUpperSize > 2) && !Point.isToTheRight(lUpper[lUpperSize - 3], lUpper[lUpperSize - 2], lUpper[lUpperSize - 1]))
			{
				// Remove the middle point of the three last
				lUpper[lUpperSize - 2] = lUpper[lUpperSize - 1];
				lUpperSize--;
			}
		}

		final Point[] lLower = new Point[n];

		lLower[0] = xSorted.get(n - 1);
		lLower[1] = xSorted.get(n - 2);

		int lLowerSize = 2;

		for (int i = n - 3; i >= 0; i--)
		{
			lLower[lLowerSize] = xSorted.get(i);
			lLowerSize++;

			while ((lLowerSize > 2) && !Point.isToTheRight(lLower[lLowerSize - 3], lLower[lLowerSize - 2], lLower[lLowerSize - 1]))
			{
				// Remove the middle point of the three last
				lLower[lLowerSize - 2] = lLower[lLowerSize - 1];
				lLowerSize--;
			}
		}

		final ArrayList<Point> result = new ArrayList<Point>();

		for (int i = 0; i < lUpperSize; i++)
		{
			result.add(lUpper[i]);
		}

		for (int i = 1; i < (lLowerSize - 1); i++)
		{
			result.add(lLower[i]);
		}

		return result;
	}
}
