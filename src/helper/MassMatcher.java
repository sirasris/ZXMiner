import java.util.ArrayList;

// match a target mass against a sorted mass list in ascending order
// return all matches within a give ppm tolerance
public final class MassMatcher
{
	private MassMatcher()
	{}

	public static ArrayList<Integer> match(Double[] sorted, double target, double ppm)
	{
		ArrayList<Integer> result = new ArrayList<Integer>();
		match(sorted, target, ppm, 0, sorted.length, result);

		return result;
	}

	// sub-routine
	private static void match(Double[] sorted, double target, double ppm, int begin, int end, ArrayList<Integer> result)
	{
		if (begin >= end) // invalid interval
			return;

		// the interval is too low or too high
		if (HelperFunctions.getppmError(target, sorted[end - 1]) < -ppm || HelperFunctions.getppmError(target, sorted[begin]) > ppm)
			return;

		if (end - begin == 1 && HelperFunctions.getAbsppmError(target, sorted[begin]) < ppm) // base case n = 1
		{
			result.add(new Integer(begin));
			return;
		}

		if (end - begin == 2) // base case n = 2
		{
			if (HelperFunctions.getAbsppmError(target, sorted[begin]) < ppm)
				result.add(new Integer(begin));
			if (HelperFunctions.getAbsppmError(target, sorted[begin + 1]) < ppm)
				result.add(new Integer(begin + 1));
			return;
		}

		int mid = (end + begin) / 2; // split the problem around midpoint
		double error = HelperFunctions.getppmError(target, sorted[mid]);

		if (error < -ppm) // too low, look at the top half
		{
			match(sorted, target, ppm, mid + 1, end, result);
			return;
		}

		if (error > ppm) // too high, look at the bottom half
		{
			match(sorted, target, ppm, begin, mid, result);
			return;
		}

		if (Math.abs(error) < ppm) // look around the middle
		{
			result.add(new Integer(mid));

			int ind = mid - 1;
			double temperror = 0;

			while (ind >= begin && temperror < ppm) // look on the lower side
			{
				temperror = match(sorted, target, ppm, ind, result);
				ind--;
			}

			ind = mid + 1;
			temperror = 0;

			while (ind < end && temperror < ppm) // look on the lower side
			{
				temperror = match(sorted, target, ppm, ind, result);
				ind++;
			}
		}
	}

	private static double match(Double[] sorted, double target, double ppm, int ind, ArrayList<Integer> result)
	{
		double error = HelperFunctions.getAbsppmError(target, sorted[ind]);

		if (error < ppm)
			result.add(new Integer(ind));

		return error;
	}
}