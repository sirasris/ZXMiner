import java.util.Comparator;

// generalized comparator method that return the indexes of the sorted data
// default to ascending order
public class IndexedComparator implements Comparator<Integer>
{
    private final Object[] reference;

    public IndexedComparator(Object[] reference)
    { this.reference = reference; }

	// initialize the index array
    public Integer[] initIndexArray()
    {
        Integer[] index = new Integer[reference.length];

        for (int i = 0; i < reference.length; i++)
            index[i] = i;

        return index;
    }

    public int compare(Integer i1, Integer i2)
    {
    	if (reference[0] instanceof Integer)
    		return ((Integer) reference[i1]).compareTo((Integer) reference[i2]);
    	if (reference[0] instanceof Double)
    		return ((Double) reference[i1]).compareTo((Double) reference[i2]);
    	if (reference[0] instanceof String)
    		return ((String) reference[i1]).compareTo((String) reference[i2]);

		// warn user for undeterminable class type
		HelperFunctions.debug("IndexedComparator", "Unknown reference array class type");
    	return 0;
    }
}