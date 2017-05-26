import java.util.Comparator;

// comparator for sorting 'IsotopeEnvelopeStruct'
public class IsotopeEnvelopeComparator implements Comparator<IsotopeEnvelopeStruct>
{
	public IsotopeEnvelopeComparator()
	{}

	public int compare(IsotopeEnvelopeStruct ie1, IsotopeEnvelopeStruct ie2)
    {
    	if (ie1.size() == ie2.size()) // same size, prefer higher charge state
    	{
    		if (ie2.chargeState == ie1.chargeState) // same charge state also, prefer lower starting mass
				return ie1.get(0) - ie2.get(0);
    		else // otherwise, prefer higher charge state
    			return ie2.chargeState - ie1.chargeState;
    	}

    	else // otherwise, prefer larger envelope
    		return ie2.size() - ie1.size();
    }
}