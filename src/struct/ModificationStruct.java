// data structure for amino acid modification
public class ModificationStruct
{
	public final int entryID;
	public final String name;
	public final double deltaMass;

	public ModificationStruct(int entryID, String name, double deltaMass)
	{
		this.entryID = entryID;
		this.name = name;
		this.deltaMass = deltaMass;
	}

	public String toString()
	{
		return name + "(" + deltaMass + ")";
	}
}