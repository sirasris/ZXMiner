import java.util.ArrayList;

// data structure for crosslinked peptides
public class CrosslinkStruct
{
	public final PeptideStruct peptideA, peptideB;

	public CrosslinkStruct(PeptideStruct peptideA, PeptideStruct peptideB)
	{
		this.peptideA = peptideA;
		this.peptideB = peptideB;
	}

	// for loop crosslink product
	public CrosslinkStruct(PeptideStruct peptideA)
	{
		this.peptideA = peptideA;
		this.peptideB = null; // not implemented yet
	}

	// have to add case where peptideB = null
	public CrosslinkStruct(ArrayList<ProteinStruct> proteins, int[] dataA, int[] dataB)
	{
		this.peptideA = new PeptideStruct(proteins, dataA);
		this.peptideB = new PeptideStruct(proteins, dataB);
	}

	public String toString()
	{
		String details = "peptideA: " + peptideA.toString();
		details += ", peptideB: " + peptideB.toString();

		return details;
	}
	
	public boolean isForward()
	{ return (peptideA.isForward() && peptideB.isForward()); }
	
	// generic crosslink report
	public String toReport()
	{ return peptideA.toReport() + "\tn/a\t" + peptideB.toReport() + "\tn/a";}
	
	// specific crosslink report
	public String toReport(int siteA, int siteB)
	{ return peptideA.toReport() + "\t" + siteA + "\t" + peptideB.toReport() + "\t" + siteB; }
}