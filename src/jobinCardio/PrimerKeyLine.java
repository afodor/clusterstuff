package jobinCardio;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.List;
import java.util.regex.Pattern;

import utils.Translate;

public class PrimerKeyLine
{
	private final int sampleIndex;
	private final String forwardKey;
	private final String reverseKey;
	private final String experiment;
	private final int experimentNum;
	private final String group;
	private final Pattern forwardPattern;
	private final Pattern reversePattern;
	
	private static final String FORWARD_PRIMER = 
			"AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT";
	
	private static final String REVERSE_PRIMER = 
			"CAAGCAGAAGACGGCATACGAGATCGGCATTCCTGCTGAACCGCTCTTCCGATCT";
	
	public int getSampleIndex()
	{
		return sampleIndex;
	}

	public String getForwardKey()
	{
		return forwardKey;
	}

	public String getReverseKey()
	{
		return reverseKey;
	}

	public String getExperiment()
	{
		return experiment;
	}

	public int getExperimentNum()
	{
		return experimentNum;
	}

	public String getGroup()
	{
		return group;
	}

	public static String getForwardPrimer()
	{
		return FORWARD_PRIMER;
	}
	
	public boolean matchesForward(String s )
	{
		s = s.substring(0, forwardKey.length() + 5);
		return forwardPattern.matcher(s).find();
	}
	
	public boolean matchesReverse(String s )
	{
		s = s.substring(0, reverseKey.length() + 5);
		return reversePattern.matcher(s).find();
	}

	private PrimerKeyLine(String s) throws Exception
	{
		String[] splits = s.split("\t");
		this.sampleIndex = Integer.parseInt(splits[0]);
		this.forwardKey = removeForwardPrimer(splits[1]);
		this.reverseKey = removeReversePrimer(splits[2]);
		this.experiment = splits[3];
		this.experimentNum = Integer.parseInt(splits[4]);
		this.group = splits[5];
		this.forwardPattern= Pattern.compile(getForwardPatternString());
		this.reversePattern = Pattern.compile(getReversePatternString());
	}
	
	private final String getReversePatternString() throws Exception
	{
		StringBuffer buff = new StringBuffer();
		
		for(int x=0; x < this.reverseKey.length(); x++)
		{
			char c = this.reverseKey.charAt(x);
			
			if( c == 'A' || c== 'C' || c == 'G' || c == 'T')
				buff.append(c);
			else if ( c == 'Y')
				buff.append("[CT]");
			else if ( c == 'R')
				buff.append("[AG]");
			else throw new Exception("No");
				
		}
		
		return buff.toString();
	}
	
	private final String getForwardPatternString() throws Exception
	{
		StringBuffer buff = new StringBuffer();
		
		for(int x=0; x < this.forwardKey.length(); x++)
		{
			char c = this.forwardKey.charAt(x);
			
			if( c == 'A' || c== 'C' || c == 'G' || c == 'T')
				buff.append(c);
			else if ( c == 'M')
				buff.append("[AC]");
			else throw new Exception("No");
				
		}
		
		return buff.toString();
	}
	
	private static String stripSpaces(String s) throws Exception
	{
		StringBuffer buff = new StringBuffer();
		
		for( int x=0; x < s.length(); x++)
			if( s.charAt(x) != ' ')
				buff.append(s.charAt(x));
		
		return buff.toString();
	}
	
	private static String removeReversePrimer(String s) throws Exception
	{
		s = stripSpaces(s);
		
		if ( ! s.startsWith(REVERSE_PRIMER))
			throw new Exception("No");
		
		return s.replace(REVERSE_PRIMER, "");
	}
	
	private static String removeForwardPrimer(String s) throws Exception
	{
		s = stripSpaces(s);
		
		if( ! s.startsWith(FORWARD_PRIMER))
			throw new Exception("No");
		
		return s.replace(FORWARD_PRIMER, "");
	}
	
	public static List<PrimerKeyLine> getList(String filepath) throws Exception
	{
		List<PrimerKeyLine> list = new ArrayList<PrimerKeyLine>();
		
		BufferedReader reader = new BufferedReader(new FileReader(new File(filepath)));
		
		reader.readLine();
		
		for(String s= reader.readLine();  s != null; s = reader.readLine())
			if( s.trim().length() > 0 )
			{
				PrimerKeyLine pkl = new PrimerKeyLine(s);
				list.add(pkl);
			}
		
		return list;
	}
	
	public static void main(String[] args) throws Exception
	{
		/*
		Pattern myPattern = Pattern.compile("[AG]ATTT[CG]");
		Matcher m = myPattern.matcher("AATTTG");
		*/
			
		/*
		List<PrimerKeyLine> aList = getList("C:\\JobinCardio\\barcode Run2 5-30-2015.txt");
		
		for(PrimerKeyLine pkl : aList)
			System.out.println(pkl.getExperiment());
			*/
		
		String aSeq = "TCATCCAGCCGCCGCGGTAATACGGAGGATCCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGAGCGTAGGCGGGTTGT" + 
"TAAGTCAGTTGTGAAAGTTTGCGGCTCAACCGTAAAATTGCAGTTGATACTGGCGACCTTGAGTGCAACAGAGGTAGGCG"+
"GAATTCGTGGTGTAGCGGTGAAATGCTTAGATATCACGAAGAACTCCGATTGCGAAGGCAGCTTACTGGATTGTAACTGA"+
"CGCTGATGCTC";
		
		System.out.println(Translate.safeReverseTranscribe(aSeq));

	}
	
}
