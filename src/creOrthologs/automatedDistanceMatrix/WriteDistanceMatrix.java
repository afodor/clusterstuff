package creOrthologs.automatedDistanceMatrix;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.util.List;
import java.util.StringTokenizer;

import parsers.FastaSequence;

public class WriteDistanceMatrix
{
	public static void main(String[] args) throws Exception
	{
		if( args.length != 2)
		{
			System.out.println("Usage inputalignment outputDistanceMatrix");
			System.exit(1);
		}
		
		List<FastaSequence> list = FastaSequence.readFastaFile(args[0]);
		
		BufferedWriter writer = new BufferedWriter(new FileWriter(new File(
			args[1]	)));
		
		writer.write("header1\theader\2diff\n");
		
		for(int x=0; x < list.size()-1; x++)
		{
			FastaSequence fsX = list.get(x);
			String seq1 = fsX.getSequence();
			String header1 = getHeader(fsX);
			
			for( int y=x+1; y < list.size(); y++)
			{
				FastaSequence fsY = list.get(y);
				
				writer.write( header1 + "\t" + getHeader(fsY) + "\t" + getNumDiffs(seq1, fsY.getSequence()) + "\n");
			}
		}
		
		writer.flush();  writer.close();
	}
	
	private static String getHeader(FastaSequence fs)
	{
		StringTokenizer st = new StringTokenizer(fs.getHeader());
		
		st.nextToken();
		
		return st.nextToken().replace(".scaffolds.fasta", "");
	}
	
	private static int getNumDiffs(String s1, String s2) throws Exception
	{
		int c=0;
		
		if( s1.length() != s2.length())
			throw new Exception("No");
		
		for( int x=0; x < s1.length(); x++)
			if( s1.charAt(x) != s2.charAt(x))
				c++;
		
		return c;
	}
}
