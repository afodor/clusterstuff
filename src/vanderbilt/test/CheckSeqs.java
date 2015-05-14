package vanderbilt.test;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.List;
import java.util.StringTokenizer;

public class CheckSeqs
{
	public static void main(String[] args) throws Exception
	{
		File smallFile = new File("/projects/afodor/vanderbilt/split16S/SW00248_all2.fasta");
		File largeFile = new File("/projects/afodor/vanderbilt/VanderbiltSequences_Dec52014/MSHR2/seqs_all2.fna");
		
		testAFile(smallFile, largeFile);
	}
	
	private static void testAFile(File smallFile, File bigFile) throws Exception
	{
		List<String> smallSeqs = new ArrayList<String>();
		
		BufferedReader reader = new BufferedReader(new FileReader(smallFile));
		
		String nameToken = new StringTokenizer(smallFile.getName(), "_").nextToken();
		
		for( String s= reader.readLine(); s != null; s= reader.readLine())
		{
			if( ! s.startsWith(">" + nameToken))
				throw new Exception("expected " + nameToken + " "+ s );
			
			smallSeqs.add(reader.readLine());
		}
		
		reader.close();
		
		reader = new BufferedReader(new FileReader(bigFile));
		
		int index =0;
		int numPassed = 0;
		for(String s= reader.readLine(); s != null; s = reader.readLine())
		{
			String aSeq = reader.readLine();
			if( s.startsWith(">" + nameToken) )
			{
				if( ! aSeq.equals(smallSeqs.get(index)))
					throw new Exception("Mismatch " + smallFile.getAbsolutePath() + " " + 
										bigFile.getAbsolutePath() + " " + index + 
										aSeq + " " + smallSeqs.get(index));
				else
					numPassed++;
				
				index++;
			}
		}
		
		if( index != smallSeqs.size())
			throw new Exception("Expected " + index + " " + smallSeqs.size());
		
		System.out.println("passed " + numPassed + " " + smallFile.getAbsolutePath() + " "+ bigFile.getAbsolutePath());
		
		reader.close();
	}
}
