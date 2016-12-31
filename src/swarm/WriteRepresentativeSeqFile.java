package swarm;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.HashMap;

import parsers.FastaSequence;

public class WriteRepresentativeSeqFile
{
	public static void main(String[] args) throws Exception
	{
		HashMap<String, FastaSequence> seqMap = 
				FastaSequence.getFirstTokenSequenceMap(
						"/nobackup/afodor_research/topeOneAtATime/mergedForSwarm.txt");
		
		BufferedReader reader = new BufferedReader(new FileReader(new File(
				"nobackup/afodor_research/topeOneAtATime/OTUPlusSeqIDs.txt")));
		
		BufferedWriter writer = new BufferedWriter(new FileWriter(
				new File("nobackup/afodor_research/topeOneAtATime/repSeqsForOTUs.txt")));
		
		for(String s = reader.readLine(); s != null; s = reader.readLine())
		{
			String id = getIdentifierWithThemMost(s);
			FastaSequence seq = seqMap.get(id);
			
			if( seq == null)
				throw new Exception("Could not find " + id);
				
			writer.write(">OTU_" + s.split("\t")[0] + " " + id + "\n");	
			writer.write(seq.getSequence() + "\n");
		}
		
		writer.flush();  writer.close();
		
		reader.close();
	}
	
	private static String getIdentifierWithThemMost(String line)
		throws Exception
	{
		String[] splits = line.split("\t");
		
		if( splits.length < 2)
			throw new Exception("Parsing error " + line);
		
		Integer max = null;
		String returnVal = null;
		
		for(int x=1; x < splits.length; x++)
		{
			Integer num = Integer.parseInt(splits[x].split("_")[1]);
			
			if( max == null || num > max)
			{
				max = num;
				returnVal = splits[x];
			}
		}
		
		return returnVal;
	}
	
}
