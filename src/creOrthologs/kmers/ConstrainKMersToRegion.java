package creOrthologs.kmers;

import java.io.File;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;

import creOrthologs.automatedDistanceMatrix.ScriptsForMultipleQueries;
import parsers.FastaSequence;

public class ConstrainKMersToRegion
{
	// outer key is the genome may; inner key is the k-mer
	private static HashMap<String, HashMap<String,Integer>> getBigMap() throws Exception
	{
		HashSet<String> constrainingSet = getConstrainingSet();
		
		
		for(String s : MakeKmers.KMER_DIR.list() )
			if( s.endsWith("_kmers.txt"))
			{
				
			}
		
		return null;
	}
	
	private static HashSet<String> getConstrainingSet() throws Exception
	{
		File queryFile =  ScriptsForMultipleQueries.writeOneExtractionFile(
				  ScriptsForMultipleQueries.INPUT_GENOME, "7000000220927533", 729729, 749719);
		
		List<FastaSequence> list = FastaSequence.readFastaFile(queryFile);
		
		if( list.size() != 1 )
			throw new Exception("No");
		
		HashMap<String, Integer> map = new HashMap<String,Integer>();
		
		String seq =list.get(0).getSequence();
		
		for( int x=0; x < seq.length()- MakeKmers.KMER_LENGTH; x++)
		{
			String sub = seq.substring(x, x +  MakeKmers.KMER_LENGTH);
			
			if( MakeKmers.isACGT(sub))
			{
				MakeKmers.addToMap(map, seq);
			}
		}
		
		MakeKmers.checkForNoReverse(map);
		
		return new HashSet<String>(map.keySet());
	}
}
