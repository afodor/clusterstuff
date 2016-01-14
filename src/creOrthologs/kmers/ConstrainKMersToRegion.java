package creOrthologs.kmers;

import java.io.File;
import java.util.HashSet;
import java.util.List;

import creOrthologs.automatedDistanceMatrix.ScriptsForMultipleQueries;
import parsers.FastaSequence;

public class ConstrainKMersToRegion
{
	private static HashSet<String> getConstrainingSet() throws Exception
	{
		File queryFile =  ScriptsForMultipleQueries.writeOneExtractionFile(
				  ScriptsForMultipleQueries.INPUT_GENOME, "7000000220927533", 729729, 749719);
		
		List<FastaSequence> list = FastaSequence.readFastaFile(queryFile);
		
		if( list.size() != 1 )
			throw new Exception("No");
		
		HashSet<String> set = new HashSet<String>();
		
		//for( int )
		return null;
	}
}
