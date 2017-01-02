package topeOneAtATime;

import parsers.FastaSequence;
import parsers.FastaSequenceOneAtATime;

public class PercentageReadsThatAreSingletons
{
	public static void main(String[] args) throws Exception
	{
		long totalCount=0;
		long singletons=0;
		
		FastaSequenceOneAtATime fsoat = new FastaSequenceOneAtATime(
			"/nobackup/afodor_research/topeOneAtATime/mergedForSwarmIncludingSingletons.txt");
		
		for(FastaSequence fs = fsoat.getNextSequence(); fs != null; fs = fsoat.getNextSequence())
		{
			Integer count = Integer.parseInt(fs.getFirstTokenOfHeader().split("_")[1]);
			
			totalCount += count;
			
			if( count == 1)
				singletons++;
		}
		
		fsoat.close();
		
		System.out.println(singletons + " " + totalCount + " " +((double)singletons)/totalCount);
	}
}
