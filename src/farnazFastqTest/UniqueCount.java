package farnazFastqTest;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.HashMap;

import parsers.FastQ;

public class UniqueCount
{
	public static void main(String[] args) throws Exception
	{
		int totalSeqs=0;
		HashMap<String, Integer> map = new HashMap<>();
		
		BufferedReader reader = new BufferedReader(new FileReader(new File("/scratch/afodor_research/datasets/test/806rcbc10_AN34_R1.fastq")));
		
		for( FastQ fq =  FastQ.readOneOrNull(reader); fq != null; fq =  FastQ.readOneOrNull(reader))
		{
			totalSeqs++;
			
			String key = fq.getSequence();
			
			Integer val =map.get(key);
			
			if( val ==null) 
				val =0;
			
			val++;
			
			map.put(key, val);
			totalSeqs++;
		}
		
		long sum =0;
		
		for( String s : map.keySet() )
			sum += map.get(s);
		
		System.out.println("total seqs "+ totalSeqs );
		System.out.println("Included seqs " + sum);
	}
}
