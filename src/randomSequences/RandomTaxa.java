package randomSequences;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.util.Random;

public class RandomTaxa
{
	private static final String OUT_DIR = "/users/afodor/testSeqs/";
	
	//private static final String OUT_DIR = "c:\\temp";
	
	private static final int SEQUENCE_LENGTH = 150;
	
	private static final float FRACTION_DOMINANT_NUCLEOTIDE = 0.5f;
	
	private static final int NUM_SEQUENCES = 100;
	
	private static final Random RANDOM = new Random();
	
	private static final char[] DNA = { 'A', 'C', 'G', 'T' } ;
	
	public static void main(String[] args) throws Exception
	{
		BufferedWriter writer = new BufferedWriter(new FileWriter(new File(OUT_DIR + File.separator + "randomSeq.txt")));
		
		for( int x=0; x < NUM_SEQUENCES; x++)
		{
			for( int y=0; y < DNA.length; y++)
			{
				writer.write( "TAXA_" +  DNA[y] + "\tseq_" + DNA[y] + "_" +  x + "\t"  );
				
				StringBuffer buff = new StringBuffer();
				
				for( int z = 0; z < SEQUENCE_LENGTH; z++)
				{
					buff.append(DNA[ getANucleotide(y) ]);
					
				}
				
				writer.write( buff.toString() + "\n" );
			}
		}
		
		writer.flush();  writer.close();
	}
	
	private static int getANucleotide(int num) throws Exception
	{
		if( RANDOM.nextFloat() <= FRACTION_DOMINANT_NUCLEOTIDE)
			return num;
		
		int returnVal = num;
		
		while(returnVal == num)
		{
			returnVal = RANDOM.nextInt(DNA.length);
		}
		
		return returnVal;
	}
}
