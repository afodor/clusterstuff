package chs_snps;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;

import parsers.FastQ;

public class BreakInto5Million
{
	public static void main(String[] args) throws Exception
	{
		writeFiles(
			new File("/projects/afodor_chs/CHS241.IonXpress_004.fastq"), "chs241");
		
		writeFiles(new File("/projects/afodor_chs/CHS242.IonXpress_005.fastq"), "chs242");
		
	}
	
	private static void writeFiles( File inFile, String prefix )
		throws Exception
	{
		long count =0;
		int index =1;
		
		BufferedReader reader = new BufferedReader(new FileReader(inFile));
		
		File outFile = new File("/projects/afodor_chs/fasta/" + prefix + "_"
				+index);
		
		if( outFile.exists())
			throw new Exception(outFile.getAbsolutePath() + " already exists");
		
		BufferedWriter writer = new BufferedWriter(new FileWriter( outFile));
		
		for( FastQ fastq = FastQ.readOneOrNull(reader); fastq != null;
					fastq = FastQ.readOneOrNull(reader))
		{
			count = count + 1;
			
			if( count % 5000000 == 0 )
			{
				index = index + 1;
				writer.flush();  writer.close();
				
				outFile = new File("/projects/afodor_chs/fasta/" + prefix + "_"
						+index);
				
				if( outFile.exists())
					throw new Exception(outFile.getAbsolutePath() + " already exists");
				
				writer = new BufferedWriter(new FileWriter( outFile));
			}
		}
		
		writer.flush();  writer.close();
	}
	
}
