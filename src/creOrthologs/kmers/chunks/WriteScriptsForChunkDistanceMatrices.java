package creOrthologs.kmers.chunks;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;

public class WriteScriptsForChunkDistanceMatrices
{
	private static final File CHUNK_KMER_SCRIPT_DIR = 
			new File( "/nobackup/afodor_research/af_broad/scripts/runKmerChunks");
	
	public static final String GENOME_PATH= 
			"/nobackup/afodor_research/af_broad/carolina/klebsiella_pneumoniae_chs_11.0.scaffolds.fasta";
	
	public static final String CONTIG = "7000000220927538";
	
	
	private static void writeOne(BufferedWriter allWriter, int startPos, int endPos, int index) 
		throws Exception
	{
		File shFile = new File(
				CHUNK_KMER_SCRIPT_DIR.getAbsolutePath() + File.separator + "run_" + index + ".txt"	);
		
		BufferedWriter aWriter = new BufferedWriter(new FileWriter(shFile));
		
		allWriter.write("qsub -q \"viper\" " +  shFile.getAbsolutePath());
		
		aWriter.write("java -cp /users/afodor/gitInstall/clusterstuff/bin "
				+ "creOrthologs.kmers.ConstrainKMersToRegion "  + 
				GENOME_PATH+ " " + CONTIG+ " " +  startPos+ " " 
				+ endPos+ "\n");
		
		allWriter.flush(); 
		
		aWriter.flush();  aWriter.close();
	
	}
	
	public static void main(String[] args) throws Exception
	{
		
		BufferedWriter allWriter  = new BufferedWriter(new FileWriter(new File(
			CHUNK_KMER_SCRIPT_DIR.getAbsolutePath() + File.separator + 
			"runAll.sh")));
		
		BufferedReader reader = new BufferedReader(new FileReader(new File(
				"/nobackup/afodor_research/af_broad/chunks/pneuOnlyChunks_0.85_0.9.txt")));
		
		reader.readLine();
		
		int index = 1;
		String[] lastSplits= null;
		for(String s = reader.readLine() ; s != null; s = reader.readLine())
		{
			String[] splits = s.split("\t");

			if( splits.length != 4)
				throw new Exception("No");
			
			if( lastSplits != null )
			{
				writeOne(allWriter, (Integer.parseInt(lastSplits[1]) + 6000), 
					(Integer.parseInt(splits[0]) -1000) , index);
				
				index++;
			
			}
			
			writeOne( allWriter,
					Integer.parseInt(splits[0]), Integer.parseInt(splits[1]) + 5000, index);
			index++;
			lastSplits = splits;
		}
		
		
		allWriter.flush();  allWriter.close();
		reader.close();
	}
}
