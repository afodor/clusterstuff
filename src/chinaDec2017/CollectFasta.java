package chinaDec2017;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.InputStreamReader;
import java.nio.file.CopyOption;
import java.nio.file.Files;
import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.zip.GZIPInputStream;

import parsers.FastQ;

public class CollectFasta
{
	public static final String FASTA_OUT = "/nobackup/afodor_research/datasets/chinaDec2017/MBiome/fastaFromRaw";
	public static final String FASTQ_OUT = "/nobackup/afodor_research/datasets/chinaDec2017/MBiome/gatheredFastq";
	
	private static final String[] DIRS_TO_SCAN = 
		{
			"/nobackup/afodor_research/datasets/chinaDec2017/MBiome/raw_1",
			"/nobackup/afodor_research/datasets/chinaDec2017/MBiome/raw_2",
			"/nobackup/afodor_research/datasets/chinaDec2017/MBiome/raw_3"
		};
	
	private static File findOneFile(File dir, String name, int read)
		throws Exception
	{
		String[] list = dir.list();
		
		File returnFile = null;
		
		for(String s : list)
		{
			if( s.endsWith(name + "_" + read +  ".fq.gz"))
			{
				if( returnFile != null)
					throw new Exception("Duplicate " + dir.getAbsolutePath() + " " +  name + " " + read);
				
				returnFile = new File(dir.getAbsolutePath() + File.separator +s );
			}
		}
		
		if( returnFile == null ) 
			throw new Exception("Could not find " +  dir.getAbsolutePath() + " " +  name + " " + read);
		
		return returnFile;
	}
	
	public static void fastqToFasta(File inFile, File outFile) throws Exception
	{
		BufferedReader reader =
				new BufferedReader(new InputStreamReader( 
						new GZIPInputStream( new FileInputStream( inFile))));
		
		BufferedWriter writer = new BufferedWriter(new FileWriter(outFile));
		
		for(FastQ fastq = FastQ.readOneOrNull(reader); fastq != null; 
				fastq = FastQ.readOneOrNull(reader))
		{
			writer.write(">" + fastq.getHeader()+ "\n");
			writer.write(fastq.getSequence() + "\n");
		}
		
		writer.flush();  writer.close();
		reader.close();
	}
	
	public static void main(String[] args) throws Exception
	{
		HashSet<File> dirs = getTopDirectories();
		
		int x=0;
		for(File f : dirs)
		{
			File forwardFile = findOneFile(f, f.getName(), 1);
			
			Files.copy(forwardFile.toPath(), 
					new File( FASTQ_OUT + File.separator + f.getName() + "_1.fastq.gz" ).toPath());
			
			//File outFasta1 = new File(FASTA_OUT + File.separator + f.getName() + "_1.fasta");
			
			//fastqToFasta(forwardFile, outFasta1);
			
			File backwardsFile = findOneFile(f, f.getName(), 2);
			

			Files.copy(backwardsFile.toPath(), 
					new File( FASTQ_OUT + File.separator + f.getName() + "_2.fastq.gz" ).toPath());
			
		
			//File outFasta2 = new File(FASTA_OUT + File.separator + f.getName() + "_2.fasta");
			
			//fastqToFasta(backwardsFile, outFasta2);
			x++;
			System.out.println(x);
			
		}
		
		System.out.println("Finished");
	}
	
	public static HashSet<File> getTopDirectories() throws Exception
	{
		HashSet<File> dirNames = new LinkedHashSet<File>();
		
		for(String dirPath : DIRS_TO_SCAN)
		{
			File topDir = new File(dirPath);
			
			String[] files = topDir.list();
			
			for( String f : files)
			{
				File aFile = new File(topDir.getAbsolutePath() + File.separator + f);
				
				if( aFile.exists() && aFile.isDirectory() )
				{
					if( dirNames.contains(aFile))
					{
						throw new Exception("Duplicate  " + f);
					}

					dirNames.add(aFile);
				}
				
			}
		}
		
		System.out.println("Finished with " +  dirNames.size() );
		return dirNames;
	}
}
