package creOrthologs;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.util.HashMap;

public class MakeBlastDBGenomeDirectories
{
	private static File SCRIPT_DIR = new File( "/projects/afodor_research/af_broad/scripts/makeBlastDBScripts");
	
	public static final int NUM_CORES = 50;
	
	public static void main(String[] args) throws Exception
	{
		BufferedWriter runAllWriter = new BufferedWriter(new FileWriter(new File( 
			SCRIPT_DIR + File.separator +  "runAll.sh"	)));
		HashMap<Integer, BufferedWriter> shellMap = new HashMap<Integer, BufferedWriter>();
		
		for( int x=1; x <= 50; x++)
		{
			File outFile = new File(SCRIPT_DIR.getAbsoluteFile() + File.separator + "run_" + x + ".sh" );
			
			shellMap.put(x, new BufferedWriter(new FileWriter(outFile)));
			runAllWriter.write("qsub -q\"viper\" " + outFile.getAbsolutePath() + "\n");
		}
		
		runAllWriter.flush();  runAllWriter.close();
		
		int index = 1;
		
		for(String d : RunBlastAll.DIRECTORIES)
		{
			File genomeDir = new File("/projects/afodor_research/af_broad" + File.separator + d);
			
			String[] list = genomeDir.list();

			for( String s : list)
			{	
				if( s.endsWith("fasta"))
				{
					File inSeqs= new File( genomeDir.getAbsolutePath() + File.separator + s);
			
					BufferedWriter writer = shellMap.get(index);
				
					writer.write("module load blast\n");
					writer.write("/apps/pkg/ncbi-blast-2.2.29+/rhel6_u5-x86_64/gnu/bin/makeblastdb -dbtype nucl " + 
								"-in " + inSeqs.getAbsolutePath() + "\n");
			
					writer.flush();
			
					index = index + 1;
			
					if ( index > NUM_CORES)
						index = 1;
				}
			}
		}

		for( BufferedWriter writer : shellMap.values())
		{
			writer.flush();  writer.close();
		}
	}
}
		
