package emilyJan2018;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.util.HashMap;

public class CreateRDPQSub
{
	private static int NUM_CORES =50;

	public static final File FASTA_DIR = new File("/nobackup/afodor_research/datasets/emilyJan2018/fastaOut");
	public static final File RDP_OUT_DIR = new File( "/nobackup/afodor_research/datasets/emilyJan2018/rdpOut");
	public static final File RDP_RUN_DIR = new File("/nobackup/afodor_research/datasets/emilyJan2018/rdpScripts");
	
	public static void main(String[] args) throws Exception
	{
		String[] files = FASTA_DIR.list();
		
		HashMap<Integer, BufferedWriter> writerMap = new HashMap<Integer, BufferedWriter>();
		
		BufferedWriter allWriter = new BufferedWriter(new FileWriter(new File(
			RDP_RUN_DIR.getAbsoluteFile() + File.separator + "runAll.sh")));
				
		for( int x=0; x < NUM_CORES; x++)
		{
			File f= new File(RDP_RUN_DIR.getAbsolutePath() + File.separator + "run_" + x + ".sh"	);
			
			BufferedWriter writer = new BufferedWriter(new FileWriter(f));
			
			allWriter.write("qsub -q \"copperhead\" -N \"CountJob" + x + "\" " + f.getAbsolutePath() +  "\n"  );
			writerMap.put(x, writer);
		}
		
		allWriter.flush();  allWriter.close();
		
		int countNum =-1;
		
		for(String s : files)
			if( s.endsWith("fasta"))
		{
			countNum++;
			
			if( countNum == NUM_CORES)
				countNum =0;
			
			File fastaFile = new File(FASTA_DIR.getAbsolutePath() + File.separator + s);
			
			File rdpOutFile = new File(RDP_OUT_DIR.getAbsolutePath() + File.separator + 
					s  + "toRDP.txt");
			
			BufferedWriter writer = writerMap.get(countNum);
		
			writer.write("#PBS -l procs=1\n");
			writer.write("java -jar /users/afodor/rdp/rdp_classifier_2.10.1/dist/classifier.jar " + 
					"-o \"" + rdpOutFile.getAbsolutePath()  + "\" -q \"" + fastaFile+ "\"\n" );
					
			writer.flush();  
		}
		
		for(BufferedWriter writer : writerMap.values())
		{
			writer.flush();  writer.close();
		}
	}
}
