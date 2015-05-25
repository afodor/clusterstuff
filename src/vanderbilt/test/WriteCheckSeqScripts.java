package vanderbilt.test;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;

public class WriteCheckSeqScripts
{
	private static final String CHECK_SEQ_SCRIPT_DIR = "/projects/afodor/vanderbilt/test/checkSeqScripts";
	private static final String SPLIT_16_DIR = "/projects/afodor/vanderbilt/split16S";
	private static final String FILE_ALL1 = "/projects/afodor/vanderbilt/VanderbiltSequences_Dec52014/MSHR1/seqs_all1.fna";
	private static final String FILE_ALL2 = "/projects/afodor/vanderbilt/VanderbiltSequences_Dec52014/MSHR2/seqs_all2-CORRECTED.fna";
	private static String CLASSPATH_DIR = "/users/afodor/gitInstall/clusterstuff/bin";
	
	public static void main(String[] args) throws Exception
	{
		BufferedWriter allWriter= new BufferedWriter(new FileWriter(new File(CHECK_SEQ_SCRIPT_DIR+ File.separator + "runAll.sh")));
		
		File splitDir =new File(SPLIT_16_DIR);
		
		for( String s : splitDir.list())
		{
			File fastaFile =  new File( splitDir.getAbsolutePath() + File.separator + s);
			
			File scriptFile = new File(CHECK_SEQ_SCRIPT_DIR+ File.separator + fastaFile.getName()+ ".sh");
			
			File bigFile = null;
			
			if( s.indexOf("all1") != -1 )
				bigFile = new File(FILE_ALL1);
			else if ( s.indexOf("all2") != -1)
				bigFile = new File(FILE_ALL2);
			else throw new Exception("Unexpected " + s);
			
			BufferedWriter writer = new BufferedWriter(new FileWriter(scriptFile));
			
			writer.write("java -cp " + CLASSPATH_DIR + " vanderbilt.test.CheckSeqs " +  fastaFile.getAbsolutePath() + " "+ 
							bigFile.getAbsolutePath() + "\n");
			
			allWriter.write("qsub -q \"viper\" "  + scriptFile.getAbsolutePath() +  "\n"  );
			
			writer.flush();  writer.close();
		}
		
		allWriter.flush();  allWriter.close();
	}
}
