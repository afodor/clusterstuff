package kw_creOrthologs;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;

public class WriteScriptsForConstrainedKmers
{
	
	private static final File ORTHOLOG_DIR = 
				new File( "/nobackup/afodor_research/kwinglee/cre/rbh/rbhOrthologs/chs11OrthogroupFastas");
	
	private static final File SCRIPT_DIR=
			new File("/nobackup/afodor_research/af_broad/orthologs/scripts");
	
	private static BufferedWriter makeNewWriter( BufferedWriter allWriter, int fileNum ) throws Exception
	{
		File aFile = new File(
				SCRIPT_DIR.getAbsolutePath() + File.separator + 
					"run_" + fileNum + ".sh");
			
		BufferedWriter aWriter = new BufferedWriter(new FileWriter(aFile));
		
		allWriter.write("qsub -q \"viper\" " + aFile.getAbsolutePath() + "\n");
		allWriter.flush();
		
		return aWriter;
			
	}
	
	public static void main(String[] args) throws Exception
	{
		String[] list = ORTHOLOG_DIR.list();
		
		int fileNum =0;
		int index = 0;
	
		BufferedWriter allWriter  = new BufferedWriter(new FileWriter(new File(
			SCRIPT_DIR.getAbsolutePath() + File.separator + 
			"runAll.sh")));
		
		BufferedWriter aWriter = makeNewWriter(allWriter, fileNum);
		
		for( String s : list  )
			if( s.endsWith(".fasta"))
			{
				File fastaFile = new File( ORTHOLOG_DIR.getAbsoluteFile() + File.separator + s);
				aWriter.write("java -cp /users/afodor/gitInstall/clusterstuff/bin "
						+ " kw_creOrthologs.ConstrainKMersToOrtholog "  + 
						fastaFile.getAbsolutePath() + "\n");
					
				aWriter.flush();
				index++;
				
				if( index == 40)
				{
					index=0;
					fileNum++;
					
					aWriter.flush(); aWriter.close();
					
					aWriter = makeNewWriter(allWriter, fileNum);
				}
			
		}
		
		aWriter.flush(); aWriter.close();
		allWriter.flush();  allWriter.close();
	}
}
