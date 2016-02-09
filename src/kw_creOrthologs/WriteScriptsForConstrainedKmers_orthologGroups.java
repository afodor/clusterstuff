/*
 * for orthologGroups results
 */
package kw_creOrthologs;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;

public class WriteScriptsForConstrainedKmers_orthologGroups
{
	
	private static final File ORTHOLOG_DIR = 
				new File("/nobackup/afodor_research/kwinglee/cre/rbh/rbhOrthologs/orthologGroupFastas/");
	
	private static final File SCRIPT_DIR=
			new File("/nobackup/afodor_research/af_broad/orthologs/scripts");
	
	private static BufferedWriter makeNewWriter( BufferedWriter allWriter, int fileNum ) throws Exception
	{
		File aFile = new File(
				SCRIPT_DIR.getAbsolutePath() + File.separator + 
					"runGroup_" + fileNum + ".sh");
			
		BufferedWriter aWriter = new BufferedWriter(new FileWriter(aFile));
		
		allWriter.write("qsub -q \"viper\" " + aFile.getAbsolutePath() + "\n");
		allWriter.flush();
		
		return aWriter;
			
	}
	
	public static void main(String[] args) throws Exception
	{
		//get list of orthogroups to analyze
		String[] list = new String[2406];
		BufferedReader groups = new BufferedReader(new FileReader(new File(
				"/nobackup/afodor_research/kwinglee/cre/rbh/rbhOrthologs/orthologGroups150.txt")));
		String line = groups.readLine();//header
		line = groups.readLine();
		int num = 0;
		while(line != null) {
			String[] sp = line.split("\t");
			list[num] = sp[0];
			num++;
			line = groups.readLine();
		}
		groups.close();
		
		int fileNum =0;
		int index = 0;
	
		BufferedWriter allWriter  = new BufferedWriter(new FileWriter(new File(
			SCRIPT_DIR.getAbsolutePath() + File.separator + 
			"runAll_orthologGroups.sh")));
		
		BufferedWriter aWriter = makeNewWriter(allWriter, fileNum);
		
		for( String s : list  ) {
				File fastaFile = new File( ORTHOLOG_DIR.getAbsoluteFile() + File.separator + s + ".fasta");
				aWriter.write("java -cp /users/kwinglee/git/clusterstuff/bin "
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
