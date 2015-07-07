package chs_snps;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.List;
import java.util.StringTokenizer;


public class WriteScriptsForComparePairsAll80CHS
{
	private static List<File> getFiles() throws Exception
	{
		List<File> list = new ArrayList<File>();
		
		//get all files in folder
		File folder = new File("/projects/afodor_research/kwinglee/cophylog_all80chs/contextCombined/");
		File[] files = folder.listFiles();
		
		//add to list if file is a context file
		
		for(int i = 0; i < files.length; i++) {
			if(files[i].getName().endsWith("context.gz")) {
				list.add(files[i]);
			}
		}
		
		
		return list;
	}
	
	private static String getCHS(File f)
	{
		String name = f.getName();
		
		return name.replace("_context.gz", "");
	}
	
	public static void main(String[] args) throws Exception
	{
		BufferedWriter allWriter = new BufferedWriter(new FileWriter(new File("/projects/afodor_research/kwinglee/cophylog_all80chs/runCompare/runAll.sh")));
		
		int index=1;
		List<File> list = getFiles();
		
		for( int x=0; x  < list.size() -1; x++)
		{
			File xFile = list.get(x);
			
			for( int y=x+1 ; y < list.size(); y++)
			{
				File yFile = list.get(y);
				
				File outFile = new File("/projects/afodor_research/kwinglee/cophylog_all80chs/compare/" + getCHS(xFile) + "_to_" 
						+ getCHS(yFile) + "_compare.txt");
				
				File outScriptFile = new File("/projects/afodor_research/kwinglee/cophylog_all80chs/runCompare/runC_" + index);
				allWriter.write("qsub -q \"Cobra_batch\" " + outScriptFile.getName() +  "\n");
				
				BufferedWriter scriptWriter = new BufferedWriter(new FileWriter(outScriptFile));
				scriptWriter.write("#PBS -l nodes=1:ppn=12\n");
				scriptWriter.write("java -cp /users/kwinglee/git/clusterstuff/bin -Xmx36g " + 
							"coPhylog.WriteSNPFile " + xFile.getAbsolutePath()  + " " + yFile.getAbsolutePath() + " " + 
							outFile.getAbsolutePath());
				
				scriptWriter.flush(); scriptWriter.close();
				
				index++;
			}
		}
		
		allWriter.flush();  allWriter.close();
	}
}
