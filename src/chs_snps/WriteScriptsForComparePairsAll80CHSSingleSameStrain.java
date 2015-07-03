/**
 * Write scripts to compare two files from the same strain to look at background rate
 * Compares the forward reads to the reverse reads for the SRR with the most reads
 * @author kwinglee
 * 6/17/15
 */

package chs_snps;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;


public class WriteScriptsForComparePairsAll80CHSSingleSameStrain
{
	public static String conversionFile = "/projects/afodor_research/kwinglee/cophylog_all80chs/BiggestSRRPerStrain.txt";//file containing the biggest read file for each strain
	public static String contextPath = "/projects/afodor_research/kwinglee/cophylog_all80chs/context/";//path to the contexts
	
	/**
	 * returns the list of files to compare, updating the HashMap to indicate which strain the file represents
	 * @param map	the map of file to strain
	 * @return		list of files to compare
	 * @throws Exception
	 */
	private static List<File> getFiles(HashMap<String, String> map) throws Exception
	{
		List<File> list = new ArrayList<File>();
		
		//get file containing the biggest read file for each strain
		BufferedReader convert = new BufferedReader(new FileReader (new File(conversionFile)));
		
		//add each file to list, updating map to indicate which strain it came from
		String line = convert.readLine();
		while(line != null) {
			String[] sp = line.split("\t");
			String file = "context"+sp[1]+"_context.gz";
			list.add(new File(contextPath + file));
			map.put(file, sp[0]);
			line = convert.readLine();
		}
		
		convert.close();
		
		return list;
	}
	
	public static void main(String[] args) throws Exception
	{
		BufferedWriter allWriter = new BufferedWriter(new FileWriter(new File("/projects/afodor_research/kwinglee/cophylog_all80chs/runCompareSingleToSelf/runAll.sh")));
		
		HashMap<String, String> map = new HashMap<String, String>();
		
		int index=1;
		List<File> list = getFiles(map);
		
		for( int x=0; x  < list.size(); x++)
		{
			File xFile = list.get(x);
			File yFile = new File(xFile.getParent()+"/"+xFile.getName().replace("_1", "_2"));
				
			File outFile = new File("/projects/afodor_research/kwinglee/cophylog_all80chs/compareSingle/" + map.get(xFile.getName()) + "_to_" 
					+ map.get(xFile.getName()) + "_compare.txt");
				
			File outScriptFile = new File("/projects/afodor_research/kwinglee/cophylog_all80chs/runCompareSingleToSelf/runC" + index);
				
			/*allWriter.write("qsub -q \"Cobra_batch\" " + outScriptFile.getName() +  "\n");
				
			BufferedWriter scriptWriter = new BufferedWriter(new FileWriter(outScriptFile));
			scriptWriter.write("#PBS -l nodes=1:ppn=12\n");
			scriptWriter.write("java -cp /users/kwinglee/git/clusterstuff/bin -mx30000m " + 
					"coPhylog.WriteSNPFile " + xFile.getAbsolutePath()  + " " + yFile.getAbsolutePath() + " " + 
					outFile.getAbsolutePath());*/
			
			allWriter.write("qsub -q \"viper_batch\" " + outScriptFile.getName() +  "\n");
			
			BufferedWriter scriptWriter = new BufferedWriter(new FileWriter(outScriptFile));
			scriptWriter.write("#PBS -l nodes=1:ppn=12\n");
			scriptWriter.write("#PBS -W x=NODESET:ONEOF:FEATURE:ib_qdr\n");
			scriptWriter.write("java -cp /users/kwinglee/git/clusterstuff/bin -Xmx64g " + 
					"coPhylog.WriteSNPFile " + xFile.getAbsolutePath()  + " " + yFile.getAbsolutePath() + " " + 
					outFile.getAbsolutePath());
				
			scriptWriter.flush(); scriptWriter.close();
				
			index++;
		}
		
		allWriter.flush();  allWriter.close();
	}
}
