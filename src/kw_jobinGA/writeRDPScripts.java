/**
 * Write scripts to run RDP on each sample
 */

package kw_jobinGA;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

public class writeRDPScripts {
	public static final String rootFolder = "/projects/afodor_research/kwinglee/jobin/ga-stool/";
	public static final String fastaFolder = rootFolder + "sequences/";
	public static final String rdpFolder = rootFolder + "rdpResults/";
	public static final String scriptFolder = rootFolder + "scripts/";

	public static void main(String[] args) throws IOException {
		//get list of fasta files to analyze
		File ffolder = new File(fastaFolder);
		File[] files = ffolder.listFiles();
		
		//writes script to launch all other scripts
		BufferedWriter allWriter = new BufferedWriter(new FileWriter(new File(scriptFolder + "runAll.sh")));
		
		for(int i = 0; i < files.length; i++) {
			String name = files[i].getName();
			
			//write script to run RDP on that file
			File scriptName = new File(scriptFolder + "runRDP_" + name);
			BufferedWriter scriptWriter = new BufferedWriter(new FileWriter(scriptName));
			scriptWriter.write("java -Xmx2g -jar ~/rdp/RDPTools/classifier.jar classify -h " + 
					rdpFolder + "hier_" + name.replace(".fasta", ".txt") + 
					" -o " + rdpFolder + "rdp_" + name.replace(".fasta", ".txt") + 
					" " + fastaFolder + name);
			
			//add script to full list
			allWriter.write("qsub -q \"Cobra_batch\" " + scriptName.getName() +  "\n");
			
			scriptWriter.close();
		}
		
		allWriter.close();
	}
}
