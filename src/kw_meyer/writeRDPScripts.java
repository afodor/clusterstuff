/**
 * Write scripts to run RDP on each sample
 */

package kw_meyer;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

public class writeRDPScripts {
	//directories for merged reads:
	//public static final String FASTA_FOLDER = "/nobackup/afodor_research/katieDataSep21/inFasta/";
	//public static final String RDP_FOLDER = "/nobackup/afodor_research/kwinglee/meyer/rdpResults/";
	//directories for merged reads with primers removed:
	public static final String FASTA_FOLDER = "/nobackup/afodor_research/kwinglee/meyer/fastas_merged_trimmed";
	public static final String RDP_FOLDER = "/nobackup/afodor_research/kwinglee/meyer/rdpResults_trimmed/";
	public static final String SCRIPT_FOLDER = "/projects/afodor_research/kwinglee/scripts/meyer/";

	public static void main(String[] args) throws IOException {
		//get list of fasta files to analyze
		File ffolder = new File(FASTA_FOLDER);
		File[] files = ffolder.listFiles();

		//writes script to launch all other scripts
		//BufferedWriter allWriter = new BufferedWriter(new FileWriter(new File(SCRIPT_FOLDER + "runAll.sh")));
		BufferedWriter allWriter = new BufferedWriter(new FileWriter(new File(SCRIPT_FOLDER + "runAllTrim.sh")));
		
		for(int i = 0; i < files.length; i++) {
			String name = files[i].getName().replace(".fna", "");
			//write script to run RDP on that file
			File scriptName = new File(SCRIPT_FOLDER + "runRDP_" + name);
			BufferedWriter scriptWriter = new BufferedWriter(new FileWriter(scriptName));
			scriptWriter.write("#PBS -l mem=4GB\n");
			scriptWriter.write("java -Xmx4g -jar ~/rdp/RDPTools/classifier.jar classify -h " + 
					RDP_FOLDER + "hier_" + name + ".txt" +
					" -o " + RDP_FOLDER + "rdp_" + name + ".txt" +
					" " + FASTA_FOLDER + name + ".fna");

			//add script to full list
			allWriter.write("qsub -q \"copperhead\" " + scriptName.getName() +  "\n");

			scriptWriter.close();
		}

		allWriter.close();
	}
}
