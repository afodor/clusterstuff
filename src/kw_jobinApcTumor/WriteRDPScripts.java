/**
 * Write scripts to run RDP on each sample
 * 12/4/15
 */

package kw_jobinApcTumor;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

public class WriteRDPScripts {
	public static final String rootFolder = "/nobackup/afodor_research/kwinglee/jobin/apcTumor/";
	public static final String fastaFolder = rootFolder + "fastas/";
	public static final String joinedFolder = rootFolder + "joinedReads/";
	public static final String rdpFolder = rootFolder + "rdpResults/";
	public static final String scriptFolder = "/projects/afodor_research/kwinglee/scripts/jobin/apcTumor/rdpScripts/";

	public static void main(String[] args) throws IOException {
		//writes script to launch all other scripts
		BufferedWriter allWriter = new BufferedWriter(new FileWriter(new File(scriptFolder + "runAll.sh")));
		
		//analyze separated reads
		File[] files = new File(fastaFolder).listFiles();
		for(int i = 0; i < files.length; i++) {
			if(!files[i].getName().contains("other")) {
				writeScripts(files[i], allWriter);
			}
		}

		//analyzed stitched reads
		File[] join = new File(joinedFolder).listFiles();
		for(File j : join) {
			if(j.getName().endsWith("join.fasta")) {
				writeScripts(j, allWriter);
			}
		}

		allWriter.close();
	}
	
	/*
	 * writes the individual script to run RDP on the given file
	 * Adds the script to allWriter
	 */
	private static void writeScripts(File file, BufferedWriter allWriter) throws IOException {
		String name = file.getName().replace(".fasta", "");
		File scriptName = new File(scriptFolder + "runRDP_" + name);
		BufferedWriter scriptWriter = new BufferedWriter(new FileWriter(scriptName));
		scriptWriter.write("#PBS -l mem=2GB\n");
		scriptWriter.write("java -Xmx2g -jar ~/rdp/RDPTools/classifier.jar classify -h " + 
				rdpFolder + "hierarch_" + name + ".txt" +
				" -o " + rdpFolder + "rdp_" + name + ".txt" +
				" " + file.getAbsolutePath());

		//add script to full list
		allWriter.write("qsub -q \"copperhead\" " + scriptName.getName() +  "\n");

		scriptWriter.close();
	}
}
