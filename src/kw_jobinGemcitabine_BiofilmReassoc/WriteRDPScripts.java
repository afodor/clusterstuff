/**
 * Write scripts to run RDP on each sample
 * 1/5/17
 */

package kw_jobinGemcitabine_BiofilmReassoc;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

public class WriteRDPScripts {
	public static final String rootFolder = "/nobackup/afodor_research/kwinglee/jobin/gemcitabine/";
	public static final String fastaFolder = rootFolder + "demultiplexedReads/";
	public static final String rdpFolder = rootFolder + "rdpResults/";
	public static final String scriptFolder = "/projects/afodor_research/kwinglee/scripts/jobin/gemcitabine/rdpScripts/";

	public static void main(String[] args) throws IOException {
		//writes script to launch all other scripts
		BufferedWriter allWriter = new BufferedWriter(new FileWriter(new File(scriptFolder + "runAll.sh")));
		
		//analyze separated reads
		File[] files = new File(fastaFolder).listFiles();
		for(int i = 0; i < files.length; i++) {
			if(!files[i].getName().contains("other") && 
					files[i].getName().endsWith(".fasta")) {
				writeScripts(files[i], allWriter);
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
		scriptWriter.write("#PBS -l mem=2GB,procs=1\n");
		scriptWriter.write("java -Xmx2g -jar ~/rdp/RDPTools/classifier.jar classify -h " + 
				rdpFolder + "hierarch_" + name + ".txt" +
				" -o " + rdpFolder + "rdp_" + name + ".txt" +
				" " + file.getAbsolutePath());

		//add script to full list
		allWriter.write("qsub -q \"copperhead\" " + scriptName.getName() +  "\n");

		scriptWriter.close();
	}
}
