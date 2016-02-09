/**
 * Writes scripts that upload read files to MG-RAST
 * Also outputs a list of files for metadata spreadsheet
 */
package kw_china;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

public class mgrastUploadFromMissingFile {

	public static void main(String[] args) throws IOException {
		String outDir = "/nobackup/afodor_research/kwinglee/china/mgrast/";

		BufferedWriter scriptAll = new BufferedWriter(new FileWriter(new File(outDir + "runAllMissing.sh")));//script to launch all other scripts

		BufferedReader missing = new BufferedReader(new FileReader(new File(outDir + "MGRAST_library_missing2816.txt")));
		String line = missing.readLine(); // header
		line = missing.readLine();
		while(line != null) {
			String[] sp = line.split("\t");
			String name = sp[5];
			//set up script to upload individual file
			File scriptName = new File(outDir + "runUploadMissing_" + name);
			BufferedWriter script = new BufferedWriter(new FileWriter(scriptName));
			script.write("curl -H \"auth: sjtZ9cv4cht46qafUmRkrUsWr\" -X POST -F " //authentication may need to change later
					+ "\"upload=@/nobackup/afodor_research/ChinaSequences/rdpResults/" + name +
					".gz\" \"http://api.metagenomics.anl.gov/1/inbox/\" ");// > " +
			//outDir + "curl_output_" + name + "\n");//don't need because output will just get print to error or output files of script
			script.close();
			//update scriptAll
			scriptAll.write("qsub -q \"Cobra_batch\" " + scriptName.getName() +  "\n");
			line = missing.readLine();
		}
		missing.close();
		scriptAll.close();
	}
}
