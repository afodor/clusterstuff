/**
 * Writes scripts that upload read files to MG-RAST
 * Also outputs a list of files for metadata spreadsheet
 */
package kw_china;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

public class mgrastUploadForWGS {

	public static void main(String[] args) throws IOException {
		String outDir = "/nobackup/afodor_research/kwinglee/china/mgrast/";
		String fastaDir = "/nobackup/afodor_research/kwinglee/china/wgs/fastas/";
		
		BufferedWriter fastaOut = new BufferedWriter(new FileWriter(new File(outDir + "wgsFastaList.txt")));//list of of fasta files
		BufferedWriter scriptAll = new BufferedWriter(new FileWriter(new File(outDir + "runAllWGS.sh")));//script to launch all other scripts
		
		File[] fastas = new File(fastaDir).listFiles();
		for(File f : fastas) {
			if(f.getName().endsWith(".fa")) {
				String name = f.getName();
				//write to list of fastas
				fastaOut.write(f.getName() + "\n");
				//set up script to upload individual file
				File scriptName = new File(outDir + "runUpload_" + name);
				BufferedWriter script = new BufferedWriter(new FileWriter(scriptName));
				script.write("curl -H \"auth: Fd3yYGeKNpv4HrBCWDSvmdHug\" -X POST -F " //authentication may need to change later
						+ "\"upload=@" + fastaDir + name +
						"\" \"http://api.metagenomics.anl.gov/1/inbox/\" ");// > " +
						//outDir + "curl_output_" + name + "\n");//don't need because output will just get print to error or output files of script
				script.close();
				//update scriptAll
				scriptAll.write("qsub -q \"Cobra\" " + scriptName.getName() +  "\n");
			}
		}
		fastaOut.close();
		scriptAll.close();
	}
}
