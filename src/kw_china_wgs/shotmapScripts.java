/**
 * generate the scripts to run shotmap on the China WGS sequences
 */
package kw_china_wgs;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

public class shotmapScripts {
	private static String inDir = "/projects/afodor_research/kwinglee/china/wgs/fastas/";//input data directory
	private static String sfamDB = "/projects/afodor_research/kwinglee/proteinDatabases/ShotmapSearchDB_Sfam_sequences";//search database
	private static String outDir = "/projects/afodor_research/kwinglee/china/wgs/shotmapSfamResults/";//output directory
	private static String scriptDir = "/projects/afodor_research/kwinglee/china/wgs/scripts/";
	private static String shotmap = "/users/kwinglee/git/shotmap/scripts/shotmap.pl";//location of shotmap script
	
	public static void main(String[] args) throws IOException {
		//list of fasta files to analyze -> shotmap recommends doing each separately
		File inFolder = new File(inDir);
		File[] fastas = inFolder.listFiles();
		
		//script to submit all other scripts
		BufferedWriter allWriter = new BufferedWriter(new FileWriter(new File(scriptDir + "runAllSfam.sh")));
		
		for(File f:fastas) {
			String name = f.getName();
			//write script to run shotmap on that file
			File script = new File(scriptDir + "runSfamShotmap_" + name.replace(".fa.gz", ""));
			BufferedWriter scriptWriter = new BufferedWriter(new FileWriter(script));
			scriptWriter.write("#PBS -l walltime=300:00:00\n");
			scriptWriter.write("perl " + shotmap + " -i " + inDir + name + " -d " + sfamDB +
					" -o " + outDir + "\n");
			
			//add script to full list
			allWriter.write("qsub -q \"viper_batch\" " + script.getName() +  "\n");
			
			scriptWriter.close();
		}
		allWriter.close();
	}

}
