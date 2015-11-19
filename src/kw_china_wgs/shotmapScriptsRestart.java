/**
 * generate the scripts to run shotmap on the China WGS sequences
 * folder ran out space so have to restart some scripts ->
 * restart only the relevant ones that crashed and delete the previous results
 */
package kw_china_wgs;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;
import java.util.HashSet;

public class shotmapScriptsRestart {
	private static String inDir = "/projects/afodor_research/kwinglee/china/wgs/fastas/";//input data directory
	private static String sfamDB = "/projects/afodor_research/kwinglee/proteinDatabases/ShotmapSearchDB_Sfam_sequences";//search database
	private static String outDir = "/projects/afodor_research/kwinglee/china/wgs/shotmapSfamResults/";//output directory
	private static String scriptDir = "/projects/afodor_research/kwinglee/china/wgs/scripts/";
	private static String shotmap = "/users/kwinglee/git/shotmap/scripts/shotmap.pl";//location of shotmap script
	
	public static void main(String[] args) throws IOException {
		HashSet<String> stillRunning = new HashSet<String>(Arrays.asList(new String[]{"29A_1", "34A_1", "18A_1", "109A_1", "118A_1", "120A_1", "21A_1"}));//files still running
		//list of fasta files to analyze -> shotmap recommends doing each separately
		File inFolder = new File(inDir);
		File[] fastas = inFolder.listFiles();
		
		//script to submit all other scripts
		BufferedWriter allWriter = new BufferedWriter(new FileWriter(new File(scriptDir + "runAllSfam.sh")));
		
		for(File f:fastas) {
			String name = f.getName();
			String id = name.replace(".fa.ga", "");
			if(name.endsWith(".fa.gz") && !stillRunning.contains(id)) { //ignore metadata file
				//write script to run shotmap on that file
				File script = new File(scriptDir + "runSfamShotmap_" + name.replace(".fa.gz", ""));
				BufferedWriter scriptWriter = new BufferedWriter(new FileWriter(script));
				scriptWriter.write("#PBS -l walltime=300:00:00\n");
				scriptWriter.write("rm -r " + outDir + name.replace(".fa.gz", "") + "\n");
				scriptWriter.write("module load R\n");
				scriptWriter.write("perl " + shotmap + " -i " + inDir + name + " -d " + sfamDB +
						" -o " + outDir + name.replace(".fa.gz", "") + " --nprocs=2\n");
				
				//add script to full list
				allWriter.write("qsub -q \"Cobra_batch\" " + script.getName() +  "\n");
				
				scriptWriter.close();
			}
		}
		allWriter.close();
	}

}
