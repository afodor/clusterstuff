/*
 * Generate scripts to run kraken
 * 12/9/16
 */
package kw_jobinMiRNA;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

public class KrakenScripts {
	public static String KRAKEN_DIR = "/nobackup/afodor_research/kwinglee/software/kraken/";
	public static String STD = KRAKEN_DIR + "krakenStandardDB";
	public static String MINI = KRAKEN_DIR + "minikraken_20141208";
	public static String DIR = "/nobackup/afodor_research/kwinglee/jobin/microRNA/";
	public static String FQDIR = DIR + "adapterFiltered/";
	public static String MINI_OUT = DIR + "miniKraken/";
	public static String STD_OUT = DIR + "stdKraken/";
	public static String SCRIPT_DIR = "/projects/afodor_research/kwinglee/scripts/jobin/microRNA/kraken/";
	
	public static void main(String[] args) throws IOException {
		for(String o : new String[]{MINI_OUT, STD_OUT}) {
			File out = new File(o);
			if(!out.exists()) {
				out.mkdirs();
			}	
		}
		
		//set up scripts to run everything
		BufferedWriter scriptMini = new BufferedWriter(new FileWriter(new File(
				SCRIPT_DIR + "minikrakenMiR")));
		scriptMini.write("#PBS -l walltime=100:00:00\n");
		scriptMini.write("#PBS -l mem=10GB\n");
		scriptMini.write("#PBS -l procs=1\n");
		BufferedWriter scriptStd = new BufferedWriter(new FileWriter(new File(
				SCRIPT_DIR + "stdkrakenMiR")));
		scriptStd.write("#PBS -l walltime=200:00:00,procs=1,mem=500GB\n");
		File[] files = new File(FQDIR).listFiles();
		for(File fa : files) {
			if(!fa.getName().endsWith(".fasta")) {
				String name = fa.getName().split("\\.")[0];
				String miniName = MINI_OUT + name + ".minikraken";
				String stdName = STD_OUT + name + ".stdkraken";
				
				//mini kraken
				scriptMini.write(KRAKEN_DIR + "kraken --preload --db " 
						+ MINI + " " + fa.getAbsolutePath() 
						+ " --threads 2 > " + miniName + "\n");
				//translate output
				scriptMini.write(KRAKEN_DIR + "kraken-translate --db "
						+ MINI + " " + miniName + " > " + miniName + "_translate\n");
				scriptMini.write(KRAKEN_DIR + "kraken-translate --mpa-format --db "
						+ MINI + " " + miniName + " > " + miniName + "_mpa\n");
				
				//standard kraken
				scriptStd.write(KRAKEN_DIR + "kraken --preload --db " 
						+ STD + " " + fa.getAbsolutePath() 
						+ " --threads 2 > " + stdName + "\n");
				//translate output
				scriptStd.write(KRAKEN_DIR + "kraken-translate --db "
						+ STD + " " + stdName + " > " + stdName + "_translate\n");
				scriptStd.write(KRAKEN_DIR + "kraken-translate --mpa-format --db "
						+ STD + " " + stdName + " > " + stdName + "_mpa\n");
			}
		}
		
		scriptMini.close();
		scriptStd.close();
	}

}
