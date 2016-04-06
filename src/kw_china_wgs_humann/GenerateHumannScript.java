/*
 * Generate script to create symbolic link for newly finished blasting genomes
 * then run humann with scons
 */
package kw_china_wgs_humann;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

public class GenerateHumannScript {
	public static String KEGG_DIR = "/nobackup/afodor_research/kwinglee/china/wgs/kegg_split_blastx_results_merge_filter/";
	public static String HUM_DIR = "/nobackup/afodor_research/kwinglee/humann-0.99/";
	
	public static void main(String[] args) throws IOException {
		BufferedWriter script = new BufferedWriter(new FileWriter(new File(
				"/nobackup/afodor_research/kwinglee/china/wgs/runHumScript")));
		//needs high memory and long time
		script.write("#PBS -l walltime=500:00:00\n");
		script.write("#PBS -l nodes=1:ppn=16\n");
		script.write("#PBS -W x=NODESET:ONEOF:FEATURE:ib_qdr2\n");
		//for each newly finished genome, create symbolic link
		script.write("cd " + KEGG_DIR + "\n");
		BufferedReader genList = new BufferedReader(new FileReader(new File(
				KEGG_DIR + "genomes_newly_finished.txt")));
		String gen = genList.readLine();
		while(gen != null) {
			script.write("cp kegg_merge_filter_human_" + gen + ".txt " +
					HUM_DIR + "input/.\n");
			gen = genList.readLine();
		}
		genList.close();
		
		//run scons
		script.write("cd " + HUM_DIR + "\n");
		script.write("module load python/2.7.2\n");
		script.write("scons\n");//scons -j 2 for multithreading
		script.close();
	}
	

}
