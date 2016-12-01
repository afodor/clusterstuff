/*
 * Download the HMP dataset
 */
package kw_machineLearning;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

public class DownloadHMP {
	public static String DIR = "/nobackup/afodor_research/kwinglee/machineLearning/hmp/";
	public static int NUMCMDS = 10;//number of downloads per script
	public static String SCRIPTDIR = "/projects/afodor_research/kwinglee/scripts/machineLearning/hmp/";
	public static String WEBSITE = "http://hmpdacc.org";
	
	public static void main(String[] args) throws IOException {
		//set up scripts
		BufferedWriter runAll = new BufferedWriter(new FileWriter(new File(
				SCRIPTDIR + "downloadAll.sh")));
		BufferedWriter script = new BufferedWriter(new FileWriter(new File(
				SCRIPTDIR + "download_0")));
		int numHours = 5 * NUMCMDS;
		script.write("#PBS -l walltime=" + Integer.toString(numHours) + ":00:00\n");
		runAll.write("qsub -q \"copperhead\" download_0\n");
		int numScripts = 1;
		int numCmds = 0;
		
		BufferedReader list = new BufferedReader(new FileReader(new File(
				DIR + "HMIWGS_healthy.csv")));//contains the list of read files
		list.readLine();//header
		for(String line = list.readLine(); line!= null; line = list.readLine()) {
			line = line.replace("\"", "");
			String[] sp = line.split(",");
			
			//put files in folder depending on body site
			String site = DIR + sp[1];
			File bodySite = new File(site);
			if(!bodySite.exists()) {
				bodySite.mkdirs();
			}
			
			//add to script
			script.write("cd " + site + "\n");
			script.write("wget " + WEBSITE + sp[2] + "\n");
			numCmds++;
			
			//check if need new script
			if(numCmds == NUMCMDS) {
				script.close();
				script = new BufferedWriter(new FileWriter(new File(
						SCRIPTDIR + "download_" + numScripts)));
				script.write("#PBS -l walltime=" + Integer.toString(numHours) + ":00:00\n");
				runAll.write("qsub -q \"copperhead\" download_" + numScripts + "\n");
				numScripts++;
				numCmds = 0;
			}
		}
		
		list.close();
		runAll.close();
		script.close();
		System.out.println("Number of scripts: " + numScripts);
	}
}
