/*
 * Download the cirrhosis dataset from ENA
 */
package kw_machineLearning;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

public class DownloadWT2Dena {
	public static String DIR = "/nobackup/afodor_research/kwinglee/machineLearning/wt2d/";
	public static int NUMCMDS = 10;//number of downloads per script
	public static String SCRIPTDIR = "/projects/afodor_research/kwinglee/scripts/machineLearning/wt2d/";

	public static void main(String[] args) throws IOException {
		//set up scripts
		BufferedWriter runAll = new BufferedWriter(new FileWriter(new File(
				SCRIPTDIR + "downloadAllENA.sh")));
		String scriptBase = "dwnldenaWT2D_";
		BufferedWriter script = new BufferedWriter(new FileWriter(new File(
				SCRIPTDIR + scriptBase + "0")));
		int numHours = 5 * NUMCMDS;
		script.write("#PBS -l walltime=" + Integer.toString(numHours) + ":00:00\n");
		script.write("cd " + DIR + "fastqs/\n");
		runAll.write("qsub -q \"copperhead\" " + scriptBase + "0\n");
		int numScripts = 1;
		int numCmds = 0;

		BufferedReader list = new BufferedReader(new FileReader(new File(
				DIR + "PRJEB1786.txt")));//contains the list of read files
		list.readLine();//header
		for(String line = list.readLine(); line!= null; line = list.readLine()) {
			String[] sp = line.split("\t");

			String fqs = sp[11];
			String[] files = fqs.split(";");

			//add to script
			script.write("wget " + files[0] + "\n");
			script.write("wget " + files[1] + "\n");
			numCmds+=2;

			//check if need new script
			if(numCmds == NUMCMDS) {
				script.close();
				script = new BufferedWriter(new FileWriter(new File(
						SCRIPTDIR + scriptBase + numScripts)));
				script.write("#PBS -l walltime=" + Integer.toString(numHours) + ":00:00\n");
				script.write("cd " + DIR + "fastqs/\n");
				runAll.write("qsub -q \"copperhead\" " + scriptBase + numScripts + "\n");
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
