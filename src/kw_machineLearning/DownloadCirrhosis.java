/*
 * Download the cirrhosis dataset
 */
package kw_machineLearning;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

public class DownloadCirrhosis {
	public static String DIR = "/nobackup/afodor_research/kwinglee/machineLearning/cirrhosis/";
	public static int NUMCMDS = 10;//number of downloads per script
	public static String SCRIPTDIR = "/projects/afodor_research/kwinglee/scripts/machineLearning/cirrhosis/";

	public static void main(String[] args) throws IOException {
		//set up scripts
		BufferedWriter runAll = new BufferedWriter(new FileWriter(new File(
				SCRIPTDIR + "downloadAll.sh")));
		BufferedWriter script = new BufferedWriter(new FileWriter(new File(
				SCRIPTDIR + "downloadCirr_0")));
		int numHours = 5 * NUMCMDS;
		script.write("#PBS -l walltime=" + Integer.toString(numHours) + ":00:00\n");
		script.write("cd " + DIR + "fastqs/\n");
		runAll.write("qsub -q \"copperhead\" downloadCirr_0\n");
		int numScripts = 1;
		int numCmds = 0;

		BufferedReader list = new BufferedReader(new FileReader(new File(
				DIR + "PRJEB6337.txt")));//contains the list of read files
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
						SCRIPTDIR + "downloadCirr_" + numScripts)));
				script.write("#PBS -l walltime=" + Integer.toString(numHours) + ":00:00\n");
				script.write("cd " + DIR + "fastqs/\n");
				runAll.write("qsub -q \"copperhead\" downloadCirr_" + numScripts + "\n");
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
