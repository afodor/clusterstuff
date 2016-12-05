/*
 * Convert the T2D dataset from SRA to fastq format using SRA toolkit
 */
package kw_machineLearning;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

public class ConvertT2DtoFastq {
	public static String DIR = "/nobackup/afodor_research/kwinglee/machineLearning/t2d/fastqs";
	public static int NUMCMDS = 10;//number of downloads per script
	public static String SCRIPTDIR = "/projects/afodor_research/kwinglee/scripts/machineLearning/t2d/";
	public static String PATH = "export PATH=$PATH:/nobackup/afodor_research/kwinglee/software/sratoolkit.2.8.0-ubuntu64/bin/";

	public static void main(String[] args) throws IOException {
		//set up scripts
		BufferedWriter runAll = new BufferedWriter(new FileWriter(new File(
				SCRIPTDIR + "convertAll.sh")));
		String scriptBase = "convertT2D_";
		BufferedWriter script = new BufferedWriter(new FileWriter(new File(
				SCRIPTDIR + scriptBase + "0")));
		int numHours = 5 * NUMCMDS;
		script.write("#PBS -l walltime=" + Integer.toString(numHours) + ":00:00\n");
		script.write("cd " + DIR + "\n");
		script.write(PATH + "\n");
		runAll.write("qsub -q \"copperhead\" " + scriptBase + "0\n");
		int numScripts = 1;
		int numCmds = 0;

		String[] list = new File(DIR).list();
		for(String file : list) {
			script.write("fastq-dump " + file + "\n");
			numCmds++;

			//check if need new script
			if(numCmds == NUMCMDS) {
				script.close();
				script = new BufferedWriter(new FileWriter(new File(
						SCRIPTDIR + scriptBase + numScripts)));
				script.write("#PBS -l walltime=" + Integer.toString(numHours) + ":00:00\n");
				script.write("cd " + DIR + "\n");
				script.write(PATH + "\n");
				runAll.write("qsub -q \"copperhead\" " + scriptBase + numScripts + "\n");
				numScripts++;
				numCmds = 0;
			}
		}

		runAll.close();
		script.close();
		System.out.println("Number of scripts: " + numScripts);
	}
}
