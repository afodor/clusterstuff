/*
 * get read lengths and database sizes
 */
package kw_jobinMiRNA;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;

public class GetLengths {
	public static String DIR = "/nobackup/afodor_research/kwinglee/jobin/microRNA/";
	
	public static void main(String[] args) throws IOException {
		//reads
		getReadLengths();
		
		//databases
		String[] dbs = new String[]{"/nobackup/afodor_research/kwinglee/mirbase_v21/mature.fa",
				"/nobackup/afodor_research/kwinglee/mirbase_v21/hairpin.fa",
				"/nobackup/afodor_research/kwinglee/piRBase_v1.0/piR_mouse_v1.0.fa"};
		BufferedWriter dbout = new BufferedWriter(new FileWriter(new File(
				DIR + "dbStats/dbSummary.txt")));
		dbout.write("database\tnumEntries\taveLength\tminLength\tmaxLength\n");
		for(String db : dbs) {
			double[] counts = getDBlengths(db);
			dbout.write(db);
			for(int i = 0; i < 4; i++) {
				dbout.write("\t" + counts[i]);
			}
			dbout.write("\n");
		}
		dbout.close();
		
	}

	/*
	 * returns the number of entries, average length, min length and max length
	 * also writes a file with all the lengths in case need distribution
	 */
	private static double[] getDBlengths(String db) throws IOException {
		String[] nameSplit = db.split("/");
		BufferedWriter out = new BufferedWriter(new FileWriter(new File(
				DIR + "dbStats/entryDist_" 
						+ nameSplit[nameSplit.length-1].replace(".fa", ""))));
		BufferedReader fa = new BufferedReader(new FileReader(new File(db)));
		String line = fa.readLine();//header
		line = fa.readLine();
		int min = line.length();
		int max = line.length();
		double sum = line.length();
		int numReads = 0;
		line = fa.readLine();
		while(line != null) {
			if(!line.startsWith(">")) {
				out.write(line.length() + "\n");
				numReads++;
				sum += line.length();
				if(min > line.length()) {
					min = line.length();
				}
				if(max < line.length()) {
					max = line.length();
				}
			}
			line = fa.readLine();
		}
		fa.close();
		out.close();
		return(new double[]{numReads, (sum/numReads), min, max});
	}

	private static void getReadLengths() throws IOException {
		String fqdir = DIR + "adapterFiltered/";
		String[] files = new File(fqdir).list();
		Arrays.sort(files);
		BufferedWriter out = new BufferedWriter(new FileWriter(new File(
				fqdir + "readLengths.txt")));
		out.write("sampleID\tnumReads\taveLength\tminLength\tmaxLength\n");
		for(String f : files) {
			if(f.endsWith(".fasta")) {
				BufferedReader fa = new BufferedReader(new FileReader(new File(
						fqdir + f)));
				String line = fa.readLine();//header
				line = fa.readLine();
				int min = line.length();
				int max = line.length();
				double sum = line.length();
				int numReads = 0;
				line = fa.readLine();
				while(line != null) {
					if(!line.startsWith(">")) {
						numReads++;
						sum += line.length();
						if(min > line.length()) {
							min = line.length();
						}
						if(max < line.length()) {
							max = line.length();
						}
					}
					line = fa.readLine();
				}
				fa.close();
				out.write(f.replace(".adapterfiltered.fasta", "") + "\t"
						+ numReads + "\t" + (sum / numReads) + "\t"
						+ min + "\t" + max + "\n");
			}
		}
		out.close();
	}

}
