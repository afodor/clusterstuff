/**
 * calculate snp distances from comparisons by counting number of lines in results from comparePairs scripts (which call WriteSNPFile)
 * each snp is present twice (forward and reverse), so divide final count by 2
 * @author kwinglee
 * 6/23/15
 */
package chs_snps;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

public class SNPDistances {
	
	/**
	 * Takes as input the path to a folder containing the comparisons, which have been named <name>_to_<name>_compare.txt
	 * Writes a table, SNPDistances.txt, in that folder, with the results
	 * @param args	path to folder
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		if(args.length != 1) {
			System.out.println("Usage: path");
			System.exit(1);
		}
		//set up input
		File folder = new File(args[0]);
		File[] files = folder.listFiles();
		
		//set up output
		if(!args[0].endsWith(File.separator)){
			args[0] += File.separator;
		}
		BufferedWriter out = new BufferedWriter(new FileWriter(new File(args[0]+"SNPDistances.txt")));
		out.write("xGenome\tyGenome\tDistance\n");
		
		//calculate snp distances for each file
		for(int i = 0; i < files.length; i++) {
			if(files[i].getName().endsWith("_compare.txt")) {
				BufferedReader br = new BufferedReader(new FileReader(files[i]));
				String line = br.readLine();
				int count = 0;
				while(line != null) {
					count++;
					line = br.readLine();
				}
				br.close();
				count--;//remove header
				count/=2;//divide by 2 because files contain forward and reverse sequence
				String[] name = files[i].getName().split("_");
				out.write(name[0]+"\t"+name[2]+"\t"+count+"\n");
				out.flush();
			}
		}
		
		out.close();
	}

}
