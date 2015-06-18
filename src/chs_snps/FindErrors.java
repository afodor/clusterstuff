/**
 * Given a path name, identify all files in the path that had an error
 * Write the file and error to a file
 * 
 * Also generate a file of the runs to repeat
 * @author kwinglee
 * 6/18/15
 */
package chs_snps;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

public class FindErrors {
	/**
	 * Takes in 2 args (the path to the list of runs and there error files and
	 * the queue to resubmit failed runs to)
	 * @param args
	 * @throws IOException
	 */
	public static void main(String[] args) throws IOException {
		if(args.length != 1) {
			System.out.println("Usage: path queue");
			System.exit(1);
		}
		File folder = new File(args[0]);
		File[] files = folder.listFiles();
		
		String path = args[0];
		if(!args[0].endsWith("/")) {
			path += "/";
		}
		
		BufferedWriter out = new BufferedWriter(new FileWriter(new File(path+"error.txt")));
		BufferedWriter sh = new BufferedWriter(new FileWriter(new File(path+"reRun.sh")));
		
		for(int i = 0; i < files.length; i++) {
			String name = files[i].getName();
			if(name.matches(".*\\.e[0-9]+$")) {
				BufferedReader br = new BufferedReader(new FileReader(files[i]));
				String line = br.readLine();
				if(line != null) {//file is not empty
					sh.write("qsub -q \"" + args[1] + " " + name.replace(".*\\.e[0-9]+$", "") + "\n");
					out.write("==========================\n"+name+":\n");
					while(line != null) {
						out.write(line+"\n");
						line = br.readLine();
					}
				}
				br.close();
			}
		}
		out.close();
		sh.close();
	}

}
