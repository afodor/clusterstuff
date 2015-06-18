/**
 * Given a path name, identify all files in the path that had an error
 * Write the file and error to a file
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
	public static void main(String[] args) throws IOException {
		if(args.length != 1) {
			System.out.println("Usage: path");
			System.exit(1);
		}
		File folder = new File(args[0]);
		File[] files = folder.listFiles();
		
		String path = args[0];
		if(!args[0].endsWith("/")) {
			path += "/";
		}
		
		BufferedWriter out = new BufferedWriter(new FileWriter(new File(path+"error.txt")));
		
		for(int i = 0; i < files.length; i++) {
			String name = files[i].getName();
			if(name.matches(".*e[0-9]+$")) {
				BufferedReader br = new BufferedReader(new FileReader(files[i]));
				String line = br.readLine();
				if(line != null) {//file is not empty
					out.write("==========================\n"+name+":\n");
					while(line != null) {
						out.write(line+"\n");
						line = br.readLine();
					}
					out.flush();
				}
				br.close();
			}
		}
		out.close();
	}

}
