/**
 * make table of how many contexts are in each comparison, and how many are shared
 * run on the large files containing all reads
 * @author kwinglee
 * 7/7/15
 */
package chs_snps;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import coPhylog.CoPhylogBinaryFileReader;
import coPhylog.ContextCount;

public class SharedContexts {
	private static String contextPath = "/projects/afodor_research/kwinglee/cophylog_all80chs/contextCombined/";//path to the contexts
	private static String outputFile = "/projects/afodor_research/kwinglee/cophylog_all80chs/compare/ContextSummaryTable.txt";
	private static int MIN_NUMBER_READS = 5;
	
	/**
	 * returns the list of files to compare, updating the HashMap to indicate which strain the file represents
	 * @param map	the map of file to strain
	 * @return		list of files to compare
	 * @throws Exception
	 */
	private static List<File> getFiles() throws Exception
	{
		List<File> list = new ArrayList<File>();
		
		//get list of files
		File folder = new File(contextPath);
		File[] files = folder.listFiles();
		
		//add each file to list
		for(int i = 0; i < files.length; i++) {
			if(files[i].getName().endsWith("_context.gz")) {
				list.add(files[i]);
			}
		}
		return list;
	}
	
	public static void main(String[] args) throws Exception {
		//get files and conversion
		List<File> list = getFiles();
		
		//set up output
		BufferedWriter out = new BufferedWriter(new FileWriter(new File(outputFile)));
		out.write("Strain1\tStrain2\tNumContexts1\tNumContexts2\tNumContextsShared\tProportionShared\n");
		
		for( int x=0; x  < list.size() -1; x++) {
			File xFile = list.get(x);
			try {
				HashMap<Long, ContextCount> map1 = 
						CoPhylogBinaryFileReader.readBinaryFileRequireMin(xFile, MIN_NUMBER_READS);
				
				int len1 = map1.keySet().size();//number of contexts in map1
				
				for( int y=x+1 ; y < list.size(); y++){
					File yFile = list.get(y);
					
					try {
						HashMap<Long, ContextCount> map2 = 
								CoPhylogBinaryFileReader.readBinaryFileRequireMin(yFile, MIN_NUMBER_READS);

						int len2 = map2.keySet().size();//number of contexts in map2
						int shared = 0;//number of shared contexts
						for( Long l : map1.keySet() ) {
							if( map2.containsKey(l)) {
								shared++;
							}
						}

						out.write(xFile.getName().replace("_context.gz", "")+"\t"+
								yFile.getName().replace("_context.gz", "")+"\t"+
								len1+"\t"+len2+"\t"+shared+"\t"+((double)shared/(len1+len2-shared))+"\n");

						out.flush(); 
						} catch(IOException e) {
						System.err.println("Ignored: "+e);
					}
				}
			} catch(IOException e) {
				System.err.println("Ignored: "+e);
			}
			
		}
		out.close();
	}

}
