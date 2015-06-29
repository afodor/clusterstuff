/**
 * make table of how many contexts are in each comparison, and how many are shared
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
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import coPhylog.CoPhylogBinaryFileReader;
import coPhylog.ContextCount;

public class SharedContextsSingleFile {
	private static String conversionFile = "/projects/afodor_research/kwinglee/cophylog_all80chs/BiggestSRRPerStrain.txt";//file containing the biggest read file for each strain
	private static String contextPath = "/projects/afodor_research/kwinglee/cophylog_all80chs/context/";//path to the contexts
	private static String outputFile = "/projects/afodor_research/kwinglee/cophylog_all80chs/compareSingle/SingleFileContextSummaryTable.txt";
	private static int MIN_NUMBER_READS = 5;
	
	/**
	 * returns the list of files to compare, updating the HashMap to indicate which strain the file represents
	 * @param map	the map of file to strain
	 * @return		list of files to compare
	 * @throws Exception
	 */
	private static List<File> getFiles(HashMap<String, String> map) throws Exception
	{
		List<File> list = new ArrayList<File>();
		
		//get file containing the biggest read file for each strain
		BufferedReader convert = new BufferedReader(new FileReader (new File(conversionFile)));
		
		//add each file to list, updating map to indicate which strain it came from
		String line = convert.readLine();
		while(line != null) {
			String[] sp = line.split("\t");
			String file = "context"+sp[1]+"_context.gz";
			list.add(new File(contextPath + file));
			map.put(file, sp[0].replaceAll("CHS0?", ""));
			line = convert.readLine();
		}
		
		convert.close();
		
		return list;
	}
	
	public static void main(String[] args) throws Exception {
		HashMap<String, String> chsmap = new HashMap<String, String>();
		
		//get files and conversion
		List<File> list = getFiles(chsmap);
		
		//set up output
		BufferedWriter out = new BufferedWriter(new FileWriter(new File(outputFile)));
		out.write("Strain1\tFile1\tStrain2\tFile2\tNumContexts1\tNumContexts2\tNumContextsShared\tProportionShared\n");
		
		for( int x=0; x  < list.size() -1; x++)
		{
			File xFile = list.get(x);
			try {
				HashMap<Long, ContextCount> map1 = 
						CoPhylogBinaryFileReader.readBinaryFileRequireMin(xFile, MIN_NUMBER_READS);
				
				int len1 = map1.keySet().size();//number of contexts in map1
				
				for( int y=x+1 ; y < list.size(); y++)
				{
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

						out.write(chsmap.get(xFile.getName())+"\t"+xFile.getName()+"\t"+
								chsmap.get(yFile.getName())+"\t"+yFile.getName()+"\t"+
								len1+"\t"+len2+"\t"+shared+"\t"+((double)shared/(len1+len2-shared))+"\n");

						out.flush(); 
						} catch(IOException e) {
						System.err.println("Ignored: "+e);
					}
				}
				//compare to self
				File yFile = new File(contextPath + xFile.getName().replace("_1", "_2"));
				HashMap<Long, ContextCount> map2 = 
						CoPhylogBinaryFileReader.readBinaryFileRequireMin(yFile, MIN_NUMBER_READS);
				
				int len2 = map2.keySet().size();//number of contexts in map2
				int shared = 0;//number of shared contexts
				for( Long l : map1.keySet() ) {
					if( map2.containsKey(l)) {
						shared++;
					}
				}
				
				out.write(chsmap.get(xFile.getName())+"\t"+xFile.getName()+"\t"+
						chsmap.get(xFile.getName())+"\t"+yFile.getName()+"\t"+
						len1+"\t"+len2+"\t"+shared+"\t"+((double)shared/(len1+len2-shared))+"\n");
				
				out.flush();
			} catch(IOException e) {
				System.err.println("Ignored: "+e);
			}
			
		}
		out.close();
	}

}
