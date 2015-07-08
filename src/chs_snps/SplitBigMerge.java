/**
 * Merge the context files from the split into 1 context file to proceed with analyses
 * @author kwinglee
 * 6/29/15
 */
package chs_snps;

import java.io.BufferedOutputStream;
import java.io.BufferedWriter;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.util.HashMap;
import java.util.zip.GZIPOutputStream;

import coPhylog.CoPhylogBinaryFileReader;
import coPhylog.ContextCount;

public class SplitBigMerge {
	public static String dir = "/projects/afodor_research/kwinglee/cophylog_all80chs/splitBig/";//path to all context files
	public static int MIN_READS = 2; //minimum number of reads to keep key; since files are split in 4, singletons will never be more than 5 (cutoff for context comparison) so can discard; note the code is >=
	public static int MIN_READS_WRITE = 2; //minimum number of reads to write key; note the code is >
	
	private static void writeBinaryFile(File outFile, HashMap<Long, ContextCount> map ) throws Exception {
		DataOutputStream out =new DataOutputStream( new BufferedOutputStream(new GZIPOutputStream(new FileOutputStream(outFile))));
		
		out.writeInt(map.size());
		
		for( Long l : map.keySet() )
		{
			out.writeLong(l);
			
			ContextCount cc = map.get(l);
			
			out.writeByte(cc.getAAsByte());
			out.writeByte(cc.getCAsByte());
			out.writeByte(cc.getGAsByte());
			out.writeByte(cc.getTAsByte());
		}
		
		out.flush();  out.close();
		
	}
	
	public static void main(String[] args) throws Exception {
		//set up log file
		int numDone =0;
		long lastTime = System.currentTimeMillis();
		File logFile = new File(dir + "MergeContexts_LOG.txt");
		BufferedWriter logWriter = new BufferedWriter(new FileWriter(logFile));
				
		logWriter.write("numDone\ttotalMemory\tfreeMemory\tmaxMemory\tfractionFree\tfractionAllocated\ttimeSinceLast\n");
				
		//list of file names
		String[] big = {"SRR1159063_1", "SRR1159216_1", "SRR1159123_1", "SRR1159223_2",
				"SRR1159345_1", "SRR1159182_2", "SRR1159182_1", "SRR1159063_2", 
				"SRR1159123_2", "SRR1159216_2", "SRR1159223_1", "SRR1159345_2"};//the 12 genomes
		String[] small = {"A", "B", "C", "D"};
		
		for(int b = 0; b < big.length; b++) {
			//set up hash
			HashMap<Long, ContextCount> map = new HashMap<Long, ContextCount>();
					
			//combine files into hash
			for(int s = 0; s < small.length; s++) {
				try {
					//read file
					HashMap<Long, ContextCount> m = CoPhylogBinaryFileReader.readBinaryFileRequireMin(new File(dir+"context"+big[b]+small[s]+"_context.gz"), MIN_READS);
								
					//merge maps
					for(Long key : m.keySet()) {
						ContextCount con = m.get(key);
						if(map.containsKey(key)) {
							ContextCount mapcon = map.get(key);
							mapcon.increaseA(con.getNumA());
							mapcon.increaseT(con.getNumT());
							mapcon.increaseC(con.getNumC());
							mapcon.increaseG(con.getNumG());
							map.put(key, mapcon);
						} else {
							map.put(key, con);
						}
					}
				} catch (Exception e) {
					System.err.println(e);
				}
							
							
				//update log
				numDone++;
				System.gc();
				
				double fractionFree= 1- (Runtime.getRuntime().totalMemory()- ((double)Runtime.getRuntime().freeMemory() ))
						/Runtime.getRuntime().totalMemory();

				double fractionAllocated = 1-  (Runtime.getRuntime().maxMemory()- ((double)Runtime.getRuntime().totalMemory() ))
						/Runtime.getRuntime().maxMemory();

				logWriter.write( numDone + "\t" + Runtime.getRuntime().totalMemory() + "\t" +
						Runtime.getRuntime().freeMemory()  + "\t" + 
						Runtime.getRuntime().maxMemory() + "\t" + fractionFree  + "\t" + 
						fractionAllocated + "\t" + (System.currentTimeMillis() - lastTime) / 1000f  + "\n"
						);

				lastTime = System.currentTimeMillis();

				logWriter.flush();
			}
					
			//write results
			writeBinaryFile(new File(dir + "context" + big[b] + "_context.gz"), map);
		}
		logWriter.flush(); logWriter.close();				
	}

}
