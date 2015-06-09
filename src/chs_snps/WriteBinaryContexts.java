package chs_snps;

import java.io.BufferedOutputStream;
import java.io.BufferedWriter;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.util.HashMap;
import java.util.zip.GZIPOutputStream;

import coPhylog.ContextCount;
import coPhylog.ContextHash;
import parsers.FastaSequence;
import parsers.FastaSequenceOneAtATime;

public class WriteBinaryContexts
{
	public static void main(String[] args) throws Exception
	{
		if( args.length != 2)
		{
			System.out.println("Usage inFile outFile");
			System.exit(1);
		}
		
		process(new File(args[0]), new File(args[1]));
	}
	
	public static void process(File inFile, File outFile ) throws Exception
	{
		int numDone =0;
		long lastTime = System.currentTimeMillis();
		File logFile = new File(outFile.getAbsoluteFile() + "_LOG.txt");
		BufferedWriter logWriter = new BufferedWriter(new FileWriter(logFile));
		
		logWriter.write("numDone\ttotalMemory\tfreeMemory\tmaxMemory\tfractionFree\tfractionAllocated\ttimeSinceLast\n");
		
		
		int contextSize = 13;
		
		if(outFile.exists())
			throw new Exception("Out file exists " + outFile.getAbsolutePath());
		
		FastaSequenceOneAtATime fsoat = new FastaSequenceOneAtATime(inFile);
		
		HashMap<Long, ContextCount> map = new HashMap<Long, ContextCount>();
		
		for(FastaSequence fs = fsoat.getNextSequence(); fs != null; fs = fsoat.getNextSequence())
		{
			ContextHash.addToHash(fs.getSequence(), map, contextSize);
			numDone++;
			
			if( numDone % 10000 == 0 )
			{
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
				
				logWriter.close();

			}
		}
		
		if(outFile.exists())
			throw new Exception("Out file exists " + outFile);
		
		writeBinaryFile(outFile, map);
		logWriter.flush(); logWriter.close();
	}
	

	private static void writeBinaryFile(File outFile, HashMap<Long, ContextCount> map ) throws Exception
	{
		DataOutputStream out =new DataOutputStream( new BufferedOutputStream(new GZIPOutputStream(new FileOutputStream(outFile))));
		
		out.writeInt(map.size());
		
		for( Long l : map.keySet() )
		{
			out.writeLong(l);
			
			ContextCount cc = map.get(l);
			
			out.writeByte( cc.getAAsByte() );
			out.writeByte( cc.getCAsByte());
			out.writeByte(cc.getGAsByte());
			out.writeByte( cc.getTAsByte());
		}
		
		out.flush();  out.close();
		
	}
	
}
