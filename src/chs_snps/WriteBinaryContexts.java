package chs_snps;

import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.HashMap;
import java.util.zip.GZIPOutputStream;

import coPhylog.ContextCount;
import coPhylog.ContextHash;
import parsers.FastaSequence;
import parsers.FastaSequenceOneAtATime;

public class WriteBinaryContexts
{
	public static void main(String[] args)
	{
		if( args.length != 2)
		{
			System.out.println("Usage inFile outFile");
			System.exit(1);
		}
	}
	
	public static void process(File inFile, File outFile ) throws Exception
	{
		int contextSize = 13;
		
		if(outFile.exists())
			throw new Exception("Out file exists ");
		
		FastaSequenceOneAtATime fsoat = new FastaSequenceOneAtATime(inFile);
		
		HashMap<Long, ContextCount> map = new HashMap<Long, ContextCount>();
		
		for(FastaSequence fs = fsoat.getNextSequence(); fs != null; fs = fsoat.getNextSequence())
		{
			ContextHash.addToHash(fs.getSequence(), map, contextSize);
		}
		
		if(outFile.exists())
			throw new Exception("Out file exists " + outFile);
		
		writeBinaryFile(outFile, map);
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
