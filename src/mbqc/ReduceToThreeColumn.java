package mbqc;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.util.HashMap;
import java.util.List;

import parsers.NewRDPNode;
import parsers.NewRDPParserFileLine;

public class ReduceToThreeColumn
{
	public static final File THREE_COLUMN_DIR =
		new File("/projects/afodor_research/mbqc/spreadsheets/threeColumn");

	private static int THRESHOLD = 50;
	
	public static void main(String[] args) throws Exception
	{
		if( args.length != 1)
		{
			System.out.println("Usage mbqc.ReduceToThreeColumn fileToGather"  );
			System.exit(1);
		}
		
		File file =new File(args[0]);
		
		if( !file.exists())
			throw new Exception("Could not find " + file.getAbsolutePath());
		
		writeForOne(file);
	}
	
	private static void writeForOne(File rdpFileToParse) throws Exception
	{
		HashMap<String, BufferedWriter> taxaWriters = new HashMap<String, BufferedWriter>();
		
		for( int x=1; x < NewRDPParserFileLine.TAXA_ARRAY.length; x++)
		{
			 BufferedWriter writer = new BufferedWriter(
					 	new FileWriter(new File(
					 THREE_COLUMN_DIR.getAbsolutePath() + File.separator 
					 + rdpFileToParse.getName().replace(".fastq.gz_rdpOut.txt.gz", "") 
					 			+  "_sparseThreeColumn_" + 
					 		NewRDPParserFileLine.TAXA_ARRAY[x] + ".txt")));
			 taxaWriters.put(NewRDPParserFileLine.TAXA_ARRAY[x], writer);
			 
		}
		
		List<NewRDPParserFileLine> list = NewRDPParserFileLine.getRdpListSingleThread(rdpFileToParse);
		
		String[] splits = rdpFileToParse.getName().split("_");
		String sample = splits[0] + "_" + splits[1] + "_" + splits[2];
		sample.replace(".fastq.gz", "");
		
		for( int x=1; x < NewRDPParserFileLine.TAXA_ARRAY.length; x++)
		{
			HashMap<String, Integer> countMap = 
					getCount(NewRDPParserFileLine.TAXA_ARRAY[x], list);
			
			BufferedWriter writer = taxaWriters.get(NewRDPParserFileLine.TAXA_ARRAY[x]);
			
			for(String key: countMap.keySet())
			{
				writer.write( sample+ "\t" + key + "\t" + countMap.get(key) + "\n");
			}
			
			writer.flush();
		}
		
		for(BufferedWriter writer : taxaWriters.values())
		{
			writer.flush();  writer.close();
		}
		
	}
	
	private static HashMap<String, Integer> getCount( String level, 
					List<NewRDPParserFileLine>  rdpList ) throws Exception
	{
		HashMap<String, Integer> map = new HashMap<String, Integer>();
		
		for( NewRDPParserFileLine rdp : rdpList )
		{
			NewRDPNode node = rdp.getTaxaMap().get(level);
			
			if( node != null && node.getScore() >= THRESHOLD)
			{
				Integer count = map.get(node.getTaxaName());
				
				if( count == null)
					count =0;
				
				count++;
				
				map.put(node.getTaxaName(), count);
			}
		}
		
		return map;
	}
	
}
