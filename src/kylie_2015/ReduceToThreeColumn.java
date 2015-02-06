package kylie_2015;

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
		new File("/projects/afodor_research/kylie_2015/threeColumn");

	private static int THRESHOLD = 50;
	
	public static void main(String[] args) throws Exception
	{
		String[] list = CreateRDPQSub.RDP_OUT_DIR.list();
		
		for(String s : list)
		{
			File file = new File(CreateRDPQSub.RDP_OUT_DIR.getAbsolutePath() + File.separator + s);
			System.out.println(s);
			writeForOne(file);
		}
	}
	
	private static void writeForOne(File rdpFileToParse) throws Exception
	{
		HashMap<String, BufferedWriter> taxaWriters = new HashMap<String, BufferedWriter>();
		
		for( int x=1; x < NewRDPParserFileLine.TAXA_ARRAY.length; x++)
		{
			 BufferedWriter writer = new BufferedWriter(
					 	new FileWriter(new File(
					 THREE_COLUMN_DIR.getAbsolutePath() + File.separator 
					 + rdpFileToParse.getName().replace("_TO_RDP", "") 
					 			+  "_sparseThreeColumn" + 
					 		NewRDPParserFileLine.TAXA_ARRAY[x] + ".txt")));
			 taxaWriters.put(NewRDPParserFileLine.TAXA_ARRAY[x], writer);
			 
		}
		
		List<NewRDPParserFileLine> list = NewRDPParserFileLine.getRdpListSingleThread(rdpFileToParse);
		
		String sample = rdpFileToParse.getName().replace("_TO_RDP.txt", "");
		
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
