package vanderbilt;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.util.HashMap;
import java.util.List;

import parsers.NewRDPNode;
import parsers.NewRDPParserFileLine;

public class WriteSparseThreeColumn
{
	private static int THRESHOLD = 50;
	
	public static void main(String[] args) throws Exception
	{
		File topDir = new File( "/projects/afodor/vanderbilt/rdpResults");
		
		HashMap<String, BufferedWriter> taxaWriters = new HashMap<String, BufferedWriter>();
		
		for( int x=1; x < NewRDPParserFileLine.TAXA_ARRAY.length; x++)
		{
			 BufferedWriter writer = new BufferedWriter(
					 	new FileWriter(new File(
					 "/projects/afodor/vanderbilt/spreadsheets/sparseThreeColumn_" + 
					 		NewRDPParserFileLine.TAXA_ARRAY[x] + ".txt")));
			 taxaWriters.put(NewRDPParserFileLine.TAXA_ARRAY[x], writer);
			 
		}
		
		for(String s : topDir.list())
		{
			System.out.println(s);
			if( s.endsWith("toRdp.txt"))
			{
				List<NewRDPParserFileLine> list = NewRDPParserFileLine.getRdpList(
					topDir.getAbsoluteFile() + File.separator + s	);
				
				for( int x=1; x < NewRDPParserFileLine.TAXA_ARRAY.length; x++)
				{
					HashMap<String, Integer> countMap = 
							getCount(NewRDPParserFileLine.TAXA_ARRAY[x], list);
					
					BufferedWriter writer = taxaWriters.get(NewRDPParserFileLine.TAXA_ARRAY[x]);
					
					for(String key: countMap.keySet())
					{
						writer.write( s.replaceAll("toRdp.txt", "") + "\t" + 
									key + "\t" + countMap.get(key) + "\n");
					}
					
					writer.flush();
				}
				
			}
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
