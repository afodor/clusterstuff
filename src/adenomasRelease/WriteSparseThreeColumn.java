package adenomasRelease;

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
	
	private static final File SPREADSHEET_DIR = 
			new File("/nobackup/afodor_research/adenomasRelease/spreadsheets");
	
	public static void main(String[] args) throws Exception
	{
		HashMap<String, BufferedWriter> taxaWriters = new HashMap<String, BufferedWriter>();
		
		for( int x=1; x < NewRDPParserFileLine.TAXA_ARRAY.length; x++)
		{
			 BufferedWriter writer = new BufferedWriter(
					 	new FileWriter(new File(
					SPREADSHEET_DIR.getAbsoluteFile() + File.separator + 
						"sparseThreeColumn_" + 
					 		NewRDPParserFileLine.TAXA_ARRAY[x] + ".txt")));
			 taxaWriters.put(NewRDPParserFileLine.TAXA_ARRAY[x], writer);
		}
		
		for(String s : CreateRDPQSub.RDP_OUT_DIR.list())
		{
			System.out.println(s);
			if( s.endsWith("toRdp.txt"))
			{
				List<NewRDPParserFileLine> list = NewRDPParserFileLine.getRdpListSingleThread(
					CreateRDPQSub.RDP_OUT_DIR.getAbsoluteFile() + File.separator + s	);
				
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
