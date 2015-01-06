package mark;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.util.HashSet;
import java.util.List;
import java.util.Map;

import parsers.NewRDPNode;
import parsers.NewRDPParserFileLine;

public class CollectRDPLines
{
	private static File RDP_OUT_DIR =  new File("F:\\markViperBack\\rdpOut");
	private static File FILE_DUMP_DIR = new File("D:\\MarkLyteTasteManuscript\\rdpAnalysis");
	
	public static void main(String[] args) throws Exception
	{
		HashSet<String> set = getRDPSet();
		
		BufferedWriter writer = new BufferedWriter(new FileWriter(FILE_DUMP_DIR.getAbsolutePath() + 
				File.separator + "rdpLines.txt"));
		
		for(String s : set)
			writer.write( s + "\n");
		
		writer.flush();  writer.close();
	}
	
	private static HashSet<String> getRDPSet() throws Exception
	{
		HashSet<String> set =new HashSet<String>();
		
		for( String s : RDP_OUT_DIR.list())
		{
			System.out.println( s);
			List<NewRDPParserFileLine> list = 
					NewRDPParserFileLine.getRdpList(RDP_OUT_DIR.getAbsolutePath() + File.separator + 
							s);
			
			for( NewRDPParserFileLine rdp : list)
			{
				StringBuffer buff = new StringBuffer();
				
				Map<String, NewRDPNode> rdpMap = rdp.getTaxaMap();
				
				for( int x=1; x < NewRDPParserFileLine.TAXA_ARRAY.length; x++)
				{
					String level = NewRDPParserFileLine.TAXA_ARRAY[x];
					
					NewRDPNode node = rdpMap.get(level);
					
					buff.append(level.charAt(0) + "__");
					
					if( node != null)
						buff.append(node.getTaxaName());
					
					if( x <  NewRDPParserFileLine.TAXA_ARRAY.length -1 )
						buff.append(";");
				}
				
				set.add(buff.toString());
			}
		}
		
		return set;
	}
}
