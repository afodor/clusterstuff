package mark.test;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.HashMap;

import mark.WriteSparseThreeColumn;
import parsers.NewRDPParserFileLine;
import parsers.OtuWrapper;

public class CheckPivotedCounts
{
	public static void main(String[] args) throws Exception
	{
		for( int x=1; x < NewRDPParserFileLine.TAXA_ARRAY.length; x++)
		{
			String level = NewRDPParserFileLine.TAXA_ARRAY[x];
			
			OtuWrapper wrapper = new OtuWrapper("/projects/afodor_research/mark/pivots/" + 
						level +"_AsColumns.txt" );
			
			for(int y=0;y  < wrapper.getSampleNames().size(); y++)
			{
				String sampleName = wrapper.getSampleNames().get(y);
				System.out.println( level + " " + sampleName);
				HashMap<String, Integer> expectedMap = getExpected(sampleName, level);
				
				for( int z=0; z < wrapper.getOtuNames().size(); z++ )
				{
					double count = wrapper.getDataPointsUnnormalized().get(y).get(z);
					
					String otuName = wrapper.getOtuNames().get(z);
					Integer expectedCount = expectedMap.get(otuName);
					
					if( count ==0 &&  expectedCount != null )
						throw new Exception("Expection 0  for " 
								+ level + " " + sampleName + "  " + otuName + " " + expectedCount);
					
					if( count != expectedCount )
						throw new Exception("Failed for " 
								+ level + " " + sampleName + "  " + otuName + " " + count + " " +  expectedCount);
					
				}	
				System.out.println("Passed " + sampleName + " " + level);
			}
		}
		
		System.out.println("global pass");
	}
	
	private static HashMap<String, Integer> getExpected(String sample, String level) throws Exception
	{
		HashMap<String, Integer> map = new HashMap<String, Integer>();
		
		BufferedReader reader = new BufferedReader(new FileReader(
				new File("/projects/afodor_research/mark/rdpOut/" + sample + "toRdp.txt")));
		
		
		for(String s = reader.readLine(); s != null; s = reader.readLine())
		{
			String[] splits = s.split("\t");
			
			for( int x=1; x < splits.length -1; x++)
			{
				if( splits[x].equals(level) )
				{
					String taxa = splits[x-1].replaceAll("\"", "");
					Integer score = (int) (Double.parseDouble(splits[x+1]) * 100.000001);
					
					if( score >= WriteSparseThreeColumn.THRESHOLD)
					{
						Integer count = map.get(taxa);
						
						if(count == null)
							count = 0;
						
						count++;
						map.put(taxa, count);
					}
				}
			}
		}
		
		reader.close();
		
		return map;
	}
	
	
}
