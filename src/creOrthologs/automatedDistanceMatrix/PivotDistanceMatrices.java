package creOrthologs.automatedDistanceMatrix;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.StringTokenizer;

public class PivotDistanceMatrices
{
	private static File distanceMatrix = new File("/nobackup/afodor_research/af_broad/distanceMatrices");
	
	public static void main(String[] args) throws Exception
	{
		List<String> columnNames = getAllColumnNames();
		writeResults(columnNames);
	}
	
	private static void writeResults(List<String> keyList ) throws Exception
	{
		BufferedWriter writer = new BufferedWriter(new FileWriter(new File(
				"/nobackup/afodor_research/af_broad/pivotedTo11.txt")));
		
		writer.write("region");
		
		for(String s : keyList)
			writer.write("\t" + s);
		
		writer.write("\n");
		
		String[] files = distanceMatrix.list();
		
		int numDone =1;
		for(String s : files)
		{
			numDone++;
			System.out.println("reading " +  numDone + " "+ files.length + " "+  s);
			if( s.startsWith("klebsiella_pneumoniae_chs_11.0.scaffolds.fasta_"))
			{
				writer.write(s.replaceAll("klebsiella_pneumoniae_chs_11.0.scaffolds.fasta_", ""));
				
				HashMap<String, Integer> innerMap = parseOne(
					new File(distanceMatrix.getAbsolutePath() + File.separator + s)	);
				
				for(String s2 : keyList)
				{
					Integer val = innerMap.get(s2);
					
					if( val == null)
						writer.write("\tNA");
					else
						writer.write("\t" + val);
				}
				
				writer.write("\n");
			}
		}
	
		writer.flush();  writer.close();
	}
	
	private static HashMap<String, Integer> parseOne(File file) throws Exception
	{
		HashMap<String, Integer> innerMap = new HashMap<String,Integer>();
		BufferedReader reader = new BufferedReader(new FileReader(file));
			
		reader.readLine();
			
		for(String s2 = reader.readLine(); s2 != null; s2 =reader.readLine())
		{
			StringTokenizer sToken = new StringTokenizer(s2);
				
			List<String> list = new ArrayList<String>();
			list.add(sToken.nextToken());
			list.add(sToken.nextToken());
			Collections.sort(list);
				
			String innerKey = list.get(0) + "@" + list.get(1);
				
			if( innerMap.containsKey(innerKey))
					throw new Exception("No");
				
				innerMap.put(innerKey, Integer.parseInt(sToken.nextToken()));
				
				if( sToken.hasMoreTokens())
					throw new Exception("No " + sToken.nextToken());
		}
			
		reader.close();
	
		return innerMap;
	}
	
	// outer key is region (e.g. 7000000220927531_3367_4032 ) 
	// inner key is genome1@genome2 ; inner value is count
	private static List<String> getAllColumnNames() throws Exception
	{
		HashSet<String> names = new HashSet<String>();
		
		String[] files = distanceMatrix.list();
		
		int numDone =1;
		for(String s : files)
		{
			numDone++;
			System.out.println("reading " +  numDone + " "+ files.length + " "+  s);
			if( s.startsWith("klebsiella_pneumoniae_chs_11.0.scaffolds.fasta_"))
			{
				HashMap<String, Integer> innerMap = parseOne(
					new File(distanceMatrix.getAbsolutePath() + File.separator + s)	);
				
				names.addAll(innerMap.keySet());
			}
		}
		
		List<String> list = new ArrayList<String>(names);
		Collections.sort(list);
		return list;
	}
}
