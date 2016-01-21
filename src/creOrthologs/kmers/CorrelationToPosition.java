package creOrthologs.kmers;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.List;
import java.util.StringTokenizer;

import utils.Pearson;
import utils.Spearman;

public class CorrelationToPosition
{
	public static final int EXPECTED_SIZE = 339;
	
	private static List<Double> reduce(List<Double> inList) throws Exception
	{
		if( inList.size() != EXPECTED_SIZE * EXPECTED_SIZE)
			throw new Exception("Unexpected size " + inList.size());
		
		List<Double> newList = new ArrayList<Double>();
		
		int index=-1;
		
		for( int x=0; x < EXPECTED_SIZE; x++)
			for(int y=0; y < EXPECTED_SIZE; y++)
			{
				index++;
				
				if( x < y)
					newList.add(inList.get(index));
			}
		
		return newList;
	}
	
	public static void main(String[] args) throws Exception
	{
		List<Double> refList = getAllPositions(GatherDistanceMatrix.GATHERED_DIR + File.separator + 
						"allDist.txt");
		refList = reduce(refList);
		
		BufferedWriter writer = new BufferedWriter(new FileWriter(new File(
			"/nobackup/afodor_research/af_broad/ConstrainedMapAllContigs.txt"	)));
		
		writer.write("contig\tstartPos\tPearson\tSpearman\n");
		
		
		for(String s : GatherDistanceMatrix.GATHERED_DIR.list())
			if( s.startsWith("klebsiella_pneumoniae_") && s.endsWith("_dist.txt"))
			{
				List<Double> otherList = getAllPositions(GatherDistanceMatrix.GATHERED_DIR + 
						File.separator + s);
				otherList = reduce(otherList);
				
				if( otherList.size() == refList.size())
				{
					String[] splits = s.split("_");
					
					writer.write( splits[4] + "\t" + splits[5] + "\t" + 
								Pearson.getPearsonR(refList, otherList) + "\t" + 
										Spearman.getSpearFromDouble(refList, otherList).getRs() + "\n");
					
					writer.flush();
				}
				else
				{
					System.out.println("Wrong size " + s +  " " + otherList.size());
				}
			}
		
		writer.flush();  writer.close();
	}
	
	private static List<Double> getAllPositions(String filepath) throws Exception
	{
		List<Double> list = new ArrayList<Double>();
		
		BufferedReader reader = new BufferedReader(new FileReader(new File(filepath)));
		
		reader.readLine();
		
		for(String s= reader.readLine(); s != null; s = reader.readLine())
		{
			StringTokenizer sToken = new StringTokenizer(s);
			
			sToken.nextToken();
			
			while( sToken.hasMoreTokens())
			{
				list.add(Double.parseDouble(sToken.nextToken()));
			}
		}
		
		return list;
	}
}
