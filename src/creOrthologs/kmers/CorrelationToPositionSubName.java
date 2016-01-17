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

public class CorrelationToPositionSubName
{
	public static final int EXPECTED_SIZE = 339;
	public static final String SUB_NAME = "pneu";
	
	public static void main(String[] args) throws Exception
	{
		List<String> allNameList = getAllListName();
		
		if( allNameList.size() != EXPECTED_SIZE)
			throw new Exception("Unexpected size " + allNameList.size());
		
		List<Double> refList = getAllPositions(GatherDistanceMatrix.GATHERED_DIR + File.separator + 
						"allDist.txt", SUB_NAME, allNameList);
		
		BufferedWriter writer = new BufferedWriter(new FileWriter(new File(
			"/nobackup/afodor_research/af_broad/initialConstrainedMap_"+  SUB_NAME +"_only.txt"	)));
		
		writer.write("startPos\tPearson\tSpearman\n");
		
		if( refList.size() != 	EXPECTED_SIZE * EXPECTED_SIZE)
			throw new Exception("Unexpected size " + refList.size());
		
		for(String s : GatherDistanceMatrix.GATHERED_DIR.list())
			if( s.startsWith("klebsiella_pneumoniae_") && s.endsWith("_dist.txt"))
			{
				List<Double> otherList = getAllPositions(GatherDistanceMatrix.GATHERED_DIR + 
						File.separator + s, SUB_NAME, allNameList);
				
				if( otherList.size() == refList.size())
				{
					String[] splits = s.split("_");
					
					writer.write(splits[5] + "\t" + 
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
	
	private static List<String> getAllListName() throws Exception
	{
		List<String> list = new ArrayList<String>();
		
		BufferedReader reader = new BufferedReader(new FileReader(new File(
				GatherDistanceMatrix.GATHERED_DIR + File.separator + 
				"allKey.txt")));
		
		// later version of this file may have a header line
		
		for(String s = reader.readLine(); s != null; s = reader.readLine())
		{
			StringTokenizer sToken = new StringTokenizer(s);
			sToken.nextToken();
			
			list.add(sToken.nextToken());
		}
		
		return list;
	}
	
	private static List<Double> getAllPositions(String filepath, String subName, List<String> nameList) throws Exception
	{
		List<Double> list = new ArrayList<Double>();
		
		BufferedReader reader = new BufferedReader(new FileReader(new File(filepath)));
		
		reader.readLine();
		
		int xPos = -1;
		
		for(String s= reader.readLine(); s != null; s = reader.readLine())
		{
			xPos++;
			int yPos = -1;
			
			StringTokenizer sToken = new StringTokenizer(s);
			
			sToken.nextToken();
			
			while( sToken.hasMoreTokens())
			{
				yPos++;
				
				if( nameList.get(xPos).indexOf(subName) != 1 && nameList.get(yPos).indexOf(subName) != -1)				
					list.add(Double.parseDouble(sToken.nextToken()));
			}
		}
		
		return list;
	}
}
