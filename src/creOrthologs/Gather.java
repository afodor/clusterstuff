package creOrthologs;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;

import parsers.HitScores;

public class Gather
{
	public static void main(String[] args) throws Exception
	{
		HashMap<String, Float> countMap = getCountMap();
		writeFile(countMap);
	}
	
	private static List<String> getVals(HashMap<String, Float> counts ,boolean first)
			throws Exception
		{
			HashSet<String> set= new HashSet<String>();
			
			for(String s : counts.keySet())
			{
				String[] splits = s.split("@");
				
				if( splits.length != 2  )
					throw new Exception("No");
				
				if( first)
					set.add(splits[0]);
				else
					set.add(splits[1]);
			}
			
			List<String> list = new ArrayList<String>(set);
			Collections.sort(list);
			return list;
		}
	
	private static void writeFile(HashMap<String, Float> counts) throws Exception
	{
		System.out.println("Writing file");
		List<String> orthologKeys = getVals(counts, true);
		
		List<String> genomes = getVals(counts, false);
		
		BufferedWriter writer = new BufferedWriter(new FileWriter(new File(
				"/projects/afodor_research/af_broad/bitScoreOrthologsAsColumns.txt")));
		
		writer.write("genome");
		
		for(String s : orthologKeys)
			writer.write("\tLine_" + s);
		
		writer.write("\n");
		
		for(String g : genomes)
		{
			writer.write(g);
			
			for(String s : orthologKeys)
			{
				Float val = counts.get(s + "@" + g);
				
				if( val == null)
					val = 0f;
				
				writer.write("\t" + val);
			}
			
			writer.write("\n");
		}
		
		writer.flush(); writer.close();
	}
	
	
	private static HashMap<String, Float> getCountMap() throws Exception
	{
		HashMap<String, Float> countMap = new HashMap<String, Float>();
		HashMap<String, HashSet<Integer>> orthoMap = DefaultTabParser.getFileLineMap();
	
		long numDone =0;
		long numFound = 0;
		
		BufferedWriter logWriter =new BufferedWriter(new FileWriter(new File( 
				"/projects/afodor_research/af_broad/logOrthoPivot.txt")));
		
		String[] innerList = MakeBlastDB.ORTHOLOG_DIR.list();
		
		for(String d : RunBlastAll.DIRECTORIES)
		{
			File genomeDir = new File("/projects/afodor_research/af_broad" + File.separator + d);
			
			String[] list = genomeDir.list();

			for( String s : list)
			{
				if( s.endsWith("fasta"))
				{
					File outSubDir = new File( RunBlastAll.BLAST_RESULTS_PATH + File.separator + s.replaceAll(".scaffolds.fasta",""));
					
					for(String s2 : innerList)
					{
						if( s2.endsWith("geneseq"))
						{
							File outFile = new File( outSubDir.getAbsolutePath()+ File.separator + 
									s.replaceAll(".scaffolds.fasta","") + "_" + d + "_to_" + 
									s2.replaceAll("geneseq.", "") + "txt.gz");
							
							if( outFile.exists())
							{
								numFound++;
								
								List<HitScores> hsList = HitScores.getAsList(outFile,true);
								
								for( HitScores hs : hsList )
								{
									String genome = hs.getQueryId();
									String target = hs.getTargetId();
									
									if( genome.indexOf("@") != -1 ) 
										throw new Exception("No " + genome);
									
									HashSet<Integer> lineIDs = orthoMap.get(target);
									
									if( lineIDs != null)
									{
										for(Integer i : lineIDs)
										{
											String key = i +"@" + genome;
											
											Float old = countMap.get(key);
											
											if( old == null || hs.getBitScore() > old)
												countMap.put(key,hs.getBitScore());
										}
									}
									else
									{
										logWriter.write("Could not find cluster " + target + "\n");
									}
								}
								
							}
							else
							{
								logWriter.write("Could not find file " + outFile.getAbsoluteFile() + "\n");
							}
								
							numDone++;
							
							if( numDone % 1000 == 0 )
							{
								logWriter.write(numDone + " " + numFound + " " +((double)numFound / numDone) + "\n");
								logWriter.flush();
							}	
						}
					}
				}
			}
			
		}
		
		logWriter.flush();  logWriter.close();
		return countMap;
	}
}
