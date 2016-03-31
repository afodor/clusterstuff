package creOrthologs.kmers;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.StringTokenizer;

import utils.Spearman;

public class WriteSpearmanFromRandom
{
	private static final File KMER_DIST_DIRECTORY = 
			new File("/nobackup/afodor_research/af_broad/randomKMerMatrices");
	
	public static final int EXPECTED_NUM_LINES = 340;
	
	public static List<String> getGenomeNames() throws Exception
	{
		String[] list = KMER_DIST_DIRECTORY.list();
		
		List<String> returnVals = null;
		
		for( String s : list)
		{
			if( s.endsWith("key.txt"))
			{
				File f = new File(KMER_DIST_DIRECTORY.getAbsolutePath() + File.separator + 
											s);
				
				if(returnVals == null)
				{
					returnVals = getNames(f);
				}
				else
				{
					if( ! returnVals.equals(getNames(f)))
						throw new Exception("No " + f.getAbsolutePath());
				}
			}
		}
		
		return returnVals;
	}
	
	public static List<String> getNames(File file) throws Exception
	{
		List<String> list = new ArrayList<String>();
		BufferedReader reader = new BufferedReader(new FileReader(file));
		
		reader.readLine();
		
		for(String s= reader.readLine() ; s != null; s = reader.readLine())
		{
			StringTokenizer sToken = new StringTokenizer(s);
			sToken.nextToken();
			list.add(new String(sToken.nextToken()));
		}
		
		reader.close();
		return list;
	}
	
	public static List<Float> getValsOrNull(File f, HashSet<Integer> include) throws Exception
	{
		List<Float> vals = new ArrayList<Float>();
		
		BufferedReader reader = new BufferedReader(new FileReader(f));
		
		reader.readLine();
		
		int numRead =0;
		for(String s=  reader.readLine(); s != null; s= reader.readLine())
		{
			if( include == null ||  include.contains(numRead))
			{
				int tokenNumber = 0;
				
				StringTokenizer sToken = new StringTokenizer(s);

				sToken.nextToken();
				
				while(sToken.hasMoreTokens())
				{
					float aVal = Float.parseFloat(sToken.nextToken());
					
					if( tokenNumber > numRead && ( include == null ||  include.contains(tokenNumber)) )
						vals.add(aVal);
					
					tokenNumber++;
				}
			}
					
			numRead++;
		}
		
		reader.close();
		
		if( numRead != EXPECTED_NUM_LINES-1)
		{
			System.out.println(numRead + " skipping " + f.getAbsolutePath());
			return null;
		}
		
		return vals;
	}
	
	public static HashSet<Integer> getIncludeIndexes( List<String> genomeNames )
	{
		HashSet<Integer> set = new HashSet<Integer>();
		
		for( int x=0; x < genomeNames.size(); x++)
			if( genomeNames.get(x).toLowerCase().indexOf("pneu") != -1 )
				set.add(x);
		
		return set;
	}
	
	public static void main(String[] args) throws Exception
	{
		List<String> genomeNames = getGenomeNames();
		HashSet<Integer> indexes = getIncludeIndexes(genomeNames);
		
		BufferedWriter writer = new BufferedWriter(
				new FileWriter("/nobackup/afodor_research/af_broad/randomSpearman.txt"));
		
		writer.write("aFile\tbFile\tdistanceAll\tdistancePneuOnly\n");
		
		String[] list = KMER_DIST_DIRECTORY.list();
		
		for(int x=0; x < list.length-1; x++)
		{
			String xName = list[x];
			
			if( xName.endsWith("dist.txt"))
			{
				List<Float> aVals = getValsOrNull(new File(
						KMER_DIST_DIRECTORY.getAbsolutePath() + File.separator + xName),null);
				
				if( aVals != null)
				{
					List<Float> aValsPneuOnly =
							getValsOrNull(new File(
									KMER_DIST_DIRECTORY.getAbsolutePath() + File.separator + xName),indexes);
					
					for( int y=x+1; y < list.length; y++)
					{
						String yName = list[y];
						
						if( yName.endsWith("dist.txt"))
						{
							List<Float> bVals = getValsOrNull(new File(
									KMER_DIST_DIRECTORY.getAbsolutePath() + File.separator + yName),null);
									
							if(bVals != null)
							{
								List<Float> bValsPneuOnly = 
										getValsOrNull(new File(
												KMER_DIST_DIRECTORY.getAbsolutePath() + File.separator + yName),indexes);
								
								writer.write(xName + "\t");
								writer.write(yName + "\t");
								writer.write(Spearman.getSpear(aVals, bVals).getRs() + "\t");
								writer.write(Spearman.getSpear(aValsPneuOnly,bValsPneuOnly).getRs() + "\n");
								writer.flush();
							}
						}
					}
				}
			}
		}
		writer.flush(); writer.close();
	}	
}