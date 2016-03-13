package creOrthologs.kmers;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.List;
import java.util.StringTokenizer;

import utils.Spearman;

public class WriteSpearmanFromRandom
{
	private static final File KMER_DIST_DIRECTORY = 
			new File("/nobackup/afodor_research/af_broad/randomKMerMatrices");
	
	public static final int EXPECTED_NUM_LINES = 340;
	
	private static List<Float> getValsOrNull(File f) throws Exception
	{
		List<Float> vals = new ArrayList<Float>();
		
		BufferedReader reader = new BufferedReader(new FileReader(f));
		
		reader.readLine();
		
		int numRead =0;
		for(String s=  reader.readLine(); s != null; s= reader.readLine())
		{
			StringTokenizer sToken = new StringTokenizer(s);

			sToken.nextToken();
			
			while(sToken.hasMoreTokens())
				vals.add(Float.parseFloat(sToken.nextToken()));
			
			numRead++;
		}
		
		reader.close();
		
		if( numRead != EXPECTED_NUM_LINES)
		{
			System.out.println(numRead + " skipping " + f.getAbsolutePath());
			return null;
		}
		
		return vals;
	}
	
	public static void main(String[] args) throws Exception
	{
		BufferedWriter writer = new BufferedWriter(
				new FileWriter("/nobackup/afodor_research/af_broad/randomSpearman.txt"));
		
		writer.write("aFile\tbFile\tdistance\n");
		
		String[] list = KMER_DIST_DIRECTORY.list();
		
		for(int x=0; x < list.length-1; x++)
		{
			String xName = list[x];
			
			if( xName.endsWith("dist"))
			{
				List<Float> aVals = getValsOrNull(new File(
						KMER_DIST_DIRECTORY.getAbsolutePath() + File.separator + xName));
				
				if( aVals != null)
				{
					for( int y=x+1; y < list.length; y++)
					{
						String yName = list[y];
						
						if( yName.endsWith("dist"))
						{
							List<Float> bVals = getValsOrNull(new File(
									KMER_DIST_DIRECTORY.getAbsolutePath() + File.separator + yName));
							
							if(bVals != null)
							{
								writer.write(xName + "\t");
								writer.write(yName + "\t");
								writer.write(Spearman.getSpear(aVals, bVals) + "\n");
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