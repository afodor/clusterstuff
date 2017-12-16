package chinaDec2017;

import java.io.File;
import java.util.HashSet;

public class CollectFasta
{
	private static final String[] DIRS_TO_SCAN = 
		{
			"/nobackup/afodor_research/datasets/chinaDec2017/MBiome/raw_1",
			"/nobackup/afodor_research/datasets/chinaDec2017/MBiome/raw_2",
			"/nobackup/afodor_research/datasets/chinaDec2017/MBiome/raw_3"
		};
	
	
	public static void main(String[] args) throws Exception
	{
		int numDuplicates =0;
		HashSet<String> dirNames = new HashSet<String>();
		
		for(String dirPath : DIRS_TO_SCAN)
		{
			File topDir = new File(dirPath);
			
			String[] files = topDir.list();
			
			for( String f : files)
			{
				File aFile = new File(topDir.getAbsolutePath() + File.separator + f);
				
				if( aFile.exists() && aFile.isDirectory() && f.toLowerCase().endsWith(".xls"))
				{
					if( dirNames.contains(f))
					{
						System.out.println("Duplicate " + f + " " + aFile.getAbsolutePath());
						numDuplicates++;
					}
				}
				
				dirNames.add(f);
			}
		}
		
		for(String s : dirNames)
		{
			System.out.println(s);
		}
		
		System.out.println("Finished with " +  dirNames.size() + " " + numDuplicates);
	}
}
