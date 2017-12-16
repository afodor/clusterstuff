package chinaDec2017;

import java.io.File;
import java.util.HashSet;
import java.util.LinkedHashSet;

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
		HashSet<File> dirs = getTopDirectories();
		
		for(File f : dirs)
		{
			File forwardFile = new File(f.getAbsolutePath() + File.separator + 
					"DMP0" + f.getName().substring(1) + "_L1_" + f.getName() + "_1.fq.gz");
			
			if( ! forwardFile.exists())
				throw new Exception("Could not find " + forwardFile);
			
			File backwardsFile = new File( f.getAbsolutePath() + File.separator + 
					"DMP0" + f.getName().substring(1) + "_L1_" + f.getName() + "_2.fq.gz" );
			
			if( ! backwardsFile.exists())
				throw new Exception("Could not find "+ backwardsFile);
		}
		
		System.out.println("Finished");
	}
	
	public static HashSet<File> getTopDirectories() throws Exception
	{
		HashSet<File> dirNames = new LinkedHashSet<File>();
		
		for(String dirPath : DIRS_TO_SCAN)
		{
			File topDir = new File(dirPath);
			
			String[] files = topDir.list();
			
			for( String f : files)
			{
				File aFile = new File(topDir.getAbsolutePath() + File.separator + f);
				
				if( aFile.exists() && aFile.isDirectory() )
				{
					if( dirNames.contains(aFile))
					{
						throw new Exception("Duplicate  " + f);
					}

					dirNames.add(aFile);
				}
				
			}
		}
		
		System.out.println("Finished with " +  dirNames.size() );
		return dirNames;
	}
}
