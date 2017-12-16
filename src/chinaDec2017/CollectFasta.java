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
	
	private static File findOneFile(File dir, String name, int read)
		throws Exception
	{
		String[] list = dir.list();
		
		File returnFile = null;
		
		for(String s : list)
		{
			if( s.endsWith(name + "_" + read +  ".fq.gz"))
			{
				if( returnFile != null)
					throw new Exception("Duplicate " + dir.getAbsolutePath() + " " +  name + " " + read);
				
				returnFile = new File(dir.getAbsolutePath() + File.separator +s );
			}
		}
		
		if( returnFile == null ) 
			throw new Exception("Could not find " +  dir.getAbsolutePath() + " " +  name + " " + read);
		
		return returnFile;
	}
	
	public static void main(String[] args) throws Exception
	{
		HashSet<File> dirs = getTopDirectories();
		
		for(File f : dirs)
		{
			File forwardFile = findOneFile(f, f.getName(), 1);
			File backwardsFile = findOneFile(f, f.getName(), 2);
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
