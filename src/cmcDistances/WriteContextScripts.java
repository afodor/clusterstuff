package cmcDistances;

import java.io.File;

public class WriteContextScripts
{
	private static File INPUT_DIRECTORY = new File(
			"/projects/afodor_chs/From Nury 5-2016");
	
	private static File OUTPUT_DIRECTORY = new File(
			"/nobackup/afodor_research/fromNury52016/contexts");
	
	public static void main(String[] args) throws Exception
	{
		String[] files = INPUT_DIRECTORY.list();
		
		for(String s : files)
		{
			if ( s.endsWith(".fastq"))
			{
				
			}
		}
	}
}
