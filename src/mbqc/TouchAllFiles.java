package mbqc;

import java.io.File;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;

public class TouchAllFiles
{
	public static final File BIOINFORMATICS_DIR = 
			new File("/projects/afodor_research/mbqc/bioinformatics_distribution");
	
	
	public static void main(String[] args) throws Exception
	{
		List<File> files= new ArrayList<File>();
		
		recurisvelyAddAllFiles(BIOINFORMATICS_DIR, files);
		
		HashSet<String> names = new HashSet<String>();
		for( File f : files)
		{
			String name = f.getParentFile().getParentFile().getName() + "_" + 
							f.getParentFile().getName() + f.getName();
			names.add( name);
			System.out.println(f.getAbsolutePath());
		}
		
		System.out.println(files.size() + " " + names.size());
	}
	
	private static void recurisvelyAddAllFiles(File startDir, List<File> files) 
				throws Exception
	{
		for(String s : startDir.list())
		{
			File f = new File(startDir.getAbsolutePath() + File.separator + s);
			
			if( f.isDirectory())
			{
				recurisvelyAddAllFiles(f, files);
			}
			else if( ! s.toLowerCase().startsWith("umatched") && s.toLowerCase().endsWith("fastq.gz") )
			{
				files.add(f);
			}
		}
		
	}
}
