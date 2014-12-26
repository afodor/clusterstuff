package mbqc;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Random;

public class CreateRDPQSub
{
	public static final File BIOINFORMATICS_DIR = 
			new File("/projects/afodor_research/mbqc/bioinformatics_distribution");
	
	public static final File FASTA_DIR = new File("/projects/afodor_research/mbqc/fastaSeqs");
	
	public static final File RDP_RUN_DIR = new File("/projects/afodor_research/mbqc/rdpRunDir");
	
	public static final File RDP_OUT_DIR = new File("/projects/afodor_research/mbqc/rdpOutDir");
	
	public static final int NUMBER_JOBS_PER_QSUB = 40;
	
	public static void main(String[] args) throws Exception
	{
		List<File> files= new ArrayList<File>();
		
		// currently adds 4,404 files
		recurisvelyAddAllFiles(BIOINFORMATICS_DIR, files);
		Collections.shuffle(files, new Random(322221));
		
		int numJobsDone =0;
		int fileIndex= 1;
		
		BufferedWriter allShWriter= new BufferedWriter(new FileWriter(new File(RDP_RUN_DIR + File.separator + 
				"runAll.sh")));
		
		File shFile = new File( RDP_RUN_DIR + File.separator + "run_" + fileIndex + ".sh");
		
		BufferedWriter shWriter = new BufferedWriter(new FileWriter(shFile));
		
		allShWriter.write("qsub -q \"viper\" " + shFile.getAbsolutePath() +  "\n"  );
		
		HashSet<String> names = new HashSet<String>();
		for( File f : files)
		{
			String name = f.getParentFile().getParentFile().getName() + "_" + 
							f.getParentFile().getName() + f.getName();
			names.add( name);
			
			File fastAFile = new File(FASTA_DIR.getAbsolutePath() + File.separator + 
										name + "_toFasta.txt"	);
			
			/*
			shWriter.write("java -cp /users/afodor/gitInstall/clusterstuff/bin " + 
					"parsers.FastQToFastA" + " " + f.getAbsolutePath() +
					" " + fastAFile + "\n"
								);
								*/
			
			File rdpFile = new File(RDP_OUT_DIR.getAbsolutePath() + File.separator + 
					name + "_rdpOut.txt"	);
			
			
			/*
			shWriter.write("java -jar /users/afodor/rdp/rdp_classifier_2.10.1/dist/classifier.jar " + 
					"-o "+ rdpFile.getAbsolutePath() + " -q " + fastAFile+ "\n" );
			
			shWriter.write("gzip " + fastAFile.getAbsolutePath() + "\n");
			shWriter.write("gzip " + rdpFile.getAbsolutePath() + "\n");
			*/
			
			shWriter.write("java -cp /users/afodor/gitInstall/clusterstuff/bin " + 
					"parsers.ReduceToThreeColumn" + " " + rdpFile.getAbsolutePath() + "\n"
								);
			
			shWriter.flush();
			
			numJobsDone++;
			
			if( numJobsDone >= NUMBER_JOBS_PER_QSUB)
			{
				shWriter.flush();  shWriter.close();
				numJobsDone = 0;
				fileIndex++;
				shFile = new File( RDP_RUN_DIR + File.separator + "run_" + fileIndex + ".sh");
				shWriter = new BufferedWriter(new FileWriter(shFile));
				allShWriter.write("qsub -q \"viper\" " + shFile.getAbsolutePath() +  "\n"  );
				allShWriter.flush();
			}
		}
		
		if( files.size() != names.size() )
			throw new Exception("Logic error");
		
		shWriter.flush();  shWriter.close();
		allShWriter.flush();  allShWriter.close();
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
