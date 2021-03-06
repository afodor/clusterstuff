package vanderbilt;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;

public class MakeKrakenScripts
{
	private static HashSet<String> names= new HashSet<String>();
	
	public static void main(String[] args) throws Exception
	{
		List<File> list = new ArrayList<File>();
		addForADirectory(list, "/projects/afodor_research/vanderbilt/rc-ns-ftp.its.unc.edu/MSHRM1/MSHRM1");
		addForADirectory(list, "/projects/afodor_research/vanderbilt/rc-ns-ftp.its.unc.edu/MSHRM2/MSHRM2");
		
		BufferedWriter writer = new BufferedWriter(new FileWriter(new File("/projects/afodor_research/vanderbilt/runKraken/runAll.sh")));
		
		for(File f : list)
		{
			writer.write("qsub -q \"viper\"  " + " -N " + f.getName().replaceAll("\"", "") + " "  +  f.getAbsolutePath()  + "\n");
		}
		
		writer.flush();  writer.close();
	}
	
	private static void addForADirectory( List<File> shFiles, String directorToAdd ) throws Exception
	{
		File topDir = new File(directorToAdd);
		
		for(String s : topDir.list())
		{
			File nextDir = new File(topDir.getAbsolutePath() + File.separator + s);
			
			String[] seqFiles= nextDir.list();
			
			if( seqFiles.length != 2)
				throw new Exception("Expecting exactly two files in " + nextDir.getAbsolutePath());
			
			String name = nextDir.getName() + "_" + topDir.getName() ;
			
			if( names.contains(name))
				throw new Exception("Duplicate name");
			
			// report will go to standard output
			File shFile = new File("/projects/afodor_research/vanderbilt/runKraken/run" + name + ".sh" );
			
			shFiles.add(shFile);
			
			BufferedWriter writer = new BufferedWriter(new FileWriter(shFile));
			
			// request 128 GB box 
			//writer.write("#PBS -l nodes=1:ppn=16\n");
			//writer.write("#PBS -W x=NODESET:ONEOF:FEATURE:ib_qdr2\n");
			
			writer.write("#/projects/afodor_research/krakenInstall/kraken --threads 15 " + 
			"--db /projects/afodor_research/krakenInstall/krakenStandardDB2 " + 
			"--output /projects/afodor_research/vanderbilt/krakenOut/" + nextDir.getName() + "_" + topDir.getName() + "_krakenData.txt " + 
			"--fastq-input --gzip-compressed " + 
			"--paired " + nextDir.getAbsolutePath() + "/" + seqFiles[0] + " " + nextDir.getAbsolutePath() + "/" + seqFiles[1] + "\n");
			
			/*
			// report goes to standard output
			writer.write("#/projects/afodor_research/krakenInstall/kraken-mpa-report " + 
					"--db /projects/afodor_research/krakenInstall/krakenStandardDB2 " + 
					 " /projects/afodor_research/vanderbilt/krakenOut/" + nextDir.getName() + "_" + topDir.getName() + "_krakenData.txt \n");
					 */
			
			writer.write("/projects/afodor_research/krakenInstall/kraken-report " + 
					"--db /projects/afodor_research/krakenInstall/krakenStandardDB2 " + 
					 " /projects/afodor_research/vanderbilt/krakenOut/" + nextDir.getName() + "_" + topDir.getName() + "_krakenData.txt " +
				      " > /projects/afodor_research/vanderbilt/krakenOut/standardReport_for_" + nextDir.getName() + "_" + topDir.getName()   + ".txt\n");
			
			writer.write("#/projects/afodor_research/krakenInstall/kraken --threads 15 " + 
					"--db /projects/afodor_research/krakenInstall/krakenHumanDB" + 
					" --output /projects/afodor_research/vanderbilt/krakenOut/" + nextDir.getName() + "_" + topDir.getName() + "_krakenHumanData.txt " + 
					" --fastq-input --gzip-compressed " + 
					"--paired " + nextDir.getAbsolutePath() + "/" + seqFiles[0] + " " + nextDir.getAbsolutePath() + "/" + seqFiles[1] + "\n");
			
			writer.flush();  writer.close();
		}
	}
}
