package vanderbilt;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;

public class MakeKrakenScriptsFor16S
{
	private static HashSet<String> names= new HashSet<String>();
	
	public static void main(String[] args) throws Exception
	{
		List<File> list = new ArrayList<File>();
		addForADirectory(list, SplitIntoSeparateFiles.outDir.getAbsolutePath());
		
		BufferedWriter writer = new BufferedWriter(new FileWriter(new File("/projects/afodor_research/vanderbilt/runKraken/runAll16S.sh")));
		
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
			File seqFile= new File(topDir.getAbsolutePath() + File.separator + s);
			String name = seqFile.getName().replaceAll(".fasta", "");
			
			if( names.contains(name))
				throw new Exception("No");
			
			names.add(name);
			
			File shFile = new File("/projects/afodor_research/vanderbilt/runKraken/run" + name + ".sh" );
			
			shFiles.add(shFile);
			
			BufferedWriter writer = new BufferedWriter(new FileWriter(shFile));
			
			// request 128 GB box 
			writer.write("#PBS -l nodes=1:ppn=16\n");
			writer.write("#PBS -W x=NODESET:ONEOF:FEATURE:ib_qdr2\n");
			
			File krakenOut = new File("/projects/afodor_research/vanderbilt/krakenOut/" + name +"_16S_to_krakenData.txt ");
			
			writer.write("/projects/afodor_research/krakenInstall/kraken --threads 15 " + 
			"--db /projects/afodor_research/krakenInstall/krakenStandardDB2 " + 
			"--output " + krakenOut.getAbsolutePath() +  " " + seqFile.getAbsolutePath() +  "\n");
			
			
			writer.write("/projects/afodor_research/krakenInstall/kraken-report " + 
					"--db /projects/afodor_research/krakenInstall/krakenStandardDB2 " + 
					   krakenOut.getAbsolutePath()  +
				      " > /projects/afodor_research/vanderbilt/krakenOut/standardReport_for_" + name+ "_to16S.txt\n");
			
			
			writer.flush();  writer.close();
		}
	}
}
