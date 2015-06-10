package coPhylog;

import java.io.File;

public class WriteAFewContextsToConsole
{
	public static void main(String[] args) throws Exception
	{
		File f = new File("/projects/afodor_chs/context/chs_241_1_context.gz");
		
		CoPhylogBinaryFileReader.dumpFirstToConsole(f, 100);
	}
}
