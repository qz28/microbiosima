

package microbiosima;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;

/*
 * @author John
 */
public class CustomFileReader {
	private String[] commands;

	public CustomFileReader(String file_path, int number_of_line) throws IOException {
		BufferedReader bf = new BufferedReader(new FileReader(file_path));
		String aLine;
		commands = new String[number_of_line];
		int i = 0;
		//FIXME: What happen if number_of_line < actual numbers of lines in the file?
		while ((aLine = bf.readLine()) != null) {
			commands[i] = aLine;
			i++;
		}
		bf.close();
	}

	/**
	 * @return the commands
	 */
	public String getCommand(int index) {
		return commands[index];
	}
	
	public String[] getCommandSpilt(int index){
		return commands[index].split("\t");
	}
	
	

}
