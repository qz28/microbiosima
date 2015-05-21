/*******************************************************************************
 *
 * Copyright (C) 2014, 2015 Qinglong Zeng, Jeet Sukumaran, Steven Wu and Allen Rodrigo
 *
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *******************************************************************************/

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
