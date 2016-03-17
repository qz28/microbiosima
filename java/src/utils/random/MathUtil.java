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
package utils.random;

import org.apache.commons.math3.random.MersenneTwister;
import org.apache.commons.math3.random.RandomDataGenerator;

public class MathUtil {

	private static RandomDataGenerator rd = new RandomDataGenerator(
			new MersenneTwister());
	
	

	public static double nextFloat(){
		return rd.nextUniform(0, 1);
	}
	public static int getNextInt(int max){
		return rd.nextInt(0, max);
	}
        public static double getNextFloat(double max){
		return rd.nextUniform(0, 1)*max;
	}
        public static void setSeed(){
            rd.reSeed();
        }
        public static void setSeed(long s){
            rd.reSeed(s);
        }
	
	
}
