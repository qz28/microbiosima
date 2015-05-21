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

import java.util.Arrays;

public class Multinomial2 {

	private double[] distribution;
	private int rangeOne;

	public Multinomial2(int range) {
		this.rangeOne = range - 1;
		distribution = new double[range];

	}

	public void updateProb(double[] probabilities) {
		double cum_prob = 1;
		for (int i = 0; i < distribution.length; i++) {
			if (probabilities[i] == 0)
				distribution[i] = 0;
			else {
				distribution[i] = probabilities[i] / cum_prob;
				cum_prob = cum_prob - probabilities[i];
			}
		}
	}

	public double[] multisample(double[] sampledNumbers, int sampleSize) {

		for (int i = 0; i < rangeOne; i++) {
			int rand = Binomial.binomial_sampling(sampleSize, distribution[i]);
			sampledNumbers[i] = rand;
			sampleSize -= rand;

			if (sampleSize == 0) {
				Arrays.fill(sampledNumbers, i + 1, sampledNumbers.length, 0);
				break;
			}
		}
		sampledNumbers[rangeOne] = sampleSize;
		return sampledNumbers;
	}
}
