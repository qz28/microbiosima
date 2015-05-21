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
/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package microbiosima;

import java.util.Random;

import utils.random.Binomial;

public class Multinomial {
//	static Random generator = new Random();
	double[] distribution;
	public int range;

	// Constructor
	Multinomial(double[] probabilities) {
		range = probabilities.length;
		distribution = new double[range];
		double cum_prob = 1;
		for (int i = 0; i < range; i++) {
			if (probabilities[i] == 0)
				distribution[i] = 0;
			else {
				distribution[i] = probabilities[i] / cum_prob;
				cum_prob = cum_prob - probabilities[i];
			}
		}
	}

	public void multisample(double[] sampled_numbers, int sample_size) {
		int N = sample_size;
		for (int i = 0; i < range - 1; i++) {
			sampled_numbers[i] = Binomial.binomial_sampling(N, distribution[i]);
			N = N - (int) sampled_numbers[i];
		}
		sampled_numbers[range - 1] = N;
	}
}
