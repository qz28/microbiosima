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

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Set;

import utils.DifferentElement;
import utils.Distance;
import utils.RandomSample;
import utils.VectorAddition;
import utils.random.MathUtil;
import utils.random.Multinomial2;

/**
 *
 * @author John
 */
public class Population {
	// FIXME: encapsulation these fields, way too many public
	private final int numberOfMicrobePerHost;
	private final int numberOfEnvironmentalSpecies;
	private final int numberOfIndividual;
	private final double environmentalFactor;
	private final double percentageOfpooledOrFixed;

	private int numberOfGeneration;
	private double[] initialEnvironment;
	private double[] microbiomeSum;
	private double g_diversity;
	private double a_diversity;
	private double b_diversity;
	private String[] alignment;

	private double[][] parentalContribution;
	private double[] environmentalContribution;
	private double[] mixedContribution;
	private double coefficient;
	private int sampleReplicates;
	private int numberOfSamples;
	private double[] previous_microbiomes;
	private Individual[] compositionOfIndividuals;
	private List<Set<Integer>> samples = new ArrayList<>();
	private List<Integer> host_index = new ArrayList<>();
	private double beta_diversity_coef;

	private Multinomial2 multiNomialDist;

	public Population(int numberOfMicrobePerHost, double[] environment,
			int noOfIndividual, double environmentalFactor,
			double pooledOrFiexd, String type, int numberOfSamples,
			int sampleReplicates0) throws IOException {
		this.numberOfMicrobePerHost = numberOfMicrobePerHost;
		this.numberOfEnvironmentalSpecies = environment.length;
		this.numberOfIndividual = noOfIndividual;
		this.environmentalFactor = environmentalFactor;
		this.percentageOfpooledOrFixed = pooledOrFiexd;// NOTE: "or" sounds like
														// a boolean variable

		this.numberOfGeneration = 0;
		this.sampleReplicates = sampleReplicates0;
		this.numberOfSamples = numberOfSamples;
		coefficient = (1 - environmentalFactor) / this.numberOfMicrobePerHost; // NOTE:
																				// change
																				// to
																				// this.
		initialEnvironment = environment;

		compositionOfIndividuals = new Individual[numberOfIndividual];
		parentalContribution = new double[noOfIndividual][numberOfEnvironmentalSpecies];// NOTE:
																						// change
																						// to
																						// numberOfIndividual
		microbiomeSum = new double[numberOfEnvironmentalSpecies];
		environmentalContribution = new double[numberOfEnvironmentalSpecies];
		mixedContribution = new double[numberOfEnvironmentalSpecies];
		previous_microbiomes = new double[numberOfEnvironmentalSpecies];

		multiNomialDist = new Multinomial2(numberOfEnvironmentalSpecies);

		beta_diversity_coef = 2 / numberOfIndividual / (numberOfIndividual - 1)
				/ numberOfEnvironmentalSpecies;
		for (int i = 0; i < numberOfIndividual; i++) {
			host_index.add(i);
		}
		for (int i = 0; i < sampleReplicates; i++) {
			samples.add(RandomSample.randomSample(host_index,
					this.numberOfSamples));
		}
		CustomFileReader fr = new CustomFileReader(type
				+ "_simulated_microbiomes.txt", numberOfIndividual);
		for (int i = 0; i < noOfIndividual; i++) {
			compositionOfIndividuals[i] = new Individual(fr.getCommandSpilt(i),
					this.numberOfMicrobePerHost, numberOfEnvironmentalSpecies);
		}
	}

	public void sumSpecies() {
		Arrays.fill(microbiomeSum, 0);
		for (Individual host : compositionOfIndividuals) {
			VectorAddition.additionOfVectorsSum(microbiomeSum, microbiomeSum,
					host.getMicrobiome());
		}
		for (int i = 0; i < numberOfEnvironmentalSpecies; i++) {
			microbiomeSum[i] = microbiomeSum[i] / numberOfIndividual
					/ numberOfMicrobePerHost;
		}
	}

	public void microbiomeSequenceAlignment() {
		alignment = new String[numberOfIndividual];// NOTE: pre declared this
													// somewhere
		int i = 0;
		for (Individual host : compositionOfIndividuals) {
			alignment[i] = host.microbial_sequences();
			i++;// NOTE: Since you use i++ here, might be better with standard
				// for loop
			// for (int j = 0; j < compositionOfIndividuals.length; j++) {
		}
	}

	public double interGenerationDistance() {
		if (getNumberOfGeneration() == 0)
			return 0;
		else {
			return Distance.getDistance(previous_microbiomes, microbiomeSum);
		}
	}

	public double environmentPopulationDistance() {
		if (getNumberOfGeneration() == 0)
			// TODO: can we merge these two??
			return Distance.getDistance(initialEnvironment, microbiomeSum);
		else
			return Distance.getDistance(environmentalContribution,
					microbiomeSum);
	}

	public void parentalInheritanceAndEnvironmentalAcquisition() {

		for (int i = 0; i < numberOfIndividual; i++) {
			System.arraycopy(compositionOfIndividuals[MathUtil
					.getNextInt(numberOfIndividual - 1)].getMicrobiome(), 0,
					parentalContribution[i], 0, numberOfEnvironmentalSpecies);
		}

		environmentalContribution = VectorAddition.additionOfVectorsMap(
				environmentalContribution, percentageOfpooledOrFixed,
				1 - percentageOfpooledOrFixed, microbiomeSum,
				initialEnvironment);

		for (int i = 0; i < numberOfIndividual; i++) {
			mixedContribution = VectorAddition.additionOfVectorsMap(
					mixedContribution, coefficient, environmentalFactor,
					parentalContribution[i], environmentalContribution);

			multiNomialDist.updateProb(mixedContribution);
			multiNomialDist.multisample(
					compositionOfIndividuals[i].getMicrobiome(),
					numberOfMicrobePerHost);

		}
	}

	public void getNextGen() {
		parentalInheritanceAndEnvironmentalAcquisition();
		for (int i = 0; i < sampleReplicates; i++) {
			samples.set(i,
					RandomSample.randomSample(host_index, numberOfSamples));
		}

		System.arraycopy(microbiomeSum, 0, previous_microbiomes, 0,
				microbiomeSum.length);

		// example of IDE overdo it
		// setNumberOfGeneration(getNumberOfGeneration() + 1);
		numberOfGeneration++;
	}

	private static double getPlogP(double abundance, int normaliser) {
		double relative_abundance = abundance / normaliser;
		double pLogP = relative_abundance * Math.log(relative_abundance);
		return pLogP;
	}

	public double alphaDiversity(boolean sampleOrNot) {
		if (sampleOrNot) {
			a_diversity = 0;
			for (Individual host : compositionOfIndividuals) {
				for (double abundance : host.getMicrobiome()) {
					if (abundance > 0 && abundance < numberOfMicrobePerHost) {
						a_diversity -= getPlogP(abundance,
								numberOfMicrobePerHost);
						// double relative_abundance = abundance
						// / numberOfMicrobePerHost;
						// relative_abundance
						// * Math.log(relative_abundance);
					}
				}
			}
			a_diversity /= numberOfIndividual;
		} else {
			// double temp_a = 0;
			for (Set<Integer> sample : samples) {
				for (Integer index : sample) {
					for (double abundance : compositionOfIndividuals[index]
							.getMicrobiome()) {
						if (abundance > 0 && abundance < numberOfMicrobePerHost) {
							a_diversity -= getPlogP(abundance,
									numberOfMicrobePerHost);
						}
					}
				}
			}
			a_diversity /= (numberOfSamples * sampleReplicates);
		}
		return a_diversity;

	}

	public double betaDiversity(boolean sampleOrNot) {
		if (sampleOrNot) {
			b_diversity = 0;
			for (int i = 1; i < numberOfIndividual; i++) {
				for (int j = 0; j < i; j++) {
					b_diversity += DifferentElement.differentElement(
							alignment[i], alignment[j]);
				}
			}

			b_diversity *= beta_diversity_coef;
			return b_diversity;
		} else {
			double temp_b = 0;
			for (Set<Integer> sample : samples) {
				String[] temp_alignment = new String[numberOfSamples];
				int m = 0;
				for (Integer index : sample) {
					temp_alignment[m] = alignment[index];
					m++;// NOTE: for loop or arrayList, track one less variable
				}
				for (int i = 1; i < numberOfSamples; i++) {
					for (int j = 0; j < i; j++) {
						temp_b += DifferentElement.differentElement(
								temp_alignment[i], temp_alignment[j]);
					}
				}
			}
			b_diversity = beta_diversity_coef * temp_b / sampleReplicates;
			return b_diversity;
		}
	}

	public double gammaDiversity(boolean sampleOrNot) {
		if (sampleOrNot) {
			g_diversity = 0;
			for (double abundance : microbiomeSum) {
				if (abundance > 0 && abundance < 1) {
					g_diversity -= getPlogP(abundance);
				}
			}
			return g_diversity;
		} else {
			double temp_g = 0;
			double[] temp_sum = new double[numberOfEnvironmentalSpecies];
			for (Set<Integer> sample : samples) {
				Arrays.fill(temp_sum, 0);
				for (Integer index : sample) {
					VectorAddition.additionOfVectorsSum(temp_sum, temp_sum,
							compositionOfIndividuals[index].getMicrobiome());
				}
				for (double abundance : temp_sum) {
					double r_abundance = abundance / numberOfMicrobePerHost
							/ numberOfSamples;
					if (r_abundance > 0 && r_abundance < 1) {
						temp_g -= getPlogP(r_abundance);
					}
				}
			}
			g_diversity = temp_g / sampleReplicates;
			return g_diversity;
		}
	}

	private static double getPlogP(double abundance) {
		return getPlogP(abundance, 1);
	}

	public String printOut() {
		// String content = Double.toString(microbiomeSum[0]);
		// for (int i = 1; i < numberOfEnvironmentalSpecies; i++) {
		// content = content + "\t" + Double.toString(microbiomeSum[i]);
		// }
		// Better with StringBuilder
		StringBuilder sb = new StringBuilder();
		for (int i = 1; i < numberOfEnvironmentalSpecies; i++) {
			sb.append(microbiomeSum[i]).append("\t");
		}
		return sb.toString().trim();
	}

	public int getNumberOfGeneration() {
		return numberOfGeneration;
	}

	public void resetGeneration() {
		numberOfGeneration = 0;

	}

}
