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

	double[] initialEnvironment;
	double[] microbiomeSum;
	double[][] parentalContribution;
	double[] environmentalContribution;
	double[] mixedContribution;
	double coefficient;
	
	private final int numberOfMicrobePerHost;
	private final int numberOfEnvironmentalSpecies;
	private final int numberOfIndividual;
	private final double environmentalFactor;
	private final double percentageOfPooledFixed;
	
	private int numberOfGeneration;
	private double g_diversity;
	private double a_diversity;
	private double b_diversity;
    private double weighted_b_diversity;
	private String[] alignment;
	
	private int sampleReplicates;
	private int numberOfSamples;
	private double[] previous_microbiomes;
	private double beta_diversity_coef;
    private double beta_diversity_coef_2;
    private double beta_diversity_coef_3;
    private double beta_diversity_coef_4;
	
    private Individual[] compositionOfIndividuals;
	private List<Set<Integer>> samples = new ArrayList<>();
	private List<Integer> host_index = new ArrayList<>();
	
	private Multinomial2 multiNomialDist;
	
	public Population(int numberOfMicrobePerHost, double[] environment, int noOfIndividual, 
			double environmentalFactor, double percentageOfPooledFixed, 
			int numberOfSamples, int sampleReplicates0) throws IOException {
		this.numberOfMicrobePerHost = numberOfMicrobePerHost;
		this.numberOfEnvironmentalSpecies = environment.length;
		this.numberOfIndividual = noOfIndividual;
		this.environmentalFactor = environmentalFactor;
		this.percentageOfPooledFixed = percentageOfPooledFixed;
		
		this.numberOfGeneration = 0;
		this.sampleReplicates = sampleReplicates0;				
		this.numberOfSamples = numberOfSamples;
		coefficient = (1 - environmentalFactor) / this.numberOfMicrobePerHost;
		initialEnvironment = environment;
		
		compositionOfIndividuals = new Individual[numberOfIndividual];
		parentalContribution = new double[noOfIndividual][numberOfEnvironmentalSpecies];
		microbiomeSum = new double[numberOfEnvironmentalSpecies];
		environmentalContribution = new double[numberOfEnvironmentalSpecies];
		mixedContribution = new double[numberOfEnvironmentalSpecies];
		previous_microbiomes = new double[numberOfEnvironmentalSpecies];
		
		multiNomialDist = new Multinomial2(numberOfEnvironmentalSpecies);
		
		alignment = new String[numberOfIndividual];
		
		beta_diversity_coef = 2.0 / numberOfIndividual
				/ (numberOfIndividual - 1) / numberOfEnvironmentalSpecies;
		beta_diversity_coef_2 = 1.0 / numberOfIndividual
				/ (numberOfIndividual - 1) / numberOfMicrobePerHost;
		beta_diversity_coef_3 = 1.0 / numberOfSamples / (numberOfSamples - 1)
				/ numberOfMicrobePerHost / sampleReplicates;
		beta_diversity_coef_4 = 2.0 / numberOfSamples / (numberOfSamples - 1)
				/ numberOfEnvironmentalSpecies / sampleReplicates;
		
		for (int i = 0; i < numberOfIndividual; i++) {
			host_index.add(i);
		}
		for (int i = 0; i < sampleReplicates; i++) {
			samples.add(RandomSample.randomSample(host_index,
					this.numberOfSamples));
		}
        multiNomialDist.updateProb(initialEnvironment);
		for (int i = 0; i < noOfIndividual; i++) {
			compositionOfIndividuals[i] = new Individual(
					multiNomialDist, this.numberOfMicrobePerHost,
					numberOfEnvironmentalSpecies);
		}
	}

	public void sumSpecies() {
		Arrays.fill(microbiomeSum, 0);
		for (Individual host : getIndividuals()) {
			VectorAddition.additionOfVectors(microbiomeSum, 1, 1,
					microbiomeSum, host.getMicrobiome());
		}
		for (int i = 0; i < numberOfEnvironmentalSpecies; i++) {
			microbiomeSum[i] = microbiomeSum[i] / numberOfIndividual
					/ numberOfMicrobePerHost;
		}
	}

	public void microbiomeSequenceAlignment() {
		
		int i = 0;
		for (Individual host : getIndividuals()) {
			alignment[i] = host.microbial_sequences();
			i++;
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

			return Distance.getDistance(initialEnvironment, microbiomeSum);
		else
			return Distance.getDistance(environmentalContribution,
					microbiomeSum);
	}



	public void parentalInheritanceAndEnvironmentalAcquisition() {
		
		for (int i = 0; i < numberOfIndividual; i++) {
			System.arraycopy(compositionOfIndividuals[MathUtil.getNextInt(numberOfIndividual-1)].getMicrobiome(),
                                0,parentalContribution[i], 0,numberOfEnvironmentalSpecies);
		}
		
		VectorAddition.additionOfVectors(environmentalContribution,
				percentageOfPooledFixed, 1 - percentageOfPooledFixed, microbiomeSum,
				initialEnvironment);

		for (int i = 0; i < numberOfIndividual; i++) {
			VectorAddition.additionOfVectors(mixedContribution, coefficient,
					environmentalFactor, parentalContribution[i],
					environmentalContribution);

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
		System.arraycopy(microbiomeSum, 0, previous_microbiomes, 0, microbiomeSum.length);
		numberOfGeneration++;
	}

	public double alphaDiversity(boolean sampleOrNot) {
		if (sampleOrNot) {
			a_diversity = 0;
			for (Individual host : getIndividuals()) {
				for (double abundance : host.getMicrobiome()) {
					if (abundance > 0 && abundance < numberOfMicrobePerHost) {
						a_diversity -= getPlogP(abundance, numberOfMicrobePerHost);
//						double relative_abundance = abundance
//								/ numberOfMicrobePerHost;
//								relative_abundance
//								* Math.log(relative_abundance);
					}
				}
			}
			a_diversity /=  numberOfIndividual;
		} else {
 			a_diversity = 0;
			for (Set<Integer> sample : samples) {
				for (Integer index : sample) {
					for (double abundance : getIndividuals()[index].getMicrobiome()) {
						if (abundance > 0
								&& abundance < numberOfMicrobePerHost) {
							a_diversity -= getPlogP(abundance, numberOfMicrobePerHost);
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
		} 
                else {
			double temp_b = 0;
			for (Set<Integer> sample : samples) {
				String[] temp_alignment = new String[numberOfSamples];
				int m = 0;
				for (Integer index : sample) {
					temp_alignment[m] = alignment[index];
					m++;
				}
				for (int i = 1; i < numberOfSamples; i++) {
					for (int j = 0; j < i; j++) {
						temp_b += DifferentElement.differentElement(
								temp_alignment[i], temp_alignment[j]);
					}
				}
			}
			b_diversity = beta_diversity_coef_4 *  temp_b;
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
					VectorAddition.additionOfVectors(temp_sum, 1, 1,
							temp_sum,
							getIndividuals()[index].getMicrobiome());
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
        
        public double BrayCurtis(boolean sampleOrNot){
            if (sampleOrNot){
                weighted_b_diversity=0;
                for (int i=1;i<numberOfIndividual;i++){
                    for (int j=0;j<i;j++){
                        for (int k=0;k<numberOfEnvironmentalSpecies;k++){
                            weighted_b_diversity+=
                                    Math.abs(getIndividuals()[i].getMicrobiome()[k]-getIndividuals()[j].getMicrobiome()[k]);
                        }
                    }
                }
                weighted_b_diversity*=beta_diversity_coef_2;
                return weighted_b_diversity;
            }
            else{
                weighted_b_diversity=0;
                for (Set<Integer> sample : samples) {
                    double[][] temp_microbiomes=new double[numberOfSamples][];
                    int m=0;
                    for (Integer index:sample){
                        temp_microbiomes[m]=getIndividuals()[index].getMicrobiome();
                        m++;
                    }
                    for (int i=1;i<numberOfSamples;i++){
                        for(int j=0;j<i;j++){
                            for (int k=0;k<numberOfEnvironmentalSpecies;k++){
                                weighted_b_diversity+=Math.abs(temp_microbiomes[i][k]-temp_microbiomes[j][k]);
                            }
                        }
                    }
                }
                weighted_b_diversity*=beta_diversity_coef_3;
                return weighted_b_diversity;
            }
        }

	public String printOut() {
		StringBuilder sb = new StringBuilder();
		for (int i = 1; i < numberOfEnvironmentalSpecies; i++) {
			sb.append(microbiomeSum[i]).append("\t");
		}
		return sb.toString().trim();
	}
    
    public Individual[] getIndividuals(){
        return compositionOfIndividuals;
    }

	public int getNumberOfGeneration() {
		return numberOfGeneration;
	}

	public void resetGeneration() {
		numberOfGeneration = 0;
		
	}

	private static double getPlogP(double abundance) {
		return getPlogP(abundance, 1);
	}

	private static double getPlogP(double abundance, int normaliser) {
		double relative_abundance = abundance / normaliser;
		double pLogP = relative_abundance * Math.log(relative_abundance);
		return pLogP;
	}
        
     


}

