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

import utils.Distance;
import utils.DiversityIndex;
import utils.RandomSample;
import utils.VectorAddition;
import utils.random.MathUtil;
import utils.random.Multinomial2;

/**
 *
 * @author John
 */
public class Population {
	
	final int numberOfMicrobePerHost;
	final int numberOfEnvironmentalSpecies;
	final int numberOfIndividual;
	final double environmentalFactor;
	final double percentageOfpooledOrFixed;
	
	double coefficient;
    double[][] parentalContribution;
	double[] environmentalContribution;
	double[] mixedContribution;
	double[] initialEnvironment;
	double[] microbiomeSum;
	private double[] previous_microbiomes;
	int[] parentalIndex;
    int[] ancestryIndex;
	
	
	private double g_diversity;
	private double a_diversity;
	private double b_diversity;
	private double beta_diversity_coef;
    private double beta_diversity_coef_2;
    private double beta_diversity_coef_3;
    private double beta_diversity_coef_4;	
    private double weighted_b_diversity;
	
	private int numberOfGeneration;
	private int sampleReplicates;
	private int numberOfSamples;	

	private Individual[] compositionOfIndividuals;

	private List<Set<Integer>> samples = new ArrayList<>();
	
	List<Integer> host_index = new ArrayList<>();
	
	Multinomial2 multiNomialDist;
	
	public Population(int numberOfMicrobePerHost, double[] environment, int noOfIndividual, 
			double environmentalFactor,
			double pooledOrFiexd, int numberOfSamples, int sampleReplicates0) throws IOException {
		this.numberOfMicrobePerHost = numberOfMicrobePerHost;
		this.numberOfEnvironmentalSpecies = environment.length;
		this.numberOfIndividual = noOfIndividual;
		this.environmentalFactor = environmentalFactor;
		this.percentageOfpooledOrFixed = pooledOrFiexd;//NOTE: "or" sounds like a boolean variable
		
		this.numberOfGeneration = 0;
		this.sampleReplicates = sampleReplicates0;				
		this.numberOfSamples = numberOfSamples;
		coefficient = (1 - environmentalFactor) / this.numberOfMicrobePerHost; //NOTE: change to this.
		initialEnvironment = environment;
		
		compositionOfIndividuals = new Individual[numberOfIndividual];
		parentalContribution = new double[noOfIndividual][numberOfEnvironmentalSpecies];//NOTE: change to numberOfIndividual
		microbiomeSum = new double[numberOfEnvironmentalSpecies];
		environmentalContribution = new double[numberOfEnvironmentalSpecies];
		mixedContribution = new double[numberOfEnvironmentalSpecies];
		previous_microbiomes = new double[numberOfEnvironmentalSpecies];
                parentalIndex=new int[numberOfIndividual];
                ancestryIndex=new int[numberOfIndividual];
		
		multiNomialDist = new Multinomial2(numberOfEnvironmentalSpecies);
		
		beta_diversity_coef = 2.0  / numberOfIndividual
				/ (numberOfIndividual - 1)/numberOfMicrobePerHost;
        beta_diversity_coef_2 = 2.0  / numberOfIndividual
				/ (numberOfIndividual - 1);
        beta_diversity_coef_3 = 2.0  / numberOfSamples
				/ (numberOfSamples - 1)
				/ sampleReplicates;
        beta_diversity_coef_4 = 2.0  / numberOfSamples
				/ (numberOfSamples - 1)
				/ sampleReplicates/numberOfMicrobePerHost;
		for (int i = 0; i < numberOfIndividual; i++) {
			host_index.add(i);
                        ancestryIndex[i]=i;
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

	public double interGenerationDistance() {
		if (getNumberOfGeneration() == 0)
			return 0;
		else {
			return Distance.getDistance(previous_microbiomes, microbiomeSum);
		}
	}

	public double environmentPopulationDistance() {
		if (getNumberOfGeneration() == 0)
			//TODO: can we merge these two??
			return Distance.getDistance(initialEnvironment, microbiomeSum);
		else
			return Distance.getDistance(environmentalContribution,
					microbiomeSum);
	}



	public void parentalInheritanceAndEnvironmentalAcquisition() {
		
        int[] oldAncestryIndex= new int[numberOfIndividual];
        System.arraycopy(ancestryIndex,0,oldAncestryIndex, 0,numberOfIndividual);
		for (int i = 0; i < numberOfIndividual; i++) {
            parentalIndex[i]=MathUtil.getNextInt(numberOfIndividual-1);
            ancestryIndex[i]=oldAncestryIndex[parentalIndex[i]];
			System.arraycopy(compositionOfIndividuals[parentalIndex[i]].getMicrobiome(),
                                0,parentalContribution[i], 0,numberOfEnvironmentalSpecies);
		}
		//NOTE: Maybe we can merge some of these environmentalContribution, mixedContribution, parentalContribution
		//These fields are not used at other places
		
		VectorAddition.additionOfVectors(environmentalContribution,
				percentageOfpooledOrFixed, 1 - percentageOfpooledOrFixed, microbiomeSum,
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
        
    public double[][] sample(int i){
            Set<Integer> sampleId=samples.get(i);
            double[][] sampleMicrobiomes=new double[numberOfSamples][];
            int m=0;
            for (Integer index:sampleId){
                 sampleMicrobiomes[m]=getIndividuals()[index].getMicrobiome();
                 m++;
            }
            return sampleMicrobiomes;
    }
        
    public double[] pooledSample(int i){
            Set<Integer> sampleId=samples.get(i);
            double[] temp_sum = new double[numberOfEnvironmentalSpecies];
            for (Integer index : sampleId) {
		VectorAddition.additionOfVectors(temp_sum, 1, 1.0/numberOfSamples,
							temp_sum,
							getIndividuals()[index].getMicrobiome());
				}
            return temp_sum;
    }

	public double alphaDiversity(boolean sampleOrNot) {
                a_diversity = 0;
		if (sampleOrNot) {
			for (Individual host : getIndividuals()) {
				a_diversity+=DiversityIndex.shannonWienerIndex(host.getMicrobiome(),numberOfMicrobePerHost);
			}
			a_diversity /=  numberOfIndividual;
		}
                else {
			for (Set<Integer> sample : samples) {
				for (Integer index : sample) {
					a_diversity+=DiversityIndex.shannonWienerIndex(getIndividuals()[index].getMicrobiome(),numberOfMicrobePerHost);
				}
			}
			a_diversity /= (numberOfSamples * sampleReplicates);
		}
		return a_diversity;

	}

	
	public double betaDiversity(boolean sampleOrNot) {
                b_diversity=0;
		if (sampleOrNot) {
                        if (numberOfIndividual>1){
			for (int i = 1; i < numberOfIndividual; i++) {
				for (int j = 0; j < i; j++) {
					b_diversity+=DiversityIndex.PiIndex(getIndividuals()[i].getMicrobiome(), getIndividuals()[j].getMicrobiome());
				}
			}
			
			b_diversity *= beta_diversity_coef_2;
		} }
                else {
                        if (numberOfSamples>1){
                        for (int index=0;index<sampleReplicates;index++){
                            double[][] temp_microbiomes=sample(index);
                            for (int i = 1; i < numberOfSamples; i++) {
                                for (int j = 0; j < i; j++) {
                                    b_diversity += DiversityIndex.PiIndex(temp_microbiomes[i], temp_microbiomes[j]);
                                }
                            }
                        }
			b_diversity *= beta_diversity_coef_3; 
		}}
                return b_diversity;
	}

	public double gammaDiversity(boolean sampleOrNot) {
                g_diversity=0;
		if (sampleOrNot) {
			g_diversity = DiversityIndex.shannonWienerIndex(microbiomeSum);
		} else {
                        for(int index=0;index<sampleReplicates;index++){
                            double[] temp_sum=pooledSample(index);
                            g_diversity+=DiversityIndex.shannonWienerIndex(temp_sum);
                        }
			g_diversity = g_diversity / sampleReplicates;	
		}
                return g_diversity;
	}
        
    public double BrayCurtis(boolean sampleOrNot){
            weighted_b_diversity=0;
            if (sampleOrNot){
                for (int i=1;i<numberOfIndividual;i++){
                    for (int j=0;j<i;j++){
                        weighted_b_diversity+= DiversityIndex.BrayCurtis(getIndividuals()[i].getMicrobiome(), getIndividuals()[j].getMicrobiome());
                    }
                }
                weighted_b_diversity*=beta_diversity_coef;
            }
            else{
                for (int index=0;index<sampleReplicates;index++){
                            double[][] temp_microbiomes=sample(index);
                            for (int i = 1; i < numberOfSamples; i++) {
                                for (int j = 0; j < i; j++) {
                                    weighted_b_diversity += DiversityIndex.BrayCurtis(temp_microbiomes[i], temp_microbiomes[j]);
                                }
                            }
                        }
			weighted_b_diversity *= beta_diversity_coef_4;
            }
            return weighted_b_diversity;
    }

	public String printOut() {
		StringBuilder sb = new StringBuilder();
		for (int i = 0; i < numberOfEnvironmentalSpecies; i++) {
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
	
    public boolean checkAncestry() {
        for (int i=1; i<numberOfIndividual;i++){
            if (ancestryIndex[i-1]!=ancestryIndex[i])
                return false;
            }
            return true;
        }
        
    public String indexPrintOut(int[] indexArray) {
		StringBuilder sb = new StringBuilder();
		for (int i = 0; i < indexArray.length; i++) {
			sb.append(indexArray[i]).append("\t");
		}
		return sb.toString().trim();
	}	

        
     


}

