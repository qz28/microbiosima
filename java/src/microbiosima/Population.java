

package microbiosima;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;
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
	//FIXME: encapsulation these fields, way too many public
	final int numberOfMicrobePerHost;
	final int numberOfEnvironmentalSpecies;
	final int numberOfIndividual;
	final double environmentalFactor;
	final double percentageOfpooledOrFixed;
	
	private int numberOfGeneration;
	double[] initialEnvironment;
	double[] microbiomeSum;
	private double g_diversity;
	private double a_diversity;
	private double b_diversity;
        private double weighted_b_diversity;
	private String[] alignment;
	
	double[][] parentalContribution;
	double[] environmentalContribution;
	double[] mixedContribution;
	double coefficient;
	private int sampleReplicates;
	private int numberOfSamples;
	private double[] previous_microbiomes;
	private Individual[] compositionOfIndividuals;
	private List<Set<Integer>> samples = new ArrayList<>();
	private List<Integer> host_index = new ArrayList<>();
	private double beta_diversity_coef;
        private double beta_diversity_coef_2;
        private double beta_diversity_coef_3;
        private double beta_diversity_coef_4;
	
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
		
		multiNomialDist = new Multinomial2(numberOfEnvironmentalSpecies);
		
		beta_diversity_coef = 2.0  / numberOfIndividual
				/ (numberOfIndividual - 1)
				/ numberOfEnvironmentalSpecies;
                beta_diversity_coef_2= 1.0  / numberOfIndividual
				/ (numberOfIndividual - 1)
				/ numberOfMicrobePerHost;
                beta_diversity_coef_3=1.0/numberOfSamples/(numberOfSamples-1)/numberOfMicrobePerHost/sampleReplicates;
                beta_diversity_coef_4 = 2.0  / numberOfSamples
				/ (numberOfSamples - 1)
				/ numberOfEnvironmentalSpecies/sampleReplicates;
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
		alignment = new String[numberOfIndividual];//NOTE: pre declared this somewhere
		int i = 0;
		for (Individual host : getIndividuals()) {
			alignment[i] = host.microbial_sequences();
			i++;//NOTE: Since you use i++ here, might be better with standard for loop
//			for (int j = 0; j < compositionOfIndividuals.length; j++) {
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
		
		for (int i = 0; i < numberOfIndividual; i++) {
			System.arraycopy(compositionOfIndividuals[MathUtil.getNextInt(numberOfIndividual-1)].getMicrobiome(),
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
		
//		for (int i = 0; i < numberOfEnvironmentalSpecies; i++) {
//			previous_microbiomes[i] = microbiomeSum[i];
//		}
		//NOTE: System.arraycopy is faster
		System.arraycopy(microbiomeSum, 0, previous_microbiomes, 0, microbiomeSum.length);

		//example of IDE overdo it
//		setNumberOfGeneration(getNumberOfGeneration() + 1);
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
					m++;//NOTE: for loop or arrayList, track one less variable
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

	private static double getPlogP(double abundance) {
		return getPlogP(abundance, 1);
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
        
     


}

