package microbiosima;

import utils.random.Multinomial2;

/**
 *
 * @author John
 */
public class Individual {

	double[] microbiome;
	private int numberEnvironmentalSpecies;
	
	
	public Individual(Multinomial2 MultDist, int nomph, int noes) {
		numberEnvironmentalSpecies = noes;
		microbiome = new double[numberEnvironmentalSpecies];
                MultDist.multisample(microbiome, nomph);
		numberMicrobePerHost = nomph;
	}

	public String microbial_sequences() {
		char[] microbiome_sequence = new char[numberEnvironmentalSpecies];
		for (int i = 0; i < numberEnvironmentalSpecies; i++) {
			if (microbiome[i] >=0.5) {
				microbiome_sequence[i] = '1';
			} else {
				microbiome_sequence[i] = '0';
			}
		}
		return new String(microbiome_sequence);
	}


	public double[] getMicrobiome() {
		return microbiome;
	}
        
	public String printOut() {
		StringBuilder sb = new StringBuilder();
		for (int i = 1; i < numberEnvironmentalSpecies; i++) {
			sb.append(microbiome[i]).append("\t");
		}
		return sb.toString().trim();
	}

	@Deprecated
	private int numberMicrobePerHost;
	
}


