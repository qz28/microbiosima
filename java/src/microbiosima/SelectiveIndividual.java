/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package microbiosima;

import utils.random.MathUtil;
import utils.random.Multinomial2;


public class SelectiveIndividual extends Individual {
    
    private double hostFitness;
    int[] microbeGenesFitness;
    int[] hostGenesFitness;
    private int numOfMGenes;
    double[] microbeFitnessRecords;
    double[] hostFitnessRecords;

    public SelectiveIndividual(Multinomial2 MultDist, int nomph, int noes, SelectiveSpeciesRegistry ssr,boolean HMS_or_TMS) {
        super(MultDist, nomph, noes);
        numOfMGenes=ssr.numOfTolGenes;
        microbeGenesFitness=new int[numOfMGenes];
        hostGenesFitness=new int[numOfMGenes];
        microbeFitnessRecords=new double[noes];
        hostFitnessRecords=new double[noes];
		if (HMS_or_TMS){
			for(int i=0;i<numOfMGenes;i++){
				microbeGenesFitness[i]=MathUtil.getNextInt(2)-1;
			}
		}else microbeGenesFitness=ssr.mfr;
		hostGenesFitness=ssr.hfr;
        ssr.getFitness(microbeFitnessRecords,hostFitnessRecords, microbeGenesFitness,hostGenesFitness);
        hostFitness=ssr.getTotalFitness(microbiome,hostFitnessRecords);
    }
    
    public double getHostFitness(double hostSelectionCoef){
        if(hostSelectionCoef==0)
            return hostFitness;
        else
            return Math.pow(hostSelectionCoef,hostFitness);
    }
    
    public double getCosTheta(){
        double innerProduct=0;
        double Mlength=0;
        double Hlength=0;
        for(int i=0;i<numOfMGenes;i++){
            innerProduct+=microbeGenesFitness[i]*hostGenesFitness[i];
            Mlength+=microbeGenesFitness[i]*microbeGenesFitness[i];
            Hlength+=hostGenesFitness[i]*hostGenesFitness[i];
        }
        return innerProduct/Math.sqrt(Hlength*Mlength);
    }
	
    public double[]  goodVersusBad(){
        double[] bacteriaContents=new double[3];
        for(int i=0;i<hostFitnessRecords.length;i++){
            if (hostFitnessRecords[i]<0)
                 bacteriaContents[0]+=getMicrobiome()[i];
            else if (hostFitnessRecords[i]==0)
                 bacteriaContents[1]+=getMicrobiome()[i];
            else
                 bacteriaContents[2]+=getMicrobiome()[i];
        }
        return bacteriaContents;
    }
  
    public String printBacteriaContents(){
       StringBuilder sb=new StringBuilder();
       double[] bacteriaContents=goodVersusBad();
       for (int i=0;i<bacteriaContents.length;i++){
            sb.append(bacteriaContents[i]).append('\t');
       }
       return sb.toString().trim();
    }
    
}
