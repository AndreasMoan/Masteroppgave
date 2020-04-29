package HGSADCwSO.implementations;

import HGSADCwSO.HGSADCwSOmain;
import HGSADCwSO.Individual;
import HGSADCwSO.ProblemData;
import HGSADCwSO.Utilities;
import HGSADCwSO.protocols.FitnessEvaluationProtocol;
import HGSADCwSO.protocols.SurvivorSelectionProtocol;

import java.util.ArrayList;
import java.util.Collections;

public class SurvivorSelectionStandard implements SurvivorSelectionProtocol {
    
    private ProblemData problemData;


    public SurvivorSelectionStandard(ProblemData problemData) {
        this.problemData = problemData;
    }


    public void selectSurvivors(ArrayList<Individual> subpopulation, ArrayList<Individual> otherSubpopulation, FitnessEvaluationProtocol fitnessEvaluationProtocol) {
        int populationSize = problemData.getHeuristicParameterInt("Population size");
        ArrayList<Individual> clones = getClones(subpopulation, fitnessEvaluationProtocol);
        //sorts the clones and the subpopulation so that the individuals with the highest biased fitness are first
        Collections.sort(clones, Collections.reverseOrder(Utilities.getBiasedFitnessComparator()));
        Collections.sort(subpopulation, Collections.reverseOrder(Utilities.getBiasedFitnessComparator()));
        while (subpopulation.size() > populationSize) {
            if (clones.size() > 0) {
                HGSADCwSOmain.removeFromSubpopulation(subpopulation, clones.remove(0), otherSubpopulation, fitnessEvaluationProtocol, true);
            }
            else {
                HGSADCwSOmain.removeFromSubpopulation(subpopulation, subpopulation.get(0), otherSubpopulation , fitnessEvaluationProtocol, true);
            }
            Collections.sort(clones, Collections.reverseOrder(Utilities.getBiasedFitnessComparator()));
            Collections.sort(subpopulation, Collections.reverseOrder(Utilities.getBiasedFitnessComparator()));
        }
    }


    public ArrayList<Individual> getClones (ArrayList<Individual> subpopulation, FitnessEvaluationProtocol fitnessEvaluationProtocol) {
        ArrayList<Individual> clones = new ArrayList<Individual>();
        for (int i = 0; i < subpopulation.size()-1; i++) {
            for (int j = i+1; j < subpopulation.size(); j++) {
                Individual individual = subpopulation.get(i);
                Individual otherIndividual = subpopulation.get(j);
                if ((fitnessEvaluationProtocol.getHammingDistance(individual, otherIndividual) == 0)
                        && !(clones.contains(individual) || clones.contains(otherIndividual))) {
                    clones.add(otherIndividual);
                }
            }
        }
        return clones;
    }
}
