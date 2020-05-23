package main.java.HGSADCwSO.protocols;

import main.java.HGSADCwSO.files.Individual;

import java.util.ArrayList;

public interface SurvivorSelectionProtocol {

    public void selectSurvivors(ArrayList<Individual> subpopulation, ArrayList<Individual> otherSubpopulation, FitnessEvaluationProtocol fitnessEvaluationProtocol);
    public ArrayList<Individual> getClones(ArrayList<Individual> subpopulation, FitnessEvaluationProtocol fitnessEvaluationProtocol);


}
