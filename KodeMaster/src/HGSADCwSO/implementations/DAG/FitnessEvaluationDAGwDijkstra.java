package HGSADCwSO.implementations.DAG;

import HGSADCwSO.Individual;
import HGSADCwSO.ProblemData;
import HGSADCwSO.Vessel;
import HGSADCwSO.Voyage;
import HGSADCwSO.protocols.FitnessEvaluationProtocol;

import java.util.ArrayList;

public class FitnessEvaluationDAGwDijkstra implements FitnessEvaluationProtocol {

    private ProblemData problemData;

    FitnessEvaluationDAGwDijkstra(ProblemData problemData){
        this.problemData = problemData;
    }


    public void evaluate(Individual individual){

        double fitness = 0;

        int multiplier = (int) (1/problemData.getHeuristicParameterDouble("Length of time period"));

        for (Vessel vessel : problemData.getVessels()){
            DAG graph = new DAG(individual.getVesselTourChromosome().get(vessel.getNumber()), vessel.getReturnDay()*24*multiplier, this.problemData);
            fitness += graph.getShortestFeasiblePathCost();
        }

        individual.setFitness(fitness);
    }
}
