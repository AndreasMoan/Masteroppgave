package HGSADCwSO.implementations.DAG;

import HGSADCwSO.*;
import HGSADCwSO.implementations.FitnessEvaluationBaseline;
import HGSADCwSO.protocols.FitnessEvaluationProtocol;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Map;


public class FitnessEvaluationDAG extends FitnessEvaluationBaseline {

    private ProblemData problemData;
    private double nCloseProp;
    protected double nEliteProp;
    private double numberOfOrders;
    private double durationViolationPenalty;
    private double capacityViolationPenalty;
    private double deadlineViolationPenalty;
    private HashMap<Individual, HashMap<Individual, Double>> hammingDistances;

    public FitnessEvaluationDAG(ProblemData problemData){
        super(problemData);
        this.problemData = problemData;
    }


    @Override
    public void evaluate(Individual individual) {

        Genotype genotype = individual.getGenotype();

        boolean feasibility = true;
        double fitness = 0;

        int multiplier = (int) problemData.getHeuristicParameterDouble("Number of time periods per hour");

        for (Vessel vessel : problemData.getVessels()){
            DAG graph = new DAG(genotype.getVesselTourChromosome().get(vessel.getNumber()), vessel.getReturnDay()*24*multiplier, this.problemData);
            fitness += graph.getShortestFeasiblePathCost();
            feasibility = graph.getFeasibility() && feasibility;

        }

        System.out.println("Fitness this iteration is equal to:          " + fitness);

        individual.setFitness(fitness);
        individual.setFeasibility(feasibility);
    }
    
}
