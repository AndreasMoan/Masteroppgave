package main.java.HGSADCwSO.implementations;

import main.java.HGSADCwSO.files.ProblemData;
import main.java.HGSADCwSO.protocols.DiversificationProtocol;

import java.util.ArrayList;

public class DiversificationAndStopping implements DiversificationProtocol {

    private int iterationsBeforeDiversify, iterationsBeforeStopping; // Number of iterations without improvement before diversifying
    private int diversificationIterationsWithoutImprovement, iterationsWithoutImprovement;
    private ArrayList<Integer> diversificationNumbers;
    private long startTime, maxTime;
    private ProblemData problemData;

    public DiversificationAndStopping(ProblemData problemData) {
        this.iterationsBeforeDiversify = problemData.getHeuristicParameterInt("Iterations before diversify");
        this.iterationsBeforeStopping = problemData.getHeuristicParameterInt("Iterations before stopping");
        this.diversificationIterationsWithoutImprovement = 0;
        this.iterationsWithoutImprovement = 0;
        this.diversificationNumbers = new ArrayList<Integer>(); //iteration numbers where diversification is performed

        startTime = System.nanoTime()/1000000000; //seconds

        if (maxTime <= 0){
            this.maxTime = Long.MAX_VALUE;
        }
        else {
            this.maxTime = problemData.getHeuristicParameterInt("Max time");
        }
    }

    public void incrementCounters(){
        diversificationIterationsWithoutImprovement++;
        iterationsWithoutImprovement++;
    }

    public boolean isDiversifyIteration(){
        return diversificationIterationsWithoutImprovement >= iterationsBeforeDiversify;
    }

    public boolean isStoppingIteration() {
        long currentTime = System.nanoTime()/1000000000;
        long elapsedTime = currentTime - startTime;
        return iterationsWithoutImprovement >= iterationsBeforeStopping || elapsedTime >= maxTime;
    }

    public void resetDiversificationCounter() {
        diversificationIterationsWithoutImprovement = 0;
    }

    private void resetStoppingCounter() { iterationsWithoutImprovement = 0; }

    public void addDiversification(int iteration) {
        diversificationNumbers.add(iteration);
    }

    public ArrayList<Integer> getDiversificationNumbers() {
        return diversificationNumbers;
    }

    public void updateIterationsSinceImprovementCounter(boolean improvingSolutionFound) {
        if (improvingSolutionFound){
            resetStoppingCounter();
            resetDiversificationCounter();
        }
        else {
            incrementCounters();
        }
    }
}