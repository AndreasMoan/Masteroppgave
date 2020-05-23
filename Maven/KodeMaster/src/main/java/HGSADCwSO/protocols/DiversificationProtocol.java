package main.java.HGSADCwSO.protocols;

import java.util.ArrayList;

public interface DiversificationProtocol {

    public void incrementCounters();

    public boolean isDiversifyIteration();

    public boolean isStoppingIteration() ;

    public void resetDiversificationCounter() ;

    public void addDiversification(int iteration) ;

    public ArrayList<Integer> getDiversificationNumbers() ;

    public void updateIterationsSinceImprovementCounter(boolean improvingSolutionFound) ;

    
}