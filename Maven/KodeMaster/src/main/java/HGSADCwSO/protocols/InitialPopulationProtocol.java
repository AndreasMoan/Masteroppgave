package main.java.HGSADCwSO.protocols;

import main.java.HGSADCwSO.files.Individual;

public interface InitialPopulationProtocol {

    public Individual createIndividual();

    public int getNumberOfConstructionHeuristicRestarts();

}
