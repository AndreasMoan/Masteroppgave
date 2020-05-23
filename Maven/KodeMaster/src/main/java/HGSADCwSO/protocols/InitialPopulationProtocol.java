package HGSADCwSO.protocols;

import HGSADCwSO.files.Individual;

public interface InitialPopulationProtocol {

    public Individual createIndividual();

    public int getNumberOfConstructionHeuristicRestarts();

}
