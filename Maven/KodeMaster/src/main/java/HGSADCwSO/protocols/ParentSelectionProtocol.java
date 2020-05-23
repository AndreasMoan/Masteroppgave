package HGSADCwSO.protocols;


import HGSADCwSO.files.Individual;

import java.util.ArrayList;

public interface ParentSelectionProtocol {

    public ArrayList<Individual> selectParents(ArrayList<Individual> population);


}
