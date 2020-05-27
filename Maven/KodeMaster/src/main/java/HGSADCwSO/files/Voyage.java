package main.java.HGSADCwSO.files;

import main.java.HGSADCwSO.implementations.DAG.Edge;
import main.java.HGSADCwSO.implementations.DAG.Node;

import java.util.ArrayList;
import java.util.Collections;

public class Voyage {

    private Vessel vessel;
    private Node best_solution_node;

    public Voyage(Vessel vessel, Node best_solution_node) {
        this.vessel = vessel;
        this.best_solution_node = best_solution_node;
    }

    public void printVoyage() {
        ArrayList<Edge> solution_edges = best_solution_node.getAllBestParentEdges();
        Collections.sort(solution_edges, Utilities.getEdgeTimeComparator());
        for (int i = 0; i < solution_edges.size(); i++) {
            Edge edge = solution_edges.get(i);
            // TODO: Print
        }
    }
}
