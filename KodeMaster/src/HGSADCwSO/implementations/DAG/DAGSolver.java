package HGSADCwSO.implementations.DAG;

import HGSADCwSO.Order;
import HGSADCwSO.ProblemData;

import java.util.*;

import static java.lang.Math.*;

public class DAGSolver {

    /*
    private double shortestPathCost;
    private double shortestFeasiblePathCost;
    private boolean feasibility;

    private ProblemData problemData;

    private Graph graph;

    private double multiplier;

    private HashMap<Integer, Order> orderByNumber;
    
    private HashMap<ArrayList<Integer>, Graph> cachedGraphs = new HashMap<>();


    public DAGSolver(ProblemData problemData) {

        this.problemData = problemData;

        this.multiplier = this.problemData.getHeuristicParameterDouble("Number of time periods per hour");

        this.orderByNumber = problemData.getOrdersByNumber();

    }

    public void solve(ArrayList<Integer> voyage, int vesselReturnTime) {
        if (voyage.size() == 0) {
            shortestPathCost = 0;
            shortestFeasiblePathCost = 0;
        }
        else {
            if (cachedGraphs.containsKey(voyage)) {
                this.graph = cachedGraphs.get(voyage);
            }
            else {
                this.graph = new Graph(voyage, problemData);
                cachedGraphs.put(voyage, graph);
            }
            doDijkstra(vesselReturnTime);
            setCost(vesselReturnTime);
            setFeasibility();
        }
    }


 */
}
