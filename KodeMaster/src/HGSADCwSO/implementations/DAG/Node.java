package HGSADCwSO.implementations.DAG;

import java.util.ArrayList;

public class Node {

    private double time;
    private int orderNumber;
    private ArrayList<Edge> parentEdges;
    private ArrayList<Edge> childEdges;

    private double bestCost;
    private Edge bestParentEdge;

    private boolean feasibility;

    public Node(double time, int orderNumber){
        this.time = time;
        this.orderNumber = orderNumber;
        this.childEdges = new ArrayList<Edge>();
        this.parentEdges = new ArrayList<Edge>();
    }

    public void setFeasibility(boolean feasibility) {
        this.feasibility = feasibility;
    }

    public boolean getFeasibility() {
        return feasibility;
    }

    public void addParentEdge(Edge edge){
        parentEdges.add(edge);
    }

    public void addChildEdge(Edge edge){
        childEdges.add(edge);
    }

    public ArrayList<Edge> getParentEdges() {
        return parentEdges;
    }

    public ArrayList<Edge> getChildEdges() {
        return childEdges;
    }

    public double getTime() {
        return time;
    }

    public int getOrderNumber() {
        return orderNumber;
    }

    public double getBestCost() {
        return bestCost;
    }

    public void setBestCost(double bestCost) {
        this.bestCost = bestCost;
    }

    public Edge getBestParentEdge() {
        return bestParentEdge;
    }

    public void setBestParentEdge(Edge bestParentEdge) {
        this.bestParentEdge = bestParentEdge;
    }
}
