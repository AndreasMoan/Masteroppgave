package main.java.HGSADCwSO.files;

public class SailingLeg {

    private double departureTime;
    private double arrivalTime;
    private double startServicingTime;
    private double endServicingTime;

    private double speed;

    private double cost;

    private int order;

    public SailingLeg(double cost, double departureTime, double arrivalTime, double startServicingTime, double endServicingTime, double speed, int orderNumber){
        this.cost = cost;
        this.departureTime = departureTime;
        this.arrivalTime = arrivalTime;
        this.startServicingTime = startServicingTime;
        this.endServicingTime = endServicingTime;
        this.speed = speed;
        this.order = orderNumber;
    }

    public void setCost(double cost) {
        this.cost = cost;
    }

    public void setArrivalTime(double arrivalTime) {
        this.arrivalTime = arrivalTime;
    }

    public void setDepartureTime(double departureTime) {
        this.departureTime = departureTime;
    }

    public void setEndServicingTime(double endServicingTime) {
        this.endServicingTime = endServicingTime;
    }

    public void setOrder(int order) {
        this.order = order;
    }

    public void setSpeed(double speed) {
        this.speed = speed;
    }

    public void setStartServicingTime(double startServicingTime) {
        this.startServicingTime = startServicingTime;
    }

    public double getCost() {
        return cost;
    }

    public double getArrivalTime() {
        return arrivalTime;
    }

    public double getDepartureTime() {
        return departureTime;
    }

    public double getEndServicingTime() {
        return endServicingTime;
    }

    public double getSpeed() {
        return speed;
    }

    public double getStartServicingTime() {
        return startServicingTime;
    }

    public int getOrder() {
        return order;
    }
}
