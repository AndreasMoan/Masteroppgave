package main.java.HGSADCwSO.implementations;

import main.java.HGSADCwSO.files.*;
import main.java.HGSADCwSO.protocols.EducationProtocol;
import main.java.HGSADCwSO.protocols.FitnessEvaluationProtocol;
//import javafx.util.Pair;

import java.lang.reflect.Array;
import java.util.*;

public class EducationStandard implements EducationProtocol {

    protected ProblemData problemData;
    protected FitnessEvaluationProtocol fitnessEvaluationProtocol;
    protected PenaltyAdjustmentProtocol penaltyAdjustmentProtocol;
    protected boolean isRepair;
    protected int penaltyMultiplier;
    private int counter = 0;
    private boolean deterministic_search = true;

    public EducationStandard(ProblemData problemData, FitnessEvaluationProtocol fitnessEvaluationProtocol, PenaltyAdjustmentProtocol penaltyAdjustmentProtocol) {
        this.problemData = problemData;
        this.fitnessEvaluationProtocol = fitnessEvaluationProtocol;
        this.penaltyAdjustmentProtocol = penaltyAdjustmentProtocol;
        this.isRepair = false;
        this.penaltyMultiplier = 1;
    }

    @Override
    public void educate(Individual individual) { //TODO - klarer ikke å finne at B&L educater nye offspring med p^EDU?. Dette må såfall legges inn i main. Kan ikke ligge inni denne funksjonen da alle initielle individer skal educates


        fitnessEvaluationProtocol.evaluate(individual);
        System.out.println("Penalized cost before education: " + individual.getPenalizedCost());

        if (individual.getVesselTourChromosome().size() > 1) {
            new_inter_voyage_improvement(individual);
        }

        fitnessEvaluationProtocol.evaluate(individual);
        System.out.println("Penalized cost after ivi: " + individual.getPenalizedCost());

        neighbourhoodSearch(individual);

        fitnessEvaluationProtocol.evaluate(individual);
        System.out.println("Penalized cost after nhs: " + individual.getPenalizedCost());

        fitnessEvaluationProtocol.evaluate(individual);
        //System.out.println("Penalized cost after education: " + individual.getPenalizedCost());
    }


    @Override
    public void repairEducate(Individual individual, int penaltyMultiplier) {
        isRepair = true;
        this.penaltyMultiplier = penaltyMultiplier;

        educate(individual);

        isRepair = false;
        this.penaltyMultiplier = 1;
    }

    protected void neighbourhoodSearch(Individual individual) {

        HashMap<Integer, ArrayList<Integer>> chromosome = individual.getVesselTourChromosome(); //first key is a vessel number, and the value is a set of orders

        for (int vessel : chromosome.keySet()) {
            ArrayList<Integer> voyage = chromosome.get(vessel);
            if (voyage.size() != 0) {
                ArrayList<Integer> improvedVoyage = getImprovedVoyage(voyage, vessel);
                chromosome.put(vessel, new ArrayList<>(improvedVoyage));
            }
        }
        fitnessEvaluationProtocol.evaluate(individual);

        counter ++;
    }
    /*
    public void neighbourizeFriendlyInstallations(Individual individual){
        HashMap<Integer, ArrayList<Integer>> vesselTourChromosome = Utilities.deepCopyVesselTour(individual.getVesselTourChromosome());
        ArrayList<Pair<Integer, Integer>> zeroDistanceOrders = new ArrayList<>();
        for (int order1 = 1; order1 <= problemData.getNumberOfOrders(); order1++){
            for (int order2 = 1; order2 <= problemData.getNumberOfOrders(); order2++){
                if (order2 > order1)
                    if (problemData.getDistancesBetweenOrderNumbersByDay().get(0).get(order1).get(order2) == 0.00){
                        zeroDistanceOrders.add(new Pair<>(order1, order2));
                    }
            }
        }
        for (Pair<Integer, Integer> pair : zeroDistanceOrders) {
            int voyageFirstOrder = 0;
            int voyageSecondOrder = 0;
            int indexFirstOrder = 0;
            int indexSecondOrder = 0;
            int firstOrder = pair.getKey();
            int secondOrder = pair.getValue();
            boolean firstOrderInChromosome = false;
            boolean secondOrderInChromosome = false;
            for (Integer voyage : vesselTourChromosome.keySet()){
                if (vesselTourChromosome.get(voyage).contains(firstOrder)){
                    firstOrderInChromosome = true;
                    voyageFirstOrder = voyage;
                    indexFirstOrder = vesselTourChromosome.get(voyage).indexOf(firstOrder);
                }
                if (vesselTourChromosome.get(voyage).contains(secondOrder)){
                    secondOrderInChromosome = true;
                    voyageSecondOrder = voyage;
                    indexSecondOrder = vesselTourChromosome.get(voyage).indexOf(secondOrder);
                }
            }
            if(firstOrderInChromosome && secondOrderInChromosome){
                int orderToMove = secondOrder;
                vesselTourChromosome.get(voyageSecondOrder).remove(indexSecondOrder);
                vesselTourChromosome.get(voyageFirstOrder).add(indexFirstOrder, orderToMove);
            }
        }
        individual.setVesselTourChromosome(vesselTourChromosome);
    }


     */

    private ArrayList<Integer> getImprovedVoyage(ArrayList<Integer> voyage, int vessel) {
        ArrayList<Integer> orders = new ArrayList<>(voyage);
        ArrayList<Integer> untreatedOrders = new ArrayList<>(voyage);

        double oldVoyagePenalizedCost = fitnessEvaluationProtocol.getPenalizedCostOfVoyage(orders, vessel);

        while (untreatedOrders.size() > 0) {
            Integer u = Utilities.pickAndRemoveRandomElementFromList(untreatedOrders);
            ArrayList<Integer> neighbours = getNeighbours(u, orders);
            while (neighbours.size() > 0) {
                Integer v = Utilities.pickAndRemoveRandomElementFromList(neighbours);
                orders = doRandomMove(u, v, orders, vessel, oldVoyagePenalizedCost);
            }
        }
        return new ArrayList<>(orders);
    }


    private ArrayList<Integer> doRandomMove(Integer u, Integer v, ArrayList<Integer> orders, int vessel, double oldVoyagePenalizedCost) {
        ArrayList<Integer> unusedMoves = new ArrayList<>(); //list of unused moves
        Move move = new Move();
        for (int i = 0; i < move.getNumberOfMoves(); i++) {
            unusedMoves.add(i + 1); //moves are 1-indexed, as in Borthen & Loennechen (2016) and Vidal et al (2012)
        }

        while (unusedMoves.size() > 0) { //try moves in random order
            int moveNumber = Utilities.pickAndRemoveRandomElementFromList(unusedMoves);
            ArrayList<Integer> newOrders = move.doMove(u, v, orders, moveNumber);
            double newVoyagePenalizedCost;

            if (!isRepair) { //Normal education

                newVoyagePenalizedCost = fitnessEvaluationProtocol.getPenalizedCostOfVoyage(newOrders, vessel);
                counter ++;
            }
            else {           //Repair education
                double durationViolationPenalty = fitnessEvaluationProtocol.getDurationViolationPenalty() * penaltyMultiplier;
                double capacityViolationPenalty = fitnessEvaluationProtocol.getCapacityViolationPenalty() * penaltyMultiplier;
                double deadlineViolationPenalty = fitnessEvaluationProtocol.getDeadlineViolationPenalty() * penaltyMultiplier;

                oldVoyagePenalizedCost = fitnessEvaluationProtocol.getPenalizedCostOfVoyage(orders, vessel, durationViolationPenalty, capacityViolationPenalty, deadlineViolationPenalty);
                counter ++;

                newVoyagePenalizedCost = fitnessEvaluationProtocol.getPenalizedCostOfVoyage(newOrders, vessel, durationViolationPenalty, capacityViolationPenalty, deadlineViolationPenalty);
                counter ++;
            }
            if (newVoyagePenalizedCost < oldVoyagePenalizedCost) {
                return newOrders;
            }
        }
        return orders;
    }

    private ArrayList<Integer> getNeighbours(Integer order, ArrayList<Integer> orders) {
        double granularityThreshold = problemData.getHeuristicParameterDouble("Granularity threshold in RI"); //share of neighbourhood
        int numberOfNeighboursAllowed = (int) ((orders.size()-1) * granularityThreshold);
        if (numberOfNeighboursAllowed < 3){numberOfNeighboursAllowed = (int) ((orders.size()-1));} //if the neighbourhood is already small, don't make it smaller

        ArrayList<Integer> neighbours = new ArrayList<>(orders); //create set of neighbours
        neighbours.remove(order); //remove this order from the set

        //get the distance from all other installations (orders) to this installation (order) as (key,value)-pairs
        ArrayList<Map.Entry<Integer, Double>> distancesByOrderNumber = new ArrayList<>(problemData.getDistancesBetweenOrderNumbersByDay().get(0).get(order).entrySet());

        //remove (key,value)-pairs that are not neighbours
        ArrayList<Map.Entry<Integer, Double>> removeList = new ArrayList<>();
        for (Map.Entry<Integer, Double> distance : distancesByOrderNumber) { //removing the order considered from the distances-ArrayList
            if (!neighbours.contains(distance.getKey())) {
                removeList.add(distance);
            }
        }
        distancesByOrderNumber.removeAll(removeList);

        //sort the (key,value)-pairs by distance, having the pairs with the highest distance first
        Collections.sort(distancesByOrderNumber, Collections.reverseOrder(Utilities.getMapEntryWithDoubleComparator()));

        //removes the neighbours with the highest distance until the correct number of neighbours is obtained
        while (neighbours.size() > numberOfNeighboursAllowed) {
            //System.out.println(orders);
            neighbours.remove(distancesByOrderNumber.remove(0).getKey());
        }
        return neighbours;
    }






    public void new_inter_voyage_improvement(Individual individual) { //TODO still has problems finding improvements

        HashMap<Integer, ArrayList<Integer>> chromosome_copy = Utilities.deepCopyVesselTour(individual.getVesselTourChromosome());
        HashMap<Integer, ArrayList<Integer>> best_chromosome = chromosome_copy;
        fitnessEvaluationProtocol.evaluate(individual);
        counter++;
        double cost_to_beat = individual.getPenalizedCost();

        // Creating all results from a move

        for (int remove_vessel_number = 0; remove_vessel_number < chromosome_copy.size(); remove_vessel_number++) {

            if (chromosome_copy.get(remove_vessel_number).size() == 0) {
                continue;
            }

            for (int remove_index = 0; remove_index < chromosome_copy.get(remove_vessel_number).size(); remove_index++) {

                int order_number = chromosome_copy.get(remove_vessel_number).get(remove_index);

                for (int insert_vessel_number = 0; insert_vessel_number < chromosome_copy.size() && insert_vessel_number != remove_vessel_number; insert_vessel_number ++) {

                    for (int insertion_index = 0; insertion_index <= chromosome_copy.get(insert_vessel_number).size(); insertion_index++) {

                        HashMap<Integer, ArrayList<Integer>> changed_chromosome_copy = Utilities.deepCopyVesselTour(chromosome_copy);

                        changed_chromosome_copy.get(remove_vessel_number).remove(remove_index);
                        changed_chromosome_copy.get(insert_vessel_number).add(insertion_index, order_number);

                        Individual temp_individual = new Individual(changed_chromosome_copy, fitnessEvaluationProtocol);

                        fitnessEvaluationProtocol.evaluate(temp_individual);
                        double candidate_cost = temp_individual.getPenalizedCost();

                        if (candidate_cost < cost_to_beat) {
                            cost_to_beat = candidate_cost;
                            best_chromosome = changed_chromosome_copy;
                        }
                    }
                }
            }
        }
        individual.setVesselTourChromosome(best_chromosome);
    }

    public void inter_voyage_improvement(Individual individual) { //TODO This one does not currently work, find out why

        // System.out.println("Before: " + individual.getVesselTourChromosome());

        HashMap<Integer, ArrayList<Integer>> chromosome_copy = Utilities.deepCopyVesselTour(individual.getVesselTourChromosome());

        double cost_to_beat = individual.getPenalizedCost();

        // double costUnchangedChromosme = fitnessEvaluationProtocol.getPenalizedCostOfVoyage();

        ArrayList<Double> cost_reductions_from_removal = new ArrayList<>();
        ArrayList<Integer> remove_indicies = new ArrayList<>();
        ArrayList<Integer> remove_vessel_numbers = new ArrayList<>();


        for (int vessel_number = 0; vessel_number < chromosome_copy.size(); vessel_number++) {

            if (chromosome_copy.get(vessel_number).size() == 0) {
                continue;
            }

            ArrayList<Integer> voyage_copy = chromosome_copy.get(vessel_number);
            double cost_of_unchanged_voyage = fitnessEvaluationProtocol.getPenalizedCostOfVoyage(voyage_copy, vessel_number);
            counter++;

            for (int index_of_order_to_be_removed = 0; index_of_order_to_be_removed < voyage_copy.size(); index_of_order_to_be_removed++) {
                ArrayList<Integer> reduced_voyage = new ArrayList<>(voyage_copy);
                reduced_voyage.remove(index_of_order_to_be_removed);
                double cost_of_reduced_voyage = fitnessEvaluationProtocol.getPenalizedCostOfVoyage(reduced_voyage, vessel_number);
                counter++;

                double cost_reduction = Math.max(0, cost_of_reduced_voyage - cost_of_unchanged_voyage);

                cost_reductions_from_removal.add(cost_reduction);
                remove_indicies.add(index_of_order_to_be_removed);
                remove_vessel_numbers.add(vessel_number);
            }
        }

        ArrayList<Double> removal_probabilities = Utilities.normalize(cost_reductions_from_removal);

        int index_of_remove_decision = Utilities.getRandomElementFromDistribution(removal_probabilities);
        int remove_decision_vessel = remove_vessel_numbers.get(index_of_remove_decision);
        int remove_decision_index = remove_indicies.get(index_of_remove_decision);
        int number_of_removed_order = chromosome_copy.get(remove_decision_vessel).get(remove_decision_index);

        HashMap<Integer, ArrayList<Integer>> reduced_chromosome_copy = Utilities.deepCopyVesselTour(chromosome_copy);
        reduced_chromosome_copy.get(remove_decision_vessel).remove(remove_decision_index);

        ArrayList<Double> total_cost_reductions_from_move = new ArrayList<>();
        ArrayList<Individual> changed_individuals = new ArrayList<>();

        boolean found_improvement = false;

        for (int second_vessel_number = 0; second_vessel_number < chromosome_copy.size(); second_vessel_number++) {

            if (chromosome_copy.get(second_vessel_number).size() == 0) {
                continue;
            }

            for (int second_index = 0; second_index <= reduced_chromosome_copy.get(second_vessel_number).size(); second_index++) {

                HashMap<Integer, ArrayList<Integer>> changed_chromosome_copy = Utilities.deepCopyVesselTour(reduced_chromosome_copy);
                changed_chromosome_copy.get(second_vessel_number).add(second_index, number_of_removed_order);

                Individual temporary_individual = new Individual(changed_chromosome_copy, fitnessEvaluationProtocol);
                fitnessEvaluationProtocol.evaluate(temporary_individual);
                counter++;
                double changed_chromosome_cost = temporary_individual.getPenalizedCost();
                double total_cost_reduction_from_move = Math.max(0, cost_to_beat - changed_chromosome_cost);

                if (total_cost_reduction_from_move > 0) {
                    found_improvement = true;
                }

                total_cost_reductions_from_move.add(total_cost_reduction_from_move);
                changed_individuals.add(temporary_individual);
            }
        }

        if (found_improvement){
            ArrayList<Double> second_probabilities = Utilities.normalize(total_cost_reductions_from_move);
            int second_decision_index = deterministic_search ? Utilities.getBestElementFromDistribution(second_probabilities) : Utilities.getRandomElementFromDistribution(second_probabilities);
            Individual improved_individual = changed_individuals.get(second_decision_index);
            HashMap<Integer, ArrayList<Integer>> improved_chromosome = Utilities.deepCopyVesselTour(improved_individual.getVesselTourChromosome());
            individual.setVesselTourChromosome(improved_chromosome);
        }
    }



/*


            double cost_of_unchanged_voyage = fitnessEvaluationProtocol.getPenalizedCostOfVoyage(chromosome.get(vessel_number), vessel_number);
            counter++;

            int index_of_order_to_be_moved = get_index_of_order_to_be_moved(chromosome.get(vessel_number), vessel_number);

            if (index_of_order_to_be_moved != -1) {

                ArrayList<Integer> reduced_voyage = new ArrayList<>(chromosome.get(vessel_number));
                int number_of_order_to_be_moved = reduced_voyage.get(index_of_order_to_be_moved);
                reduced_voyage.remove(index_of_order_to_be_moved);
                double cost_reduction_from_move = cost_of_unchanged_voyage - fitnessEvaluationProtocol.getPenalizedCostOfVoyage(chromosome.get(vessel_number), vessel_number);
                counter++;

                int[] move_destination_info = get_move_destination_info(number_of_order_to_be_moved, vessel_number, cost_reduction_from_move, chromosome);

                if (move_destination_info[0] != -1) {

                    int second_index_number = move_destination_info[0];
                    int second_vessel_number = move_destination_info[1];

                    ArrayList<Integer> incremented_voyage = new ArrayList<>(chromosome.get(second_vessel_number));
                    incremented_voyage.add(second_index_number, number_of_order_to_be_moved);

                    chromosome.put(vessel_number, reduced_voyage);
                    chromosome.put(second_vessel_number, incremented_voyage);
                }
            }
        }
        individual.setVesselTourChromosome(chromosome);
        //System.out.println("After: " + individual.getVesselTourChromosome());
    }

    private int get_index_of_order_to_be_moved(ArrayList<Integer> voyage, int vessel){

        double cost_of_unchanged_voyage = fitnessEvaluationProtocol.getPenalizedCostOfVoyage(voyage, vessel);
        counter++;

        double lowest_cost_reduction = Double.POSITIVE_INFINITY;

        ArrayList<Double> remove_cost_reductions = new ArrayList<>();


        for (int index = 0; index < voyage.size(); index ++) {

            ArrayList<Integer> reduced_voyage = new ArrayList<>(voyage);
            reduced_voyage.remove(index);

            double cost_reduction = Math.max(0, cost_of_unchanged_voyage - fitnessEvaluationProtocol.getPenalizedCostOfVoyage(reduced_voyage, vessel));
            counter++;

            remove_cost_reductions.add(cost_reduction);

            if (cost_reduction < lowest_cost_reduction) {
                lowest_cost_reduction = cost_reduction;
            }
        }
        double sum = 0;
        for (int i = 0; i < remove_cost_reductions.size(); i++){
            remove_cost_reductions.set(i, remove_cost_reductions.get(i) - lowest_cost_reduction);
            sum += remove_cost_reductions.get(i);
        }
        if (sum == 0) {
            return -1;
        }

        ArrayList<Double> probabilities = Utilities.normalize(remove_cost_reductions);

        return Utilities.getRandomElementFromDistribution(probabilities);
    }

    private int[] get_move_destination_info(int order_number, int vessel_number, double cost_reduction_from_move, HashMap<Integer, ArrayList<Integer>> chromosome) {

        HashMap<Integer, ArrayList<Integer>> unchanged_chromosome = Utilities.deepCopyVesselTour(chromosome);

        ArrayList<Double> total_cost_reductions_from_move = new ArrayList<>();
        ArrayList<Integer> move_index = new ArrayList<>();
        ArrayList<Integer> move_vessel = new ArrayList<>();

        boolean found_improvement = false;

        for (int second_vessel_number = 0; second_vessel_number < chromosome.size() && second_vessel_number != vessel_number; second_vessel_number ++) {

            ArrayList<Integer> second_unchenged_voyage = unchanged_chromosome.get(second_vessel_number);

            if (second_unchenged_voyage.size() == 0){
                continue;
            }

            double cost_of_second_unchanged_voyage = fitnessEvaluationProtocol.getPenalizedCostOfVoyage(second_unchenged_voyage, second_vessel_number);
            counter++;

            for (int second_index = 0; second_index <= second_unchenged_voyage.size(); second_index++) {

                ArrayList<Integer> second_extended_voyage = new ArrayList<>(second_unchenged_voyage);
                second_extended_voyage.add(second_index, order_number);

                double cost_increment_from_move = fitnessEvaluationProtocol.getPenalizedCostOfVoyage(second_extended_voyage, second_vessel_number) - cost_of_second_unchanged_voyage;
                counter++;

                double total_cost_reduction_from_move = Math.max(0, cost_reduction_from_move - cost_increment_from_move);

                if (total_cost_reduction_from_move > 0) {
                    found_improvement = true;
                }

                total_cost_reductions_from_move.add(total_cost_reduction_from_move);
                move_index.add(second_index);
                move_vessel.add(second_vessel_number);
            }
        }
        if (found_improvement) {
            System.out.println(total_cost_reductions_from_move);
            ArrayList<Double> second_probabilities = Utilities.normalize(total_cost_reductions_from_move);
            int chosen_move_index;
            if (deterministic_search) {
                chosen_move_index = Utilities.getBestElementFromDistribution(second_probabilities);
            }
            else {
                chosen_move_index = Utilities.getRandomElementFromDistribution(second_probabilities);
            }
            return new int[] {move_index.get(chosen_move_index), move_vessel.get(chosen_move_index)};
        }
        return new int[] {-1, -1};
    }


 */
    /*

    private void intra_voyage_improvement(Individual individual) {
        Move move = new Move();
        HashMap<Integer, ArrayList<Integer>> chromosme_copy = Utilities.deepCopyVesselTour(individual.getVesselTourChromosome());
        for (int vessel_number : chromosme_copy.keySet()) {

            ArrayList<Integer> old_voyage = chromosme_copy.get(vessel_number);

            if (old_voyage.size() == 0) {
                continue;
            }

            double penalized_cost_of_old_voyage = fitnessEvaluationProtocol.getPenalizedCostOfVoyage(old_voyage, vessel_number);

            ArrayList<Double> list_of_cost_improvements = new ArrayList<>();
            ArrayList<Integer> list_of_first_indices = new ArrayList<>();
            ArrayList<Integer> list_of_second_indices = new ArrayList<>();
            ArrayList<Integer> list_of_move_numbers = new ArrayList<>();

            for (int first_index = 0; first_index < old_voyage.size(); first_index++) {
                for (int second_index = 0; second_index < old_voyage.size() && second_index != first_index; second_index++) {
                    ArrayList<Integer> new_voyage = new ArrayList<>(old_voyage);
                    new_voyage = move.doMove(first_index, second_index, new_voyage, 1);
                    double penalized_cost_of_new_voyage = fitnessEvaluationProtocol.getPenalizedCostOfVoyage(new_voyage,vessel_number);
                    double cost_improvement = Math.max(0, penalized_cost_of_old_voyage - penalized_cost_of_new_voyage);

                    list_of_cost_improvements.add(cost_improvement);
                    list_of_first_indices.add(first_index);
                    list_of_second_indices.add(second_index);
                    list_of_move_numbers.add(1);

                }
            }

            Utilities.normalize(list_of_cost_improvements)
        }
    }

     */

    /*

    public void mergeVoyages(Individual individual) {

        Set<Integer> departingVessels = new HashSet<> (individual.getDepartingVessels());
        if (departingVessels.size() > 1){ //if there is more than one departures the next departure day
            Set<Set<Integer>> allDepartingVesselCombinations = Utilities.cartesianProduct(departingVessels);
            Integer bestVesselToKeep = -1;
            Integer bestVesselToRemove = -1;
            ArrayList<Integer> bestNewVoyage = null;
            double bestCostReduction = 0;
            for (Set<Integer> vesselPair : allDepartingVesselCombinations){
                //Randomly select which vessel to remove and which to keep. Select randomly only if both vessels runs voyages of equal size
                Integer vesselNumberToKeep = Utilities.pickAndRemoveRandomElementFromSet(vesselPair);
                Integer vesselNumberToRemove = Utilities.pickAndRemoveRandomElementFromSet(vesselPair);
                if (individual.getVesselTourChromosome().get(vesselNumberToRemove).size() < individual.getVesselTourChromosome().get(vesselNumberToKeep).size() ) { //if one of the voyages are longer than the other, keep the shortest. You may get better results from picking multiple good insertions //TODO - ANDREAS ENIG?? Eller droppe å gjøre dette? Kan være bra for diversity
                    Integer tempVesselNumberToKeep = vesselNumberToKeep;
                    vesselNumberToKeep = vesselNumberToRemove;
                    vesselNumberToRemove = tempVesselNumberToKeep;
                }

                ArrayList<Integer> voyageToMergeInto = individual.getVesselTourChromosome().get(vesselNumberToKeep); //order sequence of the voyage to keep
                ArrayList<Integer> voyageToMove = individual.getVesselTourChromosome().get(vesselNumberToRemove); //order sequence of the voyage to remove
                double currentPenalizedCost = fitnessEvaluationProtocol.getPenalizedCostOfVoyage(voyageToMergeInto, vesselNumberToKeep) + fitnessEvaluationProtocol.getPenalizedCostOfVoyage(voyageToMove, vesselNumberToRemove);
                counter += 2;

                ArrayList<Integer> newVoyage = new ArrayList<>(voyageToMergeInto);
                //insert each order in voyageToMove into voyageToMergeInto
                for (Integer order : voyageToMove) {
                    VoyageInsertion bestInsertion = getBestInsertionIntoVoyage(vesselNumberToKeep, problemData.getOrdersByNumber().get(order), newVoyage);
                    int bestPos = bestInsertion.positionInVoyageToInsertInto;
                    newVoyage.add(bestPos, order);
                }

                double newPenalizedCost = fitnessEvaluationProtocol.getPenalizedCostOfVoyage(newVoyage, vesselNumberToKeep);
                counter ++;

                double costReduction = currentPenalizedCost-newPenalizedCost;
                if (costReduction > bestCostReduction){
                    bestVesselToKeep = vesselNumberToKeep;
                    bestVesselToRemove = vesselNumberToRemove;
                    bestNewVoyage = newVoyage;
                    bestCostReduction = costReduction;
                }
            }
            //if there exist a merger resulting in reduced cost, change the tour
            if (bestCostReduction > 0){
                //copy current giant tour
                HashMap<Integer, ArrayList<Integer>> vesselTour = Utilities.deepCopyVesselTour(individual.getVesselTourChromosome());

                vesselTour.put(bestVesselToKeep, bestNewVoyage);
                vesselTour.put(bestVesselToRemove, new ArrayList<>());
                individual.setVesselTourChromosome(vesselTour);
                fitnessEvaluationProtocol.evaluate(individual);
                counter ++;
            }
        }
    }

     */


    protected VoyageInsertion getBestInsertionIntoVoyage(Integer vesselTakingOverOrder, Order order, ArrayList<Integer> voyageToMergeInto) { //TESTET OG FUNKER

        int orderToAdd = order.getNumber();
        int indexWhereNewOrderIsPlaced = 0;
        double penalizedCostBeforeOrderAdded = fitnessEvaluationProtocol.getPenalizedCostOfVoyage(voyageToMergeInto, vesselTakingOverOrder);
        counter ++;
        double bestNewPenalizedCost = Double.MAX_VALUE;

        for (int index = 0; index < voyageToMergeInto.size()+1; index++){
            ArrayList<Integer> testVoyage = new ArrayList<>(voyageToMergeInto);
            testVoyage.add(index, orderToAdd);
            double testVoyagePenalizedCost = fitnessEvaluationProtocol.getPenalizedCostOfVoyage(testVoyage, vesselTakingOverOrder);
            counter ++;
            if (testVoyagePenalizedCost < bestNewPenalizedCost){
                indexWhereNewOrderIsPlaced = index;
                bestNewPenalizedCost = testVoyagePenalizedCost;
            }
        }
        double insertionCost = bestNewPenalizedCost - penalizedCostBeforeOrderAdded;
        return new VoyageInsertion(vesselTakingOverOrder, orderToAdd, indexWhereNewOrderIsPlaced, insertionCost);
    }


    protected Integer[] getCheapestVoyageAndPositionToInsertAnOrderTo(Order order, HashMap<Integer, ArrayList<Integer>> vesselTourChromosome, Integer forbiddenVoyageToInsertInto){
        HashMap<Integer, ArrayList<Integer>> vesselTour = Utilities.deepCopyVesselTour(vesselTourChromosome);
        double bestInsertionCost = Double.MAX_VALUE;
        VoyageInsertion bestVoyageInsertion = null;

        for (Integer voyage : vesselTourChromosome.keySet()) {//if the insertion cost is the same for many voyages, then the order is put into the first voyage with this insertion cost
            if (!voyage.equals(forbiddenVoyageToInsertInto) && vesselTourChromosome.get(voyage).size() != 0) {
                VoyageInsertion currentVoyage = getBestInsertionIntoVoyage(voyage, order, vesselTour.get(voyage));
                if (currentVoyage.getInsertionCost() <= bestInsertionCost) {
                    bestInsertionCost = currentVoyage.getInsertionCost();
                    bestVoyageInsertion = currentVoyage;
                }
            }
        }
        return new Integer[]{bestVoyageInsertion.getVessel(), bestVoyageInsertion.getPositionInVoyageToInsertInto()};
    } //.getVessel() may give nullptrexception because bestVoyageInsertion will always be assigned a value. Best insertionCost is never max_value



    protected void voyageReduction(Individual individual){

        //  1.	Evaluate the penalized cost of the unchanged individual
        Individual unchangedIndividual = new Individual (Utilities.deepCopyVesselTour(individual.getVesselTourChromosome()), fitnessEvaluationProtocol);
        Individual modifiedIndividual = new Individual (Utilities.deepCopyVesselTour(individual.getVesselTourChromosome()), fitnessEvaluationProtocol);

        HashMap<Integer, ArrayList<Integer>> modifiedChromosome = new HashMap<>(modifiedIndividual.getVesselTourChromosome());

        fitnessEvaluationProtocol.evaluate(unchangedIndividual);
        counter ++;
        double unchangedIndividualPenalizedCost = unchangedIndividual.getPenalizedCost();

        //  2. Checking that we have two or more voyages departing the next day
        int countVoyagesDeparting = unchangedIndividual.getDepartingVessels().size();
        if (countVoyagesDeparting > 1) {

            //  3.	Find the shortest voyage in the individual. If there exist more than one of the shortest voyage, chose the one with highest penalized cost to be removed from the individual.
            Integer voyageToRemove = selectVoyageAndSequenceToRemoveFromChromosome(individual);

            //  4.	Remove the orders from the chromosome and put them in a separate list
            ArrayList<Integer> ordersToReallocateIntFormat = new ArrayList<>();

            while (modifiedChromosome.get(voyageToRemove).size() != 0) { //while the voyage to remove is not empty
                ordersToReallocateIntFormat.add(modifiedChromosome.get(voyageToRemove).get(0));//get the first order in the voyage to remove
                modifiedChromosome.get(voyageToRemove).remove(0); //remove the order you just places in another list
            }

            //  5. Sort the orders that are to be removed from the terminated voyage in decending due dates.

            //Convert ordersToRelocate from Integer to Order-format
            ArrayList<Order> ordersToReallocateOrderFormat = new ArrayList<>();
            for (int i : ordersToReallocateIntFormat) {
                ordersToReallocateOrderFormat.add(problemData.getOrdersByNumber().get(i));
            }

            //sort
            ArrayList<Order> ordersToReallocateOrderFormatSorted = new ArrayList<>(ordersToReallocateOrderFormat);
            ordersToReallocateOrderFormatSorted.sort(Utilities.getDeadlineComparator());

            //sort ordersToReallocateIntFormat
                /*ArrayList<Integer> ordersToReallocateIntFormatSorted = new ArrayList<Integer>();
                for (Order order : ordersToReallocateOrderFormatSorted){
                    ordersToReallocateIntFormatSorted.add(order.getNumber());
                }*/

            //  6.	For each order in the list of orders that has to be reallocated
            //      a.	Find the least cost insertion for the order into the other voyages in the individual. Assign the order to the voyage that has the cheapest insertion and delete from the ordersToReallocate-list
            for (Order order : ordersToReallocateOrderFormatSorted) {
                Integer[] cheapestVoyageAndPositionToInsertAnOrderTo = getCheapestVoyageAndPositionToInsertAnOrderTo(order, modifiedChromosome, voyageToRemove); //TODO - legg på "if getcheapestvoy.... is not null...

                Integer insertOrderInVoyageNumber = cheapestVoyageAndPositionToInsertAnOrderTo[0];
                Integer insertAtPositionInVoyage = cheapestVoyageAndPositionToInsertAnOrderTo[1];

                modifiedChromosome.get(insertOrderInVoyageNumber).add(insertAtPositionInVoyage, order.getNumber());
            }

            //  7.	Evaluate the penalized cost of the new individual
            modifiedIndividual.setVesselTourChromosome(modifiedChromosome);
            fitnessEvaluationProtocol.evaluate(modifiedIndividual);
            counter ++;
            double modifiedIndividualPenalizedCost = modifiedIndividual.getPenalizedCost();

            //  8.	If the new individual has a lower penalized cost than the old, perform the voyage reduction

            if (modifiedIndividualPenalizedCost < unchangedIndividualPenalizedCost) {

                individual.setVesselTourChromosome(modifiedChromosome);
            }
        }
    }


    protected Integer selectVoyageAndSequenceToRemoveFromChromosome(Individual individual){
        // Find the shortest voyage in the individual. If there exist more than one of the shortest voyage, chose the one with highest penalized cost to be removed from the individual.
        Individual modifiedIndividual = new Individual (Utilities.deepCopyVesselTour(individual.getVesselTourChromosome()), fitnessEvaluationProtocol);
        HashMap<Integer, ArrayList<Integer>> modifiedChromosome = new HashMap<>(modifiedIndividual.getVesselTourChromosome());

        int shortestVoyageSize = Integer.MAX_VALUE;
        HashMap<Integer, ArrayList<Integer>> shortestVoyages = new HashMap<>(); //list of the shortest voyages (of equal size), in case there are more than one
        Integer voyageToBeTerminated = 0;
        ArrayList<Integer> voyageSequenceToBeTerminated = new ArrayList<>(); //the shortest voyage with the highest penalized cost

        for (Integer voyage : modifiedChromosome.keySet()){
            if (modifiedChromosome.get(voyage).size() > 0){ //making sure the ArrayList is not empty
                int currentVoyageSize = modifiedChromosome.get(voyage).size();
                if (currentVoyageSize == shortestVoyageSize){
                    shortestVoyages.put(voyage, modifiedChromosome.get(voyage));
                    shortestVoyageSize=currentVoyageSize;
                }
                else if (currentVoyageSize < shortestVoyageSize){
                    shortestVoyages = new HashMap<>();
                    shortestVoyages.put(voyage, modifiedChromosome.get(voyage));
                    shortestVoyageSize=currentVoyageSize;
                }
            }
        }
        if (shortestVoyages.keySet().size() == 1) {
            Integer voyageToRemove = 0;
            for (Integer i : shortestVoyages.keySet()){ //there is only one key in the keyset, and hence, the value of the key will be returned
                voyageToRemove += i;
            }
            voyageToBeTerminated = voyageToRemove;
            voyageSequenceToBeTerminated = shortestVoyages.get(voyageToBeTerminated);
        }
        else if (shortestVoyages.keySet().size() > 1){ //if there are multiple shortest voyages, choose the one with largest penalized cost. If no voyages have penalized cost, choose the first vessel in the chromosome.
            double highestPenalizedCostForShortestVoyages = 0;
            Integer voyageToRemove = 0;

            for (Integer vessel : shortestVoyages.keySet()){
                double currentPenalizedCost = fitnessEvaluationProtocol.getPenalizedCostOfVoyage(shortestVoyages.get(vessel), vessel);
                counter ++;
                if(currentPenalizedCost > highestPenalizedCostForShortestVoyages){
                    highestPenalizedCostForShortestVoyages = currentPenalizedCost;
                    voyageToRemove = vessel;
                }
            }
            voyageToBeTerminated = voyageToRemove;
            voyageSequenceToBeTerminated = shortestVoyages.get(voyageToBeTerminated);
        }

        return voyageToBeTerminated;
    }

//************************* NOT USED ****************************

    protected HashMap<Integer, ArrayList<Integer>> getCopyOfVesselTourWithoutOrder(Order order, HashMap<Integer, ArrayList<Integer>> vesselTour){

        HashMap<Integer, ArrayList<Integer>> vesselTourCopy = Utilities.deepCopyVesselTour(vesselTour);

        for (ArrayList<Integer> voyage : vesselTourCopy.values()){
            voyage.remove(Integer.valueOf(order.getNumber()));
        }
        return vesselTourCopy;
    }


    /*protected HashMap<Integer, ArrayList<Integer>> getCopyOfVesselTourWithoutSpecificOrders(ArrayList<Order> orders, Individual individual){

        HashMap<Integer, ArrayList<Integer>> vesselTourCopy = Utilities.deepCopyVesselTour(individual.getVesselTourChromosome());

        for (ArrayList<Integer> voyage : vesselTourCopy.values()){
            for (Order order : orders){
                voyage.remove(Integer.valueOf(order.getNumber()));
            }
        }
        return vesselTourCopy;
    }*/
}
