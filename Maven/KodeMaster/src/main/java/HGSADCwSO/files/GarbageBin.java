package HGSADCwSO.files;

import main.java.HGSADCwSO.implementations.Move;

import java.util.ArrayList;
import java.util.HashMap;

public class GarbageBin {

    /*

    private void intra_voyage_improvement(Individual individual) {
        Move move = new Move();
        HashMap<Integer, ArrayList<Integer>> chromosme_copy = Utilities.deepCopyVesselTour(individual.getVesselTourChromosome());

        int local_counter = 0;

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
            boolean improvement = false;

            for (int first_index = 0; first_index < old_voyage.size(); first_index++) {
                for (int second_index = 0; second_index < old_voyage.size() && second_index != first_index; second_index++) {
                    for (int move_number = 1; move_number < 8; move_number++) {

                        try {


                            ArrayList<Integer> new_voyage = new ArrayList<>(old_voyage);
                            new_voyage = move.doMove(first_index, second_index, new_voyage, move_number);
                            double penalized_cost_of_new_voyage = fitnessEvaluationProtocol.getPenalizedCostOfVoyage(new_voyage,vessel_number);
                            counter++;
                            double cost_improvement = Math.max(0, penalized_cost_of_old_voyage - penalized_cost_of_new_voyage);

                            if (cost_improvement > 0) {
                                improvement = true;
                                // - System.out.println(old_voyage + " " + first_index + " " + second_index + " " + move_number + " Improvement found at: " + local_counter + " with move number: " + move_number);
                                int a = improvements_by_move_number.get(move_number-1);
                                a++;
                                improvements_by_move_number.set(move_number -1, a);
                            }
                            else {
                                // - System.out.println(old_voyage + " " + first_index + " " + second_index + " " + move_number);
                            }

                            list_of_cost_improvements.add(cost_improvement);
                            list_of_first_indices.add(first_index);
                            list_of_second_indices.add(second_index);
                            list_of_move_numbers.add(move_number);

                            local_counter ++;
                        }
                        catch (Exception e) {
                            // - System.out.println("Skipped " + move_number);
                            continue;
                        }
                    }
                }
            }

            // - System.out.println(improvements_by_move_number);

            ArrayList<Double> probabilities = Utilities.normalize(list_of_cost_improvements);

            if (probabilities.size() == 0 || !improvement) {
                continue;
            }

            // - System.out.println("p: " + probabilities);
            int index_of_picked_improvement;

            if (deterministic_improvement) {
                index_of_picked_improvement = Utilities.getBestElementFromDistribution(probabilities);
            }
            else {
                index_of_picked_improvement = Utilities.getRandomElementFromDistribution(probabilities);
            }


            // - System.out.println("size: " + list_of_cost_improvements.size() + " index: " + index_of_picked_improvement);

            int first_index = list_of_first_indices.get(index_of_picked_improvement);
            int second_index = list_of_second_indices.get(index_of_picked_improvement);
            int move_number = list_of_move_numbers.get(index_of_picked_improvement);
            ArrayList<Integer> new_tour = new ArrayList<>(chromosme_copy.get(vessel_number));

            move.doMove(first_index, second_index, new_tour, move_number);

            chromosme_copy.put(vessel_number, new_tour);
        }
        individual.setVesselTourChromosome(chromosme_copy);
    }

     */
}
