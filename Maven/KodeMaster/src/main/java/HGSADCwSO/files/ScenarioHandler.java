package main.java.HGSADCwSO.files;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;

public class ScenarioHandler {

    ArrayList<Installation> installations;
    ArrayList<ArrayList<Integer>> order_sizes = new ArrayList<>();
    ArrayList<ArrayList<Integer>> order_deadlines = new ArrayList<>();
    ArrayList<ArrayList<Integer>> vessel_return_days = new ArrayList<>();
    HashMap<Integer,Double> standard_order_size = new HashMap<>();
    HashMap<Integer, Double> order_size_variations = new HashMap<>();
    
    public ScenarioHandler(ArrayList<Installation> installations){
        this.installations = installations;
        add_all_scenarios();
        add_scenario_data();
    }
    
    private void add_all_scenarios(){
        order_sizes.add(new ArrayList<>(Arrays.asList(0, 0, 0, 0, 4, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 0, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))); //
        order_sizes.add(new ArrayList<>(Arrays.asList(3, 0, 0, 0, 0, 0, 0, 0, 4, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 4, 0, 0, 0, 0, 0, 3, 0)));
        order_sizes.add(new ArrayList<>(Arrays.asList(0, 3, 3, 0, 0, 0, 0, 0, 2, 3, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))); //
        order_sizes.add(new ArrayList<>(Arrays.asList(0, 2, 0, 0, 0, 0, 3, 0, 0, 0, 4, 0, 2, 0, 0, 2, 0, 0, 1, 0, 0, 2, 0, 0, 0, 0, 0)));
        order_sizes.add(new ArrayList<>(Arrays.asList(0, 0, 3, 4, 3, 0, 0, 0, 0, 0, 2, 3, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 5, 0, 0))); //
        order_sizes.add(new ArrayList<>(Arrays.asList(0, 0, 4, 0, 2, 0, 0, 3, 4, 0, 0, 0, 0, 5, 0, 2, 3, 0, 0, 0, 0, 3, 0, 0, 0, 2, 0)));
        order_sizes.add(new ArrayList<>(Arrays.asList(0, 0, 4, 0, 2, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 3, 0, 4, 4, 5, 3, 0, 0, 0, 0, 3))); // 6
        order_sizes.add(new ArrayList<>(Arrays.asList(0, 0, 4, 0, 2, 4, 0, 0, 0, 2, 3, 0, 0, 0, 3, 0, 0, 4, 0, 2, 0, 0, 0, 3, 3, 0, 1)));
        order_sizes.add(new ArrayList<>(Arrays.asList(2, 3, 2, 4, 0, 0, 0, 0, 4, 3, 0, 0, 0, 0, 0, 3, 0, 0, 0, 3, 3, 4, 5, 0, 0, 0, 0)));
        order_sizes.add(new ArrayList<>(Arrays.asList(2, 0, 0, 4, 0, 0, 3, 0, 0, 3, 3, 0, 2, 2, 0, 4, 1, 0, 4, 5, 0, 0, 4, 0, 0, 0, 4)));
        order_sizes.add(new ArrayList<>(Arrays.asList(0, 2, 3, 0, 0, 3, 3, 2, 1, 0, 0, 0, 0, 4, 3, 2, 5, 0, 0, 0, 0, 0, 4, 3, 3, 0, 0)));
        order_sizes.add(new ArrayList<>(Arrays.asList(2, 2, 0, 3, 0, 1, 0, 5, 0, 2, 0, 4, 3, 4, 0, 0, 0, 3, 4, 2, 1, 0, 0, 2, 0, 4, 0)));
        order_sizes.add(new ArrayList<>(Arrays.asList(0, 2, 3, 2, 3, 4, 5, 0, 0, 0, 3, 4, 5, 3, 0, 0, 0, 2, 2, 1, 0, 0, 2, 3, 0, 0, 0)));
        order_sizes.add(new ArrayList<>(Arrays.asList(0, 0, 3, 0, 0, 2, 0, 2, 3, 1, 3, 0, 0, 4, 2, 0, 4, 3, 4, 2, 3, 2, 3, 0, 3, 2, 0)));
        order_sizes.add(new ArrayList<>(Arrays.asList(3, 3, 2, 4, 5, 1, 0, 0, 0, 3, 3, 5, 5, 2, 0, 0, 0, 1, 3, 5, 0, 0, 0, 4, 5, 0, 0)));
        order_sizes.add(new ArrayList<>(Arrays.asList(3, 3, 1, 0, 0, 2, 4, 0, 2, 2, 0, 3, 0, 0, 4, 5, 0, 3, 1, 2, 3, 4, 0, 2, 3, 3, 3)));
        order_sizes.add(new ArrayList<>(Arrays.asList(3, 2, 3, 2, 3, 0, 4, 0, 3, 1, 0, 3, 3, 3, 3, 0, 2, 1, 0, 2, 2, 3, 4, 0, 3, 1, 1)));
        order_sizes.add(new ArrayList<>(Arrays.asList(3, 2, 3, 3, 2, 0, 2, 3, 4, 3, 5, 2, 2, 3, 0, 3, 3, 3, 4, 3, 2, 1, 0, 3, 4, 3, 0)));
        order_sizes.add(new ArrayList<>(Arrays.asList(3, 3, 4, 3, 0, 5, 3, 1, 3, 2, 3, 3, 3, 3, 3, 2, 4, 3, 3, 3, 3, 3, 4, 3, 5, 2, 0)));
        order_sizes.add(new ArrayList<>(Arrays.asList(3, 3, 4, 3, 3, 2, 2, 3, 3, 4, 2, 3, 2, 2, 4, 3, 3, 4, 3, 4, 1, 3, 3, 4, 3, 3, 3)));

        order_sizes.add(new ArrayList<>(Arrays.asList(0, 0, 0, 3, 2, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 3, 0, 4, 4, 5, 3, 0, 0, 0, 0, 3))); // 6
        
        order_deadlines.add(new ArrayList<>(Arrays.asList(0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)));
        order_deadlines.add(new ArrayList<>(Arrays.asList(4, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 2, 0)));
        order_deadlines.add(new ArrayList<>(Arrays.asList(0, 3, 2, 0, 0, 0, 0, 0, 2, 3, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)));
        order_deadlines.add(new ArrayList<>(Arrays.asList(0, 3, 0, 0, 0, 0, 4, 0, 0, 0, 3, 0, 3, 0, 0, 2, 0, 0, 2, 0, 0, 3, 0, 0, 0, 0, 0)));
        order_deadlines.add(new ArrayList<>(Arrays.asList(0, 0, 3, 2, 3, 0, 0, 0, 0, 0, 4, 3, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 3, 0, 0)));
        order_deadlines.add(new ArrayList<>(Arrays.asList(0, 0, 4, 0, 3, 0, 0, 2, 3, 0, 0, 0, 0, 3, 0, 2, 3, 0, 0, 0, 0, 3, 0, 0, 0, 3, 0)));
        order_deadlines.add(new ArrayList<>(Arrays.asList(0, 0, 3, 0, 2, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 3, 0, 3, 3, 2, 3, 0, 0, 0, 0, 3)));
        order_deadlines.add(new ArrayList<>(Arrays.asList(0, 0, 3, 0, 3, 2, 0, 0, 0, 3, 3, 0, 0, 0, 3, 0, 0, 2, 0, 3, 0, 0, 0, 3, 4, 0, 3)));
        order_deadlines.add(new ArrayList<>(Arrays.asList(3, 4, 3, 3, 0, 0, 0, 0, 3, 2, 0, 0, 0, 0, 0, 3, 0, 0, 0, 2, 3, 3, 4, 0, 0, 0, 0)));
        order_deadlines.add(new ArrayList<>(Arrays.asList(3, 0, 0, 3, 0, 0, 3, 0, 0, 3, 3, 0, 3, 3, 0, 3, 3, 0, 4, 3, 0, 0, 4, 0, 0, 0, 4)));
        order_deadlines.add(new ArrayList<>(Arrays.asList(0, 3, 3, 0, 0, 2, 3, 3, 2, 0, 0, 0, 0, 3, 2, 3, 4, 0, 0, 0, 0, 0, 4, 4, 4, 0, 0)));
        order_deadlines.add(new ArrayList<>(Arrays.asList(3, 3, 0, 3, 0, 3, 0, 4, 0, 3, 0, 3, 3, 3, 0, 0, 0, 4, 4, 4, 4, 0, 0, 3, 0, 4, 0)));
        order_deadlines.add(new ArrayList<>(Arrays.asList(0, 3, 3, 3, 4, 4, 4, 0, 0, 0, 4, 4, 3, 2, 0, 0, 0, 3, 3, 3, 0, 0, 4, 4, 0, 0, 0)));
        order_deadlines.add(new ArrayList<>(Arrays.asList(0, 0, 3, 0, 0, 4, 0, 3, 3, 3, 3, 0, 0, 3, 3, 0, 4, 4, 4, 4, 4, 4, 4, 0, 4, 4, 0)));
        order_deadlines.add(new ArrayList<>(Arrays.asList(3, 4, 3, 3, 3, 3, 0, 0, 0, 3, 4, 4, 4, 4, 0, 0, 0, 3, 4, 3, 0, 0, 0, 3, 4, 0, 0)));
        order_deadlines.add(new ArrayList<>(Arrays.asList(3, 3, 3, 0, 0, 3, 3, 0, 3, 3, 0, 3, 0, 0, 3, 3, 0, 3, 3, 3, 3, 3, 0, 3, 3, 3, 3)));
        order_deadlines.add(new ArrayList<>(Arrays.asList(4, 4, 4, 4, 4, 0, 3, 0, 3, 3, 0, 3, 3, 3, 3, 0, 3, 4, 0, 4, 4, 4, 3, 0, 3, 3, 3)));
        order_deadlines.add(new ArrayList<>(Arrays.asList(3, 3, 3, 3, 3, 0, 3, 3, 3, 3, 4, 4, 4, 4, 0, 4, 4, 4, 4, 4, 4, 4, 0, 3, 3, 3, 0)));
        order_deadlines.add(new ArrayList<>(Arrays.asList(3, 3, 3, 3, 0, 3, 3, 4, 4, 4, 4, 4, 4, 4, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 0)));
        order_deadlines.add(new ArrayList<>(Arrays.asList(4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4)));

        order_deadlines.add(new ArrayList<>(Arrays.asList(0, 0, 0, 3, 2, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 3, 0, 3, 3, 2, 3, 0, 0, 0, 0, 3)));

        vessel_return_days.add(new ArrayList<>(Arrays.asList(4, 0, 0, 0, 0, 0)));
        vessel_return_days.add(new ArrayList<>(Arrays.asList(3, 0, 0, 0, 0, 0)));
        vessel_return_days.add(new ArrayList<>(Arrays.asList(4, 0, 0, 0, 0, 0)));
        vessel_return_days.add(new ArrayList<>(Arrays.asList(4, 0, 0, 0, 0, 0)));
        vessel_return_days.add(new ArrayList<>(Arrays.asList(2, 0, 0, 0, 0, 0)));
        vessel_return_days.add(new ArrayList<>(Arrays.asList(3, 4, 0, 0, 0, 0)));
        vessel_return_days.add(new ArrayList<>(Arrays.asList(4, 3, 0, 0, 0, 0)));
        vessel_return_days.add(new ArrayList<>(Arrays.asList(4, 0, 0, 0, 0, 0)));
        vessel_return_days.add(new ArrayList<>(Arrays.asList(4, 3, 0, 0, 0, 0)));
        vessel_return_days.add(new ArrayList<>(Arrays.asList(3, 4, 0, 0, 0, 0)));
        vessel_return_days.add(new ArrayList<>(Arrays.asList(3, 3, 0, 0, 0, 0)));
        vessel_return_days.add(new ArrayList<>(Arrays.asList(4, 3, 3, 0, 0, 0)));
        vessel_return_days.add(new ArrayList<>(Arrays.asList(4, 4, 2, 0, 0, 0)));
        vessel_return_days.add(new ArrayList<>(Arrays.asList(4, 4, 4, 0, 0, 0)));
        vessel_return_days.add(new ArrayList<>(Arrays.asList(4, 3, 4, 0, 0, 0)));
        vessel_return_days.add(new ArrayList<>(Arrays.asList(4, 3, 3, 0, 0, 0)));
        vessel_return_days.add(new ArrayList<>(Arrays.asList(4, 3, 3, 0, 0, 0)));
        vessel_return_days.add(new ArrayList<>(Arrays.asList(3, 3, 2, 2, 0, 0)));
        vessel_return_days.add(new ArrayList<>(Arrays.asList(4, 4, 3, 2, 0, 0)));
        vessel_return_days.add(new ArrayList<>(Arrays.asList(4, 3, 3, 2, 3, 0)));

        vessel_return_days.add(new ArrayList<>(Arrays.asList(4, 3, 0, 0, 0, 0)));

    }

    public void add_scenario_data(){
        standard_order_size.put(1, 15.0);
        standard_order_size.put(2, 20.0);
        standard_order_size.put(3, 15.0);
        standard_order_size.put(4, 20.0);
        standard_order_size.put(5, 20.0);
        standard_order_size.put(6, 20.0);
        standard_order_size.put(7, 20.0);
        standard_order_size.put(8, 15.0);
        standard_order_size.put(9, 15.0);
        standard_order_size.put(10, 15.0);
        standard_order_size.put(11, 15.0);
        standard_order_size.put(12, 7.5);
        standard_order_size.put(13, 15.0);
        standard_order_size.put(14, 20.0);
        standard_order_size.put(15, 10.0);
        standard_order_size.put(16, 15.0);
        standard_order_size.put(17, 15.0);
        standard_order_size.put(18, 15.0);
        standard_order_size.put(19, 15.0);
        standard_order_size.put(20, 10.0);
        standard_order_size.put(21, 17.5);
        standard_order_size.put(22, 17.5);
        standard_order_size.put(23, 17.5);
        standard_order_size.put(24, 17.5);
        standard_order_size.put(25, 17.5);
        standard_order_size.put(26, 17.5);
        standard_order_size.put(27, 15.0);

        order_size_variations.put(1, 0.5);
        order_size_variations.put(2, 0.75);
        order_size_variations.put(3, 1.0);
        order_size_variations.put(4, 1.25);
        order_size_variations.put(5, 1.5);
    }

    public ArrayList<Order> get_orders_in_scenario(int scenario_number) {
        ArrayList<Integer> scenario_demand = order_sizes.get(scenario_number);
        ArrayList<Integer> scenario_deadlines = order_deadlines.get(scenario_number);
        ArrayList<Order> orders = new ArrayList<>();
        int counter = 1;
        for (int i = 0; i < scenario_demand.size(); i++ ){
            if (scenario_demand.get(i) != 0) {
                int order_demand = (int) Math.floor(standard_order_size.get(i + 1) * order_size_variations.get(scenario_demand.get(i)));
                int deadline_day = scenario_deadlines.get(i);
                Installation installation = this.installations.get(i + 1);
                orders.add(new Order(order_demand, 0, installation, counter, deadline_day));
                counter++;
            }
        }
        return orders;
    }

    public ArrayList<Vessel> get_vessels_in_scenario(int scenario_number) {
        ArrayList<Integer> scenario_vessels = vessel_return_days.get(scenario_number);
        ArrayList<Vessel> vessels = new ArrayList<>();
        int counter = 0;
        boolean spot_included = false;
        for (int i = 0; i < scenario_vessels.size() && !spot_included; i++) {
            if (scenario_vessels.get(i) != 0) {
                String name = "PSV_" + counter;
                int return_day = scenario_vessels.get(i);
                vessels.add(new Vessel(name, 100, 7, 14, 0, counter, return_day, false));
                counter++;
            }
            else {
                String name = "PSV_" + counter;
                int cost = 608;
                vessels.add(new Vessel(name, 100, 7, 14, cost, counter, 7, true));
                spot_included = true;
                counter++;
            }
        }
        return vessels;
    }
}
