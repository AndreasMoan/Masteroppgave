package main.java.HGSADCwSO.files;

import jdk.dynalink.beans.StaticClass;

import java.time.LocalDateTime;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;

public class Main
{
    public static void main(String[] args) throws ExecutionException
    {
        //Executor service instance
        ExecutorService executor = Executors.newFixedThreadPool(1);

        List<Callable<String>> tasksList = get_task_list();

        //1. execute tasks list using invokeAll() method
        try
        {
            List<Future<String>> results = executor.invokeAll(tasksList);

            for(Future<String> result : results) {
                System.out.println(result.get());
            }
        }
        catch (InterruptedException e1)
        {
            e1.printStackTrace();
        }

        //2. execute individual tasks using submit() method
        Future<String> result = executor.submit(callableTask);

        while(result.isDone() == false)
        {
            try
            {
                System.out.println("The method return value : " + result.get());
                break;
            }
            catch (InterruptedException | ExecutionException e)
            {
                e.printStackTrace();
            }

            //Sleep for 1 second
            try {
                Thread.sleep(1000L);
            } catch (InterruptedException e) {
                e.printStackTrace();
            }
        }

        //Shut down the executor service
        executor.shutdownNow();
    }


    private static List<Callable<String>> get_task_list() {
        List<Callable<Boolean> tasksList = new ArrayList<>();
        List<Integer> scenario_numbers = Arrays.asList(3, 6, 9, 12);
        List<String> education_rates = Arrays.asList("0", "0,1", "0,2", "0,5", "1");
        List<String> repair_rates = Arrays.asList("0,25", "0,5", "0,75");
        List<String> move_rates = Arrays.asList("0,2", "0,3", "0,5", "1");
        List<String> min_pop_sizes = Arrays.asList("25");
        List<String> generation_sizes = Arrays.asList("75");
        for (int scenario_number : scenario_numbers) {
            for (String education_rate : education_rates){
                for (String repair_rate : repair_rates) {
                    for (String move_rate : move_rates) {
                        for (String min_pop_size : min_pop_sizes) {
                            for (String generation_size : generation_sizes) {
                                Callable<Void> task = () -> {
                                    HGSADCwSOmain run = new HGSADCwSOmain(scenario_number, education_rate, repair_rate, move_rate, min_pop_size, generation_size);
                                    return run.fullEvolutionaryRun();
                                }
                            }
                        }
                    }
                }
            }
        }

    }
}