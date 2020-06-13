package main.java.HGSADCwSO.files;

import jdk.dynalink.beans.StaticClass;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
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

class Main
{
    /*
    public static void main(String[] args) throws ExecutionException, IOException {
        //Executor service instance
        HGSADCwSOmain run = new HGSADCwSOmain(1, 19);
        System.out.println(run.fullEvolutionaryRun());
        System.out.println("Done");
    }
    */

    public static void main(String[] args) throws ExecutionException
    {
        //Executor service instance
        ExecutorService executor = Executors.newFixedThreadPool(1);

        List<Callable<String>> tasksList = get_task_list();

        //1. execute tasks list using invokeAll() method
        try {
            List<Future<String>> results = executor.invokeAll(tasksList);

            File file;
            String filename = "datastack.txt";
            try {
                file = new File(filename);
                if (file.createNewFile()) {
                    System.out.println("File created: " + file.getName());
                } else {
                    System.out.println("File already exists.");
                }
                FileWriter fileWriter = new FileWriter(filename);
                for (Future<String> result : results) {
                    System.out.println(result.get());
                    fileWriter.write(result.get());
                    fileWriter.flush();
                }
                fileWriter.close();
            } catch (IOException e) {
                System.out.println("An error occurred.");
                e.printStackTrace();
            }

        }
        catch (InterruptedException e1) {
            e1.printStackTrace();
        }

        //2. execute individual tasks using submit() method

        //Shut down the executor service
        executor.shutdownNow();
        System.out.println("Done");
    }


    private static List<Callable<String>> get_task_list() {
        List<Callable<String>> tasks_list = new ArrayList<>();
        List<Integer> scenario_numbers = Arrays.asList(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19);
        List<String> discos = Arrays.asList("4");//"15", "25", "35", "50");

        int run_id = 0;
        for (int scenario_number : scenario_numbers) {
            for (String disco : discos) {
                int id = run_id;
                Callable<String> task = () -> {
                    String s = "";
                    try {
                        HGSADCwSOmain run = new HGSADCwSOmain(id, scenario_number, disco);
                        s += run.fullEvolutionaryRun();
                    } catch (Exception e) {
                        e.printStackTrace();
                        s += id + "|error|error";
                    }
                    return s + "| " + scenario_number + "|" + disco + "\n";
                };
                tasks_list.add(task);
                run_id++;
            }
        }
        return tasks_list;
    }
}