package spiderAlgorithm;

import java.util.Random;

import geneticAlgorithm.GeneticAlgorithm;

public class SpiderAlgorithm {
	
	public static void run(String[] args)
	{			
		/* create random generator and verify number of args*/
		long seed = System.currentTimeMillis();//1429463964577L;
		boolean updateSeed = GeneticAlgorithm.verifyNumberArgs(args);
		if (updateSeed){ seed = Long.parseLong(args[3]); }
		Random generator = new Random();
		generator.setSeed(seed);
		
		/* read dataset */
		String datasetFilename = args[0];//"datasets//iris.data";
		String columnsFilename = args[1];//"datasets//iris_cols.txt";
		int numberPoints = GeneticAlgorithm.countLinesFile(datasetFilename);
		int numberColumns = GeneticAlgorithm.countLinesFile(columnsFilename);
		int[] columns = new int[numberColumns];
		double[][] dataset = new double[numberPoints][];

		GeneticAlgorithm.readColumns(columnsFilename, columns);
		GeneticAlgorithm.readDataset(datasetFilename, dataset, columns, false);
		
		/* Set parameters for SSA */
		int numberClusters = Integer.parseInt(args[2]);//3; 
		int numberGenerations = 100;//Set other values
		int populationSize = 100;
		
		/* SSO Algorithm */
		Population pop = new Population(generator);
		pop.generateInitialPopulation(populationSize, numberClusters, dataset);
		pop.calculateWeightPopulation(dataset);
		
		for(int i=2; i<=numberGenerations; i++){
			pop.femaleCooperativeOperator();
			pop.maleCooperativeOperator(dataset);
			pop.matingOperator(dataset);
			pop.replacement();
			pop.calculateWeightPopulation(dataset);
		}
		
		//Show results
		boolean showAllResults = true;
		if (showAllResults){
			System.out.printf("Dataset number of points: %d\n\n", dataset.length);
			pop.showParameters(numberGenerations);
			pop.showBestSpider();        
			pop.showClusterBestSpider(dataset); 
		}
		else{ /*show just best metric*/
			pop.showMetricBestSpider();		
		}
	}

}
