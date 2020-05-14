package spiderAlgorithm;

import java.util.Locale;
import java.util.Random;
import java.util.List;
import java.util.ArrayList;

public class Population {
	private int 				numberFemales;
	private int 				numberMales;
	private int 				numClusters;
	private int 				pointDimension;
	private double 				radiusMating;
	private Spider[]			spiders;
	private List<Spider>		offspring;
	private double 				medianWeight;
	private int 				indexBest;//index of the best result
	private int 				indexWorst;
	private Random 				generator;//random generator of numbers
	
	public Population(Random pGenerator){
		numberFemales = 0;
		numberMales = 0;
		numClusters = 0;
		pointDimension = 0;
		radiusMating = 0;
		spiders = null;
		offspring = null;
		medianWeight = 0;
		indexBest = 0;//Let the first be the best
		indexWorst = 0;//Let the first be the worst
		generator = pGenerator;
	}
	
	/** Generate initial population: males and females, and calculate radius **/
	public void generateInitialPopulation(int populationSize, 
			int pNumberClusters, double[][] dataset){
		
		//*Random generator = new Random();
		numberFemales = 
				(int)Math.floor((0.9 - generator.nextDouble()*0.25)*populationSize);
		numberMales = populationSize - numberFemales;
		numClusters = pNumberClusters;
		pointDimension = dataset[0].length;
		int pos = -1;//invalid position
		spiders = new Spider[populationSize];
		
		//generate FEMALES by choosing points randomly from dataset
		for(int i=0; i<numberFemales; i++){
			spiders[i] = new Spider(0, numClusters, dataset, generator);//female=0
			for(int j=0; j<numClusters; j++){
				//random integer in range [0 , dataset.length-1]
				pos = generator.nextInt(dataset.length);
				
				spiders[i].getClusterCenters()[j] = new ClusterCenter(pointDimension);
				for(int k=0; k<pointDimension; k++){
					spiders[i].getClusterCenters()[j].getPoint()[k] = dataset[pos][k];	
				}	
			}
		}
		
		//generate MALES by choosing points randomly from dataset
		for(int i=numberFemales; i<populationSize; i++){
			spiders[i] = new Spider(1, numClusters, dataset, generator);//male=1 
			for(int j=0; j<numClusters; j++){
				//random integer in range [0 , dataset.length-1]
				pos = generator.nextInt(dataset.length);
				
				spiders[i].getClusterCenters()[j] = new ClusterCenter(pointDimension);
				for(int k=0; k<pointDimension; k++){
					spiders[i].getClusterCenters()[j].getPoint()[k] = dataset[pos][k];	
				}	
			}
		}
		
		//calculate the radius of mating
		double d = maximumDistance(dataset);
		radiusMating = d / 2; 
		
	}
	
	private double maximumDistance(double[][] dataset){
		double maximumDistance = 0;
		double distance = 0;
		for(int i=0; i<dataset.length-1; i++){
			for(int j=i+1; j<dataset.length; j++){
				distance = ClusterCenter.calculateDistance(i, j, dataset);
				if (distance > maximumDistance)
					maximumDistance = distance;
			}
		}
		
		return maximumDistance;
	}
	
	/** Calculate the weight of every spider in the population **/
	public void calculateWeightPopulation(double[][] dataset){
		//Evaluate the fitness of the population 
		evaluateFitnessPopulation(dataset);
		//Calculate the best and the worst individual of pop.
		calculateBestWorstPopulation();
		//Calculate weight of population
		for(int i=0; i<spiders.length; i++){
			double w = (spiders[indexWorst].getFitness() - spiders[i].getFitness()) 
				/ (spiders[indexWorst].getFitness() - spiders[indexBest].getFitness());
			
			spiders[i].setWeight(w);
		}	
	}
	
	private void evaluateFitnessPopulation(double[][] dataset){
		for(int i=0; i<spiders.length; i++){
			//keeping the best spider the same
			if (i == indexBest){
				continue;
			}
			
			spiders[i].evaluateFitness(dataset);
		}
	}
	
	private void calculateBestWorstPopulation(){
		double minFitness = spiders[0].getFitness();
		double maxFitness = spiders[0].getFitness();
		
		for(int i=0; i<spiders.length; i++){
			if (spiders[i].getFitness() < minFitness){
				minFitness = spiders[i].getFitness();
				indexBest = i;
			}
			
			if (spiders[i].getFitness() > maxFitness){
				maxFitness = spiders[i].getFitness();
				indexWorst = i;
			}
		}
	}
	
	
	/** Apply the female cooperative operator **/
	public void femaleCooperativeOperator(){
		//*Random generator = new Random();
		double thresholdPF = 0.5;// Test with other values
		double w = 0, d = 0, d2 = 0;
		double vibci = 0;
		double vibbi = 0;
		double rm = 0, alpha = 0, beta = 0, delta = 0, rand = 0;
		
		for(int i=0; i<numberFemales; i++){
			//keeping the best spider the same
			if (i == indexBest){
				continue;
			}
			
			//Calculate Vibci: Vib. of ind. nearest with higher weight compared to i
			int index = nearestWithHigherWeightTo(i);
			w = spiders[index].getWeight();
			d = spiders[i].calculateDistance(spiders[index]);
			d2 = Math.pow(d, 2);
			vibci = w / Math.pow(Math.E, d2);
			
			//Calculate Vibbi: vibration of individual with best fitness
			w = spiders[indexBest].getWeight(); 
			d = spiders[i].calculateDistance(spiders[indexBest]);
			d2 = Math.pow(d, 2);
			vibbi = w / Math.pow(Math.E, d2);
			
			rm = generator.nextDouble();
			alpha = generator.nextDouble();
			beta = generator.nextDouble();
			delta = generator.nextDouble();
			rand = generator.nextDouble();
			ClusterCenter[] clusCenters;
			double cons, signal;
			
			//Define if the movement is: attraction or repulsion
			if (rm < thresholdPF){
				signal = 1; // --> 1=sum
			}
			else{
				signal = -1; // --> -1=subtraction
			}
			
			//Sum expression with alpha
			clusCenters = new ClusterCenter[numClusters];
			Spider.diffSpiders(spiders[index], spiders[i], clusCenters);
			cons = alpha * vibci;
			Spider.mulClusterCentersByConstant(clusCenters, cons);
			spiders[i].sumSpider(clusCenters, signal);
		
			//Sum expression with beta
			clusCenters = new ClusterCenter[numClusters]; 
			Spider.diffSpiders(spiders[indexBest], spiders[i], clusCenters);
			cons = beta * vibbi;
			Spider.mulClusterCentersByConstant(clusCenters, cons);
			spiders[i].sumSpider(clusCenters, signal); 

			//Sum expression with delta
			clusCenters = new ClusterCenter[numClusters];
			for(int j=0; j<numClusters; j++){
				clusCenters[j] = new ClusterCenter(pointDimension);
				for(int k=0; k<pointDimension; k++){
					clusCenters[j].getPoint()[k] = delta * (rand - 0.5);
				}
			}
			spiders[i].sumSpider(clusCenters, signal);	
		}
	}
	
	/** Return the index of the nearest individual with higher weight
	 * compared to individual with index i **/
	private int nearestWithHigherWeightTo(int i){
		int index = 0;
		double nearestDistance = spiders[index].calculateDistance(spiders[i]);
		
		for(int k=0; k<spiders.length; k++){
			if (spiders[k].getWeight() > spiders[i].getWeight()){
				double newDistance = spiders[k].calculateDistance(spiders[i]);
				if (newDistance < nearestDistance){
					nearestDistance = newDistance;
					index = k;
				}
			}
		}
		return index;
	}
	
	/** Apply the male cooperative operator **/
	public void maleCooperativeOperator(double[][] dataset){
		//Calculate the median weight of male population
		calculateMedianWeightOfMales(); 
		//Calculate the male spider with weighted mean
		Spider spiderWeightedMean = new Spider(1, numClusters, dataset, generator);
		calculateMaleSpiderWeightedMean(spiderWeightedMean);
		
		//*Random generator = new Random();
		double w = 0, d = 0, d2 = 0;
		double vibfi = 0;
		double alpha = 0, delta = 0, rand = 0, cons = 0;
		
		for(int i=numberFemales; i<spiders.length; i++){//males
			//keeping the best spider the same
			if (i == indexBest){
				continue;
			}
			
			//calculate vibfi: vibration of nearest female
			int index = nearestFemaleTo(i);
			w = spiders[index].getWeight();
			d = spiders[i].calculateDistance(spiders[index]);
			d2 = Math.pow(d, 2);
			vibfi = w / Math.pow(Math.E, d2);
			
			//Define if movement is attraction to females or to the mean
			ClusterCenter[] clusCenters;
			alpha = generator.nextDouble();
			delta = generator.nextDouble();
			rand = generator.nextDouble();
			if (spiders[i].getWeight() > medianWeight){//male is dominant(D)
				//Sum expression with alpha
				clusCenters = new ClusterCenter[numClusters];
				Spider.diffSpiders(spiders[index], spiders[i], clusCenters);
				cons = alpha * vibfi;
				Spider.mulClusterCentersByConstant(clusCenters, cons);
				spiders[i].sumSpider(clusCenters, 1);//sum = 1
				
				//Sum expression with delta
				clusCenters = new ClusterCenter[numClusters];
				for(int j=0; j<numClusters; j++){
					clusCenters[j] = new ClusterCenter(pointDimension);
					for(int k=0; k<pointDimension; k++){
						clusCenters[j].getPoint()[k] = delta * (rand - 0.5);
					}
				}
				spiders[i].sumSpider(clusCenters, 1);//sum = 1
				
			}
			else{//male is not dominant(ND)
				//Sum expression with alpha
				clusCenters = new ClusterCenter[numClusters];
				Spider.diffSpiders(spiderWeightedMean, spiders[i], clusCenters);
				Spider.mulClusterCentersByConstant(clusCenters, alpha);
				spiders[i].sumSpider(clusCenters, 1);//sum = 1
			}
		}
	}
	
	private int nearestFemaleTo(int i){
		int index = 0;
		double nearestDistance = spiders[index].calculateDistance(spiders[i]);
		
		for(int k=0; k<numberFemales; k++){//females
			double newDistance = spiders[k].calculateDistance(spiders[i]);
			if (newDistance < nearestDistance){
				nearestDistance = newDistance;
				index = k;
			}
		}
		return index;
	}
	
	/* calculate median weight of males using variant of quicksort */
	private void calculateMedianWeightOfMales(){
		//Copy the weight of males into males array
		double[] males = new double[numberMales];
		int k=0;
		for(int i=numberFemales; i<spiders.length; i++){
			males[k] = spiders[i].getWeight();
			k++;
		}
		
		//Calculate the median position
		int medianPos = (0 + numberMales-1) / 2;
		
		//Calculate the median value		
		int begin, end, p, r, i;
		double pivot;
		
		begin = 0;
		end = numberMales-1;
		
		while(true){
			p = begin;
			r = end;
			pivot = males[r];
			i = p-1;
			for(int j=p; j<r; j++){
				if (males[j] < pivot){
					i++;
					//swap males[i] and males[j]
					double tmp = males[i];
					males[i] = males[j];
					males[j] = tmp;
				}
			}
			//swap males[i+1] and males[r]
			double tmp = males[i+1];
			males[i+1] = males[r];
			males[r] = tmp;
			
			if (i+1 == medianPos){
				break;
			}
			else{
				if (medianPos > i+1){
					//search right
					begin = i+2;
				}
				else{
					//search left
					end = i;
				}			
			}
		}
		
		medianWeight = males[medianPos];
	}
	
	private void calculateMaleSpiderWeightedMean(Spider spi){
		//calculate total weight of males
		double totalWeight = 0;
		for(int i=numberFemales; i<spiders.length; i++){
			totalWeight += spiders[i].getWeight();
		}
		
		//calculate spiders multiplied by their weights
		ClusterCenter[] spidersWeights = new ClusterCenter[numClusters];
		for(int i=0; i<numClusters; i++){
			spidersWeights[i] = new ClusterCenter(pointDimension);
			for(int j=0; j<pointDimension; j++){
				spidersWeights[i].getPoint()[j] = 0;
			}
		}
		
		for(int i=numberFemales; i<spiders.length; i++){
			for(int j=0; j<numClusters; j++){
				for(int k=0; k<pointDimension; k++){
					spidersWeights[j].getPoint()[k] += 
							spiders[i].getClusterCenters()[j].getPoint()[k] * spiders[i].getWeight();
				}
			}
		}
		
		//calculate the weighted mean
		for(int j=0; j<numClusters; j++){
			spi.getClusterCenters()[j] = new ClusterCenter(pointDimension);
			for(int k=0; k<pointDimension; k++){
				spi.getClusterCenters()[j].getPoint()[k] = 
							spidersWeights[j].getPoint()[k] / totalWeight;
				
			}
		}
	}
	
	/** Apply the mating operator **/
	public void matingOperator(double[][] dataset)
	{
		offspring = new ArrayList<Spider>();//create a new offspring
		//Begin mating
		for(int i=numberFemales; i<spiders.length; i++){//males
			
			if (spiders[i].getWeight() > medianWeight){//male is dominant 
				//Calculate females in the radius of male "i"
				List<Integer> matingGroup = new ArrayList<Integer>();//indexes  
				matingGroup.add(i);//add male index as first element
				for(int j=0; j<numberFemales; j++){//females
					double distance = spiders[i].calculateDistance(spiders[j]);
					
					if (distance < radiusMating){
						matingGroup.add(j);//add female index
					}
				}

				if (matingGroup.size() > 1){//do mating

					//Create mating roulette
					double[] matingRoulette = new double[matingGroup.size()];//females +  1 male
					createMatingRoulette(matingRoulette, matingGroup);
					//Create the new spider using mating roulette
					Spider spi = new Spider(0, numClusters, dataset, generator);
					//*Random generator = new Random();
					for(int j=0; j<numClusters; j++){
						double rand = generator.nextDouble();// 0.0 <= rand < 1.0
						//Go through the mating roulette
						for(int k=0; k<matingRoulette.length; k++){
							if (rand < matingRoulette[k]){
								//Copy cluster "j"
								spi.getClusterCenters()[j] = new ClusterCenter(pointDimension);
								for(int h=0; h<pointDimension; h++){
									spi.getClusterCenters()[j].getPoint()[h] = 
											spiders[matingGroup.get(k)].getClusterCenters()[j].getPoint()[h];	
								}
								break;
							}
						}
					}
					//Calculate fitness of new spider and put into offspring
					spi.evaluateFitness(dataset); 
					offspring.add(spi);
				}
			}
		}//end-mating
		
	}
	
	/** Create Mating Roulette for a male spider **/
	private void createMatingRoulette(double[] matingRoulette, List<Integer> matingGroup){
		
		//Sum fitness of mating spiders
		double total =  0;
		for(int i=0; i<matingGroup.size(); i++){
			total += spiders[matingGroup.get(i)].getFitness();
		}
		//Calculate values of the roulette
		matingRoulette[0] = spiders[matingGroup.get(0)].getFitness() / total;
		for(int i=1; i<matingGroup.size(); i++){
			matingRoulette[i] = matingRoulette[i-1] + 
					spiders[matingGroup.get(i)].getFitness() / total;
			
		}
		
	}
	
	/** Replace offspring into spiders **/
	public void replacement()
	{
		//Create replacement roulette for all spiders, giving more prob. to worst spi.
		double[] replacementRoulette = new double[spiders.length];
		createReplacementRoulette(replacementRoulette);
		//Replace worst spider by offspring by comparing its fitness
		//*Random generator = new Random();
		for(int i=0; i<offspring.size(); i++){
			double rand = generator.nextDouble();// 0.0 <= rand < 1.0
			//Go through the replacement roulette
			for(int j=0; j<replacementRoulette.length; j++){
				if (rand < replacementRoulette[j]){
					//replace spider "j" if it is worst than offspring "i"
					if (offspring.get(i).getFitness() < spiders[j].getFitness()){
						for(int k=0; k<numClusters; k++){
							for(int h=0; h<pointDimension; h++){
								spiders[j].getClusterCenters()[k].getPoint()[h] = 
										offspring.get(i).getClusterCenters()[k].getPoint()[h];
								
							}
						}
						spiders[j].setFitness(offspring.get(i).getFitness());
						spiders[j].setWeight(0);
						for(int k=0; k<offspring.get(i).getDatasetClusters().length; k++){
							spiders[j].getDatasetClusters()[k] = 
									offspring.get(i).getDatasetClusters()[k];
						}
					}
					
					break;
				}
			}
		}
		
	}
	
	//using a factor
	private void createReplacementRoulette(double[] replacementRoulette){
		//Sum fitness of all spiders
		double factor = 1;  
		double total=0;
		for(int i=0; i<spiders.length;i++){
			if (spiders[i].getWeight() > medianWeight)
				total += spiders[i].getFitness();
			else
				total += spiders[i].getFitness() * factor;
		}
		
		//calculate values of the roulette
		replacementRoulette[0] = spiders[0].getFitness() / total;
		for(int i=1; i<spiders.length; i++){
			if (spiders[i].getWeight() > medianWeight)
				replacementRoulette[i] = replacementRoulette[i-1] + 
						spiders[i].getFitness() / total;
			else
				replacementRoulette[i] = replacementRoulette[i-1] + 
						(spiders[i].getFitness() * factor) / total;
		}
		
	}
	
	
	/** Show best fitness **/
	public void showBestFitness(int numberGeneration){
		System.out.printf("Generation: %d, Best Fitness: %f, Worst Fitness: %f \n",
				numberGeneration, spiders[indexBest].getFitness(), spiders[indexWorst].getFitness());
	}
	
	/** Show Population **/
	public void showPopulation(){
		System.out.printf("\n\n ----> POPULATION\n");
		for(int i=0; i<spiders.length; i++){
			System.out.printf("%d: [",i);
			for(int j=0; j<numClusters; j++){
				System.out.printf("(");
				for(int k=0; k<pointDimension; k++){
					//System.out.printf("%f; ",
					//		spiders[i].getClusterCenters()[j].getPoint()[k]);
				}
				System.out.printf("),");
			}
			System.out.printf("] Fitness: %f \n", spiders[i].getFitness());
		}
		System.out.printf("*Pop.\n");
	}
	
	/** Show some parameters **/
	public void showParameters(int numberGenerations){
		System.out.printf("----------------------------\n");
		System.out.printf("----S. Spider Algorithm-----\n");
		System.out.printf("----------------------------\n");
		System.out.printf("Number of Generations: \t%d\n", numberGenerations);
		System.out.printf("Population Size: \t%d\n", spiders.length);
		System.out.printf("Number of females: \t%d\n", numberFemales);
		System.out.printf("Number of males: \t%d\n", numberMales);
	}
	
	/** Show best spider **/
	public void showBestSpider(){
		System.out.printf("----------------------------\n");
		System.out.printf("Best Spider: \n");
		System.out.printf("[\n");
		for(int j=0; j<numClusters; j++){
			System.out.printf("(");
			for(int k=0; k<pointDimension; k++){
				System.out.printf("%f; ",
						spiders[indexBest].getClusterCenters()[j].getPoint()[k]);
			}
			System.out.printf("),\n");
		}
		System.out.printf("] \nFitness: %f \n", spiders[indexBest].getFitness());
		
	}
	
	/** Show cluster generated by the best spider **/
	public void showClusterBestSpider(double[][] dataset){
		System.out.printf("----------------------------\n");
		System.out.printf("Cluster Generated by Best Spider: \n");
		for(int i=0; i<numClusters; i++){
			System.out.printf("Cluster %d: ", i);
			for(int j=0; j<spiders[indexBest].getDatasetClusters().length; j++){
				if (i == spiders[indexBest].getDatasetClusters()[j]){
					System.out.printf("%d, ", j);
				}
			}
			System.out.printf("\n -->Center %d: (",i);
			for(int j=0; j<pointDimension; j++){
				System.out.printf("%f; ",
						spiders[indexBest].getClusterCenters()[i].getPoint()[j]);
			}
			System.out.printf(") \n\n");
		}
	}
	
	public void showMetricBestSpider(){
		System.out.printf(Locale.US,"%f",spiders[indexBest].getFitness()); /* US decimal separator is "." */ 
	}
	
}
