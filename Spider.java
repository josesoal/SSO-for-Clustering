package spiderAlgorithm;

import java.util.Random;

public class Spider {
	private int 				type; //female=0, male=1 
	private ClusterCenter[] 	centers;
	private int[] 				datasetClusters;
	private double 				fitness;
	private double 				weight;
	private Random 				generator;//random generator of numbers
	
	public Spider(int pType, int length, double[][] dataset, Random pGenerator){
		type = pType;
		centers = new ClusterCenter[length];
		datasetClusters = new int[dataset.length];
		fitness = Integer.MAX_VALUE;
		weight = 0;
		generator = pGenerator;
	}
	
	public int getType(){
		return type;
	}
	
	public ClusterCenter[] getClusterCenters(){
		return centers;
	}
	
	public int[] getDatasetClusters(){
		return datasetClusters;
	}
	
	public double getFitness(){
		return fitness;
	}
	
	public void setFitness(double pFitness){
		fitness = pFitness;
	}
	
	public double getWeight(){
		return weight;
	}
	
	public void setWeight(double pWeight){
		weight = pWeight; 
	}
	
	/*Methods*/
	
	public void evaluateFitness(double[][] dataset){
		//build clusters according to its centers
		buildClusters(dataset);
		//compute new cluster centers
		computeNewClusterCenters(dataset);
		//calculate clustering metric
		fitness = calculateMetric(dataset);
	}
	
	public double calculateMetric(double[][] dataset){
		double distances = 0;
		for(int i=0; i<dataset.length; i++){
			int index = datasetClusters[i];
			distances += centers[index].calculateDistance(dataset[i]);
		}
		return distances;
	}
	
	private void buildClusters(double[][] dataset){
		double[] distances = null;
		int numClusters = centers.length;
		
        //build the dataset clusters		
		for(int i=0; i<dataset.length; i++){
			distances = new double[numClusters];
			for(int j=0; j<numClusters; j++){
				//euclidean distance from a dataset point "i"
				//to each cluster center
				distances[j] = centers[j].calculateDistance(dataset[i]);
			}
			
			//determine minimum distance
			double min = distances[0];
			int clusterIndex = 0;
			for(int j=1; j<numClusters; j++){
				if (distances[j] < min){
					min = distances[j];
					clusterIndex = j;
				}
			}
			
			//assign cluster index of the minimum distance
			//to the dataset point "i" 
			datasetClusters[i] = clusterIndex;
			
		}
		
		/*return datasetClusters;*/
	}
	
	private void computeNewClusterCenters(double[][] dataset){
		//*Random generator = new Random();
		//calculate total sum for each cluster
		int numClusters = centers.length;
		int[] numPointsForCluster = new int[numClusters];
		int pointDimension = dataset[0].length;
		double[][] sum = new double[numClusters][pointDimension];
		
		for(int i=0; i<dataset.length; i++){
			int index = datasetClusters[i];
			for(int j=0; j<pointDimension; j++){
				sum[index][j] = sum[index][j] + dataset[i][j];
			}
			numPointsForCluster[index]++;
		}
		
		//replace mean points as new cluster centers
		double newValue;
		for(int i=0; i<numClusters; i++){
			for(int j=0; j<pointDimension; j++){
				if (numPointsForCluster[i] == 0){
					//when all points of cluster have the same value
					//choose new value from random point at dataset
					int pos = generator.nextInt(dataset.length);
					newValue = dataset[pos][j];
				}
				else{
					//calculate mean point
					newValue = sum[i][j] / numPointsForCluster[i];
				}
				centers[i].getPoint()[j] = newValue;
			}
		}	
	}

	/** Calculate the dist. between "current" spider and spider "spi"
	 * Note: **/
	public double calculateDistance(Spider spi){
		double totalDistance = 0;
		for(int i=0; i<centers.length; i++){		
			totalDistance += 
				centers[i].calculateDistance(spi.centers[i].getPoint());
		}
		
		return totalDistance;
	}
	
	/** Difference of the cluster centers of two spiders **/
	public static void diffSpiders(Spider spi1, Spider spi2, ClusterCenter[] c){
		
		for (int i=0; i<spi1.centers.length; i++){
			c[i] = new ClusterCenter(spi1.centers[0].getPoint().length);
			
			for (int j=0; j<spi1.centers[0].getPoint().length; j++){
				c[i].getPoint()[j] = spi1.centers[i].getPoint()[j] - 
										spi2.centers[i].getPoint()[j];
			}
		}
		
	}
	
	/** Sum the cluster centers of a spider with other cluster centers **/
	public void sumSpider(ClusterCenter[] c, double signal){
		for (int i=0; i<centers.length; i++){
			for (int j=0; j<centers[0].getPoint().length; j++){
				centers[i].getPoint()[j] += signal * c[i].getPoint()[j]; 
				//System.out.printf("\n p %f\n", c[i].getPoint()[j]);	//--					 
			}
		}
	}
	
	/** Multiply cluster centers by a constant **/
	public static void mulClusterCentersByConstant(ClusterCenter[] c,  double cons){
		for(int i=0; i<c.length; i++){
			for(int j=0; j<c[0].getPoint().length; j++){
				c[i].getPoint()[j] *= cons;
			}
		}
	}
	
	
	public void showSpider(){
		System.out.printf("[");
		for(int j=0; j<centers.length; j++){
			System.out.printf("(");
			for(int k=0; k<centers[0].getPoint().length; k++){
				System.out.printf("%f; ",
					centers[j].getPoint()[k]);
			}
			System.out.printf("),");
		}
		System.out.printf("] Fitness: %f \n", getFitness());
	}
}
