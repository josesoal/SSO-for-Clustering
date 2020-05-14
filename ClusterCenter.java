package spiderAlgorithm;

public class ClusterCenter {
	private double[] point;//N-dimensional point
	
	public ClusterCenter(int length){
		point = new double[length];//N = length
	}
	
	public double[] getPoint(){
		return point;
	}
	
	/*Methods*/
	public double calculateDistance(double[] pPoint){
		int pointDimension = point.length;
		double totalSum = 0;
		for(int i=0; i<pointDimension; i++){
			totalSum += (point[i]-pPoint[i])*(point[i]-pPoint[i]);
		}
		
		return (double)Math.sqrt(totalSum);
	}
	
	public static double calculateDistance(int i, int j, double[][] dataset){
		int pointDimension = dataset[0].length;//  point.length;
		double totalSum = 0;
		for(int k=0; k<pointDimension; k++){
			totalSum += (dataset[i][k]-dataset[j][k])*(dataset[i][k]-dataset[j][k]);
		}
		
		return (double)Math.sqrt(totalSum);
	}

}
