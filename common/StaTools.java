/*
 * Terrier - Terabyte Retriever 
 * Webpage: http://ir.dcs.gla.ac.uk/terrier 
 * Contact: terrier{a.}dcs.gla.ac.uk
 * University of Glasgow - Department of Computing Science
 * http://www.gla.ac.uk/
 * 
 * The contents of this file are subject to the Mozilla Public License
 * Version 1.1 (the "License"); you may not use this file except in
 * compliance with the License. You may obtain a copy of the License at
 * http://www.mozilla.org/MPL/
 *
 * Software distributed under the License is distributed on an "AS IS"
 * basis, WITHOUT WARRANTY OF ANY KIND, either express or implied. See
 * the License for the specific language governing rights and limitations
 * under the License.
 *
 * The Original Code is StaTools.java.
 *
 * The Original Code is Copyright (C) 2004-2009 the University of Glasgow.
 * All Rights Reserved.
 *
 * Contributor(s):
 *   Ben He <ben{a.}dcs.gla.ac.uk>
 */
package common;

import java.util.*;
/**
 * This class implements a series of basic statistical functions.
 */
public class StaTools {
	/**
	 * This method provides the contract for implementing the Stirling formula for the power series.
	 * @param n The parameter of the Stirling formula.
	 * @param m The parameter of the Stirling formula.
	 * @return the approximation of the power series
	 */
	//    public static double stirlingPower(double n, double m) {
	//        double dif = n - m;
	//        return (m + 0.5d) * Idf.log(n / m) + dif * Idf.log(n);
	//    }

	/**
	 * This method returns the standard error of the mean for an array of data.
	 * @param data The sampled data.
	 * @return The standard error of the mean.
	 */
	public static double stdErrorOfTheMean(double[] data){
		return standardDeviation(data) / Math.sqrt(data.length);
	}


	/**
	 * The sum of an array of integers.
	 * @param data The integers.
	 * @return The sum.
	 */
	public static int sum(int[] data){
		int sum = 0;
		for (int i = 0; i < data.length; i++)
			sum+=data[i];
		return sum;
	}
	public static double sum(double[] data){
		int sum = 0;
		for (int i = 0; i < data.length; i++)
			sum+=data[i];
		return sum;
	}   
	public static HashMap<String, Double> uniformarray(HashMap<String, Double> topTermWeightVecotr) {
		double temp = 0;
		HashMap<String, Double> newVector = new HashMap<String, Double>();
		for (String term : topTermWeightVecotr.keySet()) {
			temp += topTermWeightVecotr.get(term);
		}

		for (String term : topTermWeightVecotr.keySet()) {
			if(temp!=0){
				newVector.put(term, topTermWeightVecotr.get(term)/temp);
			}else{
				newVector.put(term, 0d);
			}
		}
		return newVector;
	}    

	public static HashMap<String, Double> softmax(HashMap<String, Double> topTermWeightVecotr) {
		double temp = 0;
		HashMap<String, Double> newVector = new HashMap<String, Double>();
		for (String term : topTermWeightVecotr.keySet()) {
			temp += Math.exp(topTermWeightVecotr.get(term));
		}

		for (String term : topTermWeightVecotr.keySet()) {
			if(temp!=0){
				newVector.put(term, Math.exp(topTermWeightVecotr.get(term))/temp);
			}else{
				newVector.put(term, 0d);
			}
		}
		return newVector;
	}  

	public static HashMap<String, Double> Divide(HashMap<String, Double> topTermWeightVecotr, double temp) {
		HashMap<String, Double> newVector = new HashMap<String, Double>();

		for (String term : topTermWeightVecotr.keySet()) {
			newVector.put(term, topTermWeightVecotr.get(term)/temp);
		}
		return newVector;
	} 

	public static HashMap<String, Double> MAXNormalizedarray(HashMap<String, Double> topTermWeightVecotr) {
		double maxV = Double.MIN_VALUE;
		HashMap<String, Double> newVector = new HashMap<String, Double>();
		for (String term : topTermWeightVecotr.keySet()) {
			if (maxV<topTermWeightVecotr.get(term)){
				maxV = topTermWeightVecotr.get(term);
			}
		}

		for (String term : topTermWeightVecotr.keySet()) {
			if(maxV!=0){
				newVector.put(term, topTermWeightVecotr.get(term)/maxV);
			}else{
				newVector.put(term, 0d);
			}
		}
		return newVector;
	} 

	public static HashMap<String, Double> Normalizedarray(HashMap<String, Double> topTermWeightVecotr) {
		double maxV = Double.MIN_VALUE;
		double minV = Double.MAX_VALUE;
		HashMap<String, Double> newVector = new HashMap<String, Double>();
		for (String term : topTermWeightVecotr.keySet()) {
			if (maxV<topTermWeightVecotr.get(term)){
				maxV = topTermWeightVecotr.get(term);
			}
			if (minV>topTermWeightVecotr.get(term)){
				minV = topTermWeightVecotr.get(term);
			}
		}

		for (String term : topTermWeightVecotr.keySet()) {
			if(maxV!=minV){
				newVector.put(term, (topTermWeightVecotr.get(term)-minV)/(maxV-minV));
			}else{
				newVector.put(term, 0d);
			}
		}
		return newVector;
	} 

	public static double[] uniformarray(double[] data) {
		double temp = sum(data);
		if(temp!=0){
			for (int i=0; i<data.length; i++)
				data[i]=data[i]/temp;
		}
		return data;
	}    



	public static double[] softmax(double[] data) {
		double temp = 0d;
		for (int i=0; i<data.length; i++){
			temp += Math.exp(data[i]);
		}
		if(temp!=0){
			for (int i=0; i<data.length; i++)
				data[i] = Math.exp(data[i])/temp;
		}
		return data;
	}  


	/**
	 * The mean of an array of double values.
	 * @param data The double values.
	 * @return The mean.
	 */
	public static double mean(double[] data) {
		double mean = 0d;
		for (int i=0; i<data.length; i++)
			mean+=data[i];
		mean/=data.length;
		return mean;
	}

	/**
	 * The mean of a sub-array of an array of double values.
	 * @param data The array of double values.
	 * @param start The starting index of the sub-array.
	 * @param length The length of the sub-array.
	 * @param ascending Is the starting index the left (true) or 
	 * right (false) end of the sub-array?
	 * @return The mean of the sub-array.
	 */
	public static double mean(double[] data, int start, int length, boolean ascending) {
		double mean = 0d;
		if (ascending)
			for (int i = start; i < length; i++)
				mean += data[i];
		else
			for (int i = 0; i < length; i++)
				mean += data[start - i];
		mean /= length;
		return mean;
	}

	/**
	 * The mean of an array of integers.
	 * @param data The array of integers.
	 * @return The mean.
	 */
	public static double mean(int[] data) {
		double mean = 0d;
		for (int i=0; i<data.length; i++)
			mean+=data[i];
		mean/=data.length;
		return mean;
	}

	/**
	 * The median of an array of double values.
	 * @param data The array of double values.
	 * @return The median.
	 */
	public static double median(double[] data) {
		double[] copy = (double[])data.clone();
		Arrays.sort(copy);
		return data[(copy.length-1)/2];
	}
	/**
	 * The standard deviation of an array of double values.
	 * @param data The array of double values.
	 * @return The standrad deviation.
	 */
	public static double standardDeviation(double[] data) {	
		return Math.sqrt(variance(data));
	}
	/**
	 * The variance of an array of double values. 
	 * @param data The array of double values.
	 * @return The variance.
	 */
	public static double variance(double[] data) {
		double var = 0d;
		int n = data.length;
		double mean =mean(data);
		for (int i=0; i<n; i++)
			var+=(data[i]-mean)*(data[i]-mean);
		var /= n;

		return var;
	}
	public static double max(double[] data) {
		double amax=0;
		if (data.length>0){
			amax = data[0];
			int ldata = data.length;
			for(int i = 0; i < ldata; i++){
				amax = Math.max(amax,data[i]);
			}
		}
		return amax;
	}
	public static double min(double[] data) {
		double amin;
		amin=10000;
		if (data.length>0){
			amin = data[0];
			int ldata = data.length;
			for(int i = 0; i < ldata; i++){
				amin = Math.min(amin,data[i]);
			}
		}
		return amin;
	}

	public static double[] normalizedarray(double[] data){
		int data_length = data.length;
		double[] normalizedarray = new double[data_length];
		double dmin = min(data);
		double dmax = max(data);

		for (int i=0; i<data_length; i++){
			if (dmin == dmax || Double.isNaN(dmin) || Double.isNaN(dmax)){
				normalizedarray[i]=0;
			}else{
				normalizedarray[i] = (data[i]-dmin)/(dmax-dmin); 
			}
		}

		return normalizedarray;
	}
	/**
	 * The add of two arrays of double values with same length. 
	 * @param data1,data2 The two arrays of double values
	 * @return The mean.
	 */
	public static double[] add(double[] data1, double[] data2 ) {
		double []add=new double[data1.length];
		for (int i=0; i<data1.length; i++)
			add[i]=data1[i]+data2[i];
		return add;
	}
	public static double[] add(double[] data1, double[] data2, double alpha, double beta) {
		double []add=new double[data1.length];
		for (int i=0; i<data1.length; i++)
			add[i]=alpha*data1[i]+beta*data2[i];
		return add;
	}

}
