package srp.dbgraph;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Enumeration;
import java.util.HashMap;
import java.util.HashSet;


import srp.haplotypes.Haplotype;
import srp.haplotypes.HaplotypeModel;

public class DeBruijnGraph {
//	private OpenIntIntHashMap
	private HashMap<Integer, String> allNodes;
	private HashMap<Integer, Integer> allLength;
	private HashMap<Integer, Integer> allEdges;
	private int maxNode;
	private ArrayList<CompatibleSet> allSets;
	
	public DeBruijnGraph() {
		allNodes = new HashMap<Integer, String>();
		allEdges = new HashMap<Integer, Integer>();
		maxNode = 0;
	}
	
	public void addNode(int node, String seqs) {
		allNodes.put(node, seqs);
		if(node> maxNode){
			maxNode = node;
		}
		allLength.put(node, seqs.length());
	}

	public void addEdge(int node, int node2) {
		allEdges.put(node, node2);
		
	}
	
	public void postPrecossDeBruijnGraph(){
		
		
	}
	
	public void test(){
		for (int i = 0; i < maxNode; i++) {
			System.out.println(i +"\t"+ allEdges.get(i) +"\t"+ allNodes.get(i));
			
		}
	}
	public double calculateLikelihood(HaplotypeModel haplotypeModel){
		
		
		double likelihood = 0;
		
//		#print "In compute set \n";
//		my($llik) = 0; my($min) = 0; my($max) = -9**9**9; 
		double min = 0;
		double Scale = haplotypeModel.getHaplotypeCount() * haplotypeModel.getHaplotypeLength();
		
		
		double[][] d_hashTable = computeDHashTable(haplotypeModel);
		double[][] length_diff = new double[maxNode][maxNode];
		int totalCount = 0;
		for (CompatibleSet k1CSet : allSets) {
			int k1 = k1CSet.getNode();
			ArrayList<Integer> cNodeList = k1CSet.getCNodeList();
			for (Integer k2 : cNodeList) {
				totalCount += k1CSet.getCNodeCount(k2);

				length_diff[k1][k2] = k1CSet.getCNodeDepth(k2) - k1CSet.getDepth() + allLength.get(k1)+ allLength.get(k2);
			}
				
		}
//}{$v2} = $v2_depth - $v1_depth + $length_ver{$v2} + $length_ver{$v1};
		for (CompatibleSet k1CSet : allSets) {
			int k1 = k1CSet.getNode();
			ArrayList<Integer> cNodeList = k1CSet.getCNodeList();
			for (Integer ik2 : cNodeList) {
				int k2 = ik2;
				if(d_hashTable[k1][k2]>0){
					
					double temp0 = d_hashTable[k1][k2]/(Scale-length_diff[k1][k2]);
					double temp1 = totalCount - k1CSet.getCNodeCount(k2);
					double temp2 = 1 - temp0;
					double val = temp1*Math.log(temp2) + k1CSet.getCNodeCount(k2) * Math.log(temp0);
					likelihood += val;
				}else				{
					likelihood += min;
				}
			}
		}
//		for $k1( keys %compatible_set ) 
//		{
//		#	print "Processing $k1 compatible set \n";
//			for $k2(keys %{$compatible_set{$k1}}) 
//			{
//				if(exists($d_hashtable{$k1}{$k2}))
//				{
//					double temp0 = $d_hashtable{$k1}{$k2}/($Scale-$length_diff{$k1}{$k2});
//					double temp1 = $total_count - $compatible_set{$k1}{$k2};
//					double temp2 = 1 - ($d_hashtable{$k1}{$k2}/($Scale-$length_diff{$k1}{$k2}));
//					double val = temp1*Math.log(temp2) + $compatible_set{$k1}{$k2} * Math.log(temp0);
//					$llik += $val;
//				}else				{
//					likelihood += min;
//				}
//			}
//		}

		
		return likelihood;
	}
/*	
	sub compute_paired_likelihood
	{
		# Compute the likelihood of a set of paths based on the paired k-mers that support it 
		# Need to know d_ijk, 
		# Break the node paths into compatible sets 
		# See if some compatible sets are shared across two or more haplotypes 
		# %set_paths contains the set of haplotypes for which likelihood is to be computed... 
		%d_hashtable = ();
		$num_haps = scalar(keys %set_paths);
		$genome_length = 1200;
	#	print "Number of haplotypes is $num_haps \n";
	#	print "Computing d-hasttable \n";
		compute_d_hashtable();
	#	print "D hash table done,Starting on likelihood \n";
		$like = compute_set_likelihood_using_d();
	#	print "Likelihood of set $like \n";
		return $like;
	}

	sub compute_set_likelihood_using_d
	{
		#print "In compute set \n";
		my($llik) = 0; my($min) = 0; my($max) = -9**9**9; 
		$Scale = $num_haps*$genome_length; 
		for $k1( keys %compatible_set ) 
		{
		#	print "Processing $k1 compatible set \n";
			for $k2(keys %{$compatible_set{$k1}}) 
			{
				if(exists($d_hashtable{$k1}{$k2}))
				{
					$temp0 = $d_hashtable{$k1}{$k2}/($Scale-$length_diff{$k1}{$k2});
					$temp1 = $total_count - $compatible_set{$k1}{$k2};
					$temp2 = 1 - ($d_hashtable{$k1}{$k2}/($Scale-$length_diff{$k1}{$k2}));
					$val = $temp1*log($temp2) + $compatible_set{$k1}{$k2} * log($temp0);
					$llik += $val;
				}else 
				{
					#compatible set contains a haplotype which is not supported by the set of haplotypes given, reject the set of haplotypes 	
	#				print "Set of haplotypes not compatible based on data $k1 $k2 $compatible_set{$k1}{$k2} \n";
					return $min;
				}
			}
		}
		#print "Max $max and Min $min $t1 $t2 $compatible_set{$t1}{$t2}\n";
		#print "Likelihood is $llik \n";
		return $llik;
	}
	*/

	private double[][] computeDHashTable(HaplotypeModel haplotypeModel) {
		
		double[][] d_hashTable = new double[maxNode][maxNode];
		
		for (int i = 0; i < haplotypeModel.getHaplotypeCount(); i++) {
			Haplotype haplotype = haplotypeModel.getHaplotype(i);
			ArrayList<Integer> path = new ArrayList<Integer>();
			
			HashSet<Integer> tempSet = new HashSet<>();
			for (Integer k1 : path) {
				tempSet.add(k1);
				//TODO: redo this part, need create PATH class
//				my(%temp_set) = ();
//				for $k(@one_path) 
//				{	
//					$temp_set{$k} = 1;
//				}
			}
			
			for (Integer k1 : path) {
				CompatibleSet compatibleSet = allSets.get(k1);	
				int nodeIndex = compatibleSet.getNode();
				if (tempSet.contains(nodeIndex)){
					d_hashTable[k1][nodeIndex]++;
				}
			}
		}
//			for $n1( keys %set_paths ) 
//			{
//				#pick a path in the haplotype set
//				my(@one_path) = @{$set_paths{$n1}};
//				#print "Operating on $n1 ".scalar(@one_path)." nodes in it \n";
//				my(%temp_set) = ();
//				for $k(@one_path) 
//				{	
//					$temp_set{$k} = 1;
//				}
//				for (my($i)=0;$i<scalar(@one_path);$i++)
//				{
//					for $k(keys %{$compatible_set{$one_path[$i]}})
//					{
//						if(exists($temp_set{$k})) 
//						{
//							$d_hashtable{$one_path[$i]}{$k}++;
//						}
//					}
//				}
//			}
//		}
		return null;
	}
}
