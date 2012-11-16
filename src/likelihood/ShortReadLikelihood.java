package likelihood;

import java.util.ArrayList;
import org.apache.commons.math3.util.ArithmeticUtils;
import org.w3c.dom.Document;
import org.w3c.dom.Element;

import dr.evolution.alignment.SimpleAlignment;
import dr.evolution.sequence.Sequence;
import dr.evolution.sequence.Sequences;
import dr.evolution.tree.Tree;
import dr.inference.model.AbstractModelLikelihood;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.inference.model.Variable.ChangeType;



public class ShortReadLikelihood extends AbstractModelLikelihood {

    public static final String SHORT_READ_LIKELIHOOD = "ShortReadLikelihood";
	public static final String NAME = SHORT_READ_LIKELIHOOD;
	public static final double ERROR_RATE = 0.0107;
	public static final double C = 1e-200;
	public static final double LOG_C = Math.log(C);
	
	ArrayList<String> shortRead = new ArrayList<>();
	ArrayList<String> haplotypes = new ArrayList<>();
	

    @Override
	public Element createElement(Document d) {
        throw new RuntimeException("Not implemented yet!");
    }

    public ShortReadLikelihood() {

        super(SHORT_READ_LIKELIHOOD);
//
//        this.trialsParameter = trialsParameter;
//        this.proportionParameter = proportionParameter;
//        addVariable(trialsParameter);
//        addVariable(proportionParameter);
//        this.counts = counts;

    }

//	public ShortReadLikelihood(String name){
//		
//	}

	public ShortReadLikelihood(Sequences reads, SimpleAlignment alignment){
		super(SHORT_READ_LIKELIHOOD);
		for (int i = 0; i < alignment.getSequenceCount(); i++) {
			this.haplotypes.add(alignment.getAlignedSequenceString(i));
		}
		for (int i = 0; i < reads.getSequenceCount(); i++) {
			this.shortRead.add(reads.getSequence(i).getSequenceString());
		}

	}

	
	public ShortReadLikelihood(ArrayList<String> shortRead, ArrayList<String> haplotypes){
		super(SHORT_READ_LIKELIHOOD);
		this.shortRead = shortRead;
		this.haplotypes = haplotypes;
	}

	@Override
	public double getLogLikelihood(){

		double logLikelihood = 0;
		LikelihoodScaler liS = new LikelihoodScaler(LOG_C);
		for (String reads : shortRead) {
			int srLength = reads.length();
			double plambda = srLength * ERROR_RATE;
			double logPlambda = Math.log(plambda);

			for (String h : haplotypes) {
				int hLength = h.length();
				int noComb = hLength - srLength + 1;

				for (int j = 0; j < noComb; j++) {
					String subH = h.substring(j, j + srLength);
					int dist = LikelihoodUtils.hamDist(reads, subH);
					double logProb = -plambda + dist * logPlambda - ArithmeticUtils.factorialLog(dist);
					liS.scaleLogProb(logProb);
				}

			
			}	
			logLikelihood += liS.getLogLikelihood();
			liS.reset(LOG_C);
		}

	    return(logLikelihood);
	}
	
	/**
	 * @param shortRead the shortRead to set
	 */
	public void setShortRead(ArrayList<String> shortRead) {
		this.shortRead = shortRead;
	}

	/**
	 * @param haplotypes the haplotypes to set
	 */
	public void setHaplotypes(ArrayList<String> haplotypes) {
		this.haplotypes = haplotypes;
	}

	@Override
	public Model getModel() {
		return this;
		
	}



	@Override
	public void makeDirty() {
		// TODO Auto-generated method stub
		
	}

	@Override
	protected void handleModelChangedEvent(Model model, Object object, int index) {
		// TODO Auto-generated method stub
		
	}

	@Override
	protected void handleVariableChangedEvent(Variable variable, int index,
			ChangeType type) {
		// TODO Auto-generated method stub
		
	}

	@Override
	protected void storeState() {
		// TODO Auto-generated method stub
		
	}

	@Override
	protected void restoreState() {
		// TODO Auto-generated method stub
		
	}

	@Override
	protected void acceptState() {
		// TODO Auto-generated method stub
		
	}
	
}
/*
 * 


def hamdist(str1, str2):
    """Count the # of differences between equal length strings str1 and str2"""
    diffs = 0
    for ch1, ch2 in zip(str1, str2):
            if ch1 != ch2:
                    diffs += 1
    return diffs



def calLikelihood2(short_reads, haplotypes):
    """
    with scaling method to aviod underflow
    """
    likelihood = 0
    likelihood2 = 0
    C = 1e-200
    logC = math.log(C)
    for s, srp in enumerate(short_reads):
#        print s, "========="
        srp_length = len(srp)
        plambda = srp_length * ERROR_RATE
        log_lambda = math.log(plambda)
#        pdist = poisson(plambda)

#        pdist = poisson(srp_length * 0.99)
#        print plambda
        prob = 0
        prob2 = 0
        for h in haplotypes:

            h_length = len(h)
            no_comp = h_length - srp_length + 1

            for i in range(no_comp):
                sub_h = h[i:(i + srp_length)]
                dist = hamdist(sub_h, srp)
#                print -plambda, dist * log_lambda , math.log(math.factorial(dist))
                logProb = -plambda + dist * log_lambda - math.log(math.factorial(dist))
#                print "log:", dist, logProb, pdist.logpmf(dist)
                prob += avoid_underflow(logProb, logC)
#                prob2 += pdist.pmf(dist)
#                print Decimal(pdist.pmf(dist))
#print Decimal(1) / Decimal(7)
#            print prob * C, prob2, prob * C - prob2

#            print "prob: ", prob
#        print prob != 0
#        print math.log(1e-323)
#        if prob2 != 0:
#            print "m:", prob, math.log(prob)
#        print prob, C, prob * C,
        likelihood += (math.log(prob) + logC)
#        likelihood2 += math.log(prob2)
#        print math.log(prob) + logC, likelihood, likelihood2
#        else:
#            likelihood2 = float("-inf")
#            print "hl: ", prob2, likelihood2
#    print haplotypes
#    print likelihood, likelihood2
    return(likelihood)





def avoid_underflow(logA, logC):
    """
    A/C = e^log(A/C) = e^(logA - logC) = e^B
    """

    expB = math.exp(logA - logC)
    return expB

##########################################################################
"""
Created on Aug 21, 2012

@author: Steven Wu
"""
from scipy.stats.distributions import poisson
import math

import sys
print(sys.version)
from Bio import SeqIO
import random

import numpy.random
import os
from abi import path_utils, sequence_formatter
import sys
from Bio.Align.Applications._Clustalw import ClustalwCommandline

cwd = os.getcwd()
data_dir = path_utils.get_data_dir()
print cwd, data_dir

DNA = "ACGT"
ERROR_RATE = (1.07 / 100.0)
LOG_ERROR_RATE = math.log(ERROR_RATE)
ONE_MINUS_ERROR = 1 - ERROR_RATE

def _calLikelihood(short_reads, haplotypes):
    """
    @deprecated: 
    WILL have underflow
    """
    likelihood = 0
    for s, srp in enumerate(short_reads):
#        print s, "========="
        srp_length = len(srp)
        pdist = poisson(srp_length * 0.0075)
#        pdist = poisson(srp_length * 0.99)
#        print srp_length * 0.0075
        prob = 0
        for h in haplotypes:
            h_length = len(h)
            no_comp = h_length - srp_length + 1

            for i in range(no_comp):
                sub_h = h[i:(i + srp_length)]
                dist = hamdist(sub_h, srp)
#                if dist < 200:
#                    print dist, pdist.pmf(dist), sub_h
#                    print dist, pdist.pmf(dist), srp, "\n"
#                print  ("%d %d" % (dist, pdist.pmf(dist)))
                prob += pdist.pmf(dist)
#                print Decimal(pdist.pmf(dist))
#print Decimal(1) / Decimal(7)


#            print "prob: ", prob
#        print prob != 0
#        print math.log(1e-323)
        if prob != 0:
#            print "m:", prob, math.log(prob)
            likelihood += math.log(prob)
#            print likelihood
        else:
            likelihood = float("-inf")
#            print "hl1: ", prob, likelihood
#    print haplotypes
#    print likelihood
    return(likelihood)




def _calLikelihood3(short_reads, haplotypes):
    """
    @deprecated: 
    with some sort of threshold to avoid underflow
    """
    likelihood = 0
    for s, srp in enumerate(short_reads):
#        print s, "========="
        srp_length = len(srp)
        pdist = poisson(srp_length * 0.0075)


        half_length = math.floor(srp_length / 2)

        if half_length > 200:
            half_length = 200
        min_pdist = pdist.pmf(half_length)

#        pdist = poisson(srp_length * 0.99)
#        print srp_length * 0.0075
        prob = 0
        for h in haplotypes:
            h_length = len(h)
            no_comp = h_length - srp_length + 1

            for i in range(no_comp):
                sub_h = h[i:(i + srp_length)]
                dist = hamdist(sub_h, srp)
                if dist > 200:
                    prob += min_pdist
                else:
                    prob += pdist.pmf(dist)
#                print dist, pdist.pmf(dist), sub_h
#                print dist, pdist.pmf(dist), srp, "\n"


#            print "prob: ", prob
#        print prob != 0
#        print math.log(1e-323)
#        if prob != 0:
#            print "m:", prob, math.log(prob)
        likelihood += math.log(prob)
#            print likelihood
#        else:
#            likelihood = float("-inf")
#            print "hl3: ", prob, likelihood
#    print haplotypes
#    print likelihood
    return(likelihood)


def _test_likelihood_function_freq():
    #    short_reads = ["AACCGGTTACGT", "ACTGG"]
#    haplotypes = ["CGATGTGTTCTTGGAATCACCTACCCTAGTGAC",
##                  "AAAAAACACACACAACCGGTTACGTGTGTGTGT",
#                  "AAAAAACACACACTGCATTGGCCAAGTGTGTGT"]

#    shortread_file = "/home/sw167/Postdoc/Project_ABI/data/metaSim/testData/jt10-454.250.fna"
#    haplotype_file = "/home/sw167/Postdoc/Project_ABI/data/shorah_result/20H_10kR/jt20-454_10k_cor.popl"
#    haplotype_file = "/home/sw167/Postdoc/Project_ABI/data/metaSim/testData/jt10.fasta"

    shortread_file = "/home/sw167/Postdoc/Project_ABI/data/metaSim/testData/jt3-454.30.fna"
#    shortread_file = "/home/sw167/Postdoc/Project_ABI/data/metaSim/testData/jt3-454.l200.fna"
    haplotype_file = "/home/sw167/Postdoc/Project_ABI/data/metaSim/testData/jt3.fasta"

    shortread_list = list(SeqIO.parse(shortread_file, "fasta"))
    haplotype_list = list(SeqIO.parse(haplotype_file, "fasta"))

    short_reads = []
    haplotypes = []

    for s in shortread_list:
        short_reads.append(str(s.seq))
    for h in haplotype_list:
        haplotypes.append(str(h.seq))

    hfc = dict({'A':0, 'C':0, 'G':0, 'T':0})
    hap_fre_count = []
    srp_fre_count = []
    lens = len(haplotypes[0])

    for i in range(lens):
        hap_fre_count.append(hfc.copy())
        srp_fre_count.append(hfc.copy())

    add_prop = 1.0 / len(haplotypes)
    for h in haplotypes:
        for l in range(len(h)):
            hap_fre_count[l][h[l]] += add_prop
#            for k in hap_fre_count.keys():


    print hap_fre_count
    print srp_fre_count
    max_score = 0
    temp_input = "/home/sw167/Postdoc/Project_ABI/data/alpha.fasta"
    needle_out = "/home/sw167/Postdoc/Project_ABI/data/needle.fasta"

    temp_template = ""
    for i, seq in enumerate(haplotype_list):
        print seq.id, seq.seq
        temp_template += ">%s\n%s\n" % (seq.id, str(seq.seq))

    temp_file = open(temp_input, "w")
#    temp_file.write(temp_template)
    for i, s in enumerate(short_reads):
        temp_file = open(temp_input, "w")
        temp_file.write(">test%d\n%s\n" % (i, s))

        temp_file.close()


    clustal_cline = ClustalwCommandline("/home/sw167/Postdoc/Project_ABI/data/clustalw2",
                                        profile2=temp_input, PROFILE1=haplotype_file,
                                        OUTFILE="/home/sw167/Postdoc/Project_ABI/data/alpha_zzz",
                                        SEQUENCES=" ",
                                        ALIGN=" ",
#                                        PROFILE=" "
                                        )


    print i, clustal_cline
    clustal_cline()
#
#        needle_handle = open(needle_out + str(i), "rU")
#        for line in needle_handle:
#            if line.startswith("# Identity"):
#                score = line[line.index("(") + 1:line.index("%)")]
#                if(score > max_score):
#                    max_score = score

    print max_score
#Identity:     900/1000 (90.0%)
    print max_score



def _test_likelihood_function_allmatches():
#    short_reads = ["AACCGGTTACGT", "ACTGG"]
#    haplotypes = ["CGATGTGTTCTTGGAATCACCTACCCTAGTGAC",
##                  "AAAAAACACACACAACCGGTTACGTGTGTGTGT",
#                  "AAAAAACACACACTGCATTGGCCAAGTGTGTGT"]

#    shortread_file = "/home/sw167/Postdoc/Project_ABI/data/metaSim/testData/jt10-454.250.fna"
#    haplotype_file = "/home/sw167/Postdoc/Project_ABI/data/shorah_result/20H_10kR/jt20-454_10k_cor.popl"
#    haplotype_file = "/home/sw167/Postdoc/Project_ABI/data/metaSim/testData/jt10.fasta"

    shortread_file = "/home/sw167/Postdoc/Project_ABI/data/metaSim/testData/jt3-454.30.fna"
#    shortread_file = "/home/sw167/Postdoc/Project_ABI/data/metaSim/testData/jt3-454.l200.fna"
    haplotype_file = "/home/sw167/Postdoc/Project_ABI/data/metaSim/testData/jt3.fasta"

    shortread_list = list(SeqIO.parse(shortread_file, "fasta"))
    haplotype_list = list(SeqIO.parse(haplotype_file, "fasta"))

    short_reads = []
    haplotypes = []
    for s in shortread_list:
        short_reads.append(str(s.seq))
    for h in haplotype_list:
        haplotypes.append(str(h.seq))

#    short_reads = ["AGAG"*100, "A"*400]
#    haplotypes = ["TTCC"*250,
#                  "AGAG"*250,
#                  "AGTG"*250,
#                  "T"*1000,
#                  "A"*1000]

    print len(haplotypes)
    print len(short_reads), short_reads[0]
    print len(short_reads[0])

#    logC = math.log(1e-200)
#    testA = [1e-250, 2e-251, 2.5e-251, 3.2e-252, 2.2e-251, 1.5e-250, 3e-260, 4e-270 ]
#    R = 0
#    R2 = 0
#    for A in testA:
#
#        logA = math.log(A)
#        expB = avoid_underflow(logA, logC)
#        R += expB
#        R2 += A
##        print R, expB
#    print R * 1e-200, R2, math.log(R2), logA

#    base_likelihood = _calLikelihood(short_reads, haplotypes)
#    print(base_likelihood)

#    base_likelihood = _calLikelihood(short_reads, haplotypes)
#    print(base_likelihood)
#    base_likelihood = calLikelihood2(short_reads, haplotypes)
#    print(base_likelihood)
    base_likelihood = _calLikelihood3(short_reads, haplotypes)
    print(base_likelihood)

    likelihood_mutate = []
    for i in range(100):
        haplotypes_2 = mutates_list(haplotypes, 0.010)
#        print _calLikelihood(short_reads, haplotypes_2), \
#            calLikelihood2(short_reads, haplotypes_2), \
#        print _calLikelihood3(short_reads, haplotypes_2)

        likelihood_mutate.append(_calLikelihood3(short_reads, haplotypes_2))
        print(likelihood_mutate[i])

    result_file = open("/home/sw167/Postdoc/Project_ABI/data/likelihood3S30R0.01E.result", "w")
    result_file.write(str(base_likelihood) + "\n")
    result_file.write(str(likelihood_mutate))
    result_file.close()



def simulatePointMutation(haplotype_file, outfile, error_rate):

    haplotypes = fasta_file_to_list(haplotype_file)

    seq_len = len(haplotypes[0])
    out = ""
    noReadPerSeq = 50
    for i, s in enumerate(haplotypes):
        for j in range(noReadPerSeq):
            l = int(numpy.random.normal(400, 20))
            start = random.randint(0, seq_len - l)
            s_temp = mutate_single(s[start:(start + l)], error_rate)
#            s_temp = mutate_single(s[start:(start + l)], 0)
            out += ">" + str(i * noReadPerSeq + j) + "_" + str(start) + "\n" + s_temp + "\n"

    print outfile
    out_handle = open(outfile, "w")
    out_handle.write(out)
    out_handle.close()


def test_likelihood_function_allmatches2(haplotype_file, srp_file, error_rate):

    haplotypes = fasta_file_to_list(haplotype_file)
    short_reads = fasta_file_to_list(srp_file)

#    short_reads = short_reads[0:5]
#    haplotypes = haplotypes[0:3]
    print len(haplotypes), len(haplotypes[0]), len(short_reads), "E:", error_rate

    base_likelihood = calLikelihood2(short_reads, haplotypes)
    print("%f%s" % (base_likelihood , "\n"))

    result_file = open(data_dir + "likelihoodAllMatches" + str(error_rate) + "E.result", "w")
    result_file.write(str(base_likelihood) + "\n")

    likelihood_mutate = []
    for i in range(100):
        haplotypes_mutate = mutates_list(haplotypes, error_rate)

        likelihood_mutate.append(calLikelihood2(short_reads, haplotypes_mutate))
        print(likelihood_mutate[i])
        result_file.write(str(likelihood_mutate[i]) + "\t")
        result_file.flush()


#    result_file.write(str(likelihood_mutate))
    result_file.close()



def cal_likelihood_freq(all_prop, all_freq):


    likelihood = 0.0
    seq_lens = len(all_prop)
    for i in range(seq_lens):
        prop = all_prop[i]
        freq = all_freq[i]
#        base_likelihood += cal_likelihood_freq(hap_fre_prop[i], srp_fre_count[i])

        bin_cat = []
        non_zero_cat = []
        bin_freq = 0
        total_freq = 0
        for k, v in prop.items():
            total_freq += freq[k]
            if v == 0.0:
                bin_cat.append(k)
                bin_freq += freq[k]
            else:
                non_zero_cat.append(k)


        for k in non_zero_cat:
            prop[k] *= ONE_MINUS_ERROR
            likelihood += (freq[k] * math.log(prop[k]))
            likelihood -= math.log(math.factorial(freq[k]))
    #        print "p:", prop, likelihood

        likelihood += (bin_freq * LOG_ERROR_RATE)
        likelihood -= math.log(math.factorial(bin_freq))
#    print likelihood
    ## TODO: the factorial parts, 
#    for k in prop.keys():
#        if not bin_cat.count(k):
#
#            likelihood += math.pow(prop[k], freq[k])
#        print likelihood, likelihood / 1000
#    print freq, bin_cat, bin_freq
#    print sum(freq)
        likelihood += math.log(math.factorial(total_freq))
#    print likelihood
    return likelihood

def test_likelihood_function_freq2(haplotype_file, srp_file, error_rate):

    haplotypes = fasta_file_to_list(haplotype_file)
#    srps = fasta_file_to_list(srp_file)



    hfc = dict({'A':0, 'C':0, 'G':0, 'T':0})
    hap_fre_prop = []
    srp_fre_count = []
    seq_lens = len(haplotypes[0])

    for i in range(seq_lens):
        hap_fre_prop.append(hfc.copy())
        srp_fre_count.append(hfc.copy())

    for h in haplotypes:
        for l in range(len(h)):
            hap_fre_prop[l][h[l]] += 1
#            hap_fre_prop[l][h[l]] = math.fsum([add_prop, hap_fre_prop[l][h[l]] ])
#            print hap_fre_prop[l][h[l]], hap_fre_prop[l]
#            for k in hap_fre_prop.keys():
    no_hap_float = float(len(haplotypes))
    for h in range(len(hap_fre_prop)):
        for k in hfc.keys():
            hap_fre_prop[h][k] /= no_hap_float

#    print hap_fre_prop

    srp_dict = SeqIO.to_dict(SeqIO.parse(srp_file, "fasta"))

#    srp_dict = dict({
#                    "0_0":"AACCAACC",
#                    "1_1": "ACCAACCT",
#                    "2_5":    "GGGGGG"
#                    })
    for k, v in srp_dict.items():
        index = k.find("_") + 1
        start_pos = int(k[index:len(k)])
        seq = v.seq
        for s in seq:
            srp_fre_count[start_pos][s] += 1
            start_pos += 1
#        print start_pos - len(seq) , start_pos, srp_fre_count[(start_pos - len(seq)) : start_pos]
#        break
#        print s.id, s.seq
#        seq_list.append(str(s.seq))
#    base_likelihood = 0
#    for i in range(seq_lens):
    base_likelihood = cal_likelihood_freq(hap_fre_prop, srp_fre_count)


    #### permutation
#    print("%f%s" % (base_likelihood , "\n"))

    result_file = open(data_dir + "likelihoodFreqs" + str(error_rate) + "E.result", "w")
    result_file.write(str(base_likelihood) + "\n")

    likelihood_mutate = []
    for p in range(100):
        haplotypes_mutate = mutates_list(haplotypes, error_rate)

        hfc = dict({'A':0, 'C':0, 'G':0, 'T':0})
        hap_fre_prop = []
        for i in range(seq_lens):
            hap_fre_prop.append(hfc.copy())

        for h in haplotypes_mutate:
            for l in range(len(h)):
                hap_fre_prop[l][h[l]] += 1
        no_hap_float = float(len(haplotypes_mutate))
        for h in range(len(hap_fre_prop)):
            for k in hfc.keys():
                hap_fre_prop[h][k] /= no_hap_float



        likelihood_mutate.append(cal_likelihood_freq(hap_fre_prop, srp_fre_count))
        print(likelihood_mutate[p])
        result_file.write(str(likelihood_mutate[p]) + "\t")
        result_file.flush()


#    result_file.write(str(likelihood_mutate))
    result_file.write("")
    result_file.close()
#    short_reads = ["AGAG"*100, "A"*400]
#    haplotypes = ["TTCC"*250,
#                  "AGAG"*250,


#                  "AGTG"*250,
#                  "T"*1000,
#                  "A"*1000]

############################################################
#
#class likelihood(object):
#    """
#    classdocs
#    """
##    test_likelihood_function_allmatches()
##    test_likelihood_function_freq()
#
#    def __init__(self):
#        """
#        Constructor
#        """

if __name__ == '__main__':

    #old data
#    srp_file = "/home/sw167/Postdoc/Project_ABI/data/metaSim/testData/jt3-454.30.fna"
#    haplotype_file = "/home/sw167/Postdoc/Project_ABI/data/metaSim/testData/jt3.fasta"

    srp_file = data_dir + "srp_likelihood.fasta"
    haplotype_file = data_dir + "jt_likelihood.fasta"
#    print sys.argv, sys.argv[0], sys.argv[1]
#    simulatePointMutation(haplotype_file, "/home/sw167/Postdoc/Project_ABI/data/srp_likelihood.fasta", ERROR_RATE)

#    if len(sys.argv) is 2:
#        param = sys.argv[1]
#        if param is "0":
#            test_likelihood_function_allmatches2(haplotype_file, srp_file, 0.001)
#        if param is "1":
#            test_likelihood_function_allmatches2(haplotype_file, srp_file, 0.005)
#        if param is "2":
#            test_likelihood_function_allmatches2(haplotype_file, srp_file, ERROR_RATE)
#        if param is "3":
#            test_likelihood_function_allmatches2(haplotype_file, srp_file, 0.01)
#        if param is "4":
#            test_likelihood_function_allmatches2(haplotype_file, srp_file, 0.05)
#    #    test_likelihood_function_freq()
#    else:
#        test_likelihood_function_freq2(haplotype_file, srp_file, 0.001)

#        test_likelihood_function_freq2(haplotype_file, srp_file, 0.005)
#        test_likelihood_function_freq2(haplotype_file, srp_file, ERROR_RATE)
#        test_likelihood_function_freq2(haplotype_file, srp_file, 0.01)
#        test_likelihood_function_freq2(haplotype_file, srp_file, 0.05)
#        test_likelihood_function_allmatches2(haplotype_file, srp_file, 0.05)
    pass




    srp_file = data_dir + "srp_likelihood.fasta"
    haplotype_file = data_dir + "jt_likelihood.fasta"

#    simulatePointMutation(haplotype_file, "/home/sw167/Postdoc/Project_ABI/data/srp_likelihood.fasta", ERROR_RATE)
    ref_infile = "/home/sw167/Postdoc/Project_ABI/YT_data/p1.fasta"
    short_read_outfile = "/home/sw167/Postdoc/Project_ABI/YT_data/p1_pointMutation.fasta"


    ref_infile = "/home/sw167/Postdoc/Project_A2BI_temp/data/Stage0/haplotype.fasta"
    consensus_file = "/home/sw167/Postdoc/Project_A2BI_temp/data/Stage0/haplotype.consensus"
    short_read_outfile = "/home/sw167/Postdoc/Project_A2BI_temp/data/Stage0/srp.fasta"
    simulatePointMutation(ref_infile, short_read_outfile, ERROR_RATE)
    sequence_formatter.create_consensus_seq(ref_infile, consensus_file)

*/