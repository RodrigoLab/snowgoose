package likelihood;

public class LikelihoodUtils {

	public LikelihoodUtils(){
		
	}


	public static int hamDist(String s1, String s2) {

		if (s1.length() != s2.length()) {
			System.err.println("Different length\t"+  s1.length() +"\t"+ s2.length()) ;
			System.exit(-1);
		}
		int count = 0;
		for (int i = 0; i < s1.length(); i++) {
			if (s1.charAt(i) != s2.charAt(i) ){
				count++;
			}
		}
		return count;
	}

}
