import java.io.*;
import java.lang.Math;
import java.util.Arrays;
import java.util.List;
import java.util.ArrayList;

public class oap {

    public static String[] ReadFasta(String filename) throws Exception {
        File file = new File(filename);
        BufferedReader br = new BufferedReader(new FileReader(file));
        String[] dnas = new String[]{"", ""};
        String st;
        int index = -1;
        while ((st = br.readLine()) != null) {
            if (st.charAt(0) == '>') {
                index++;
            } else {
                dnas[index] += st;
            }
        }
        return dnas;
    }

    public static int[] GetMaxIndex(int[] candidates) {
        // returns int[]{index of max, max}
        int[] ToBeReturned = new int[2];
        int[] candidates_copy = new int[candidates.length];
        for (int i = 0; i < candidates.length; i++) {
            candidates_copy[i] = candidates[i];
        }
        Arrays.sort(candidates_copy);
        int MaxScore = candidates_copy[candidates.length - 1];
        for (int i = 0; i < candidates.length; i++) {
            if (candidates[i] == MaxScore) {
                ToBeReturned[0] = i;
                ToBeReturned[1] = MaxScore;
                break;
            }
        }
        return ToBeReturned;
    }

    public static boolean CheckScore(String s, String t, int OptimalScore, int MatchScore, int SubstitutionPenalty, int GapPenalty) {
        /*
        s and t are two reconstructed strings from alignment,
        OptimalScore is the best alignment score
        returns a boolean that indicates whether s and t are indeed optimal alignment strings
        */
        if (s.length() != t.length()) {
            return false;
        }
        int score = 0;

        for (int k = 0; k < s.length(); k++) {
            if (s.charAt(k) == t.charAt(k)) {
                score += MatchScore;
            } else {
                if (s.charAt(k) == '-' | t.charAt(k) == '-') {
                    score -= GapPenalty;
                } else {
                    score -= SubstitutionPenalty;
                }
            }
        }
        return score == OptimalScore;

    }

    public static void OverLapAlignmentLinear(String s, String t, int MatchScore, int SubstitutionPenalty, int GapPenalty, boolean WriteToFile) throws Exception{
        int[][] S = new int[s.length()+1][t.length()+1];
        
        //initialization
        
        // iteration
        for (int i = 1; i < s.length()+1; i++) {
            for (int j = 1; j < t.length()+1; j++) {
                int OneMatch = (s.charAt(i-1) == t.charAt(j-1) ? 1:0) * MatchScore - (s.charAt(i-1) != t.charAt(j-1) ? 1:0) * SubstitutionPenalty;
                int[] candidates = new int[] {S[i-1][j-1] + OneMatch, S[i-1][j] - GapPenalty, S[i][j-1] - GapPenalty};
                Arrays.sort(candidates);
                S[i][j] = candidates[2];
            }
        }

        int[] LastRow = new int[t.length()+1];

        for (int j = 0; j < t.length()+1; j++) {
            LastRow[j] = S[s.length()][j]; 
        } 

        int[] MaxIndex = GetMaxIndex(LastRow);

        int best = MaxIndex[1];
        System.out.println(best);

        // traceback

        String s_reconstruct = "";
        String t_reconstruct = "";

        int i;
        int j;
        int max_i;
        int max_j;

        max_j = MaxIndex[0];
        j = max_j;
        i = s.length();

        while (j > 0) {
            int OneMatch = (s.charAt(i-1) == t.charAt(j-1) ? 1:0) * MatchScore - (s.charAt(i-1) != t.charAt(j-1) ? 1:0) * SubstitutionPenalty;

            if (S[i][j] == S[i-1][j-1] + OneMatch) {
                s_reconstruct = String.valueOf(s.charAt(i-1)) + s_reconstruct;
                t_reconstruct = String.valueOf(t.charAt(j-1)) + t_reconstruct;
                i-=1;
                j-=1;
            } else {
                if (S[i][j] == S[i-1][j] - GapPenalty) {
                    s_reconstruct = String.valueOf(s.charAt(i-1)) + s_reconstruct;
                    t_reconstruct = "-" + t_reconstruct;
                    i-=1;
                } else {
                    t_reconstruct = String.valueOf(t.charAt(j-1)) + t_reconstruct;
                    s_reconstruct = "-" + s_reconstruct;
                    j-=1;
                }
            }
  
        }

        if (!CheckScore(s_reconstruct, t_reconstruct, best, MatchScore, SubstitutionPenalty, GapPenalty)) {
            System.out.println("Something wrong with tracing");
            return;
        }
        
        if (WriteToFile) {
            File myObj = new File("oap_rslt.txt");
            FileWriter myWriter = new FileWriter("oap_rslt.txt");
            myWriter.write(String.valueOf(best));
            myWriter.write("\n");
            myWriter.write(s_reconstruct);
            myWriter.write("\n");
            myWriter.write(t_reconstruct);
            myWriter.close(); // to prevent some part of the test got flushed: https://stackoverflow.com/questions/13426142/bufferedwriter-not-writing-everything-to-its-output-file
        } else {
            System.out.println(s_reconstruct);
            System.out.println(t_reconstruct);
        }

    }

    public static void main(String[] args) throws Exception {
        String[] dnas = new String[] {"CTAAGGGATTCCGGTAATTAGACAG", "ATAGACCATATGTCAGTGACTGTGTAA"};
        OverLapAlignmentLinear(dnas[0], dnas[1], 1, 2, 2, false); // test case
        dnas = ReadFasta("rosalind_oap.txt");
        long time = System.nanoTime();
        OverLapAlignmentLinear(dnas[0], dnas[1], 1, 2, 2, true); 
        System.out.println("Alignment Time: " + (System.nanoTime() - time) / 1e9 + " s"); // benchmark time
    }

}
